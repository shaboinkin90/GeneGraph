#!/usr/bin/env python3
##################################################################
#
# A small project I made to let me visualize the rather opaque view of genetics from undergrad biology textbooks. 
#
# Queries NCBI for information on a list of genes supplied. 
# Creates a database storing all information queried.
# Creates a gml file connecting nodes of genes to nodes their 
# functions and processes within a cellular component
#
# Ideas:
#   Incorporate other sources of information
#   Filter on coding and non-coding genes
#   Filter on transcription factors, enzymes, structural/non-catalytic proteins, specific protein domains within a single gene
#   Show what tissues in a specific cellular component a specific gene is found to be expressed in
#   Show relationships with protein domains and their cellular process or function
#   Show miRNA interactions with protein coding genes (https://mirtarbase.cuhk.edu.cn/~miRTarBase/miRTarBase_2022/php/index.php)
#   Expand the usefulness of the database
#   Get information on medications that are known to have some influence on gene expression
#   Get information on externally caused interactions where research shows gene expression is influenced
#   
# Author: Daniel Kulas
#
##################################################################

import argparse
import csv
import json
import networkx as nx
import os
import sqlalchemy as sa
import zlib

from Bio import Entrez
from collections import defaultdict
from collections import deque
from copy import deepcopy
from sqlalchemy.orm import sessionmaker
from time import sleep
from tqdm import tqdm
from typing import List

# TODO: I could probably remove this
from EntrezGeneratedClasses import EntrezGeneElement

Base = sa.orm.declarative_base()

class GeneDatabase:
    def __init__(self, engine: sa.Engine, session: sessionmaker):
        self.engine = engine
        self.session = session


def open_database() -> GeneDatabase:
    engine = sa.create_engine('sqlite:///genes.db', echo=False)
    Base.metadata.create_all(engine)
    SessionMaker = sessionmaker(bind=engine)
    session = SessionMaker()
    return GeneDatabase(engine, session)


class TsvGeneInfo:
    def __init__(self, gene_id: str, gene_symbol: str, gene_desc: str, tax_name: str, gene_type: str, gene_group_id: str = None, gene_group_method: str = None):
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol
        self.gene_desc = gene_desc
        self.tax_name = tax_name
        self.gene_type = gene_type
        # These appear to be optional
        self.gene_group_id = gene_group_id
        self.gene_group_method = gene_group_method


# Database Table
class GeneEntry(Base):
    __tablename__ = "gene_table"

    gene_id = sa.Column(sa.Integer, primary_key=True)
    gene_symbol = sa.Column(sa.String)
    gene_response = sa.Column(sa.String)

    def __init__(self, gene_id, gene_symbol, gene_response):
        self.gene_id = gene_id
        self.gene_symbol = gene_symbol
        self.gene_response = gene_response
        self.ontology: dict[str, List[Ontology]] = None


class GeneDictionary:
    def __init__(self, species: str):
        self.species = species
        self.genes: dict[str, GeneEntry] = {}

class Node:
    def __init__(self, gene):
        self.gene = gene
        self.children: List[GeneEntry] = []

    def add_child(self, node):
        self.children.append(node)


def parse_tsv_gene_list(file_path: str) -> List[TsvGeneInfo]:
    gene_list: List[TsvGeneInfo] = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if not row[0].startswith('NCBI'):
                # default to None for missing fields - Arabidopsis thaliana is missing `Gene Group Identifier` and `Gene Group Method`, for example
                row += [None] * (7 - len(row))
                gene_list.append(TsvGeneInfo(*row))
    return gene_list


def compress_string(string: str):
    bytes = string.encode('utf-8')
    compressed = zlib.compress(bytes, level=zlib.Z_BEST_COMPRESSION)
    return compressed


def decompress_to_string(bytes: bytes):
    data = zlib.decompress(bytes)
    return data.decode('utf-8')


def fetch_ncbi_gene(gene_ids: List[str]):
    try:
        # NCBI allows at most 10 requets per second with API key, 3 per second otherwises.
        if Entrez.api_key == None:
            sleep(.334)
        else:
            sleep(.1)

        handle = Entrez.efetch(db="gene", id=",".join(gene_ids), rettype="xml")
        # return type is xml, Entrez turns it into a dictionary. 
        gene_records = Entrez.read(handle)
        handle.close()
        return gene_records
    except Exception as e:
        print("Failed to process genes")
        print(f"Exception: {e}")
        return None
    

def extract_gene_symbol_from_response(gene_response: dict) -> str:
    try:
        # Preferred, but not all have this set
        entrezgene_properties = gene_response.get('Entrezgene_properties', [{}])[0]
        gene_commentary_properties = entrezgene_properties.get('Gene-commentary_properties', [{}])[0]
        symbol = gene_commentary_properties.get('Gene-commentary_text')
    except IndexError:
        symbol = None
    
    if symbol == None:
        entrezgene_gene = gene_response.get('Entrezgene_gene', {})
        gene_ref = entrezgene_gene.get('Gene-ref', {})
        symbol = gene_ref.get('Gene-ref_locus', {})
        assert symbol != None, "No symbol in the NCBI response?"
    return symbol


def extract_gene_id_from_response(gene_response: dict) -> str:
    entrezgene_track = gene_response.get('Entrezgene_track-info', {})
    gene_track = entrezgene_track.get('Gene-track', {})
    id = gene_track.get('Gene-track_geneid')
    assert id != None, "No gene ID in the NCBI response?"
    return id


# Downloads NCBI data for each gene in the gene list file and saves data into the database
# If the gene already exists in the table, it does not query NCBI unless `redownload` is set to `True`
def download_genes(gene_list_file: str, redownload: bool, database: GeneDatabase):
    gene_list = parse_tsv_gene_list(gene_list_file)
    assert len(gene_list) > 0

    # Batch up requests to NCBI
    chunk_size = 100
    print("Downloading genes ...")
    for i in tqdm(range(0, len(gene_list), chunk_size)):
        chunk = gene_list[i : i + chunk_size]
        gene_ids = [gene.gene_id for gene in chunk]
        unprocessed_genes = []
        for id in gene_ids:
            exists = does_gene_exist_in_database(id, database)
            if not exists or redownload:
                unprocessed_genes.append(id)

        if len(unprocessed_genes) == 0:
            continue

        ncbi_gene_responses = fetch_ncbi_gene(gene_ids=unprocessed_genes)
        if ncbi_gene_responses == None:
            continue

        for response in ncbi_gene_responses:
            gene_id = extract_gene_id_from_response(response)
            gene_symbol = extract_gene_symbol_from_response(response)
            if does_gene_exist_in_database(gene_id, database):
                update_gene_in_database(gene_id, response, database)
            else:
                save_gene_to_database(gene_id, gene_symbol, response, database)

# Adds genes from the gene_list_file using data stored in the database, 
# and places them in a GeneDictionary object: key = gene id, value = GeneEntry
# All Homo sapien genes results in about 55GB of RAM.
def generate_gene_dictionary(gene_list_file: str, database: GeneDatabase) -> GeneDictionary:
    gene_list = parse_tsv_gene_list(gene_list_file)
    assert(len(gene_list) > 0)

    gene_dict = GeneDictionary(species=gene_list[0].tax_name)

    print("Building gene dictionary ...")
    for gene in tqdm(gene_list[:1]):
        db_gene = get_gene_entry_in_database('gene_id', gene.gene_id, database)
        if db_gene:
            gene_dict.genes[gene.gene_id] = db_gene

    return gene_dict


def does_gene_exist_in_database(gene_id: str, database: GeneDatabase) -> bool:
    query = database.session.query(GeneEntry).filter(GeneEntry.gene_id == gene_id)
    (exists,) = database.session.query(query.exists()).one()
    return exists


def get_gene_entry_in_database(column: str, value: object, database: GeneDatabase) -> GeneEntry:
    try:
        entry = database.session.query(GeneEntry).filter(getattr(GeneEntry, column) == value).one()
        entry_copy = deepcopy(entry)
        entry_copy.gene_response = json.loads(decompress_to_string(entry_copy.gene_response))
        return entry_copy
    except sa.exc.NoResultFound:
        print(f"No entry found for {column} - {value}")
        return None
    except Exception as e:
        print(f"Problem querying for {column} - {value}.")
        print(e)
        return None


def save_gene_to_database(gene_id: str, gene_symbol: str, gene_response: dict, database: GeneDatabase):
    # Compress the response when saving in the database, otherwise you're looking at 15GB+ of text for Homo sapiens genes.
    compressed_gene_response = compress_string(json.dumps(gene_response))
    entry = GeneEntry(gene_id, gene_symbol, compressed_gene_response)
    database.session.add(entry)
    database.session.commit()


def update_gene_in_database(gene_id: str, gene_response: dict, database: GeneDatabase):
    gene = get_gene_entry_in_database('gene_id', gene_id, database)
    compressed_gene_response = compress_string(json.dumps(gene_response))
    gene.gene_response = compressed_gene_response
    database.session.commit()


class OntologyData:
    def __init__(self, pre_text, anchor, post_text):
        self.pre_text = pre_text
        self.anchor = anchor
        self.post_text = post_text

class Ontology:
    def __init__(self, components=List[OntologyData], functions=List[OntologyData], processes=List[OntologyData]):
        self.components = components
        self.functions = functions
        self.processes = processes

def get_ontology_data(comment: dict, match: str) -> List[OntologyData]:

    get_fields = lambda comment, match: \
       comment['Gene-commentary_source'][0].get(match)
    
    data: List[OntologyData] = []

    for ontology_comment in comment.get('Gene-commentary_comment', []):
        pre_text = get_fields(ontology_comment, 'Other-source_pre-text')
        anchor = get_fields(ontology_comment, 'Other-source_anchor')
        post_text = get_fields(ontology_comment, 'Other-source_post-text') 
        data.append(OntologyData(pre_text, anchor, post_text))

    return data


def get_ontology(properties: dict):
    for prop in properties:
        if prop.get('Gene-commentary_heading') == "GeneOntology":
            prop_comments = prop.get('Gene-commentary_comment', [])
            for comment in prop_comments:
                label = comment['Gene-commentary_label']
                if label == 'Component':
                    components = get_ontology_data(comment, 'Component')
                elif label == 'Function':
                   functions = get_ontology_data(comment, 'Function')
                elif label == 'Process':
                    processes = get_ontology_data(comment, 'Process')
    return Ontology(components, functions, processes)


def add_edge_to_node(graph: nx.Graph, new_node: str, node_desc: str, node_color: str, from_node: str) -> bool:
    if not graph.has_node(new_node):
        graph.add_node(new_node, color=node_color)
    if graph.has_edge(from_node, new_node):
        return False
    graph.add_edge(from_node, new_node, label=node_desc)
    return True


def generate_cell_component_graph_data(genes:List[GeneEntry], separate_graphs=False):
    assert(genes != None and len(genes) >= 1)
    print("Construction cell component graph")

    component_tree = {}
    # Build trees where roots are cellular components and sub nodes are genes that exist in those components
    component_tree = defaultdict(list)

    print("Preprocess...")
    for gene in tqdm(genes):
        properties = gene.gene_response.get('Entrezgene_properties', [])
        # FIXME: move all additional data gathering from the original response after querying from the table
        # Probably save in database for quicker access
        gene.ontology = get_ontology(properties)
        for component in gene.ontology.components:
            # anchor is the cellular component, check if there's already a root in the tree
            if component.anchor in component_tree:
                # Add gene to the cell component root
                # A single gene can have duplicate "cell component" fields, but these just have a different value for post-text evidence.
                if gene not in component_tree[component.anchor]:
                    component_tree[component.anchor].append(gene)
            else:
                component_tree[component.anchor].append(gene)

    ordered_component_tree = dict(sorted(component_tree.items(), key=lambda item: len(item[1]), reverse=True))
    
    results_dir = "component_graph_results"
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)

    graph: nx.graph = None
    if not separate_graphs:
        graph = nx.Graph()

    node_color_component = "yellow"
    node_color_gene = "blue"
    node_color_function = "purple"
    node_color_process = "orange"
    
    number_of_genes_file = open(f"{results_dir}/number_of_genes_per_component.txt", mode='w')
    number_of_genes_file.write("Cellular component | # of genes\n")
    print("Constructing..")
    for node in tqdm(ordered_component_tree.items()):
        if separate_graphs:
            graph = nx.Graph()

        # The component
        component_root_node = node[0]
        graph.add_node(component_root_node, color=node_color_component)

        # The list of genes at this component
        genes: List[GeneEntry] = node[1]

        number_of_genes_file.write(f"{component_root_node} | {len(genes)}\n")

        # child = gene
        for gene in genes:
            # Genes can have different cellular components. To avoid genes forming edges to different component, prepend the component to the gene symbol string 
            if separate_graphs:
                gene_name_node = gene.gene_symbol
            else:
                gene_name_node = f"{component_root_node}_{gene.gene_symbol}"
            
            add_edge_to_node(graph=graph, new_node=gene_name_node, node_desc="", node_color=node_color_gene, from_node=component_root_node)
            
            for function in gene.ontology.functions:
                if separate_graphs:
                    func_node_name = function.anchor
                else:
                    func_node_name = f"{component_root_node}_{function.anchor}"
                add_edge_to_node(graph=graph, new_node=func_node_name, node_desc=function.pre_text, node_color=node_color_function, from_node=gene_name_node)

            for process in gene.ontology.processes:
                if separate_graphs:
                    process_node_name = process.anchor
                else:
                    process_node_name = f"{component_root_node}_{process.anchor}"
                add_edge_to_node(graph=graph, new_node=process_node_name, node_desc=process.pre_text, node_color=node_color_process, from_node=gene_name_node)

        if separate_graphs:
            component_root_node = component_root_node.replace('/','_').replace(':','_')
            filename = f"{component_root_node}_testcomponent_graph.gml"
            nx.write_gml(graph, f"{results_dir}/{filename}")

    number_of_genes_file.close()

    if not separate_graphs:
        nx.write_gml(graph, f"{results_dir}/testcomponent_graph.gml")


def get_interacting_genes(gene: GeneEntry) -> List[str]:
    interacting_gene_ids = []
    comments = gene.gene_response.get('Entrezgene_comments', [])
    for comment in comments:
        if comment.get('Gene-commentary_heading', {}) == 'Interactions':
            for interaction_comment in comment.get('Gene-commentary_comment', []):
                for sub_comment in interaction_comment.get('Gene-commentary_comment', []):
                    for source in sub_comment.get('Gene-commentary_source', []):
                        if source.get('Other-source_src', {}).get('Dbtag', {}).get('Dbtag_db', {}) == 'GeneID':
                            gene_id = source.get('Other-source_src', {}).get('Dbtag', {}).get('Dbtag_tag', {}).get('Object-id', {}).get('Object-id_id', {})
                            interacting_gene_ids.append(gene_id)

    return interacting_gene_ids


# TODO: Double check this implementation 
def generate_gene_interaction_graph(gene_dictionary: GeneDictionary, database: GeneDatabase):
    graph = nx.Graph()
    original = deepcopy(gene_dictionary)
    for item in original.genes.items():
            
        gene_entry = item[1]
        root = Node(gene=gene_entry)
        queue = deque([root])
        
        root_sym = gene_entry.gene_symbol
        curr_sym = root_sym

        while queue:
            current_node = queue.popleft()
            interacting_gene_ids = get_interacting_genes(current_node.gene)
            curr_sym = extract_gene_symbol_from_response(current_node.gene.gene_response)

            if not graph.has_node(curr_sym):
                graph.add_node(curr_sym, color='green')
            # This takes so long if you try to link up all interactions..
            for gene_id in interacting_gene_ids[:3]:
                interacting_gene = None
                if not gene_dictionary.genes.get(gene_id):
                    interacting_gene = get_gene_entry_in_database(column='gene_id', value=gene_id, database=database)
                    if interacting_gene == None:
                        print(f'gene id {gene_id} not in database..skipping')
                        continue
                    # Cache to avoid requerying if another gene interacts with this one when the gene_dictionary does not contain all genes
                    gene_dictionary.genes[interacting_gene.gene_id] = interacting_gene
                else:
                    interacting_gene = gene_dictionary.genes.get(gene_id)
                print(f"{interacting_gene.gene_symbol} interacts with {curr_sym}")
                added_node = add_edge_to_node(graph=graph, new_node=interacting_gene.gene_symbol, node_desc="", node_color='white', from_node=curr_sym)
                if not added_node:
                    continue
                child_node = Node(interacting_gene)
                current_node.add_child(child_node)
                queue.append(child_node)

    nx.write_gml(graph, "interactions_graph.gml")


def main():
    parser = argparse.ArgumentParser(description="""
    This script does reads in a tsv formatted list of genes for a species as supplied by NCBI and queries NCBI for general information about the gene,
    saves the results in a local database, then creates a node-graph showing relationships between genes, the cellular component 
    it is found to be expressed in, their biochemcial functions, and what biochemical processes it is involved in.
        
    This script allows you supply your own NCBI API key. NCBI allows more requests per second with the key (10 verses 3 without).

    For NCBI Entrez library information, see: https://www.ncbi.nlm.nih.gov/books/NBK25497/

    NCBI API keys are access via environment variables:
        NCBI: os.getenv("NCBI_API_KEY")
    Change the "NCBI_API_KEY" string to whatever you already have if required.

    The Entrez utility also requires you to specify your email address. 
    The script accesses your email through this environment variable
        os.getenv("NCBI_EMAIL_ADDRESS")

    A SQLite database will be created to cache the results queried from NCBI so as to not require requerying. 
    If a gene is not found in the database, then the script will fetch for it.
    
    As an example to acquire an input tsv file from NCBI, see: https://www.ncbi.nlm.nih.gov/datasets/gene/taxon/9606/
    In the table, check the box next to "Gene ID", then hit the download button and select "Download Table". Supply this file to the script.
    """)
    parser.add_argument("--tsv_gene_list_file", type=str, required=True, help="Supply the list of genes as downloaded from NCBI you wish to build a graph of - in NCBI TSV format")
    parser.add_argument("--redownload", action='store_true', help="Set to redownload the response from NCBI and resave in the database")
    parser.add_argument("--skip_download", action='store_true', help="Set to bypass downloading from NCBI and attempt to query directly from the database.")
    parser.add_argument("--separate_graphs", action='store_true', help="Set to create individual gml files for each cellular component instead of having all cellular component graphs in one gml")
    parser.add_argument("--option", type=str, required=True, help="Specify what the script will do. Supported options: `interaction` or `component_graph")
    args = parser.parse_args()

    ncbi_api_key = os.getenv("NCBI_API_KEY")
    email_address = os.getenv("NCBI_EMAIL_ADDRESS")

    if ncbi_api_key == None:
        print("No NCBI api key supplied. Requests will be limited to 3 per second.")
    else:
        Entrez.api_key = ncbi_api_key

    if email_address == None:
        print("No email for NCBI supplied. Will not issue a query.")
        exit()
    else:
        Entrez.email = email_address
        Entrez.tool = f"Gene Graph"
    
    if not os.path.exists(args.tsv_gene_list_file):
        print(f"{args.tsv_gene_list_file} not found")
        exit()
    
    gene_database = open_database()
    if not args.skip_download:
        download_genes(gene_list_file=args.tsv_gene_list_file, redownload=args.redownload, database=gene_database)

    gene_dictionary = generate_gene_dictionary(args.tsv_gene_list_file, gene_database)

    if args.option == 'interactions':
        # I don't like the implementation of this. The interactions can also list proteins instead of a specific gene.
        # Plus I don't think you can actually infer anything from this as is other than "gee that's a lot"
        generate_gene_interaction_graph(gene_dictionary, gene_database)
    
    if args.option == 'component_graph':
        generate_cell_component_graph_data(gene_dictionary.genes.values(), separate_graphs=args.separate_graphs)

if __name__ == '__main__':
    main()
