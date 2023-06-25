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
#   Utilize interactions:
            # <Gene-commentary>
            #     <Gene-commentary_type value="generif">18</Gene-commentary_type>
            #     <Gene-commentary_heading>Interactions</Gene-commentary_heading> 
            #     <Gene-commentary_comment>
# Author: Daniel Kulas
#
##################################################################

import argparse
import concurrent.futures
import csv
import json
import networkx as nx
import os
import pickle
import sqlalchemy as sa
import zlib

from Bio import Entrez
from time import sleep
from typing import List

from EntrezGeneratedClasses import EntrezGeneElement

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


def parse_tsv_gene_list(file_path) -> List[TsvGeneInfo]:
    gene_list = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            if not row[0].startswith('NCBI'):
                # default to None for missing fields - Arabidopsis thaliana is missing `Gene Group Identifier` and `Gene Group Method`, for example
                row += [None] * (7 - len(row))
                gene_list.append(TsvGeneInfo(*row))
    return gene_list


def compress_string(string):
    bytes = string.encode('utf-8')
    compressed = zlib.compress(bytes, level=zlib.Z_BEST_COMPRESSION)
    return compressed


def decompress_to_string(bytes):
    data = zlib.decompress(bytes)
    return data.decode('utf-8')


def fetch_ncbi_gene(gene_ids):
    print(f"Downloading genes")
    try:
        # NCBI allows at most 10 requets per second with API key, 3 per second otherwises.
        if Entrez.api_key == None:
            sleep(.334)
        else:
            sleep(.1)

        handle = Entrez.efetch(db="gene", id=",".join(gene_ids), rettype="xml")
        gene_records = Entrez.read(handle)
        handle.close()
        return gene_records
    except Exception as e:
        print("Failed to process genes")
        print(f"Exception: {e}")
        return None
    

def process_genes(gene_list_file, redownload, database) -> List[EntrezGeneElement]:
    gene_list = parse_tsv_gene_list(gene_list_file)
    assert(len(gene_list) > 0)

    species = gene_list[0].tax_name
    if os.path.exists(f"{species}_gene_cache.pkl") and not redownload:
        return depicklefy(f"{species}_gene_cache.pkl")
    
    processed_genes = []
    
    # Batch up requests to NCBI, though unnessary if it's already saved in the database
    chunk_size = 100
    for i in range(0, len(gene_list), chunk_size):
        # Uncomment and change the number to start where you left off at if you're downloading *EVERYTHING*
        # I could probably make a new .tsv file leaving only what is left to be queried since there's no real ordering in the database, and use that 
        # to pick up where it left off at, or make a temporary file listing just the gene IDs remaining. 
        # That's a more elegant than me manually changing the index. But changing the number is less effort and time is not of the essence, so, change the number 
        # idx = i + 130480
        idx = i

        chunk = gene_list[idx : idx + chunk_size]
        gene_ids = [gene.gene_id for gene in chunk]
        unprocessed_genes = []
        for id in gene_ids:
            db_gene = get_gene_in_database(gene_id=id, database=database)
            if not db_gene or redownload:
                unprocessed_genes.append(id)
            else:
                decompressed_response = decompress_to_string(db_gene.gene_response)
                obj = json.loads(decompressed_response)
                processed_genes.append(obj)

        if len(unprocessed_genes) != 0:
            ncbi_gene_responses = fetch_ncbi_gene(gene_ids=unprocessed_genes)
            if ncbi_gene_responses == None:
                continue

            for response in ncbi_gene_responses:
                gene = EntrezGeneElement.from_dict(response)
                # Compress the response when saving in the database, otherwise you're looking at 15GB+ of text.
                compressed_gene = compress_string(json.dumps(response))
                save_gene_to_database(gene=gene, ncbi_raw_response=compressed_gene, database=database)
                processed_genes.append(response)

        print(f"Processed {idx + chunk_size} genes out of {len(gene_list)}")

    # Three FIXME's here that would be nice to address if I were to expand on this:
    #   1. The memory usage is pretty ridiclious. On my system, I'm using over 80GB of RAM. 
    #   2. Converting the dictionary into it's objects takes quite a long time.
    #   3. It'd be faster to parse the raw dictionaries instead of shoving the data into objects. 
    #      I did not simply due to the large number of nested entries and lists within the response.
    #      Example: 
    #           genes[x]['Entrezgene_properties'][2]['Gene-commentary_comment'][1]['Gene-commentary_comment'][0]['Gene-commentary_source'][0]['Other-source_anchor']
    #      Or maybe keep the json response compressed until I need to start generating the graph.

    print("Converting..")
    genes = None
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(EntrezGeneElement.from_dict, gene) for gene in processed_genes]
        genes = [future.result() for future in concurrent.futures.as_completed(futures) if future.result()]

    # Cache results
    picklefy(genes, f"{species}_gene_cache.pkl")
    return genes

def picklefy(data, file_name):
    with open(file_name, "wb") as file:
        print(f"picklefying {file_name}")
        pickle.dump(data, file)


def depicklefy(file_name):
    with open(file_name, "rb") as file:
        print(f"depicklefying {file_name}")
        return pickle.load(file)

def open_database():
    engine = sa.create_engine('sqlite:///gene.db', echo=False)
    metadata = sa.MetaData()

    gene_table = None

    ins = sa.inspect(engine)
    if not ins.dialect.has_table(engine.connect(), 'gene_table'):
        gene_table = sa.Table('gene_table', metadata,
            sa.Column('gene_id', sa.Integer, primary_key=True),
            sa.Column('gene_symbol', sa.String),
            sa.Column('gene_name', sa.String),          # NCBI response -> Entrezgene_properties -> Gene-commentary_properties -> "Official Full name" -> Gene-commentary_text 
            sa.Column('gene_summary', sa.String),       # NCBI response, if any
            sa.Column('gene_response', sa.String))      # Stores the complete NCBI response 
        gene_table.create(engine)
    else:
        gene_table = sa.Table("gene_table", metadata, autoload_with=engine)

    return {'engine': engine, \
            'gene_table': gene_table}


def get_gene_name_and_symbol(gene: EntrezGeneElement):
    gene_name = None
    gene_symbol = None

    if gene.entrezgene_properties:
        for property in gene.entrezgene_properties:
            if property.gene_commentary_label == "Nomenclature":
                for sub_property in property.gene_commentary_properties:
                    if sub_property.gene_commentary_label == "Official Symbol":
                        gene_symbol = sub_property.gene_commentary_text
                    if sub_property.gene_commentary_label == "Official Full Name":
                        gene_name = sub_property.gene_commentary_text

    # Some genes don't have the "Nomenclature" field, so grab this information in another location
    if gene_symbol == None:
        gene_symbol = gene.entrezgene_gene.gene_ref.gene_ref_locus
    if gene_name == None:
        gene_name = gene.entrezgene_gene.gene_ref.gene_ref_desc
    if gene_name == None:
        gene_name = gene_symbol

    return {'name': gene_name, 'symbol': gene_symbol}


def get_gene_in_database(gene_id, database):
    engine = database['engine']
    gene_table = database['gene_table']
    with engine.connect() as conn:
        query = sa.select(gene_table).where(gene_table.c.gene_id == gene_id)
        result = conn.execute(query)
        return result.fetchone()


def save_gene_to_database(gene, ncbi_raw_response, database):
    engine = database['engine']
    gene_table = database['gene_table']

    gene_id = gene.entrezgene_track_info.gene_track.gene_track_geneid
    result = get_gene_name_and_symbol(gene)
    gene_name = result['name']
    gene_symbol = result['symbol']
    gene_summary = gene.entrezgene_summary

    with engine.connect() as conn:
        query = sa.select(gene_table).where(gene_table.c.gene_id == gene_id)
        result = conn.execute(query)
        row = result.fetchone()
        entry = None
        if row:
            entry = sa.update(gene_table).where(gene_table.c.gene_id == gene_id).values(gene_symbol=gene_symbol, \
                                                        gene_name=gene_name, gene_summary=gene_summary, gene_response=ncbi_raw_response)
        else:
            entry = sa.insert(gene_table).values(gene_id=gene_id, gene_symbol=gene_symbol, \
                                                        gene_name=gene_name, gene_summary=gene_summary, gene_response=ncbi_raw_response)
        conn.execute(entry)
        conn.commit()


def get_geneontology_data_pair(gene, node_type):
    try:
        commentaries =  [c for c in gene.entrezgene_properties              if c.gene_commentary_heading == 'GeneOntology']
        comments =      [c for c in commentaries[0].gene_commentary_comment if c.gene_commentary_label == node_type]
        new_comments =  [c for comment in comments for c in comment.gene_commentary_comment if c.gene_commentary_source]
        node = new_comments[0].gene_commentary_source[0].other_source_anchor
        desc = new_comments[0].gene_commentary_source[0].other_source_pre_text.replace('_', ' ')
        return {"node": node, "description": desc}
    except:
        return None
    

def add_geneontology_to_graph(graph, gene_name_node, property, node_color, prepend=None):
    nodes_added = []
    for comment in property.gene_commentary_comment:
        if len(comment.gene_commentary_source) > 1:
            print("more than 1 gene_commentary_source?")
            breakpoint()

        node_name = comment.gene_commentary_source[0].other_source_anchor
        desc = comment.gene_commentary_source[0].other_source_pre_text.replace('_', ' ')

        if prepend:
            node_name = f"{prepend}_{node_name}"
        add_edge_to_node(graph, node_name, gene_name_node, desc, node_color)
        nodes_added.append(node_name)

    return nodes_added
    

def add_edge_to_node(graph, node, node_name, node_desc, node_color):
    if not graph.has_node(node):
        graph.add_node(node, color=node_color)
    graph.add_edge(node_name, node, label=node_desc)


def generate_graph_data(genes=List[EntrezGeneElement], separate_graphs=False, root_node_type="Component"):
    assert(genes != None and len(genes) >= 1)
    print("Making graph")
    # Build a dictionary where each key is a unique cellular component, and values are genes that 
    # have this cellular component listed in the GeneOntology section of the NCBI data 
    unordered_trees = {}
    for gene in genes:
        if gene.entrezgene_properties == None or len(gene.entrezgene_properties) == 1:
            continue

        gene_name = gene.entrezgene_gene.gene_ref.gene_ref_locus
        if gene_name == None:
            continue

        # Input data for GeneOntology section lists out: "Component", "Function", "Process"
        data_pair = get_geneontology_data_pair(gene, node_type=root_node_type)
        if data_pair == None:
            continue

        node = data_pair['node']
        if node in unordered_trees:
            unordered_trees[node].append(gene)
        else:
            unordered_trees[node] = [gene]
        
    # Not required for the graph but the seeing the number of genes in a specific component is interesting
    component_tree = dict(sorted(unordered_trees.items(), key=lambda item: len(item[1]), reverse=True))

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

    # I could probably generate various statistics here from the rest of the NCBI data, though I'm not sure what exactly right now.

    number_of_genes_file = open("component_graph_results/number_of_genes_per_component.txt", mode='w')
    number_of_genes_file.write("Cellular component - # of genes\n")
    for node in component_tree.items():
        if separate_graphs:
            graph = nx.Graph()

        root_node = node[0]
        graph.add_node(root_node, color=node_color_component)

        children = node[1]
        number_of_genes_file.write(f"{root_node} - {len(children)}\n")
        for child in children:
            gene_name_node = child.entrezgene_gene.gene_ref.gene_ref_locus

            graph.add_node(f"{gene_name_node}", color=node_color_gene)
            graph.add_edge(root_node, f"{gene_name_node}")

            for commentary in child.entrezgene_properties:
                if commentary.gene_commentary_heading == 'GeneOntology':
                    for comment in commentary.gene_commentary_comment:
                        if comment.gene_commentary_label == 'Function':
                             add_geneontology_to_graph(graph=graph, gene_name_node=f"{gene_name_node}", property=comment, node_color=node_color_function, prepend=root_node)
                        elif comment.gene_commentary_label == 'Process':
                            add_geneontology_to_graph(graph=graph, gene_name_node=f"{gene_name_node}", property=comment, node_color=node_color_process, prepend=root_node)

        if separate_graphs:
            root_node = root_node.replace('/','_')
            root_node = root_node.replace(':','_')
            filename = f"{root_node}_component_graph.gml"
            nx.write_gml(graph, f"{results_dir}/{filename}")

    number_of_genes_file.close()

    if not separate_graphs:
        nx.write_gml(graph, f"{results_dir}/component_graph.gml")

    return component_tree


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
    parser.add_argument("--separate_graphs", action='store_true', help="Set to create individual gml files for each cellular component instead of having all cellular component graphs in one gml")
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

    genes = process_genes(gene_list_file=args.tsv_gene_list_file, redownload=args.redownload, database=gene_database)
    generate_graph_data(genes=genes, separate_graphs=args.separate_graphs)


if __name__ == '__main__':
    main()
