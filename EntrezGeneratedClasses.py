# Generated from https://app.quicktype.io/
# Language: python3.6
# classes only turned off,
# transform property names to be pythonic enabled
# Make all properties optional All "other" properties enabled
# Based on ttn.json 

from typing import Optional, Any, List, TypeVar, Type, cast, Callable


T = TypeVar("T")


def from_str(x: Any) -> str:
    assert isinstance(x, str)
    return x


def from_none(x: Any) -> Any:
    assert x is None
    return x


def from_union(fs, x):
    for f in fs:
        try:
            return f(x)
        except:
            pass
    assert False


def to_class(c: Type[T], x: Any) -> dict:
    assert isinstance(x, c)
    return cast(Any, x).to_dict()


def from_list(f: Callable[[Any], T], x: Any) -> List[T]:
    assert isinstance(x, list)
    return [f(y) for y in x]


class SeqID:
    seq_id_gi: Optional[str]

    def __init__(self, seq_id_gi: Optional[str]) -> None:
        self.seq_id_gi = seq_id_gi

    @staticmethod
    def from_dict(obj: Any) -> 'SeqID':
        assert isinstance(obj, dict)
        seq_id_gi = from_union([from_str, from_none], obj.get("Seq-id_gi"))
        return SeqID(seq_id_gi)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_id_gi is not None:
            result["Seq-id_gi"] = from_union([from_str, from_none], self.seq_id_gi)
        return result


class PurpleSeq:
    seq_id: Optional[SeqID]

    def __init__(self, seq_id: Optional[SeqID]) -> None:
        self.seq_id = seq_id

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleSeq':
        assert isinstance(obj, dict)
        seq_id = from_union([SeqID.from_dict, from_none], obj.get("Seq-id"))
        return PurpleSeq(seq_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_id is not None:
            result["Seq-id"] = from_union([lambda x: to_class(SeqID, x), from_none], self.seq_id)
        return result


class SeqIntervalStrand:
    na_strand: Optional[str]

    def __init__(self, na_strand: Optional[str]) -> None:
        self.na_strand = na_strand

    @staticmethod
    def from_dict(obj: Any) -> 'SeqIntervalStrand':
        assert isinstance(obj, dict)
        na_strand = from_union([from_str, from_none], obj.get("Na-strand"))
        return SeqIntervalStrand(na_strand)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.na_strand is not None:
            result["Na-strand"] = from_union([from_str, from_none], self.na_strand)
        return result


class PurpleSeqInterval:
    seq_interval_from: Optional[str]
    seq_interval_to: Optional[str]
    seq_interval_strand: Optional[SeqIntervalStrand]
    seq_interval_id: Optional[PurpleSeq]

    def __init__(self, seq_interval_from: Optional[str], seq_interval_to: Optional[str], seq_interval_strand: Optional[SeqIntervalStrand], seq_interval_id: Optional[PurpleSeq]) -> None:
        self.seq_interval_from = seq_interval_from
        self.seq_interval_to = seq_interval_to
        self.seq_interval_strand = seq_interval_strand
        self.seq_interval_id = seq_interval_id

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleSeqInterval':
        assert isinstance(obj, dict)
        seq_interval_from = from_union([from_str, from_none], obj.get("Seq-interval_from"))
        seq_interval_to = from_union([from_str, from_none], obj.get("Seq-interval_to"))
        seq_interval_strand = from_union([SeqIntervalStrand.from_dict, from_none], obj.get("Seq-interval_strand"))
        seq_interval_id = from_union([PurpleSeq.from_dict, from_none], obj.get("Seq-interval_id"))
        return PurpleSeqInterval(seq_interval_from, seq_interval_to, seq_interval_strand, seq_interval_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_interval_from is not None:
            result["Seq-interval_from"] = from_union([from_str, from_none], self.seq_interval_from)
        if self.seq_interval_to is not None:
            result["Seq-interval_to"] = from_union([from_str, from_none], self.seq_interval_to)
        if self.seq_interval_strand is not None:
            result["Seq-interval_strand"] = from_union([lambda x: to_class(SeqIntervalStrand, x), from_none], self.seq_interval_strand)
        if self.seq_interval_id is not None:
            result["Seq-interval_id"] = from_union([lambda x: to_class(PurpleSeq, x), from_none], self.seq_interval_id)
        return result


class SeqLOCMixSeqLOCInt:
    seq_interval: Optional[PurpleSeqInterval]

    def __init__(self, seq_interval: Optional[PurpleSeqInterval]) -> None:
        self.seq_interval = seq_interval

    @staticmethod
    def from_dict(obj: Any) -> 'SeqLOCMixSeqLOCInt':
        assert isinstance(obj, dict)
        seq_interval = from_union([PurpleSeqInterval.from_dict, from_none], obj.get("Seq-interval"))
        return SeqLOCMixSeqLOCInt(seq_interval)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_interval is not None:
            result["Seq-interval"] = from_union([lambda x: to_class(PurpleSeqInterval, x), from_none], self.seq_interval)
        return result


class GeneCommentarySeq:
    seq_loc_int: Optional[SeqLOCMixSeqLOCInt]

    def __init__(self, seq_loc_int: Optional[SeqLOCMixSeqLOCInt]) -> None:
        self.seq_loc_int = seq_loc_int

    @staticmethod
    def from_dict(obj: Any) -> 'GeneCommentarySeq':
        assert isinstance(obj, dict)
        seq_loc_int = from_union([SeqLOCMixSeqLOCInt.from_dict, from_none], obj.get("Seq-loc_int"))
        return GeneCommentarySeq(seq_loc_int)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_loc_int is not None:
            result["Seq-loc_int"] = from_union([lambda x: to_class(SeqLOCMixSeqLOCInt, x), from_none], self.seq_loc_int)
        return result


class TentacledGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_seqs: Optional[List[GeneCommentarySeq]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_label: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_seqs: Optional[List[GeneCommentarySeq]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_seqs = gene_commentary_seqs

    @staticmethod
    def from_dict(obj: Any) -> 'TentacledGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_seqs = from_union([lambda x: from_list(GeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        return TentacledGeneCommentaryComment(gene_commentary_type, gene_commentary_label, gene_commentary_accession, gene_commentary_version, gene_commentary_seqs)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        return result


class FluffyGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_comment: Optional[List[TentacledGeneCommentaryComment]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_comment: Optional[List[TentacledGeneCommentaryComment]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_comment = gene_commentary_comment

    @staticmethod
    def from_dict(obj: Any) -> 'FluffyGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_comment = from_union([lambda x: from_list(TentacledGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        return FluffyGeneCommentaryComment(gene_commentary_type, gene_commentary_heading, gene_commentary_accession, gene_commentary_version, gene_commentary_comment)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(TentacledGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        return result


class PubPmid:
    pub_med_id: Optional[str]

    def __init__(self, pub_med_id: Optional[str]) -> None:
        self.pub_med_id = pub_med_id

    @staticmethod
    def from_dict(obj: Any) -> 'PubPmid':
        assert isinstance(obj, dict)
        pub_med_id = from_union([from_str, from_none], obj.get("PubMedId"))
        return PubPmid(pub_med_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.pub_med_id is not None:
            result["PubMedId"] = from_union([from_str, from_none], self.pub_med_id)
        return result


class GeneCommentaryRef:
    pub_pmid: Optional[PubPmid]

    def __init__(self, pub_pmid: Optional[PubPmid]) -> None:
        self.pub_pmid = pub_pmid

    @staticmethod
    def from_dict(obj: Any) -> 'GeneCommentaryRef':
        assert isinstance(obj, dict)
        pub_pmid = from_union([PubPmid.from_dict, from_none], obj.get("Pub_pmid"))
        return GeneCommentaryRef(pub_pmid)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.pub_pmid is not None:
            result["Pub_pmid"] = from_union([lambda x: to_class(PubPmid, x), from_none], self.pub_pmid)
        return result


class FluffySeq:
    seq_id: Optional[SeqID]

    def __init__(self, seq_id: Optional[SeqID]) -> None:
        self.seq_id = seq_id

    @staticmethod
    def from_dict(obj: Any) -> 'FluffySeq':
        assert isinstance(obj, dict)
        seq_id = from_union([SeqID.from_dict, from_none], obj.get("Seq-id"))
        return FluffySeq(seq_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_id is not None:
            result["Seq-id"] = from_union([lambda x: to_class(SeqID, x), from_none], self.seq_id)
        return result


class FluffySeqInterval:
    seq_interval_from: Optional[str]
    seq_interval_to: Optional[str]
    seq_interval_strand: Optional[SeqIntervalStrand]
    seq_interval_id: Optional[FluffySeq]

    def __init__(self, seq_interval_from: Optional[str], seq_interval_to: Optional[str], seq_interval_strand: Optional[SeqIntervalStrand], seq_interval_id: Optional[FluffySeq]) -> None:
        self.seq_interval_from = seq_interval_from
        self.seq_interval_to = seq_interval_to
        self.seq_interval_strand = seq_interval_strand
        self.seq_interval_id = seq_interval_id

    @staticmethod
    def from_dict(obj: Any) -> 'FluffySeqInterval':
        assert isinstance(obj, dict)
        seq_interval_from = from_union([from_str, from_none], obj.get("Seq-interval_from"))
        seq_interval_to = from_union([from_str, from_none], obj.get("Seq-interval_to"))
        seq_interval_strand = from_union([SeqIntervalStrand.from_dict, from_none], obj.get("Seq-interval_strand"))
        seq_interval_id = from_union([FluffySeq.from_dict, from_none], obj.get("Seq-interval_id"))
        return FluffySeqInterval(seq_interval_from, seq_interval_to, seq_interval_strand, seq_interval_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_interval_from is not None:
            result["Seq-interval_from"] = from_union([from_str, from_none], self.seq_interval_from)
        if self.seq_interval_to is not None:
            result["Seq-interval_to"] = from_union([from_str, from_none], self.seq_interval_to)
        if self.seq_interval_strand is not None:
            result["Seq-interval_strand"] = from_union([lambda x: to_class(SeqIntervalStrand, x), from_none], self.seq_interval_strand)
        if self.seq_interval_id is not None:
            result["Seq-interval_id"] = from_union([lambda x: to_class(FluffySeq, x), from_none], self.seq_interval_id)
        return result


class PurpleSeqLOCInt:
    seq_interval: Optional[FluffySeqInterval]

    def __init__(self, seq_interval: Optional[FluffySeqInterval]) -> None:
        self.seq_interval = seq_interval

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleSeqLOCInt':
        assert isinstance(obj, dict)
        seq_interval = from_union([FluffySeqInterval.from_dict, from_none], obj.get("Seq-interval"))
        return PurpleSeqLOCInt(seq_interval)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_interval is not None:
            result["Seq-interval"] = from_union([lambda x: to_class(FluffySeqInterval, x), from_none], self.seq_interval)
        return result


class PurpleGeneCommentarySeq:
    seq_loc_int: Optional[PurpleSeqLOCInt]

    def __init__(self, seq_loc_int: Optional[PurpleSeqLOCInt]) -> None:
        self.seq_loc_int = seq_loc_int

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleGeneCommentarySeq':
        assert isinstance(obj, dict)
        seq_loc_int = from_union([PurpleSeqLOCInt.from_dict, from_none], obj.get("Seq-loc_int"))
        return PurpleGeneCommentarySeq(seq_loc_int)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_loc_int is not None:
            result["Seq-loc_int"] = from_union([lambda x: to_class(PurpleSeqLOCInt, x), from_none], self.seq_loc_int)
        return result


class PurpleObjectID:
    object_id_id: Optional[str]

    def __init__(self, object_id_id: Optional[str]) -> None:
        self.object_id_id = object_id_id

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleObjectID':
        assert isinstance(obj, dict)
        object_id_id = from_union([from_str, from_none], obj.get("Object-id_id"))
        return PurpleObjectID(object_id_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.object_id_id is not None:
            result["Object-id_id"] = from_union([from_str, from_none], self.object_id_id)
        return result


class OrgRefDBDbtagTag:
    object_id: Optional[PurpleObjectID]

    def __init__(self, object_id: Optional[PurpleObjectID]) -> None:
        self.object_id = object_id

    @staticmethod
    def from_dict(obj: Any) -> 'OrgRefDBDbtagTag':
        assert isinstance(obj, dict)
        object_id = from_union([PurpleObjectID.from_dict, from_none], obj.get("Object-id"))
        return OrgRefDBDbtagTag(object_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.object_id is not None:
            result["Object-id"] = from_union([lambda x: to_class(PurpleObjectID, x), from_none], self.object_id)
        return result


class OrgRefDB:
    dbtag_db: Optional[str]
    dbtag_tag: Optional[OrgRefDBDbtagTag]

    def __init__(self, dbtag_db: Optional[str], dbtag_tag: Optional[OrgRefDBDbtagTag]) -> None:
        self.dbtag_db = dbtag_db
        self.dbtag_tag = dbtag_tag

    @staticmethod
    def from_dict(obj: Any) -> 'OrgRefDB':
        assert isinstance(obj, dict)
        dbtag_db = from_union([from_str, from_none], obj.get("Dbtag_db"))
        dbtag_tag = from_union([OrgRefDBDbtagTag.from_dict, from_none], obj.get("Dbtag_tag"))
        return OrgRefDB(dbtag_db, dbtag_tag)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.dbtag_db is not None:
            result["Dbtag_db"] = from_union([from_str, from_none], self.dbtag_db)
        if self.dbtag_tag is not None:
            result["Dbtag_tag"] = from_union([lambda x: to_class(OrgRefDBDbtagTag, x), from_none], self.dbtag_tag)
        return result


class PurpleOtherSourceSrc:
    dbtag: Optional[OrgRefDB]

    def __init__(self, dbtag: Optional[OrgRefDB]) -> None:
        self.dbtag = dbtag

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleOtherSourceSrc':
        assert isinstance(obj, dict)
        dbtag = from_union([OrgRefDB.from_dict, from_none], obj.get("Dbtag"))
        return PurpleOtherSourceSrc(dbtag)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.dbtag is not None:
            result["Dbtag"] = from_union([lambda x: to_class(OrgRefDB, x), from_none], self.dbtag)
        return result


class PurpleGeneCommentarySource:
    other_source_anchor: Optional[str]
    other_source_url: Optional[str]
    other_source_src: Optional[PurpleOtherSourceSrc]

    def __init__(self, other_source_anchor: Optional[str], other_source_url: Optional[str], other_source_src: Optional[PurpleOtherSourceSrc]) -> None:
        self.other_source_anchor = other_source_anchor
        self.other_source_url = other_source_url
        self.other_source_src = other_source_src

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        other_source_url = from_union([from_str, from_none], obj.get("Other-source_url"))
        other_source_src = from_union([PurpleOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        return PurpleGeneCommentarySource(other_source_anchor, other_source_url, other_source_src)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        if self.other_source_url is not None:
            result["Other-source_url"] = from_union([from_str, from_none], self.other_source_url)
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(PurpleOtherSourceSrc, x), from_none], self.other_source_src)
        return result


class PurpleGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_refs: Optional[List[GeneCommentaryRef]]
    gene_commentary_source: Optional[List[PurpleGeneCommentarySource]]
    gene_commentary_label: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_seqs: Optional[List[PurpleGeneCommentarySeq]]
    gene_commentary_text: Optional[str]
    gene_commentary_comment: Optional[List[FluffyGeneCommentaryComment]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_refs: Optional[List[GeneCommentaryRef]], gene_commentary_source: Optional[List[PurpleGeneCommentarySource]], gene_commentary_label: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_seqs: Optional[List[PurpleGeneCommentarySeq]], gene_commentary_text: Optional[str], gene_commentary_comment: Optional[List[FluffyGeneCommentaryComment]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_refs = gene_commentary_refs
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_seqs = gene_commentary_seqs
        self.gene_commentary_text = gene_commentary_text
        self.gene_commentary_comment = gene_commentary_comment

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_refs = from_union([lambda x: from_list(GeneCommentaryRef.from_dict, x), from_none], obj.get("Gene-commentary_refs"))
        gene_commentary_source = from_union([lambda x: from_list(PurpleGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_seqs = from_union([lambda x: from_list(PurpleGeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        gene_commentary_comment = from_union([lambda x: from_list(FluffyGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        return PurpleGeneCommentaryComment(gene_commentary_type, gene_commentary_heading, gene_commentary_refs, gene_commentary_source, gene_commentary_label, gene_commentary_accession, gene_commentary_version, gene_commentary_seqs, gene_commentary_text, gene_commentary_comment)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_refs is not None:
            result["Gene-commentary_refs"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentaryRef, x), x), from_none], self.gene_commentary_refs)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(PurpleGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(PurpleGeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(FluffyGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        return result


class FluffyDateStd:
    date_std_year: Optional[str]
    date_std_month: Optional[str]
    date_std_day: Optional[str]
    date_std_hour: Optional[str]
    date_std_minute: Optional[str]
    date_std_second: Optional[str]

    def __init__(self, date_std_year: Optional[str], date_std_month: Optional[str], date_std_day: Optional[str], date_std_hour: Optional[str], date_std_minute: Optional[str], date_std_second: Optional[str]) -> None:
        self.date_std_year = date_std_year
        self.date_std_month = date_std_month
        self.date_std_day = date_std_day
        self.date_std_hour = date_std_hour
        self.date_std_minute = date_std_minute
        self.date_std_second = date_std_second

    @staticmethod
    def from_dict(obj: Any) -> 'FluffyDateStd':
        assert isinstance(obj, dict)
        date_std_year = from_union([from_str, from_none], obj.get("Date-std_year"))
        date_std_month = from_union([from_str, from_none], obj.get("Date-std_month"))
        date_std_day = from_union([from_str, from_none], obj.get("Date-std_day"))
        date_std_hour = from_union([from_str, from_none], obj.get("Date-std_hour"))
        date_std_minute = from_union([from_str, from_none], obj.get("Date-std_minute"))
        date_std_second = from_union([from_str, from_none], obj.get("Date-std_second"))
        return FluffyDateStd(date_std_year, date_std_month, date_std_day, date_std_hour, date_std_minute, date_std_second)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.date_std_year is not None:
            result["Date-std_year"] = from_union([from_str, from_none], self.date_std_year)
        if self.date_std_month is not None:
            result["Date-std_month"] = from_union([from_str, from_none], self.date_std_month)
        if self.date_std_day is not None:
            result["Date-std_day"] = from_union([from_str, from_none], self.date_std_day)
        if self.date_std_hour is not None:
            result["Date-std_hour"] = from_union([from_str, from_none], self.date_std_hour)
        if self.date_std_minute is not None:
            result["Date-std_minute"] = from_union([from_str, from_none], self.date_std_minute)
        if self.date_std_second is not None:
            result["Date-std_second"] = from_union([from_str, from_none], self.date_std_second)
        return result


class PurpleDateStd:
    date_std: Optional[FluffyDateStd]

    def __init__(self, date_std: Optional[FluffyDateStd]) -> None:
        self.date_std = date_std

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleDateStd':
        assert isinstance(obj, dict)
        date_std = from_union([FluffyDateStd.from_dict, from_none], obj.get("Date-std"))
        return PurpleDateStd(date_std)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.date_std is not None:
            result["Date-std"] = from_union([lambda x: to_class(FluffyDateStd, x), from_none], self.date_std)
        return result


class GeneCommentaryCreateDateDate:
    date_std: Optional[PurpleDateStd]

    def __init__(self, date_std: Optional[PurpleDateStd]) -> None:
        self.date_std = date_std

    @staticmethod
    def from_dict(obj: Any) -> 'GeneCommentaryCreateDateDate':
        assert isinstance(obj, dict)
        date_std = from_union([PurpleDateStd.from_dict, from_none], obj.get("Date_std"))
        return GeneCommentaryCreateDateDate(date_std)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.date_std is not None:
            result["Date_std"] = from_union([lambda x: to_class(PurpleDateStd, x), from_none], self.date_std)
        return result


class GeneAteDate:
    date: Optional[GeneCommentaryCreateDateDate]

    def __init__(self, date: Optional[GeneCommentaryCreateDateDate]) -> None:
        self.date = date

    @staticmethod
    def from_dict(obj: Any) -> 'GeneAteDate':
        assert isinstance(obj, dict)
        date = from_union([GeneCommentaryCreateDateDate.from_dict, from_none], obj.get("Date"))
        return GeneAteDate(date)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.date is not None:
            result["Date"] = from_union([lambda x: to_class(GeneCommentaryCreateDateDate, x), from_none], self.date)
        return result


class FluffyObjectID:
    object_id_str: Optional[str]

    def __init__(self, object_id_str: Optional[str]) -> None:
        self.object_id_str = object_id_str

    @staticmethod
    def from_dict(obj: Any) -> 'FluffyObjectID':
        assert isinstance(obj, dict)
        object_id_str = from_union([from_str, from_none], obj.get("Object-id_str"))
        return FluffyObjectID(object_id_str)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.object_id_str is not None:
            result["Object-id_str"] = from_union([from_str, from_none], self.object_id_str)
        return result


class PurpleDbtagTag:
    object_id: Optional[FluffyObjectID]

    def __init__(self, object_id: Optional[FluffyObjectID]) -> None:
        self.object_id = object_id

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleDbtagTag':
        assert isinstance(obj, dict)
        object_id = from_union([FluffyObjectID.from_dict, from_none], obj.get("Object-id"))
        return PurpleDbtagTag(object_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.object_id is not None:
            result["Object-id"] = from_union([lambda x: to_class(FluffyObjectID, x), from_none], self.object_id)
        return result


class Dbtag:
    dbtag_db: Optional[str]
    dbtag_tag: Optional[PurpleDbtagTag]

    def __init__(self, dbtag_db: Optional[str], dbtag_tag: Optional[PurpleDbtagTag]) -> None:
        self.dbtag_db = dbtag_db
        self.dbtag_tag = dbtag_tag

    @staticmethod
    def from_dict(obj: Any) -> 'Dbtag':
        assert isinstance(obj, dict)
        dbtag_db = from_union([from_str, from_none], obj.get("Dbtag_db"))
        dbtag_tag = from_union([PurpleDbtagTag.from_dict, from_none], obj.get("Dbtag_tag"))
        return Dbtag(dbtag_db, dbtag_tag)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.dbtag_db is not None:
            result["Dbtag_db"] = from_union([from_str, from_none], self.dbtag_db)
        if self.dbtag_tag is not None:
            result["Dbtag_tag"] = from_union([lambda x: to_class(PurpleDbtagTag, x), from_none], self.dbtag_tag)
        return result


class FluffyOtherSourceSrc:
    dbtag: Optional[Dbtag]

    def __init__(self, dbtag: Optional[Dbtag]) -> None:
        self.dbtag = dbtag

    @staticmethod
    def from_dict(obj: Any) -> 'FluffyOtherSourceSrc':
        assert isinstance(obj, dict)
        dbtag = from_union([Dbtag.from_dict, from_none], obj.get("Dbtag"))
        return FluffyOtherSourceSrc(dbtag)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.dbtag is not None:
            result["Dbtag"] = from_union([lambda x: to_class(Dbtag, x), from_none], self.dbtag)
        return result


class EntrezgeneCommentGeneCommentarySource:
    other_source_src: Optional[FluffyOtherSourceSrc]
    other_source_anchor: Optional[str]

    def __init__(self, other_source_src: Optional[FluffyOtherSourceSrc], other_source_anchor: Optional[str]) -> None:
        self.other_source_src = other_source_src
        self.other_source_anchor = other_source_anchor

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneCommentGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_src = from_union([FluffyOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        return EntrezgeneCommentGeneCommentarySource(other_source_src, other_source_anchor)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(FluffyOtherSourceSrc, x), from_none], self.other_source_src)
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        return result


class StickyGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_text: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_source: Optional[List[EntrezgeneCommentGeneCommentarySource]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_label: Optional[str], gene_commentary_text: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_source: Optional[List[EntrezgeneCommentGeneCommentarySource]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_text = gene_commentary_text
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_source = gene_commentary_source

    @staticmethod
    def from_dict(obj: Any) -> 'StickyGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_source = from_union([lambda x: from_list(EntrezgeneCommentGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        return StickyGeneCommentaryComment(gene_commentary_type, gene_commentary_label, gene_commentary_text, gene_commentary_heading, gene_commentary_source)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneCommentGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        return result


class HilariousGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_text: Optional[str]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_text: Optional[str]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_text = gene_commentary_text

    @staticmethod
    def from_dict(obj: Any) -> 'HilariousGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        return HilariousGeneCommentaryComment(gene_commentary_type, gene_commentary_text)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        return result


class TentacledObjectID:
    object_id_str: Optional[str]
    object_id_id: Optional[str]

    def __init__(self, object_id_str: Optional[str], object_id_id: Optional[str]) -> None:
        self.object_id_str = object_id_str
        self.object_id_id = object_id_id

    @staticmethod
    def from_dict(obj: Any) -> 'TentacledObjectID':
        assert isinstance(obj, dict)
        object_id_str = from_union([from_str, from_none], obj.get("Object-id_str"))
        object_id_id = from_union([from_str, from_none], obj.get("Object-id_id"))
        return TentacledObjectID(object_id_str, object_id_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.object_id_str is not None:
            result["Object-id_str"] = from_union([from_str, from_none], self.object_id_str)
        if self.object_id_id is not None:
            result["Object-id_id"] = from_union([from_str, from_none], self.object_id_id)
        return result


class EntrezgeneUniqueKeyDbtagTag:
    object_id: Optional[TentacledObjectID]

    def __init__(self, object_id: Optional[TentacledObjectID]) -> None:
        self.object_id = object_id

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneUniqueKeyDbtagTag':
        assert isinstance(obj, dict)
        object_id = from_union([TentacledObjectID.from_dict, from_none], obj.get("Object-id"))
        return EntrezgeneUniqueKeyDbtagTag(object_id)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.object_id is not None:
            result["Object-id"] = from_union([lambda x: to_class(TentacledObjectID, x), from_none], self.object_id)
        return result


class EntrezgeneUniqueKey:
    dbtag_db: Optional[str]
    dbtag_tag: Optional[EntrezgeneUniqueKeyDbtagTag]

    def __init__(self, dbtag_db: Optional[str], dbtag_tag: Optional[EntrezgeneUniqueKeyDbtagTag]) -> None:
        self.dbtag_db = dbtag_db
        self.dbtag_tag = dbtag_tag

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneUniqueKey':
        assert isinstance(obj, dict)
        dbtag_db = from_union([from_str, from_none], obj.get("Dbtag_db"))
        dbtag_tag = from_union([EntrezgeneUniqueKeyDbtagTag.from_dict, from_none], obj.get("Dbtag_tag"))
        return EntrezgeneUniqueKey(dbtag_db, dbtag_tag)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.dbtag_db is not None:
            result["Dbtag_db"] = from_union([from_str, from_none], self.dbtag_db)
        if self.dbtag_tag is not None:
            result["Dbtag_tag"] = from_union([lambda x: to_class(EntrezgeneUniqueKeyDbtagTag, x), from_none], self.dbtag_tag)
        return result


class TentacledOtherSourceSrc:
    dbtag: Optional[EntrezgeneUniqueKey]

    def __init__(self, dbtag: Optional[EntrezgeneUniqueKey]) -> None:
        self.dbtag = dbtag

    @staticmethod
    def from_dict(obj: Any) -> 'TentacledOtherSourceSrc':
        assert isinstance(obj, dict)
        dbtag = from_union([EntrezgeneUniqueKey.from_dict, from_none], obj.get("Dbtag"))
        return TentacledOtherSourceSrc(dbtag)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.dbtag is not None:
            result["Dbtag"] = from_union([lambda x: to_class(EntrezgeneUniqueKey, x), from_none], self.dbtag)
        return result


class FluffyGeneCommentarySource:
    other_source_src: Optional[TentacledOtherSourceSrc]
    other_source_anchor: Optional[str]

    def __init__(self, other_source_src: Optional[TentacledOtherSourceSrc], other_source_anchor: Optional[str]) -> None:
        self.other_source_src = other_source_src
        self.other_source_anchor = other_source_anchor

    @staticmethod
    def from_dict(obj: Any) -> 'FluffyGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_src = from_union([TentacledOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        return FluffyGeneCommentarySource(other_source_src, other_source_anchor)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(TentacledOtherSourceSrc, x), from_none], self.other_source_src)
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        return result


class IndecentGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_source: Optional[List[FluffyGeneCommentarySource]]
    gene_commentary_comment: Optional[List[HilariousGeneCommentaryComment]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_source: Optional[List[FluffyGeneCommentarySource]], gene_commentary_comment: Optional[List[HilariousGeneCommentaryComment]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_comment = gene_commentary_comment

    @staticmethod
    def from_dict(obj: Any) -> 'IndecentGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_source = from_union([lambda x: from_list(FluffyGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_comment = from_union([lambda x: from_list(HilariousGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        return IndecentGeneCommentaryComment(gene_commentary_type, gene_commentary_source, gene_commentary_comment)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(FluffyGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(HilariousGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        return result


class TentacledGeneCommentarySource:
    other_source_pre_text: Optional[str]
    other_source_anchor: Optional[str]
    other_source_src: Optional[FluffyOtherSourceSrc]

    def __init__(self, other_source_pre_text: Optional[str], other_source_anchor: Optional[str], other_source_src: Optional[FluffyOtherSourceSrc]) -> None:
        self.other_source_pre_text = other_source_pre_text
        self.other_source_anchor = other_source_anchor
        self.other_source_src = other_source_src

    @staticmethod
    def from_dict(obj: Any) -> 'TentacledGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_pre_text = from_union([from_str, from_none], obj.get("Other-source_pre-text"))
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        other_source_src = from_union([FluffyOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        return TentacledGeneCommentarySource(other_source_pre_text, other_source_anchor, other_source_src)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_pre_text is not None:
            result["Other-source_pre-text"] = from_union([from_str, from_none], self.other_source_pre_text)
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(FluffyOtherSourceSrc, x), from_none], self.other_source_src)
        return result


class IndigoGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_source: Optional[List[TentacledGeneCommentarySource]]
    gene_commentary_comment: Optional[List[IndecentGeneCommentaryComment]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_source: Optional[List[TentacledGeneCommentarySource]], gene_commentary_comment: Optional[List[IndecentGeneCommentaryComment]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_comment = gene_commentary_comment

    @staticmethod
    def from_dict(obj: Any) -> 'IndigoGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_source = from_union([lambda x: from_list(TentacledGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_comment = from_union([lambda x: from_list(IndecentGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        return IndigoGeneCommentaryComment(gene_commentary_type, gene_commentary_heading, gene_commentary_source, gene_commentary_comment)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(TentacledGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(IndecentGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        return result


class FluffyGeneCommentarySeq:
    seq_loc_whole: Optional[PurpleSeq]

    def __init__(self, seq_loc_whole: Optional[PurpleSeq]) -> None:
        self.seq_loc_whole = seq_loc_whole

    @staticmethod
    def from_dict(obj: Any) -> 'FluffyGeneCommentarySeq':
        assert isinstance(obj, dict)
        seq_loc_whole = from_union([PurpleSeq.from_dict, from_none], obj.get("Seq-loc_whole"))
        return FluffyGeneCommentarySeq(seq_loc_whole)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_loc_whole is not None:
            result["Seq-loc_whole"] = from_union([lambda x: to_class(PurpleSeq, x), from_none], self.seq_loc_whole)
        return result


class StickyGeneCommentarySource:
    other_source_src: Optional[TentacledOtherSourceSrc]
    other_source_anchor: Optional[str]
    other_source_post_text: Optional[str]

    def __init__(self, other_source_src: Optional[TentacledOtherSourceSrc], other_source_anchor: Optional[str], other_source_post_text: Optional[str]) -> None:
        self.other_source_src = other_source_src
        self.other_source_anchor = other_source_anchor
        self.other_source_post_text = other_source_post_text

    @staticmethod
    def from_dict(obj: Any) -> 'StickyGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_src = from_union([TentacledOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        other_source_post_text = from_union([from_str, from_none], obj.get("Other-source_post-text"))
        return StickyGeneCommentarySource(other_source_src, other_source_anchor, other_source_post_text)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(TentacledOtherSourceSrc, x), from_none], self.other_source_src)
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        if self.other_source_post_text is not None:
            result["Other-source_post-text"] = from_union([from_str, from_none], self.other_source_post_text)
        return result


class PurpleGeneCommentaryProduct:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_source: Optional[List[StickyGeneCommentarySource]]
    gene_commentary_seqs: Optional[List[FluffyGeneCommentarySeq]]
    gene_commentary_comment: Optional[List[IndigoGeneCommentaryComment]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_source: Optional[List[StickyGeneCommentarySource]], gene_commentary_seqs: Optional[List[FluffyGeneCommentarySeq]], gene_commentary_comment: Optional[List[IndigoGeneCommentaryComment]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_seqs = gene_commentary_seqs
        self.gene_commentary_comment = gene_commentary_comment

    @staticmethod
    def from_dict(obj: Any) -> 'PurpleGeneCommentaryProduct':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_source = from_union([lambda x: from_list(StickyGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_seqs = from_union([lambda x: from_list(FluffyGeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        gene_commentary_comment = from_union([lambda x: from_list(IndigoGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        return PurpleGeneCommentaryProduct(gene_commentary_type, gene_commentary_heading, gene_commentary_accession, gene_commentary_version, gene_commentary_source, gene_commentary_seqs, gene_commentary_comment)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(StickyGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(FluffyGeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(IndigoGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        return result


class TentacledGeneCommentarySeq:
    seq_loc_whole: Optional[PurpleSeq]
    seq_loc_int: Optional[SeqLOCMixSeqLOCInt]

    def __init__(self, seq_loc_whole: Optional[PurpleSeq], seq_loc_int: Optional[SeqLOCMixSeqLOCInt]) -> None:
        self.seq_loc_whole = seq_loc_whole
        self.seq_loc_int = seq_loc_int

    @staticmethod
    def from_dict(obj: Any) -> 'TentacledGeneCommentarySeq':
        assert isinstance(obj, dict)
        """
            BEGIN DK
        """
        if 'Seq-loc_pnt' in obj \
            or 'Seg-loc_bond' in obj \
                or 'Seg-loc_null' in obj \
                    or 'Seg-loc_int' in obj \
                        or 'Seg-loc_empty' in obj \
                            or 'Seq-Loc_feat' in obj \
                                or 'Seg-loc_equiv' in obj \
                                    or 'Seq-loc_packed-int' in obj \
                                        or 'Seq-loc_mix' in obj:
            print("Unhandled Seq-loc_*")
            #breakpoint()
            # obj['Seq-loc_mix']['Seq-loc-mix'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
            # obj['Seq-loc_mix']['Seq-loc-mix'][1]['Seq-loc_int']['Seq-interval']['Seq-interval_from']
        """
            END DK
        """
        seq_loc_whole = from_union([PurpleSeq.from_dict, from_none], obj.get("Seq-loc_whole"))
        seq_loc_int = from_union([SeqLOCMixSeqLOCInt.from_dict, from_none], obj.get("Seq-loc_int"))
        return TentacledGeneCommentarySeq(seq_loc_whole, seq_loc_int)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_loc_whole is not None:
            result["Seq-loc_whole"] = from_union([lambda x: to_class(PurpleSeq, x), from_none], self.seq_loc_whole)
        if self.seq_loc_int is not None:
            result["Seq-loc_int"] = from_union([lambda x: to_class(SeqLOCMixSeqLOCInt, x), from_none], self.seq_loc_int)
        return result


class GeneCommentaryCommentGeneCommentaryProduct:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_source: Optional[List[FluffyGeneCommentarySource]]
    gene_commentary_seqs: Optional[List[TentacledGeneCommentarySeq]]
    gene_commentary_products: Optional[List[PurpleGeneCommentaryProduct]]
    gene_commentary_comment: Optional[List[StickyGeneCommentaryComment]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_source: Optional[List[FluffyGeneCommentarySource]], gene_commentary_seqs: Optional[List[TentacledGeneCommentarySeq]], gene_commentary_products: Optional[List[PurpleGeneCommentaryProduct]], gene_commentary_comment: Optional[List[StickyGeneCommentaryComment]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_seqs = gene_commentary_seqs
        self.gene_commentary_products = gene_commentary_products
        self.gene_commentary_comment = gene_commentary_comment

    @staticmethod
    def from_dict(obj: Any) -> 'GeneCommentaryCommentGeneCommentaryProduct':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_source = from_union([lambda x: from_list(FluffyGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_seqs = from_union([lambda x: from_list(TentacledGeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        gene_commentary_products = from_union([lambda x: from_list(PurpleGeneCommentaryProduct.from_dict, x), from_none], obj.get("Gene-commentary_products"))
        gene_commentary_comment = from_union([lambda x: from_list(StickyGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        return GeneCommentaryCommentGeneCommentaryProduct(gene_commentary_type, gene_commentary_heading, gene_commentary_accession, gene_commentary_version, gene_commentary_source, gene_commentary_seqs, gene_commentary_products, gene_commentary_comment)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(FluffyGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(TentacledGeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        if self.gene_commentary_products is not None:
            result["Gene-commentary_products"] = from_union([lambda x: from_list(lambda x: to_class(PurpleGeneCommentaryProduct, x), x), from_none], self.gene_commentary_products)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(StickyGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        return result


class IndigoGeneCommentarySource:
    other_source_src: Optional[TentacledOtherSourceSrc]
    other_source_pre_text: Optional[str]
    other_source_anchor: Optional[str]
    other_source_url: Optional[str]
    other_source_post_text: Optional[str]

    def __init__(self, other_source_src: Optional[TentacledOtherSourceSrc], other_source_pre_text: Optional[str], other_source_anchor: Optional[str], other_source_url: Optional[str], other_source_post_text: Optional[str]) -> None:
        self.other_source_src = other_source_src
        self.other_source_pre_text = other_source_pre_text
        self.other_source_anchor = other_source_anchor
        self.other_source_url = other_source_url
        self.other_source_post_text = other_source_post_text

    @staticmethod
    def from_dict(obj: Any) -> 'IndigoGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_src = from_union([TentacledOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        other_source_pre_text = from_union([from_str, from_none], obj.get("Other-source_pre-text"))
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        other_source_url = from_union([from_str, from_none], obj.get("Other-source_url"))
        other_source_post_text = from_union([from_str, from_none], obj.get("Other-source_post-text"))
        return IndigoGeneCommentarySource(other_source_src, other_source_pre_text, other_source_anchor, other_source_url, other_source_post_text)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(TentacledOtherSourceSrc, x), from_none], self.other_source_src)
        if self.other_source_pre_text is not None:
            result["Other-source_pre-text"] = from_union([from_str, from_none], self.other_source_pre_text)
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        if self.other_source_url is not None:
            result["Other-source_url"] = from_union([from_str, from_none], self.other_source_url)
        if self.other_source_post_text is not None:
            result["Other-source_post-text"] = from_union([from_str, from_none], self.other_source_post_text)
        return result


class XtraProperty:
    xtra_terms_tag: Optional[str]
    xtra_terms_value: Optional[str]

    def __init__(self, xtra_terms_tag: Optional[str], xtra_terms_value: Optional[str]) -> None:
        self.xtra_terms_tag = xtra_terms_tag
        self.xtra_terms_value = xtra_terms_value

    @staticmethod
    def from_dict(obj: Any) -> 'XtraProperty':
        assert isinstance(obj, dict)
        xtra_terms_tag = from_union([from_str, from_none], obj.get("Xtra-Terms_tag"))
        xtra_terms_value = from_union([from_str, from_none], obj.get("Xtra-Terms_value"))
        return XtraProperty(xtra_terms_tag, xtra_terms_value)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.xtra_terms_tag is not None:
            result["Xtra-Terms_tag"] = from_union([from_str, from_none], self.xtra_terms_tag)
        if self.xtra_terms_value is not None:
            result["Xtra-Terms_value"] = from_union([from_str, from_none], self.xtra_terms_value)
        return result


class EntrezgeneCommentGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_source: Optional[List[IndigoGeneCommentarySource]]
    gene_commentary_heading: Optional[str]
    gene_commentary_comment: Optional[List[PurpleGeneCommentaryComment]]
    gene_commentary_products: Optional[List[GeneCommentaryCommentGeneCommentaryProduct]]
    gene_commentary_text: Optional[str]
    gene_commentary_xtra_properties: Optional[List[XtraProperty]]
    gene_commentary_refs: Optional[List[GeneCommentaryRef]]
    gene_commentary_create_date: Optional[GeneAteDate]
    gene_commentary_update_date: Optional[GeneAteDate]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_label: Optional[str], gene_commentary_source: Optional[List[IndigoGeneCommentarySource]], gene_commentary_heading: Optional[str], gene_commentary_comment: Optional[List[PurpleGeneCommentaryComment]], gene_commentary_products: Optional[List[GeneCommentaryCommentGeneCommentaryProduct]], gene_commentary_text: Optional[str], gene_commentary_xtra_properties: Optional[List[XtraProperty]], gene_commentary_refs: Optional[List[GeneCommentaryRef]], gene_commentary_create_date: Optional[GeneAteDate], gene_commentary_update_date: Optional[GeneAteDate]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_comment = gene_commentary_comment
        self.gene_commentary_products = gene_commentary_products
        self.gene_commentary_text = gene_commentary_text
        self.gene_commentary_xtra_properties = gene_commentary_xtra_properties
        self.gene_commentary_refs = gene_commentary_refs
        self.gene_commentary_create_date = gene_commentary_create_date
        self.gene_commentary_update_date = gene_commentary_update_date

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneCommentGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_source = from_union([lambda x: from_list(IndigoGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_comment = from_union([lambda x: from_list(PurpleGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        gene_commentary_products = from_union([lambda x: from_list(GeneCommentaryCommentGeneCommentaryProduct.from_dict, x), from_none], obj.get("Gene-commentary_products"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        gene_commentary_xtra_properties = from_union([lambda x: from_list(XtraProperty.from_dict, x), from_none], obj.get("Gene-commentary_xtra-properties"))
        gene_commentary_refs = from_union([lambda x: from_list(GeneCommentaryRef.from_dict, x), from_none], obj.get("Gene-commentary_refs"))
        gene_commentary_create_date = from_union([GeneAteDate.from_dict, from_none], obj.get("Gene-commentary_create-date"))
        gene_commentary_update_date = from_union([GeneAteDate.from_dict, from_none], obj.get("Gene-commentary_update-date"))
        return EntrezgeneCommentGeneCommentaryComment(gene_commentary_type, gene_commentary_label, gene_commentary_source, gene_commentary_heading, gene_commentary_comment, gene_commentary_products, gene_commentary_text, gene_commentary_xtra_properties, gene_commentary_refs, gene_commentary_create_date, gene_commentary_update_date)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(IndigoGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(PurpleGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        if self.gene_commentary_products is not None:
            result["Gene-commentary_products"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentaryCommentGeneCommentaryProduct, x), x), from_none], self.gene_commentary_products)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        if self.gene_commentary_xtra_properties is not None:
            result["Gene-commentary_xtra-properties"] = from_union([lambda x: from_list(lambda x: to_class(XtraProperty, x), x), from_none], self.gene_commentary_xtra_properties)
        if self.gene_commentary_refs is not None:
            result["Gene-commentary_refs"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentaryRef, x), x), from_none], self.gene_commentary_refs)
        if self.gene_commentary_create_date is not None:
            result["Gene-commentary_create-date"] = from_union([lambda x: to_class(GeneAteDate, x), from_none], self.gene_commentary_create_date)
        if self.gene_commentary_update_date is not None:
            result["Gene-commentary_update-date"] = from_union([lambda x: to_class(GeneAteDate, x), from_none], self.gene_commentary_update_date)
        return result


class StickyGeneCommentarySeq:
    seq_loc_whole: Optional[FluffySeq]

    def __init__(self, seq_loc_whole: Optional[FluffySeq]) -> None:
        self.seq_loc_whole = seq_loc_whole

    @staticmethod
    def from_dict(obj: Any) -> 'StickyGeneCommentarySeq':
        assert isinstance(obj, dict)
        seq_loc_whole = from_union([FluffySeq.from_dict, from_none], obj.get("Seq-loc_whole"))
        return StickyGeneCommentarySeq(seq_loc_whole)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_loc_whole is not None:
            result["Seq-loc_whole"] = from_union([lambda x: to_class(FluffySeq, x), from_none], self.seq_loc_whole)
        return result


class IndecentGeneCommentarySource:
    other_source_src: Optional[TentacledOtherSourceSrc]
    other_source_anchor: Optional[str]
    other_source_url: Optional[str]

    def __init__(self, other_source_src: Optional[TentacledOtherSourceSrc], other_source_anchor: Optional[str], other_source_url: Optional[str]) -> None:
        self.other_source_src = other_source_src
        self.other_source_anchor = other_source_anchor
        self.other_source_url = other_source_url

    @staticmethod
    def from_dict(obj: Any) -> 'IndecentGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_src = from_union([TentacledOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        other_source_url = from_union([from_str, from_none], obj.get("Other-source_url"))
        return IndecentGeneCommentarySource(other_source_src, other_source_anchor, other_source_url)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(TentacledOtherSourceSrc, x), from_none], self.other_source_src)
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        if self.other_source_url is not None:
            result["Other-source_url"] = from_union([from_str, from_none], self.other_source_url)
        return result


class FluffyGeneCommentaryProduct:
    gene_commentary_type: Optional[str]
    gene_commentary_text: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_source: Optional[List[IndecentGeneCommentarySource]]
    gene_commentary_seqs: Optional[List[StickyGeneCommentarySeq]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_text: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_source: Optional[List[IndecentGeneCommentarySource]], gene_commentary_seqs: Optional[List[StickyGeneCommentarySeq]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_text = gene_commentary_text
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_seqs = gene_commentary_seqs

    @staticmethod
    def from_dict(obj: Any) -> 'FluffyGeneCommentaryProduct':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_source = from_union([lambda x: from_list(IndecentGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_seqs = from_union([lambda x: from_list(StickyGeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        return FluffyGeneCommentaryProduct(gene_commentary_type, gene_commentary_text, gene_commentary_accession, gene_commentary_version, gene_commentary_source, gene_commentary_seqs)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(IndecentGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(StickyGeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        return result


class IndigoGeneCommentarySeq:
    seq_loc_int: Optional[PurpleSeqLOCInt]
    seq_loc_whole: Optional[FluffySeq]

    def __init__(self, seq_loc_int: Optional[PurpleSeqLOCInt], seq_loc_whole: Optional[FluffySeq]) -> None:
        self.seq_loc_int = seq_loc_int
        self.seq_loc_whole = seq_loc_whole

    @staticmethod
    def from_dict(obj: Any) -> 'IndigoGeneCommentarySeq':
        assert isinstance(obj, dict)
        seq_loc_int = from_union([PurpleSeqLOCInt.from_dict, from_none], obj.get("Seq-loc_int"))
        seq_loc_whole = from_union([FluffySeq.from_dict, from_none], obj.get("Seq-loc_whole"))
        return IndigoGeneCommentarySeq(seq_loc_int, seq_loc_whole)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_loc_int is not None:
            result["Seq-loc_int"] = from_union([lambda x: to_class(PurpleSeqLOCInt, x), from_none], self.seq_loc_int)
        if self.seq_loc_whole is not None:
            result["Seq-loc_whole"] = from_union([lambda x: to_class(FluffySeq, x), from_none], self.seq_loc_whole)
        return result


class HilariousGeneCommentarySource:
    other_source_src: Optional[PurpleOtherSourceSrc]
    other_source_anchor: Optional[str]

    def __init__(self, other_source_src: Optional[PurpleOtherSourceSrc], other_source_anchor: Optional[str]) -> None:
        self.other_source_src = other_source_src
        self.other_source_anchor = other_source_anchor

    @staticmethod
    def from_dict(obj: Any) -> 'HilariousGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_src = from_union([PurpleOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        return HilariousGeneCommentarySource(other_source_src, other_source_anchor)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(PurpleOtherSourceSrc, x), from_none], self.other_source_src)
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        return result


class EntrezgeneCommentGeneCommentaryProduct:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_source: Optional[List[HilariousGeneCommentarySource]]
    gene_commentary_seqs: Optional[List[IndigoGeneCommentarySeq]]
    gene_commentary_products: Optional[List[FluffyGeneCommentaryProduct]]
    gene_commentary_text: Optional[str]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_source: Optional[List[HilariousGeneCommentarySource]], gene_commentary_seqs: Optional[List[IndigoGeneCommentarySeq]], gene_commentary_products: Optional[List[FluffyGeneCommentaryProduct]], gene_commentary_text: Optional[str]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_seqs = gene_commentary_seqs
        self.gene_commentary_products = gene_commentary_products
        self.gene_commentary_text = gene_commentary_text

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneCommentGeneCommentaryProduct':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_source = from_union([lambda x: from_list(HilariousGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_seqs = from_union([lambda x: from_list(IndigoGeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        gene_commentary_products = from_union([lambda x: from_list(FluffyGeneCommentaryProduct.from_dict, x), from_none], obj.get("Gene-commentary_products"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        return EntrezgeneCommentGeneCommentaryProduct(gene_commentary_type, gene_commentary_heading, gene_commentary_accession, gene_commentary_version, gene_commentary_source, gene_commentary_seqs, gene_commentary_products, gene_commentary_text)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(HilariousGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(IndigoGeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        if self.gene_commentary_products is not None:
            result["Gene-commentary_products"] = from_union([lambda x: from_list(lambda x: to_class(FluffyGeneCommentaryProduct, x), x), from_none], self.gene_commentary_products)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        return result


class EntrezgeneComment:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_comment: Optional[List[EntrezgeneCommentGeneCommentaryComment]]
    gene_commentary_refs: Optional[List[GeneCommentaryRef]]
    gene_commentary_products: Optional[List[EntrezgeneCommentGeneCommentaryProduct]]
    gene_commentary_text: Optional[str]
    gene_commentary_create_date: Optional[GeneAteDate]
    gene_commentary_update_date: Optional[GeneAteDate]
    gene_commentary_source: Optional[List[EntrezgeneCommentGeneCommentarySource]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_label: Optional[str], gene_commentary_comment: Optional[List[EntrezgeneCommentGeneCommentaryComment]], gene_commentary_refs: Optional[List[GeneCommentaryRef]], gene_commentary_products: Optional[List[EntrezgeneCommentGeneCommentaryProduct]], gene_commentary_text: Optional[str], gene_commentary_create_date: Optional[GeneAteDate], gene_commentary_update_date: Optional[GeneAteDate], gene_commentary_source: Optional[List[EntrezgeneCommentGeneCommentarySource]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_comment = gene_commentary_comment
        self.gene_commentary_refs = gene_commentary_refs
        self.gene_commentary_products = gene_commentary_products
        self.gene_commentary_text = gene_commentary_text
        self.gene_commentary_create_date = gene_commentary_create_date
        self.gene_commentary_update_date = gene_commentary_update_date
        self.gene_commentary_source = gene_commentary_source

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_comment = from_union([lambda x: from_list(EntrezgeneCommentGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        gene_commentary_refs = from_union([lambda x: from_list(GeneCommentaryRef.from_dict, x), from_none], obj.get("Gene-commentary_refs"))
        gene_commentary_products = from_union([lambda x: from_list(EntrezgeneCommentGeneCommentaryProduct.from_dict, x), from_none], obj.get("Gene-commentary_products"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        gene_commentary_create_date = from_union([GeneAteDate.from_dict, from_none], obj.get("Gene-commentary_create-date"))
        gene_commentary_update_date = from_union([GeneAteDate.from_dict, from_none], obj.get("Gene-commentary_update-date"))
        gene_commentary_source = from_union([lambda x: from_list(EntrezgeneCommentGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        return EntrezgeneComment(gene_commentary_type, gene_commentary_heading, gene_commentary_label, gene_commentary_comment, gene_commentary_refs, gene_commentary_products, gene_commentary_text, gene_commentary_create_date, gene_commentary_update_date, gene_commentary_source)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneCommentGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        if self.gene_commentary_refs is not None:
            result["Gene-commentary_refs"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentaryRef, x), x), from_none], self.gene_commentary_refs)
        if self.gene_commentary_products is not None:
            result["Gene-commentary_products"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneCommentGeneCommentaryProduct, x), x), from_none], self.gene_commentary_products)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        if self.gene_commentary_create_date is not None:
            result["Gene-commentary_create-date"] = from_union([lambda x: to_class(GeneAteDate, x), from_none], self.gene_commentary_create_date)
        if self.gene_commentary_update_date is not None:
            result["Gene-commentary_update-date"] = from_union([lambda x: to_class(GeneAteDate, x), from_none], self.gene_commentary_update_date)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneCommentGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        return result


class GeneRef:
    gene_ref_locus: Optional[str]
    gene_ref_desc: Optional[str]
    gene_ref_maploc: Optional[str]
    gene_ref_db: Optional[List[EntrezgeneUniqueKey]]
    gene_ref_syn: Optional[List[str]]

    def __init__(self, gene_ref_locus: Optional[str], gene_ref_desc: Optional[str], gene_ref_maploc: Optional[str], gene_ref_db: Optional[List[EntrezgeneUniqueKey]], gene_ref_syn: Optional[List[str]]) -> None:
        self.gene_ref_locus = gene_ref_locus
        self.gene_ref_desc = gene_ref_desc
        self.gene_ref_maploc = gene_ref_maploc
        self.gene_ref_db = gene_ref_db
        self.gene_ref_syn = gene_ref_syn

    @staticmethod
    def from_dict(obj: Any) -> 'GeneRef':
        assert isinstance(obj, dict)
        gene_ref_locus = from_union([from_str, from_none], obj.get("Gene-ref_locus"))
        gene_ref_desc = from_union([from_str, from_none], obj.get("Gene-ref_desc"))
        gene_ref_maploc = from_union([from_str, from_none], obj.get("Gene-ref_maploc"))
        gene_ref_db = from_union([lambda x: from_list(EntrezgeneUniqueKey.from_dict, x), from_none], obj.get("Gene-ref_db"))
        gene_ref_syn = from_union([lambda x: from_list(from_str, x), from_none], obj.get("Gene-ref_syn"))
        return GeneRef(gene_ref_locus, gene_ref_desc, gene_ref_maploc, gene_ref_db, gene_ref_syn)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_ref_locus is not None:
            result["Gene-ref_locus"] = from_union([from_str, from_none], self.gene_ref_locus)
        if self.gene_ref_desc is not None:
            result["Gene-ref_desc"] = from_union([from_str, from_none], self.gene_ref_desc)
        if self.gene_ref_maploc is not None:
            result["Gene-ref_maploc"] = from_union([from_str, from_none], self.gene_ref_maploc)
        if self.gene_ref_db is not None:
            result["Gene-ref_db"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneUniqueKey, x), x), from_none], self.gene_ref_db)
        if self.gene_ref_syn is not None:
            result["Gene-ref_syn"] = from_union([lambda x: from_list(from_str, x), from_none], self.gene_ref_syn)
        return result


class EntrezgeneGene:
    gene_ref: Optional[GeneRef]

    def __init__(self, gene_ref: Optional[GeneRef]) -> None:
        self.gene_ref = gene_ref

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneGene':
        assert isinstance(obj, dict)
        gene_ref = from_union([GeneRef.from_dict, from_none], obj.get("Gene-ref"))
        return EntrezgeneGene(gene_ref)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_ref is not None:
            result["Gene-ref"] = from_union([lambda x: to_class(GeneRef, x), from_none], self.gene_ref)
        return result


class GeneSource:
    gene_source_src: Optional[str]
    gene_source_src_int: Optional[str]
    gene_source_src_str2: Optional[str]

    def __init__(self, gene_source_src: Optional[str], gene_source_src_int: Optional[str], gene_source_src_str2: Optional[str]) -> None:
        self.gene_source_src = gene_source_src
        self.gene_source_src_int = gene_source_src_int
        self.gene_source_src_str2 = gene_source_src_str2

    @staticmethod
    def from_dict(obj: Any) -> 'GeneSource':
        assert isinstance(obj, dict)
        gene_source_src = from_union([from_str, from_none], obj.get("Gene-source_src"))
        gene_source_src_int = from_union([from_str, from_none], obj.get("Gene-source_src-int"))
        gene_source_src_str2 = from_union([from_str, from_none], obj.get("Gene-source_src-str2"))
        return GeneSource(gene_source_src, gene_source_src_int, gene_source_src_str2)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_source_src is not None:
            result["Gene-source_src"] = from_union([from_str, from_none], self.gene_source_src)
        if self.gene_source_src_int is not None:
            result["Gene-source_src-int"] = from_union([from_str, from_none], self.gene_source_src_int)
        if self.gene_source_src_str2 is not None:
            result["Gene-source_src-str2"] = from_union([from_str, from_none], self.gene_source_src_str2)
        return result


class EntrezgeneGeneSource:
    gene_source: Optional[GeneSource]

    def __init__(self, gene_source: Optional[GeneSource]) -> None:
        self.gene_source = gene_source

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneGeneSource':
        assert isinstance(obj, dict)
        gene_source = from_union([GeneSource.from_dict, from_none], obj.get("Gene-source"))
        return EntrezgeneGeneSource(gene_source)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_source is not None:
            result["Gene-source"] = from_union([lambda x: to_class(GeneSource, x), from_none], self.gene_source)
        return result


class MapsMethod:
    maps_method_map_type: Optional[str]

    def __init__(self, maps_method_map_type: Optional[str]) -> None:
        self.maps_method_map_type = maps_method_map_type

    @staticmethod
    def from_dict(obj: Any) -> 'MapsMethod':
        assert isinstance(obj, dict)
        maps_method_map_type = from_union([from_str, from_none], obj.get("Maps_method_map-type"))
        return MapsMethod(maps_method_map_type)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.maps_method_map_type is not None:
            result["Maps_method_map-type"] = from_union([from_str, from_none], self.maps_method_map_type)
        return result


class EntrezgeneLocation:
    maps_display_str: Optional[str]
    maps_method: Optional[MapsMethod]

    def __init__(self, maps_display_str: Optional[str], maps_method: Optional[MapsMethod]) -> None:
        self.maps_display_str = maps_display_str
        self.maps_method = maps_method

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneLocation':
        assert isinstance(obj, dict)
        maps_display_str = from_union([from_str, from_none], obj.get("Maps_display-str"))
        maps_method = from_union([MapsMethod.from_dict, from_none], obj.get("Maps_method"))
        return EntrezgeneLocation(maps_display_str, maps_method)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.maps_display_str is not None:
            result["Maps_display-str"] = from_union([from_str, from_none], self.maps_display_str)
        if self.maps_method is not None:
            result["Maps_method"] = from_union([lambda x: to_class(MapsMethod, x), from_none], self.maps_method)
        return result


class SeqLOCMix:
    seq_loc_mix: Optional[List[GeneCommentarySeq]]

    def __init__(self, seq_loc_mix: Optional[List[GeneCommentarySeq]]) -> None:
        self.seq_loc_mix = seq_loc_mix

    @staticmethod
    def from_dict(obj: Any) -> 'SeqLOCMix':
        assert isinstance(obj, dict)
        seq_loc_mix = from_union([lambda x: from_list(GeneCommentarySeq.from_dict, x), from_none], obj.get("Seq-loc-mix"))
        return SeqLOCMix(seq_loc_mix)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_loc_mix is not None:
            result["Seq-loc-mix"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentarySeq, x), x), from_none], self.seq_loc_mix)
        return result


class GeneCommentaryGenomicCoord:
    seq_loc_mix: Optional[SeqLOCMix]

    def __init__(self, seq_loc_mix: Optional[SeqLOCMix]) -> None:
        self.seq_loc_mix = seq_loc_mix

    @staticmethod
    def from_dict(obj: Any) -> 'GeneCommentaryGenomicCoord':
        assert isinstance(obj, dict)
        seq_loc_mix = from_union([SeqLOCMix.from_dict, from_none], obj.get("Seq-loc_mix"))
        return GeneCommentaryGenomicCoord(seq_loc_mix)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.seq_loc_mix is not None:
            result["Seq-loc_mix"] = from_union([lambda x: to_class(SeqLOCMix, x), from_none], self.seq_loc_mix)
        return result


class TentacledGeneCommentaryProduct:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_genomic_coords: Optional[List[GeneCommentaryGenomicCoord]]
    gene_commentary_seqs: Optional[List[FluffyGeneCommentarySeq]]
    gene_commentary_source: Optional[List[EntrezgeneCommentGeneCommentarySource]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_label: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_genomic_coords: Optional[List[GeneCommentaryGenomicCoord]], gene_commentary_seqs: Optional[List[FluffyGeneCommentarySeq]], gene_commentary_source: Optional[List[EntrezgeneCommentGeneCommentarySource]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_genomic_coords = gene_commentary_genomic_coords
        self.gene_commentary_seqs = gene_commentary_seqs
        self.gene_commentary_source = gene_commentary_source

    @staticmethod
    def from_dict(obj: Any) -> 'TentacledGeneCommentaryProduct':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_genomic_coords = from_union([lambda x: from_list(GeneCommentaryGenomicCoord.from_dict, x), from_none], obj.get("Gene-commentary_genomic-coords"))
        gene_commentary_seqs = from_union([lambda x: from_list(FluffyGeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        gene_commentary_source = from_union([lambda x: from_list(EntrezgeneCommentGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        return TentacledGeneCommentaryProduct(gene_commentary_type, gene_commentary_heading, gene_commentary_label, gene_commentary_accession, gene_commentary_version, gene_commentary_genomic_coords, gene_commentary_seqs, gene_commentary_source)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_genomic_coords is not None:
            result["Gene-commentary_genomic-coords"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentaryGenomicCoord, x), x), from_none], self.gene_commentary_genomic_coords)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(FluffyGeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneCommentGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        return result


class EntrezgeneLocusGeneCommentaryProduct:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_genomic_coords: Optional[List[GeneCommentaryGenomicCoord]]
    gene_commentary_seqs: Optional[List[FluffyGeneCommentarySeq]]
    gene_commentary_products: Optional[List[TentacledGeneCommentaryProduct]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_label: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_genomic_coords: Optional[List[GeneCommentaryGenomicCoord]], gene_commentary_seqs: Optional[List[FluffyGeneCommentarySeq]], gene_commentary_products: Optional[List[TentacledGeneCommentaryProduct]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_genomic_coords = gene_commentary_genomic_coords
        self.gene_commentary_seqs = gene_commentary_seqs
        self.gene_commentary_products = gene_commentary_products

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneLocusGeneCommentaryProduct':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_genomic_coords = from_union([lambda x: from_list(GeneCommentaryGenomicCoord.from_dict, x), from_none], obj.get("Gene-commentary_genomic-coords"))
        gene_commentary_seqs = from_union([lambda x: from_list(FluffyGeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        gene_commentary_products = from_union([lambda x: from_list(TentacledGeneCommentaryProduct.from_dict, x), from_none], obj.get("Gene-commentary_products"))
        return EntrezgeneLocusGeneCommentaryProduct(gene_commentary_type, gene_commentary_heading, gene_commentary_label, gene_commentary_accession, gene_commentary_version, gene_commentary_genomic_coords, gene_commentary_seqs, gene_commentary_products)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_genomic_coords is not None:
            result["Gene-commentary_genomic-coords"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentaryGenomicCoord, x), x), from_none], self.gene_commentary_genomic_coords)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(FluffyGeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        if self.gene_commentary_products is not None:
            result["Gene-commentary_products"] = from_union([lambda x: from_list(lambda x: to_class(TentacledGeneCommentaryProduct, x), x), from_none], self.gene_commentary_products)
        return result


class EntrezgeneLocus:
    gene_commentary_type: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_accession: Optional[str]
    gene_commentary_version: Optional[str]
    gene_commentary_seqs: Optional[List[GeneCommentarySeq]]
    gene_commentary_products: Optional[List[EntrezgeneLocusGeneCommentaryProduct]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_label: Optional[str], gene_commentary_accession: Optional[str], gene_commentary_version: Optional[str], gene_commentary_seqs: Optional[List[GeneCommentarySeq]], gene_commentary_products: Optional[List[EntrezgeneLocusGeneCommentaryProduct]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_accession = gene_commentary_accession
        self.gene_commentary_version = gene_commentary_version
        self.gene_commentary_seqs = gene_commentary_seqs
        self.gene_commentary_products = gene_commentary_products

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneLocus':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_accession = from_union([from_str, from_none], obj.get("Gene-commentary_accession"))
        gene_commentary_version = from_union([from_str, from_none], obj.get("Gene-commentary_version"))
        gene_commentary_seqs = from_union([lambda x: from_list(GeneCommentarySeq.from_dict, x), from_none], obj.get("Gene-commentary_seqs"))
        gene_commentary_products = from_union([lambda x: from_list(EntrezgeneLocusGeneCommentaryProduct.from_dict, x), from_none], obj.get("Gene-commentary_products"))
        return EntrezgeneLocus(gene_commentary_type, gene_commentary_heading, gene_commentary_label, gene_commentary_accession, gene_commentary_version, gene_commentary_seqs, gene_commentary_products)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_accession is not None:
            result["Gene-commentary_accession"] = from_union([from_str, from_none], self.gene_commentary_accession)
        if self.gene_commentary_version is not None:
            result["Gene-commentary_version"] = from_union([from_str, from_none], self.gene_commentary_version)
        if self.gene_commentary_seqs is not None:
            result["Gene-commentary_seqs"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentarySeq, x), x), from_none], self.gene_commentary_seqs)
        if self.gene_commentary_products is not None:
            result["Gene-commentary_products"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneLocusGeneCommentaryProduct, x), x), from_none], self.gene_commentary_products)
        return result


class AmbitiousGeneCommentarySource:
    other_source_src: Optional[PurpleOtherSourceSrc]
    other_source_pre_text: Optional[str]
    other_source_anchor: Optional[str]
    other_source_post_text: Optional[str]

    def __init__(self, other_source_src: Optional[PurpleOtherSourceSrc], other_source_pre_text: Optional[str], other_source_anchor: Optional[str], other_source_post_text: Optional[str]) -> None:
        self.other_source_src = other_source_src
        self.other_source_pre_text = other_source_pre_text
        self.other_source_anchor = other_source_anchor
        self.other_source_post_text = other_source_post_text

    @staticmethod
    def from_dict(obj: Any) -> 'AmbitiousGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_src = from_union([PurpleOtherSourceSrc.from_dict, from_none], obj.get("Other-source_src"))
        other_source_pre_text = from_union([from_str, from_none], obj.get("Other-source_pre-text"))
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        other_source_post_text = from_union([from_str, from_none], obj.get("Other-source_post-text"))
        return AmbitiousGeneCommentarySource(other_source_src, other_source_pre_text, other_source_anchor, other_source_post_text)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_src is not None:
            result["Other-source_src"] = from_union([lambda x: to_class(PurpleOtherSourceSrc, x), from_none], self.other_source_src)
        if self.other_source_pre_text is not None:
            result["Other-source_pre-text"] = from_union([from_str, from_none], self.other_source_pre_text)
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        if self.other_source_post_text is not None:
            result["Other-source_post-text"] = from_union([from_str, from_none], self.other_source_post_text)
        return result


class AmbitiousGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_refs: Optional[List[GeneCommentaryRef]]
    gene_commentary_source: Optional[List[AmbitiousGeneCommentarySource]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_refs: Optional[List[GeneCommentaryRef]], gene_commentary_source: Optional[List[AmbitiousGeneCommentarySource]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_refs = gene_commentary_refs
        self.gene_commentary_source = gene_commentary_source

    @staticmethod
    def from_dict(obj: Any) -> 'AmbitiousGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_refs = from_union([lambda x: from_list(GeneCommentaryRef.from_dict, x), from_none], obj.get("Gene-commentary_refs"))
        gene_commentary_source = from_union([lambda x: from_list(AmbitiousGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        return AmbitiousGeneCommentaryComment(gene_commentary_type, gene_commentary_refs, gene_commentary_source)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_refs is not None:
            result["Gene-commentary_refs"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentaryRef, x), x), from_none], self.gene_commentary_refs)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(AmbitiousGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        return result


class EntrezgenePropertyGeneCommentaryComment:
    gene_commentary_type: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_comment: Optional[List[AmbitiousGeneCommentaryComment]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_label: Optional[str], gene_commentary_comment: Optional[List[AmbitiousGeneCommentaryComment]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_comment = gene_commentary_comment

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgenePropertyGeneCommentaryComment':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_comment = from_union([lambda x: from_list(AmbitiousGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        return EntrezgenePropertyGeneCommentaryComment(gene_commentary_type, gene_commentary_label, gene_commentary_comment)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(AmbitiousGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        return result


class GeneCommentaryProperty:
    gene_commentary_type: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_text: Optional[str]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_label: Optional[str], gene_commentary_text: Optional[str]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_text = gene_commentary_text

    @staticmethod
    def from_dict(obj: Any) -> 'GeneCommentaryProperty':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        return GeneCommentaryProperty(gene_commentary_type, gene_commentary_label, gene_commentary_text)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        return result


class EntrezgenePropertyGeneCommentarySource:
    other_source_anchor: Optional[str]
    other_source_pre_text: Optional[str]
    other_source_url: Optional[str]

    def __init__(self, other_source_anchor: Optional[str], other_source_pre_text: Optional[str], other_source_url: Optional[str]) -> None:
        self.other_source_anchor = other_source_anchor
        self.other_source_pre_text = other_source_pre_text
        self.other_source_url = other_source_url

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgenePropertyGeneCommentarySource':
        assert isinstance(obj, dict)
        other_source_anchor = from_union([from_str, from_none], obj.get("Other-source_anchor"))
        other_source_pre_text = from_union([from_str, from_none], obj.get("Other-source_pre-text"))
        other_source_url = from_union([from_str, from_none], obj.get("Other-source_url"))
        return EntrezgenePropertyGeneCommentarySource(other_source_anchor, other_source_pre_text, other_source_url)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.other_source_anchor is not None:
            result["Other-source_anchor"] = from_union([from_str, from_none], self.other_source_anchor)
        if self.other_source_pre_text is not None:
            result["Other-source_pre-text"] = from_union([from_str, from_none], self.other_source_pre_text)
        if self.other_source_url is not None:
            result["Other-source_url"] = from_union([from_str, from_none], self.other_source_url)
        return result


class EntrezgeneProperty:
    gene_commentary_type: Optional[str]
    gene_commentary_label: Optional[str]
    gene_commentary_source: Optional[List[EntrezgenePropertyGeneCommentarySource]]
    gene_commentary_properties: Optional[List[GeneCommentaryProperty]]
    gene_commentary_text: Optional[str]
    gene_commentary_heading: Optional[str]
    gene_commentary_comment: Optional[List[EntrezgenePropertyGeneCommentaryComment]]

    def __init__(self, gene_commentary_type: Optional[str], gene_commentary_label: Optional[str], gene_commentary_source: Optional[List[EntrezgenePropertyGeneCommentarySource]], gene_commentary_properties: Optional[List[GeneCommentaryProperty]], gene_commentary_text: Optional[str], gene_commentary_heading: Optional[str], gene_commentary_comment: Optional[List[EntrezgenePropertyGeneCommentaryComment]]) -> None:
        self.gene_commentary_type = gene_commentary_type
        self.gene_commentary_label = gene_commentary_label
        self.gene_commentary_source = gene_commentary_source
        self.gene_commentary_properties = gene_commentary_properties
        self.gene_commentary_text = gene_commentary_text
        self.gene_commentary_heading = gene_commentary_heading
        self.gene_commentary_comment = gene_commentary_comment

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneProperty':
        assert isinstance(obj, dict)
        gene_commentary_type = from_union([from_str, from_none], obj.get("Gene-commentary_type"))
        gene_commentary_label = from_union([from_str, from_none], obj.get("Gene-commentary_label"))
        gene_commentary_source = from_union([lambda x: from_list(EntrezgenePropertyGeneCommentarySource.from_dict, x), from_none], obj.get("Gene-commentary_source"))
        gene_commentary_properties = from_union([lambda x: from_list(GeneCommentaryProperty.from_dict, x), from_none], obj.get("Gene-commentary_properties"))
        gene_commentary_text = from_union([from_str, from_none], obj.get("Gene-commentary_text"))
        gene_commentary_heading = from_union([from_str, from_none], obj.get("Gene-commentary_heading"))
        gene_commentary_comment = from_union([lambda x: from_list(EntrezgenePropertyGeneCommentaryComment.from_dict, x), from_none], obj.get("Gene-commentary_comment"))
        return EntrezgeneProperty(gene_commentary_type, gene_commentary_label, gene_commentary_source, gene_commentary_properties, gene_commentary_text, gene_commentary_heading, gene_commentary_comment)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_commentary_type is not None:
            result["Gene-commentary_type"] = from_union([from_str, from_none], self.gene_commentary_type)
        if self.gene_commentary_label is not None:
            result["Gene-commentary_label"] = from_union([from_str, from_none], self.gene_commentary_label)
        if self.gene_commentary_source is not None:
            result["Gene-commentary_source"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgenePropertyGeneCommentarySource, x), x), from_none], self.gene_commentary_source)
        if self.gene_commentary_properties is not None:
            result["Gene-commentary_properties"] = from_union([lambda x: from_list(lambda x: to_class(GeneCommentaryProperty, x), x), from_none], self.gene_commentary_properties)
        if self.gene_commentary_text is not None:
            result["Gene-commentary_text"] = from_union([from_str, from_none], self.gene_commentary_text)
        if self.gene_commentary_heading is not None:
            result["Gene-commentary_heading"] = from_union([from_str, from_none], self.gene_commentary_heading)
        if self.gene_commentary_comment is not None:
            result["Gene-commentary_comment"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgenePropertyGeneCommentaryComment, x), x), from_none], self.gene_commentary_comment)
        return result


class ProtRef:
    prot_ref_name: Optional[List[str]]
    prot_ref_desc: Optional[str]

    def __init__(self, prot_ref_name: Optional[List[str]], prot_ref_desc: Optional[str]) -> None:
        self.prot_ref_name = prot_ref_name
        self.prot_ref_desc = prot_ref_desc

    @staticmethod
    def from_dict(obj: Any) -> 'ProtRef':
        assert isinstance(obj, dict)
        prot_ref_name = from_union([lambda x: from_list(from_str, x), from_none], obj.get("Prot-ref_name"))
        prot_ref_desc = from_union([from_str, from_none], obj.get("Prot-ref_desc"))
        return ProtRef(prot_ref_name, prot_ref_desc)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.prot_ref_name is not None:
            result["Prot-ref_name"] = from_union([lambda x: from_list(from_str, x), from_none], self.prot_ref_name)
        if self.prot_ref_desc is not None:
            result["Prot-ref_desc"] = from_union([from_str, from_none], self.prot_ref_desc)
        return result


class EntrezgeneProt:
    prot_ref: Optional[ProtRef]

    def __init__(self, prot_ref: Optional[ProtRef]) -> None:
        self.prot_ref = prot_ref

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneProt':
        assert isinstance(obj, dict)
        prot_ref = from_union([ProtRef.from_dict, from_none], obj.get("Prot-ref"))
        return EntrezgeneProt(prot_ref)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.prot_ref is not None:
            result["Prot-ref"] = from_union([lambda x: to_class(ProtRef, x), from_none], self.prot_ref)
        return result


class BinomialOrgName:
    binomial_org_name_genus: Optional[str]
    binomial_org_name_species: Optional[str]

    def __init__(self, binomial_org_name_genus: Optional[str], binomial_org_name_species: Optional[str]) -> None:
        self.binomial_org_name_genus = binomial_org_name_genus
        self.binomial_org_name_species = binomial_org_name_species

    @staticmethod
    def from_dict(obj: Any) -> 'BinomialOrgName':
        assert isinstance(obj, dict)
        binomial_org_name_genus = from_union([from_str, from_none], obj.get("BinomialOrgName_genus"))
        binomial_org_name_species = from_union([from_str, from_none], obj.get("BinomialOrgName_species"))
        return BinomialOrgName(binomial_org_name_genus, binomial_org_name_species)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.binomial_org_name_genus is not None:
            result["BinomialOrgName_genus"] = from_union([from_str, from_none], self.binomial_org_name_genus)
        if self.binomial_org_name_species is not None:
            result["BinomialOrgName_species"] = from_union([from_str, from_none], self.binomial_org_name_species)
        return result


class OrgNameNameBinomial:
    binomial_org_name: Optional[BinomialOrgName]

    def __init__(self, binomial_org_name: Optional[BinomialOrgName]) -> None:
        self.binomial_org_name = binomial_org_name

    @staticmethod
    def from_dict(obj: Any) -> 'OrgNameNameBinomial':
        assert isinstance(obj, dict)
        binomial_org_name = from_union([BinomialOrgName.from_dict, from_none], obj.get("BinomialOrgName"))
        return OrgNameNameBinomial(binomial_org_name)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.binomial_org_name is not None:
            result["BinomialOrgName"] = from_union([lambda x: to_class(BinomialOrgName, x), from_none], self.binomial_org_name)
        return result


class OrgNameName:
    org_name_name_binomial: Optional[OrgNameNameBinomial]

    def __init__(self, org_name_name_binomial: Optional[OrgNameNameBinomial]) -> None:
        self.org_name_name_binomial = org_name_name_binomial

    @staticmethod
    def from_dict(obj: Any) -> 'OrgNameName':
        assert isinstance(obj, dict)
        org_name_name_binomial = from_union([OrgNameNameBinomial.from_dict, from_none], obj.get("OrgName_name_binomial"))
        return OrgNameName(org_name_name_binomial)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.org_name_name_binomial is not None:
            result["OrgName_name_binomial"] = from_union([lambda x: to_class(OrgNameNameBinomial, x), from_none], self.org_name_name_binomial)
        return result


class OrgName:
    org_name_name: Optional[OrgNameName]
    org_name_attrib: Optional[str]
    org_name_lineage: Optional[str]
    org_name_gcode: Optional[str]
    org_name_mgcode: Optional[str]
    org_name_div: Optional[str]

    def __init__(self, org_name_name: Optional[OrgNameName], org_name_attrib: Optional[str], org_name_lineage: Optional[str], org_name_gcode: Optional[str], org_name_mgcode: Optional[str], org_name_div: Optional[str]) -> None:
        self.org_name_name = org_name_name
        self.org_name_attrib = org_name_attrib
        self.org_name_lineage = org_name_lineage
        self.org_name_gcode = org_name_gcode
        self.org_name_mgcode = org_name_mgcode
        self.org_name_div = org_name_div

    @staticmethod
    def from_dict(obj: Any) -> 'OrgName':
        assert isinstance(obj, dict)
        org_name_name = from_union([OrgNameName.from_dict, from_none], obj.get("OrgName_name"))
        org_name_attrib = from_union([from_str, from_none], obj.get("OrgName_attrib"))
        org_name_lineage = from_union([from_str, from_none], obj.get("OrgName_lineage"))
        org_name_gcode = from_union([from_str, from_none], obj.get("OrgName_gcode"))
        org_name_mgcode = from_union([from_str, from_none], obj.get("OrgName_mgcode"))
        org_name_div = from_union([from_str, from_none], obj.get("OrgName_div"))
        return OrgName(org_name_name, org_name_attrib, org_name_lineage, org_name_gcode, org_name_mgcode, org_name_div)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.org_name_name is not None:
            result["OrgName_name"] = from_union([lambda x: to_class(OrgNameName, x), from_none], self.org_name_name)
        if self.org_name_attrib is not None:
            result["OrgName_attrib"] = from_union([from_str, from_none], self.org_name_attrib)
        if self.org_name_lineage is not None:
            result["OrgName_lineage"] = from_union([from_str, from_none], self.org_name_lineage)
        if self.org_name_gcode is not None:
            result["OrgName_gcode"] = from_union([from_str, from_none], self.org_name_gcode)
        if self.org_name_mgcode is not None:
            result["OrgName_mgcode"] = from_union([from_str, from_none], self.org_name_mgcode)
        if self.org_name_div is not None:
            result["OrgName_div"] = from_union([from_str, from_none], self.org_name_div)
        return result


class OrgRefOrgname:
    org_name: Optional[OrgName]

    def __init__(self, org_name: Optional[OrgName]) -> None:
        self.org_name = org_name

    @staticmethod
    def from_dict(obj: Any) -> 'OrgRefOrgname':
        assert isinstance(obj, dict)
        org_name = from_union([OrgName.from_dict, from_none], obj.get("OrgName"))
        return OrgRefOrgname(org_name)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.org_name is not None:
            result["OrgName"] = from_union([lambda x: to_class(OrgName, x), from_none], self.org_name)
        return result


class OrgRef:
    org_ref_taxname: Optional[str]
    org_ref_common: Optional[str]
    org_ref_db: Optional[List[OrgRefDB]]
    org_ref_orgname: Optional[OrgRefOrgname]

    def __init__(self, org_ref_taxname: Optional[str], org_ref_common: Optional[str], org_ref_db: Optional[List[OrgRefDB]], org_ref_orgname: Optional[OrgRefOrgname]) -> None:
        self.org_ref_taxname = org_ref_taxname
        self.org_ref_common = org_ref_common
        self.org_ref_db = org_ref_db
        self.org_ref_orgname = org_ref_orgname

    @staticmethod
    def from_dict(obj: Any) -> 'OrgRef':
        assert isinstance(obj, dict)
        org_ref_taxname = from_union([from_str, from_none], obj.get("Org-ref_taxname"))
        org_ref_common = from_union([from_str, from_none], obj.get("Org-ref_common"))
        org_ref_db = from_union([lambda x: from_list(OrgRefDB.from_dict, x), from_none], obj.get("Org-ref_db"))
        org_ref_orgname = from_union([OrgRefOrgname.from_dict, from_none], obj.get("Org-ref_orgname"))
        return OrgRef(org_ref_taxname, org_ref_common, org_ref_db, org_ref_orgname)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.org_ref_taxname is not None:
            result["Org-ref_taxname"] = from_union([from_str, from_none], self.org_ref_taxname)
        if self.org_ref_common is not None:
            result["Org-ref_common"] = from_union([from_str, from_none], self.org_ref_common)
        if self.org_ref_db is not None:
            result["Org-ref_db"] = from_union([lambda x: from_list(lambda x: to_class(OrgRefDB, x), x), from_none], self.org_ref_db)
        if self.org_ref_orgname is not None:
            result["Org-ref_orgname"] = from_union([lambda x: to_class(OrgRefOrgname, x), from_none], self.org_ref_orgname)
        return result


class BioSourceOrg:
    org_ref: Optional[OrgRef]

    def __init__(self, org_ref: Optional[OrgRef]) -> None:
        self.org_ref = org_ref

    @staticmethod
    def from_dict(obj: Any) -> 'BioSourceOrg':
        assert isinstance(obj, dict)
        org_ref = from_union([OrgRef.from_dict, from_none], obj.get("Org-ref"))
        return BioSourceOrg(org_ref)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.org_ref is not None:
            result["Org-ref"] = from_union([lambda x: to_class(OrgRef, x), from_none], self.org_ref)
        return result


class BioSourceSubtype:
    sub_source_subtype: Optional[str]
    sub_source_name: Optional[str]

    def __init__(self, sub_source_subtype: Optional[str], sub_source_name: Optional[str]) -> None:
        self.sub_source_subtype = sub_source_subtype
        self.sub_source_name = sub_source_name

    @staticmethod
    def from_dict(obj: Any) -> 'BioSourceSubtype':
        assert isinstance(obj, dict)
        sub_source_subtype = from_union([from_str, from_none], obj.get("SubSource_subtype"))
        sub_source_name = from_union([from_str, from_none], obj.get("SubSource_name"))
        return BioSourceSubtype(sub_source_subtype, sub_source_name)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.sub_source_subtype is not None:
            result["SubSource_subtype"] = from_union([from_str, from_none], self.sub_source_subtype)
        if self.sub_source_name is not None:
            result["SubSource_name"] = from_union([from_str, from_none], self.sub_source_name)
        return result


class BioSource:
    bio_source_genome: Optional[str]
    bio_source_origin: Optional[str]
    bio_source_org: Optional[BioSourceOrg]
    bio_source_subtype: Optional[List[BioSourceSubtype]]

    def __init__(self, bio_source_genome: Optional[str], bio_source_origin: Optional[str], bio_source_org: Optional[BioSourceOrg], bio_source_subtype: Optional[List[BioSourceSubtype]]) -> None:
        self.bio_source_genome = bio_source_genome
        self.bio_source_origin = bio_source_origin
        self.bio_source_org = bio_source_org
        self.bio_source_subtype = bio_source_subtype

    @staticmethod
    def from_dict(obj: Any) -> 'BioSource':
        assert isinstance(obj, dict)
        bio_source_genome = from_union([from_str, from_none], obj.get("BioSource_genome"))
        bio_source_origin = from_union([from_str, from_none], obj.get("BioSource_origin"))
        bio_source_org = from_union([BioSourceOrg.from_dict, from_none], obj.get("BioSource_org"))
        bio_source_subtype = from_union([lambda x: from_list(BioSourceSubtype.from_dict, x), from_none], obj.get("BioSource_subtype"))
        return BioSource(bio_source_genome, bio_source_origin, bio_source_org, bio_source_subtype)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.bio_source_genome is not None:
            result["BioSource_genome"] = from_union([from_str, from_none], self.bio_source_genome)
        if self.bio_source_origin is not None:
            result["BioSource_origin"] = from_union([from_str, from_none], self.bio_source_origin)
        if self.bio_source_org is not None:
            result["BioSource_org"] = from_union([lambda x: to_class(BioSourceOrg, x), from_none], self.bio_source_org)
        if self.bio_source_subtype is not None:
            result["BioSource_subtype"] = from_union([lambda x: from_list(lambda x: to_class(BioSourceSubtype, x), x), from_none], self.bio_source_subtype)
        return result


class EntrezgeneSource:
    bio_source: Optional[BioSource]

    def __init__(self, bio_source: Optional[BioSource]) -> None:
        self.bio_source = bio_source

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneSource':
        assert isinstance(obj, dict)
        bio_source = from_union([BioSource.from_dict, from_none], obj.get("BioSource"))
        return EntrezgeneSource(bio_source)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.bio_source is not None:
            result["BioSource"] = from_union([lambda x: to_class(BioSource, x), from_none], self.bio_source)
        return result


class StickyDateStd:
    date_std_year: Optional[str]
    date_std_month: Optional[str]
    date_std_day: Optional[str]

    def __init__(self, date_std_year: Optional[str], date_std_month: Optional[str], date_std_day: Optional[str]) -> None:
        self.date_std_year = date_std_year
        self.date_std_month = date_std_month
        self.date_std_day = date_std_day

    @staticmethod
    def from_dict(obj: Any) -> 'StickyDateStd':
        assert isinstance(obj, dict)
        date_std_year = from_union([from_str, from_none], obj.get("Date-std_year"))
        date_std_month = from_union([from_str, from_none], obj.get("Date-std_month"))
        date_std_day = from_union([from_str, from_none], obj.get("Date-std_day"))
        return StickyDateStd(date_std_year, date_std_month, date_std_day)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.date_std_year is not None:
            result["Date-std_year"] = from_union([from_str, from_none], self.date_std_year)
        if self.date_std_month is not None:
            result["Date-std_month"] = from_union([from_str, from_none], self.date_std_month)
        if self.date_std_day is not None:
            result["Date-std_day"] = from_union([from_str, from_none], self.date_std_day)
        return result


class TentacledDateStd:
    date_std: Optional[StickyDateStd]

    def __init__(self, date_std: Optional[StickyDateStd]) -> None:
        self.date_std = date_std

    @staticmethod
    def from_dict(obj: Any) -> 'TentacledDateStd':
        assert isinstance(obj, dict)
        date_std = from_union([StickyDateStd.from_dict, from_none], obj.get("Date-std"))
        return TentacledDateStd(date_std)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.date_std is not None:
            result["Date-std"] = from_union([lambda x: to_class(StickyDateStd, x), from_none], self.date_std)
        return result


class GeneTrackCreateDateDate:
    date_std: Optional[TentacledDateStd]

    def __init__(self, date_std: Optional[TentacledDateStd]) -> None:
        self.date_std = date_std

    @staticmethod
    def from_dict(obj: Any) -> 'GeneTrackCreateDateDate':
        assert isinstance(obj, dict)
        date_std = from_union([TentacledDateStd.from_dict, from_none], obj.get("Date_std"))
        return GeneTrackCreateDateDate(date_std)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.date_std is not None:
            result["Date_std"] = from_union([lambda x: to_class(TentacledDateStd, x), from_none], self.date_std)
        return result


class GeneTrackCreateDate:
    date: Optional[GeneTrackCreateDateDate]

    def __init__(self, date: Optional[GeneTrackCreateDateDate]) -> None:
        self.date = date

    @staticmethod
    def from_dict(obj: Any) -> 'GeneTrackCreateDate':
        assert isinstance(obj, dict)
        date = from_union([GeneTrackCreateDateDate.from_dict, from_none], obj.get("Date"))
        return GeneTrackCreateDate(date)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.date is not None:
            result["Date"] = from_union([lambda x: to_class(GeneTrackCreateDateDate, x), from_none], self.date)
        return result


class GeneTrack:
    gene_track_geneid: Optional[str]
    gene_track_status: Optional[str]
    gene_track_create_date: Optional[GeneTrackCreateDate]
    gene_track_update_date: Optional[GeneAteDate]

    def __init__(self, gene_track_geneid: Optional[str], gene_track_status: Optional[str], gene_track_create_date: Optional[GeneTrackCreateDate], gene_track_update_date: Optional[GeneAteDate]) -> None:
        self.gene_track_geneid = gene_track_geneid
        self.gene_track_status = gene_track_status
        self.gene_track_create_date = gene_track_create_date
        self.gene_track_update_date = gene_track_update_date

    @staticmethod
    def from_dict(obj: Any) -> 'GeneTrack':
        assert isinstance(obj, dict)
        gene_track_geneid = from_union([from_str, from_none], obj.get("Gene-track_geneid"))
        gene_track_status = from_union([from_str, from_none], obj.get("Gene-track_status"))
        gene_track_create_date = from_union([GeneTrackCreateDate.from_dict, from_none], obj.get("Gene-track_create-date"))
        gene_track_update_date = from_union([GeneAteDate.from_dict, from_none], obj.get("Gene-track_update-date"))
        return GeneTrack(gene_track_geneid, gene_track_status, gene_track_create_date, gene_track_update_date)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_track_geneid is not None:
            result["Gene-track_geneid"] = from_union([from_str, from_none], self.gene_track_geneid)
        if self.gene_track_status is not None:
            result["Gene-track_status"] = from_union([from_str, from_none], self.gene_track_status)
        if self.gene_track_create_date is not None:
            result["Gene-track_create-date"] = from_union([lambda x: to_class(GeneTrackCreateDate, x), from_none], self.gene_track_create_date)
        if self.gene_track_update_date is not None:
            result["Gene-track_update-date"] = from_union([lambda x: to_class(GeneAteDate, x), from_none], self.gene_track_update_date)
        return result


class EntrezgeneTrackInfo:
    gene_track: Optional[GeneTrack]

    def __init__(self, gene_track: Optional[GeneTrack]) -> None:
        self.gene_track = gene_track

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezgeneTrackInfo':
        assert isinstance(obj, dict)
        gene_track = from_union([GeneTrack.from_dict, from_none], obj.get("Gene-track"))
        return EntrezgeneTrackInfo(gene_track)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.gene_track is not None:
            result["Gene-track"] = from_union([lambda x: to_class(GeneTrack, x), from_none], self.gene_track)
        return result


class EntrezGeneElement:
    entrezgene_track_info: Optional[EntrezgeneTrackInfo]
    entrezgene_type: Optional[str]
    entrezgene_source: Optional[EntrezgeneSource]
    entrezgene_gene: Optional[EntrezgeneGene]
    entrezgene_prot: Optional[EntrezgeneProt]
    entrezgene_summary: Optional[str]
    entrezgene_location: Optional[List[EntrezgeneLocation]]
    entrezgene_gene_source: Optional[EntrezgeneGeneSource]
    entrezgene_locus: Optional[List[EntrezgeneLocus]]
    entrezgene_properties: Optional[List[EntrezgeneProperty]]
    entrezgene_comments: Optional[List[EntrezgeneComment]]
    entrezgene_unique_keys: Optional[List[EntrezgeneUniqueKey]]
    entrezgene_xtra_index_terms: Optional[List[str]]
    entrezgene_xtra_properties: Optional[List[XtraProperty]]

    def __init__(self, entrezgene_track_info: Optional[EntrezgeneTrackInfo], entrezgene_type: Optional[str], entrezgene_source: Optional[EntrezgeneSource], entrezgene_gene: Optional[EntrezgeneGene], entrezgene_prot: Optional[EntrezgeneProt], entrezgene_summary: Optional[str], entrezgene_location: Optional[List[EntrezgeneLocation]], entrezgene_gene_source: Optional[EntrezgeneGeneSource], entrezgene_locus: Optional[List[EntrezgeneLocus]], entrezgene_properties: Optional[List[EntrezgeneProperty]], entrezgene_comments: Optional[List[EntrezgeneComment]], entrezgene_unique_keys: Optional[List[EntrezgeneUniqueKey]], entrezgene_xtra_index_terms: Optional[List[str]], entrezgene_xtra_properties: Optional[List[XtraProperty]]) -> None:
        self.entrezgene_track_info = entrezgene_track_info
        self.entrezgene_type = entrezgene_type
        self.entrezgene_source = entrezgene_source
        self.entrezgene_gene = entrezgene_gene
        self.entrezgene_prot = entrezgene_prot
        self.entrezgene_summary = entrezgene_summary
        self.entrezgene_location = entrezgene_location
        self.entrezgene_gene_source = entrezgene_gene_source
        self.entrezgene_locus = entrezgene_locus
        self.entrezgene_properties = entrezgene_properties
        self.entrezgene_comments = entrezgene_comments
        self.entrezgene_unique_keys = entrezgene_unique_keys
        self.entrezgene_xtra_index_terms = entrezgene_xtra_index_terms
        self.entrezgene_xtra_properties = entrezgene_xtra_properties

    @staticmethod
    def from_dict(obj: Any) -> 'EntrezGeneElement':
        assert isinstance(obj, dict)
        entrezgene_track_info = from_union([EntrezgeneTrackInfo.from_dict, from_none], obj.get("Entrezgene_track-info"))
        entrezgene_type = from_union([from_str, from_none], obj.get("Entrezgene_type"))
        entrezgene_source = from_union([EntrezgeneSource.from_dict, from_none], obj.get("Entrezgene_source"))
        entrezgene_gene = from_union([EntrezgeneGene.from_dict, from_none], obj.get("Entrezgene_gene"))
        entrezgene_prot = from_union([EntrezgeneProt.from_dict, from_none], obj.get("Entrezgene_prot"))
        entrezgene_summary = from_union([from_str, from_none], obj.get("Entrezgene_summary"))
        entrezgene_location = from_union([lambda x: from_list(EntrezgeneLocation.from_dict, x), from_none], obj.get("Entrezgene_location"))
        entrezgene_gene_source = from_union([EntrezgeneGeneSource.from_dict, from_none], obj.get("Entrezgene_gene-source"))
        entrezgene_locus = from_union([lambda x: from_list(EntrezgeneLocus.from_dict, x), from_none], obj.get("Entrezgene_locus"))
        entrezgene_properties = from_union([lambda x: from_list(EntrezgeneProperty.from_dict, x), from_none], obj.get("Entrezgene_properties"))
        entrezgene_comments = from_union([lambda x: from_list(EntrezgeneComment.from_dict, x), from_none], obj.get("Entrezgene_comments"))
        entrezgene_unique_keys = from_union([lambda x: from_list(EntrezgeneUniqueKey.from_dict, x), from_none], obj.get("Entrezgene_unique-keys"))
        entrezgene_xtra_index_terms = from_union([lambda x: from_list(from_str, x), from_none], obj.get("Entrezgene_xtra-index-terms"))
        entrezgene_xtra_properties = from_union([lambda x: from_list(XtraProperty.from_dict, x), from_none], obj.get("Entrezgene_xtra-properties"))
        return EntrezGeneElement(entrezgene_track_info, entrezgene_type, entrezgene_source, entrezgene_gene, entrezgene_prot, entrezgene_summary, entrezgene_location, entrezgene_gene_source, entrezgene_locus, entrezgene_properties, entrezgene_comments, entrezgene_unique_keys, entrezgene_xtra_index_terms, entrezgene_xtra_properties)

    def to_dict(self) -> dict:
        result: dict = {}
        if self.entrezgene_track_info is not None:
            result["Entrezgene_track-info"] = from_union([lambda x: to_class(EntrezgeneTrackInfo, x), from_none], self.entrezgene_track_info)
        if self.entrezgene_type is not None:
            result["Entrezgene_type"] = from_union([from_str, from_none], self.entrezgene_type)
        if self.entrezgene_source is not None:
            result["Entrezgene_source"] = from_union([lambda x: to_class(EntrezgeneSource, x), from_none], self.entrezgene_source)
        if self.entrezgene_gene is not None:
            result["Entrezgene_gene"] = from_union([lambda x: to_class(EntrezgeneGene, x), from_none], self.entrezgene_gene)
        if self.entrezgene_prot is not None:
            result["Entrezgene_prot"] = from_union([lambda x: to_class(EntrezgeneProt, x), from_none], self.entrezgene_prot)
        if self.entrezgene_summary is not None:
            result["Entrezgene_summary"] = from_union([from_str, from_none], self.entrezgene_summary)
        if self.entrezgene_location is not None:
            result["Entrezgene_location"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneLocation, x), x), from_none], self.entrezgene_location)
        if self.entrezgene_gene_source is not None:
            result["Entrezgene_gene-source"] = from_union([lambda x: to_class(EntrezgeneGeneSource, x), from_none], self.entrezgene_gene_source)
        if self.entrezgene_locus is not None:
            result["Entrezgene_locus"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneLocus, x), x), from_none], self.entrezgene_locus)
        if self.entrezgene_properties is not None:
            result["Entrezgene_properties"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneProperty, x), x), from_none], self.entrezgene_properties)
        if self.entrezgene_comments is not None:
            result["Entrezgene_comments"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneComment, x), x), from_none], self.entrezgene_comments)
        if self.entrezgene_unique_keys is not None:
            result["Entrezgene_unique-keys"] = from_union([lambda x: from_list(lambda x: to_class(EntrezgeneUniqueKey, x), x), from_none], self.entrezgene_unique_keys)
        if self.entrezgene_xtra_index_terms is not None:
            result["Entrezgene_xtra-index-terms"] = from_union([lambda x: from_list(from_str, x), from_none], self.entrezgene_xtra_index_terms)
        if self.entrezgene_xtra_properties is not None:
            result["Entrezgene_xtra-properties"] = from_union([lambda x: from_list(lambda x: to_class(XtraProperty, x), x), from_none], self.entrezgene_xtra_properties)
        return result


def entrez_gene_from_dict(s: Any) -> List[EntrezGeneElement]:
    return from_list(EntrezGeneElement.from_dict, s)


def entrez_gene_to_dict(x: List[EntrezGeneElement]) -> Any:
    return from_list(lambda x: to_class(EntrezGeneElement, x), x)
