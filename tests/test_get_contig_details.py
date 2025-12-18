from check_hp import get_contig_details
import gffutils
from types import SimpleNamespace

# Mock Classes

class MockFastaSeq:
    def __init__(self, seq):
        self.seq = seq

    def __len__(self):
        return len(self.seq)
    
    def __str__(self):
        return self.seq

class MockFasta:
    def __init__(self, sequences):
        self.sequences = sequences
    
    def __getitem__(self, contig):
        return MockFastaSeq(self.sequences[contig])

class MockFeature:
    def __init__(self, start, end, strand, featuretype):
        self.start = start
        self.end = end
        self.strand = strand
        self.featuretype = featuretype

class MockDB:
    def __init__(self, contig_features=None, contig_bounds=None):
        self.contig_features = contig_features or {}
        self.contig_bounds = contig_bounds or {}
    
    def __getitem__(self, contig):
        if contig not in self.contig_bounds:
            raise gffutils.exceptions.FeatureNotFoundError(contig)
        return SimpleNamespace(start=self.contig_bounds[contig][0], end=self.contig_bounds[contig][1])
    
    def region(self, seqid, featuretype):
        for feature in self.contig_features.get(seqid):
             if feature.featuretype == featuretype:
                 yield feature


# Tests

def test_get_contig_details_with_contig_in_db():
    fasta = MockFasta({"contig1": "ATGC" * 25})

    features = [
        MockFeature(5, 20, "+", "CDS"),
        MockFeature(30, 50, "-", "CDS"),
        MockFeature(60, 80, "+", "gene")
    ]

    db = MockDB(
        contig_features={"contig1": features},
        contig_bounds={"contig1": (1, 150)},
    )

    result = get_contig_details(db, "contig1", fasta, featuretype="CDS")

    assert result["start"] == 1
    assert result["end"] == 150
    assert len(db.contig_features["contig1"]) == 3
    assert len(result["features"]) == 2
    assert all(feature.featuretype == "CDS" for feature in result["features"])
    assert str(result["fasta"]) == "ATGC" * 25

def test_get_contig_details_without_contig_in_db():
    fasta = MockFasta({"contig1": "ATGC" * 25})

    db = MockDB(
        contig_features={"contig1": []},
    )

    result = get_contig_details(db, "contig1", fasta, featuretype="CDS")

    assert result["start"] == 1
    assert result["end"] == 100

def test_get_contig_details_returns_empty_feature_list_if_no_matching_featuretype():
    fasta = MockFasta({"contig1": "ATGC" * 25})

    features = [
        MockFeature(1, 10, "+", "gene"),
        MockFeature(15, 30, "-", "gene"),
    ]

    db = MockDB(
        contig_features={"contig1": features},
        contig_bounds={"contig1": (1, 100)},
    )

    result = get_contig_details(db, "contig1", fasta, featuretype="CDS")

    assert result["features"] == []
