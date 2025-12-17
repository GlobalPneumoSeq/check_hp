from check_hp import has_promoter

class MockFeature:
    def __init__(self, start, end, strand):
        self.start = start
        self.end = end
        self.strand = strand

def test_has_promoter_plus_strand_positive_100bp():
    feature = MockFeature(start=201, end=300, strand="+")
    adjacent_features = {"upstream_end": 100, "downstream_start": 301}
    promoter_length = 100
    assert has_promoter(feature, adjacent_features, promoter_length) == True

def test_has_promoter_plus_strand_positive_300bp():
    feature = MockFeature(start=401, end=300, strand="+")
    adjacent_features = {"upstream_end": 100, "downstream_start": 301}
    promoter_length = 100
    assert has_promoter(feature, adjacent_features, promoter_length) == True

def test_has_promoter_plus_strand_negative_99bp():
    feature = MockFeature(start=201, end=300, strand="+")
    adjacent_features = {"upstream_end": 101, "downstream_start": 301}
    promoter_length = 100
    assert has_promoter(feature, adjacent_features, promoter_length) == False

def test_has_promoter_minus_strand_positive_100bp():
    feature = MockFeature(start=201, end=300, strand="-")
    adjacent_features = {"upstream_end": 200, "downstream_start": 401}
    promoter_length = 100
    assert has_promoter(feature, adjacent_features, promoter_length) == True

def test_has_promoter_minus_strand_positive_300bp():
    feature = MockFeature(start=201, end=300, strand="-")
    adjacent_features = {"upstream_end": 200, "downstream_start": 601}
    promoter_length = 100
    assert has_promoter(feature, adjacent_features, promoter_length) == True

def test_has_promoter_minus_strand_negative_99bp():
    feature = MockFeature(start=201, end=300, strand="-")
    adjacent_features = {"upstream_end": 200, "downstream_start": 400}
    promoter_length = 100
    assert has_promoter(feature, adjacent_features, promoter_length) == False
