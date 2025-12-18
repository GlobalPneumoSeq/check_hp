from check_hp import get_gap
import pytest


# Mock Classes

class MockFeature:
    def __init__(self, start, end, strand):
        self.start = start
        self.end = end
        self.strand = strand


# Tests 

def test_get_gap_plus_strand_positive_gap():
    feature = MockFeature(start=201, end=300, strand="+")
    adjacent_features = {"upstream_end": 100, "downstream_start": 301}
    assert get_gap(feature, adjacent_features) == 100

def test_get_gap_plus_strand_negative_gap():
    feature = MockFeature(start=201, end=300, strand="+")
    adjacent_features = {"upstream_end": 250, "downstream_start": 301}
    assert get_gap(feature, adjacent_features) == -50

def test_get_gap_plus_strand_no_gap():
    feature = MockFeature(start=201, end=300, strand="+")
    adjacent_features = {"upstream_end": 200, "downstream_start": 301}
    assert get_gap(feature, adjacent_features) == 0

def test_get_gap_minus_strand_positive_gap():
    feature = MockFeature(start=201, end=300, strand="-")
    adjacent_features = {"upstream_end": 200, "downstream_start": 401}
    assert get_gap(feature, adjacent_features) == 100

def test_get_gap_minus_strand_negative_gap():
    feature = MockFeature(start=201, end=300, strand="-")
    adjacent_features = {"upstream_end": 200, "downstream_start": 251}
    assert get_gap(feature, adjacent_features) == -50

def test_get_gap_minus_strand_no_gap():
    feature = MockFeature(start=201, end=300, strand="-")
    adjacent_features = {"upstream_end": 200, "downstream_start": 301}
    assert get_gap(feature, adjacent_features) == 0

def test_get_gap_invalid_strand():
    feature = MockFeature(start=100, end=150, strand="x")
    adjacent_features = {"upstream_end": 50, "downstream_start": 200}
    with pytest.raises(ValueError):
        get_gap(feature, adjacent_features)
