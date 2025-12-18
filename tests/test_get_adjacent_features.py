from check_hp import get_adjacent_features

# Mock Classes

class MockFeature:
    def __init__(self, start, end, strand):
        self.start = start
        self.end = end
        self.strand = strand


# Tests

def test_get_adjacent_features_middle():
    feature1 = MockFeature(1, 50, "+")
    feature2 = MockFeature(51, 100, "-")
    feature3 = MockFeature(101, 150, "+")
    contig_details = {"features": [feature1, feature2, feature3], "start": 1, "end": 200}

    adjacent_features = get_adjacent_features(1, contig_details)
    assert adjacent_features["upstream_feature"] == feature1
    assert adjacent_features["upstream_end"] == 50
    assert adjacent_features["downstream_feature"] == feature3
    assert adjacent_features["downstream_start"] == 101

def test_get_adjacent_features_first():
    feature1 = MockFeature(1, 50, "+")
    feature2 = MockFeature(51, 100, "-")
    contig_details = {"features": [feature1, feature2], "start": 1, "end": 200}

    adjacent_features = get_adjacent_features(0, contig_details)
    assert adjacent_features["upstream_feature"] is None
    assert adjacent_features["upstream_end"] == 0
    assert adjacent_features["downstream_feature"] == feature2
    assert adjacent_features["downstream_start"] == 51

def test_get_adjacent_features_last():
    feature1 = MockFeature(1, 50, "+")
    feature2 = MockFeature(51, 100, "-")
    contig_details = {"features": [feature1, feature2], "start": 1, "end": 200}

    adjacent_features = get_adjacent_features(1, contig_details)
    assert adjacent_features["upstream_feature"] == feature1
    assert adjacent_features["upstream_end"] == 50
    assert adjacent_features["downstream_feature"] is None
    assert adjacent_features["downstream_start"] == 201
