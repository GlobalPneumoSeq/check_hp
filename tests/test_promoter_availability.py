from intervaltree import IntervalTree
import pytest

from check_hp import promoter_feature_overlapped, promoter_within_contig, validate_feature_strand


# Mock Classes

class MockFeature:
    def __init__(self, start, end, strand, feature_id="mock_id"):
        self.start = start
        self.end = end
        self.strand = strand
        self.id = feature_id


def make_features_intervals(*features):
    # GFF3 coordinates are inclusive; IntervalTree intervals are half-open.
    return IntervalTree.from_tuples(
        (feature.start, feature.end + 1) for feature in features
    )


# Tests

@pytest.mark.parametrize(
    "testing_feature, blocking_feature, expected",
    [
        # blocking_feature on either strands which overlaps with promoter should return True otherwise False
        # promoter of the testing feature on plus strand at 101 - 200
        (MockFeature(201, 300, "+"), MockFeature(90, 100, "+"), False),
        (MockFeature(201, 300, "+"), MockFeature(90, 100, "-"), False),
        (MockFeature(201, 300, "+"), MockFeature(150, 160, "+"), True),
        (MockFeature(201, 300, "+"), MockFeature(150, 160, "-"), True),
        (MockFeature(201, 300, "+"), MockFeature(301, 400, "+"), False),
        (MockFeature(201, 300, "+"), MockFeature(301, 400, "-"), False),
        # promoter of the testing feature on minus strand at 301 - 400
        (MockFeature(201, 300, "-"), MockFeature(401, 410, "+"), False),
        (MockFeature(201, 300, "-"), MockFeature(401, 410, "-"), False),
        (MockFeature(201, 300, "-"), MockFeature(350, 360, "+"), True),
        (MockFeature(201, 300, "-"), MockFeature(350, 360, "-"), True),
        (MockFeature(201, 300, "-"), MockFeature(190, 200, "+"), False),
        (MockFeature(201, 300, "-"), MockFeature(190, 200, "-"), False),
        
        # The promoter's first and final base must be included in the overlap query
        # promoter of the testing feature on plus strand at 101 - 200
        (MockFeature(201, 300, "+"), MockFeature(101, 101, "+"), True),
        (MockFeature(201, 300, "+"), MockFeature(100, 100, "+"), False),
        (MockFeature(201, 300, "+"), MockFeature(200, 200, "+"), True),
        (MockFeature(201, 300, "+"), MockFeature(201, 201, "+"), False),
        # promoter of the testing feature on minus strand at 301 - 400
        (MockFeature(201, 300, "-"), MockFeature(301, 301, "+"), True),
        (MockFeature(201, 300, "-"), MockFeature(300, 300, "+"), False),
        (MockFeature(201, 300, "-"), MockFeature(400, 400, "+"), True),
        (MockFeature(201, 300, "-"), MockFeature(401, 401, "+"), False),
    ],
)
def test_promoter_feature_overlapped(testing_feature, blocking_feature, expected):
    features_intervals = make_features_intervals(testing_feature, blocking_feature)
    assert promoter_feature_overlapped(testing_feature, features_intervals, promoter_length=100) is expected


@pytest.mark.parametrize(
    "feature, contig_details, expected",
    [
        (MockFeature(101, 200, "+"), {"start": 1, "end": 500}, True),
        (MockFeature(100, 200, "+"), {"start": 1, "end": 500}, False),
        (MockFeature(201, 300, "-"), {"start": 1, "end": 400}, True),
        (MockFeature(201, 300, "-"), {"start": 1, "end": 399}, False),
    ],
)
def test_promoter_within_contig(feature, contig_details, expected):
    assert promoter_within_contig(feature, contig_details, promoter_length=100) is expected


@pytest.mark.parametrize("strand", ["+", "-"])
def test_validate_feature_strand_accepts_valid_strands(strand):
    validate_feature_strand(MockFeature(10, 20, strand), "contig_1")


def test_validate_feature_strand_reports_feature_and_contig():
    feature = MockFeature(10, 20, ".", feature_id="test_id")
    with pytest.raises(ValueError, match=r"Malformed GFF3: expected strand value of \+ or -, got \. for feature test_id on contig contig_1"):
        validate_feature_strand(feature, "contig_1")
