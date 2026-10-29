from check_hp import get_homopolymers


# Mock Classes

class SequenceSlice:
    def __init__(self, sequence):
        self.sequence = sequence

    def __str__(self):
        return self.sequence

    @property
    def reverse(self):
        return SequenceSlice(self.sequence[::-1])

    @property
    def complement(self):
        return SequenceSlice(self.sequence.translate(str.maketrans("ACGTUacgtu", "TGCAAtgcaa")))


class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence

    def __getitem__(self, item):
        return SequenceSlice(self.sequence[item])


class MockFeature:
    def __init__(self, start, end, strand):
        self.start = start
        self.end = end
        self.strand = strand


# Tests

def test_get_homopolymers_plus_strand_reports_correctly():
    # Positive strand, promoter occupies positions 11 - 20
    result = get_homopolymers(MockFeature(21, 30, "+"), Sequence("ATCGATCGATAAAAATCGATCGATCGATCG"), 10, "A", 3)

    assert result == [{
        "repeats": 5,
        "repeat_start_pos": 11,
        "repeat_end_pos": 15,
        "distance_from_cds": 5,
    }]


def test_get_homopolymers_minus_strand_reports_correctly():
    # Negative strand, promoter occupies positions 21 - 30, AAAA on the sequence becomes TTTT
    result = get_homopolymers(MockFeature(1, 20, "-"), Sequence("CGATCGATCGATCGATCGATAAAATTCGAT"), 10, "T", 4)

    assert result == [{
        "repeats": 4,
        "repeat_start_pos": 21,
        "repeat_end_pos": 24,
        "distance_from_cds": 0,
    }]


def test_get_homopolymers_excludes_runs_shorter_than_threshold():
    # Positive strand, promoter occupies positions 11 - 20
    result = get_homopolymers(MockFeature(21, 30, "+"), Sequence("ATCGATCGATAAACGTCGATCGATCGATCG"), 10, "A", 4)

    assert result == []


def test_get_homopolymers_is_case_insensitive():
    # Positive strand, promoter occupies positions 11 - 20
    result = get_homopolymers(MockFeature(21, 30, "+"), Sequence("ATCGATCGATAaaAATCGATCGATCGATCG"), 10, "A", 4)

    assert result == [{
        "repeats": 5,
        "repeat_start_pos": 11,
        "repeat_end_pos": 15,
        "distance_from_cds": 5,
    }]


def test_get_homopolymers_finds_multiple_runs_at_promoter_edges():
    # Positive strand, promoter occupies positions 11 - 20, 2 four-base runs touch each edge of that promoter
    sequence = Sequence("ATCGATCGATAAAATAAAAACGATCGATCG")

    result = get_homopolymers(MockFeature(21, 30, "+"), sequence, 10, "A", 4)

    assert result == [
        {
            "repeats": 4,
            "repeat_start_pos": 11,
            "repeat_end_pos": 14,
            "distance_from_cds": 6,
        },
        {
            "repeats": 5,
            "repeat_start_pos": 16,
            "repeat_end_pos": 20,
            "distance_from_cds": 0,
        }
    ]
