import pandas as pd

from check_hp import get_results_dir, save_results


# Tests

def test_save_results_writes_ordered_headers_and_empty_csvs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    results_dir = get_results_dir("/path/to/sample.gff3")
    save_results([], [], results_dir)

    promoters = pd.read_csv(results_dir / "promoters.csv")
    homopolymers = pd.read_csv(results_dir / "homopolymers.csv")
    
    assert list(promoters.columns) == ["contig", "strand", "cds_id", "cds_start_pos", "cds_end_pos", "gene", "protein", "polyA", "polyT", "polyC", "polyG"]
    assert list(homopolymers.columns) == ["contig", "strand", "cds_id", "cds_start_pos", "cds_end_pos", "gene", "protein", "repeated_base", "repeats", "repeat_start_pos", "repeat_end_pos", "distance_from_cds"]
    assert promoters.empty and homopolymers.empty


def test_save_results_sorts_promoters_naturally(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rows = [
        {"contig": "contig_10", "strand": "+", "cds_id": "gene_10", "cds_start_pos": 1, "cds_end_pos": 2, "gene": "", "protein": "", "polyA": 0, "polyT": 0, "polyC": 0, "polyG": 0},
        {"contig": "contig_2", "strand": "+", "cds_id": "gene_2", "cds_start_pos": 1, "cds_end_pos": 2, "gene": "", "protein": "", "polyA": 0, "polyT": 0, "polyC": 0, "polyG": 0},
        {"contig": "contig_2", "strand": "+", "cds_id": "gene_10", "cds_start_pos": 3, "cds_end_pos": 4, "gene": "", "protein": "", "polyA": 0, "polyT": 0, "polyC": 0, "polyG": 0}
    ]

    results_dir = get_results_dir("sample.gff3")
    save_results(rows, [], results_dir)

    promoters = pd.read_csv(results_dir / "promoters.csv")
    assert list(zip(promoters["contig"], promoters["cds_id"])) == [
        ("contig_2", "gene_2"),
        ("contig_2", "gene_10"),
        ("contig_10", "gene_10"),
    ]


def test_save_results_sorts_homopolymers_naturally(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    details = [
        {"contig": "contig_10", "strand": "+", "cds_id": "gene_10", "cds_start_pos": 1, "cds_end_pos": 2, "gene": "", "protein": "", "repeated_base": "A", "repeats": 8, "repeat_start_pos": 3, "repeat_end_pos": 10, "distance_from_cds": 0},
        {"contig": "contig_2", "strand": "+", "cds_id": "gene_2", "cds_start_pos": 1, "cds_end_pos": 2, "gene": "", "protein": "", "repeated_base": "T", "repeats": 9, "repeat_start_pos": 3, "repeat_end_pos": 11, "distance_from_cds": 0},
        {"contig": "contig_2", "strand": "+", "cds_id": "gene_2", "cds_start_pos": 1, "cds_end_pos": 2, "gene": "", "protein": "", "repeated_base": "A", "repeats": 8, "repeat_start_pos": 11, "repeat_end_pos": 18, "distance_from_cds": 10},
        {"contig": "contig_2", "strand": "+", "cds_id": "gene_2", "cds_start_pos": 1, "cds_end_pos": 2, "gene": "", "protein": "", "repeated_base": "A", "repeats": 8, "repeat_start_pos": 19, "repeat_end_pos": 26, "distance_from_cds": 2},
    ]

    results_dir = get_results_dir("sample.gff3")
    save_results([], details, results_dir)

    homopolymers = pd.read_csv(results_dir / "homopolymers.csv")
    assert list(zip(
        homopolymers["contig"],
        homopolymers["cds_id"],
        homopolymers["repeated_base"],
        homopolymers["distance_from_cds"],
    )) == [
        ("contig_2", "gene_2", "A", 2),
        ("contig_2", "gene_2", "A", 10),
        ("contig_2", "gene_2", "T", 0),
        ("contig_10", "gene_10", "A", 0),
    ]
