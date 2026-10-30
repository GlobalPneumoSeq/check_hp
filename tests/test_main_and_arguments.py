import pandas as pd
import pytest

from check_hp import main, parse_arguments


# Tests

def test_main_scans_valid_promoter_and_writes_both_reports(tmp_path, monkeypatch):
    input_dir = tmp_path / "inputs"
    input_dir.mkdir()
    gff_path = input_dir / "testfile_abc.gff3"
    monkeypatch.chdir(tmp_path)
    gff_path.write_text(
        "##gff-version 3\n"
        "contig_1	Bakta	region	1	40	.	+	.	ID=contig_1;Name=contig_1\n"
        "contig_1	Bakta	CDS	21	30	.	+	0	ID=cds1;gene=test;Name=protein\n"
        "##FASTA\n"
        ">contig_1\n"
        "ATCGATCGAT" + "A" * 8 + "C" * 2 + "G" * 10 + "T" * 10 + "\n"
    )

    main(str(gff_path), promoter_length=10, repeats=8)

    results_dir = tmp_path / "testfile_abc_results"
    promoters = pd.read_csv(results_dir / "promoters.csv")
    homopolymers = pd.read_csv(results_dir / "homopolymers.csv")
    assert promoters[["cds_id", "polyA", "polyT", "polyC", "polyG"]].to_dict("records") == [
        {"cds_id": "cds1", "polyA": 1, "polyT": 0, "polyC": 0, "polyG": 0}
    ]
    assert homopolymers[["repeated_base", "repeats", "repeat_start_pos", "repeat_end_pos", "distance_from_cds"]].to_dict("records") == [
        {"repeated_base": "A", "repeats": 8, "repeat_start_pos": 11, "repeat_end_pos": 18, "distance_from_cds": 2}
    ]


def test_main_scans_minus_strand_and_excludes_overlapping_and_out_of_bounds_promoters(tmp_path, monkeypatch):
    gff_path = tmp_path / "test.gff3"
    monkeypatch.chdir(tmp_path)
    gff_path.write_text(
        "##gff-version 3\n"
        "contig_minus	Bakta	region	1	40	.	+	.	ID=contig_minus;Name=contig_minus\n"
        "contig_minus	Bakta	CDS	11	20	.	-	0	ID=minus_cds;gene=test;Name=minus protein\n"
        "contig_overlap	Bakta	region	1	40	.	+	.	ID=contig_overlap;Name=contig_overlap\n"
        "contig_overlap	Bakta	CDS	12	13	.	-	0	ID=blocker;gene=test;Name=blocker\n"
        "contig_overlap	Bakta	CDS	21	30	.	+	0	ID=overlapped_cds;gene=test;Name=overlapped\n"
        "contig_boundary	Bakta	region	1	40	.	+	.	ID=contig_boundary;Name=contig_boundary\n"
        "contig_boundary	Bakta	CDS	5	10	.	+	0	ID=boundary_cds;gene=test;Name=boundary\n"
        "##FASTA\n"
        ">contig_minus\n"
        + "ATCG" * 5 + "C" * 8 + "A" * 2 + "G" * 10 + "\n"
        ">contig_overlap\n"
        + "ATCGATCGAT" + "A" * 10 + "G" * 20 + "\n"
        ">contig_boundary\n"
        + "A" * 40 + "\n"
    )

    main(str(gff_path), promoter_length=10, repeats=8)

    results_dir = tmp_path / "test_results"
    summary = pd.read_csv(results_dir / "promoters.csv")
    details = pd.read_csv(results_dir / "homopolymers.csv")
    assert summary["cds_id"].tolist() == ["minus_cds"]
    assert summary.loc[0, "polyG"] == 1
    assert details[["repeated_base", "repeat_start_pos", "repeat_end_pos", "distance_from_cds"]].to_dict("records") == [
        {"repeated_base": "G", "repeat_start_pos": 21, "repeat_end_pos": 28, "distance_from_cds": 0}
    ]


def test_main_rejects_cds_with_invalid_strand(tmp_path, monkeypatch):
    gff_path = tmp_path / "invalid_strand.gff3"
    monkeypatch.chdir(tmp_path)
    gff_path.write_text(
        "##gff-version 3\n"
        "contig_1	Bakta	region	1	30	.	+	.	ID=contig_1;Name=contig_1\n"
        "contig_1	Bakta	CDS	11	20	.	.	0	ID=cds1\n"
        "##FASTA\n"
        ">contig_1\n"
        + "A" * 30
    )

    with pytest.raises(ValueError, match=r"expected strand value of \+ or -, got \. for feature cds1 on contig contig_1"):
        main(str(gff_path), promoter_length=5, repeats=4)


def test_parse_arguments_uses_defaults(tmp_path, monkeypatch):
    gff_path = tmp_path / "test.gff3"
    gff_path.write_text("##gff-version 3\n")
    monkeypatch.setattr("sys.argv", ["check_hp.py", str(gff_path)])

    args = parse_arguments()

    assert (args.length, args.repeats) == (100, 8)


def test_parse_arguments_accepts_options_and_rejects_invalid_length(tmp_path, monkeypatch):
    gff_path = tmp_path / "test.gff3"
    gff_path.write_text("##gff-version 3\n")
    monkeypatch.setattr("sys.argv", ["check_hp.py", str(gff_path), "--length", "150", "--repeats", "10"])
    args = parse_arguments()
    assert (args.gff, args.length, args.repeats) == (str(gff_path), 150, 10)

    monkeypatch.setattr("sys.argv", ["check_hp.py", str(gff_path), "--length", "1"])
    with pytest.raises(SystemExit):
        parse_arguments()


@pytest.mark.parametrize(
    "argument, value",
    [
        ("gff", "missing.gff3"),
        ("gff", "wrong-header.gff3"),
        ("--length", "abc"),
        ("--length", "1"),
        ("--repeats", "abc"),
        ("--repeats", "1"),
    ],
)
def test_parse_arguments_rejects_invalid_values(tmp_path, monkeypatch, argument, value):
    gff_path = tmp_path / "input.gff3"
    gff_path.write_text("##gff-version 3\n")
    args = ["check_hp.py"]
    if argument == "gff":
        path = tmp_path / value
        if value == "wrong-header.gff3":
            path.write_text("not a GFF header\n")
        args.append(str(path))
    else:
        args.extend([str(gff_path), argument, value])
    monkeypatch.setattr("sys.argv", args)

    with pytest.raises(SystemExit):
        parse_arguments()
