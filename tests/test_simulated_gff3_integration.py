from pathlib import Path
import pandas as pd

from check_hp import main


# Tests

def test_check_hp_detects_expected_homopolymers_in_all_fixtures(tmp_path, monkeypatch):
    fixtures_dir = Path(__file__).resolve().parent / "fixtures"
    sort_columns = ["contig", "cds_id", "repeated_base", "distance_from_cds"]

    gff_files = fixtures_dir.glob("SIM*.gff3")
    expected_df = pd.read_csv(fixtures_dir / "expected_homopolymers.csv")

    monkeypatch.chdir(tmp_path)
    for gff_path in gff_files:
        main(str(gff_path), promoter_length=100, repeats=8)

        results_dir = tmp_path / f"{gff_path.stem}_results"
        detected = pd.read_csv(results_dir / "homopolymers.csv")
        expected_fixture_df = expected_df[expected_df["file"] == gff_path.name].drop(columns="file")

        pd.testing.assert_frame_equal(
            detected.sort_values(by=sort_columns, ignore_index=True),
            expected_fixture_df.sort_values(by=sort_columns, ignore_index=True),
            obj=f"Detected homopolymers in {gff_path.name}",
        )
