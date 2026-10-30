from check_hp import get_results_dir


# Tests

def test_get_results_dir_returns_subdirectory_in_current_working_directory(tmp_path, monkeypatch):
    working_dir = tmp_path / "working"
    input_dir = tmp_path / "inputs"
    working_dir.mkdir()
    input_dir.mkdir()
    monkeypatch.chdir(working_dir)

    results_dir = get_results_dir(input_dir / "random_name.gff3")

    assert results_dir.is_dir()
    assert results_dir.resolve() == working_dir / "random_name_results"
