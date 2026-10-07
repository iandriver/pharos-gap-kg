"""The registry is the list of inputs; these checks keep it usable as one."""
import pytest

from pharos_gap_kg import registry


def test_keys_are_unique_and_urls_are_https():
    keys = [s.key for s in registry._S]
    assert len(keys) == len(set(keys))
    assert all(s.url.startswith("https://") for s in registry._S)


def test_every_source_names_its_fetch_command_and_size():
    for s in registry._S:
        assert s.bytes > 0 and s.key in s.fetch_command and len(s.published) == 10


def test_a_missing_input_names_how_to_fetch_it(tmp_path, monkeypatch):
    monkeypatch.setattr(registry.paths, "RAW", tmp_path)
    with pytest.raises(SystemExit) as e:
        registry.require("pharos317")
    assert "python -m pharos_gap_kg.fetch pharos317" in str(e.value)


def test_a_truncated_input_is_not_accepted(tmp_path, monkeypatch):
    monkeypatch.setattr(registry.paths, "RAW", tmp_path)
    (tmp_path / "pharos317.sql.gz").write_bytes(b"short")
    with pytest.raises(SystemExit) as e:
        registry.require("pharos317")
    assert "interrupted or replaced" in str(e.value)


def test_every_summary_the_documents_quote_is_tracked():
    """`docs/data/` holds the figures the findings cite and the frozen split. A broad ignore pattern
    once kept every file there out of the history while the commit messages said they were in it."""
    import subprocess
    from pharos_gap_kg import paths
    if paths.ROOT is None or not (paths.ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    ignored = subprocess.run(["git", "ls-files", "--others", "--ignored", "--exclude-standard", "docs"],
                             cwd=paths.ROOT, capture_output=True, text=True, check=True).stdout.split()
    assert not ignored, f"ignored by .gitignore and so never committed: {ignored}"
