"""Import interruption and identity checks preserve the previous readable state."""
import json
from pathlib import Path
import pandas as pd
import pytest
from pharos_gap_kg import importstore


def table():
    return pd.DataFrame({'target': ['UniProtKB:P1'], 'trait': ['MONDO:1'], 'study': ['study']})


def test_first_import_and_repeat_are_verified_and_idempotent(tmp_path):
    r = importstore.commit(tmp_path, 'test', 'study', table(), {'doi': '10.1234/x'}, {'source': 'sha'})
    before = (tmp_path / 'manifest.json').read_bytes()
    again = importstore.commit(tmp_path, 'test', 'study', table(), {'doi': '10.1234/x'}, {'source': 'sha'})
    assert not r['already_present'] and again['already_present']
    assert (tmp_path / 'manifest.json').read_bytes() == before
    receipt = json.loads(Path(r['receipt']).read_text())
    assert receipt['recipe'] == {'source': 'sha'}
    assert importstore.digest(Path(r['file']).read_bytes()) == receipt['output_sha256']


def test_changed_or_legacy_study_key_is_refused_without_overwriting(tmp_path):
    importstore.commit(tmp_path, 'test', 'study', table(), {}, {'source': 'one'})
    old = (tmp_path / 'manifest.json').read_bytes()
    with pytest.raises(ValueError, match='already exists'):
        importstore.commit(tmp_path, 'test', 'study', table(), {}, {'source': 'two'})
    assert (tmp_path / 'manifest.json').read_bytes() == old
    manifest = json.loads(old)
    manifest['studies']['legacy'] = {'doi': '10.1234/legacy'}
    (tmp_path / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='unrecorded recipe'):
        importstore.commit(tmp_path, 'test', 'legacy', table(), {}, {})


def test_interrupted_commit_leaves_old_manifest_and_table_intact(tmp_path, monkeypatch):
    first = importstore.commit(tmp_path, 'test', 'first', table(), {}, {})
    before = (tmp_path / 'manifest.json').read_bytes()
    saved = Path(first['file']).read_bytes()
    def interrupted(*args):
        raise OSError('simulated interrupted manifest commit')
    monkeypatch.setattr(importstore, 'atomic_manifest', interrupted)
    with pytest.raises(OSError):
        importstore.commit(tmp_path, 'test', 'next', table(), {}, {'source': 'next'})
    assert (tmp_path / 'manifest.json').read_bytes() == before
    assert Path(first['file']).read_bytes() == saved


def test_tampered_table_or_receipt_is_not_an_idempotent_success(tmp_path):
    result = importstore.commit(tmp_path, 'test', 'study', table(), {}, {})
    Path(result['file']).write_bytes(b'changed')
    with pytest.raises(ValueError, match='table differs'):
        importstore.commit(tmp_path, 'test', 'study', table(), {}, {})
    Path(result['receipt']).write_bytes(b'changed')
    with pytest.raises(ValueError, match='receipt differs'):
        importstore.commit(tmp_path, 'test', 'study', table(), {}, {})


def test_other_release_is_refused(tmp_path):
    importstore.commit(tmp_path, 'test', 'study', table(), {}, {})
    with pytest.raises(ValueError, match='another release'):
        importstore.commit(tmp_path, 'other', 'new', table(), {}, {})


def test_same_recipe_cannot_hide_changed_output(tmp_path):
    importstore.commit(tmp_path, 'test', 'study', table(), {}, {})
    before = (tmp_path / 'manifest.json').read_bytes()
    changed = table()
    changed.loc[0, 'trait'] = 'MONDO:2'
    with pytest.raises(ValueError, match='different table'):
        importstore.commit(tmp_path, 'test', 'study', changed, {}, {})
    assert (tmp_path / 'manifest.json').read_bytes() == before


def test_concurrent_import_lock_refuses_a_second_writer(tmp_path):
    import fcntl
    with (tmp_path / '.import.lock').open('a+b') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(ValueError, match='another import'):
            importstore.commit(tmp_path, 'test', 'study', table(), {}, {})
    assert not (tmp_path / 'manifest.json').exists()


def test_tool_join_writes_a_bound_receipt_and_repeats_safely(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from pharos_gap_kg import add
    bundle = tmp_path / 'bundle' / 'test'
    bundle.mkdir(parents=True)
    (bundle / 'manifest.json').write_text('{"release":"test"}')
    source = tmp_path / 'source.tsv'
    source.write_text('gene\nAPOE\n')
    S = SimpleNamespace(release='test', G=SimpleNamespace(dir=bundle), lists='test lists',
                        refs=SimpleNamespace(**{k: {} for k in ('uniprot', 'gene', 'target', 'vocab', 'xref')}),
                        sym=SimpleNamespace(map={}, ids={}))
    monkeypatch.setattr(add, 'traits', lambda *args: [])
    monkeypatch.setattr(add, 'check2', lambda *args: {'accepted': True})
    joined = pd.DataFrame({'target': ['UniProtKB:P1'], 'trait': ['MONDO:1']})
    monkeypatch.setattr(add, 'joinable', lambda *args: joined)
    about = add.STUDIES['A']['about']
    first = add.tool_join(S, str(source), about, source.name, 1, {'map': True})
    assert first['joined'] and not first['already_present']
    receipt = json.loads(Path(first['receipt']).read_text())
    assert receipt['recipe']['source_sha256'] == importstore.digest(source.read_bytes())
    assert receipt['recipe']['mapping'] == {'map': True}
    assert receipt['recipe']['bundle_manifest_sha256']
    assert receipt['recipe']['resolver_sha256']
    again = add.tool_join(S, str(source), about, source.name, 1, {'map': True})
    assert again['joined'] and again['already_present']
    source.write_text('gene\nOTHER\n')
    refused = add.tool_join(S, str(source), about, source.name, 1, {'map': True})
    assert not refused['joined'] and 'already exists' in refused['problems'][0]


def test_replay_uses_saved_mapping_without_writing_or_calling_a_model(tmp_path, monkeypatch):
    from types import SimpleNamespace
    from pharos_gap_kg import add
    bundle = tmp_path / 'bundle'
    bundle.mkdir()
    (bundle / 'manifest.json').write_text('{}')
    source = tmp_path / 'source.tsv'
    source.write_text('gene\nAPOE\n')
    recipe = {'source_sha256': importstore.digest(source.read_bytes()),
              'bundle_manifest_sha256': importstore.digest((bundle / 'manifest.json').read_bytes()),
              'about': {}, 'mapping': {'map': True}, 'sheet': source.name, 'header_row': 1,
              'code_sha256': {}}
    receipt = tmp_path / 'receipt.json'
    receipt.write_text(json.dumps({'recipe': recipe, 'output_sha256': 'old'}))
    monkeypatch.setattr(add, 'tool_check', lambda *args: {'accepted': True, 'problems': []})
    S = SimpleNamespace(G=SimpleNamespace(dir=bundle))
    before = {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()}
    result = add.replay_import(S, str(receipt), str(source))
    assert result['accepted'] and not result['replay']['writes']
    assert {p: p.read_bytes() for p in tmp_path.rglob('*') if p.is_file()} == before
    source.write_text('changed')
    refused = add.replay_import(S, str(receipt), str(source))
    assert not refused['accepted'] and 'source' in refused['problems'][0]
