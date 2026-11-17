"""Commit an MCP import without replacing a readable manifest or prior study.

Tables and receipts are immutable, content-addressed files. The manifest is the
only commit point. An interruption before it changes may leave an unreferenced
file, but the old manifest remains readable. Existing studies require explicit
future update/version policy; this writer never silently replaces them.
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import re
from pathlib import Path
import tempfile


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def encoded(value: dict) -> bytes:
    return (json.dumps(value, indent=1, sort_keys=True, allow_nan=False) + '\n').encode()


def atomic_manifest(path: Path, data: bytes) -> None:
    fd, name = tempfile.mkstemp(prefix='.manifest-', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


def immutable(path: Path, data: bytes) -> None:
    """Publish a whole file with no overwrite, even if another writer races."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(prefix='.import-', dir=path.parent)
    try:
        with os.fdopen(fd, 'wb') as f:
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        try:
            os.link(name, path)
        except FileExistsError:
            if path.read_bytes() != data:
                raise ValueError(f'conflicting immutable import file: {path.name}')
    finally:
        Path(name).unlink(missing_ok=True)


def commit(directory: Path, release: str, key: str, table, metadata: dict, recipe: dict) -> dict:
    """Return file/receipt names, or raise ValueError without replacing a study.

    A matching recipe can be repeated without writing. Repeated imports verify
    both receipt and table bytes. The lock serialises writers using this API;
    legacy writers that do not acquire it must not run concurrently.
    """
    import io
    if re.fullmatch(r"[a-z0-9_]+", key) is None:
        raise ValueError('study key must contain only lowercase letters, digits and underscores')
    output = io.BytesIO()
    table.to_parquet(output, index=False)
    data = output.getvalue()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / '.import.lock').open('a+b') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise ValueError('another import is committing; retry after it completes') from None
        manifest = directory / 'manifest.json'
        old = json.loads(manifest.read_text()) if manifest.exists() else {'release': release, 'files': {}, 'studies': {}}
        if not isinstance(old, dict) or not isinstance(old.get('files'), dict) or not isinstance(old.get('studies'), dict):
            raise ValueError('existing manifest has no valid files/studies tables; inspect it before importing')
        if old.get('release') != release:
            raise ValueError('import directory belongs to another release')
        identity = digest(encoded({"release": release, "key": key, "recipe": recipe}))
        previous = old['studies'].get(key)
        if previous is not None:
            if previous.get('import_identity') != identity:
                raise ValueError('this study key already exists with another or unrecorded recipe; no files replaced')
            receipt_name = previous['receipt_file']
            receipt_bytes = (directory / receipt_name).read_bytes()
            if digest(receipt_bytes) != previous['receipt_sha256']:
                raise ValueError('existing import receipt differs from its manifest')
            receipt = json.loads(receipt_bytes)
            if receipt['recipe'] != recipe:
                raise ValueError('existing import recipe differs from the request')
            filename = receipt['output_file']
            if digest((directory / filename).read_bytes()) != old['files'].get(filename):
                raise ValueError('existing import table differs from its manifest')
            if digest(data) != receipt['output_sha256']:
                raise ValueError('the same recipe produced a different table; no files replaced')
            return {'file': str(directory / filename), 'receipt': str(directory / receipt_name), 'already_present': True}
        filename = f'{key}-{digest(data)}.parquet'
        receipt = {'recipe': recipe, 'import_identity': identity, 'output_file': filename, 'output_sha256': digest(data)}
        receipt_bytes = encoded(receipt)
        receipt_name = f'receipts/{identity}.json'
        immutable(directory / filename, data)
        immutable(directory / receipt_name, receipt_bytes)
        new = json.loads(json.dumps(old))
        new['files'][filename] = digest(data)
        new['studies'][key] = {**metadata, 'import_identity': identity, 'receipt_file': receipt_name,
                               'receipt_sha256': digest(receipt_bytes)}
        atomic_manifest(manifest, encoded(new))
        return {'file': str(directory / filename), 'receipt': str(directory / receipt_name), 'already_present': False}
