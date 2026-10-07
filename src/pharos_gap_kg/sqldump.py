r"""Read a published MySQL dump as a stream: which tables it holds, and chosen tables as typed parquet.

The Pharos snapshots are `mysqldump` files, about 2 GB compressed. Restoring one needs a MySQL server
and tens of GB of disk; reading it as text needs neither, and leaves a lineage that is one digest long:
every table written here cites the SHA-256 of the file it was read from.

HOW A DUMP IS READ. A table is a `CREATE TABLE` block, whose column lines give names and types, then
`INSERT INTO ... VALUES (..),(..),...;` lines of about 1 MB each. A row is found by a pattern that
knows a quoted string may hold parentheses, commas and escaped quotes; splitting on "),(" would cut a
row in two wherever a description contains that text. The rows of a chunk are then parsed as
single-quoted CSV by Arrow with the column types the table declares. A value that does not fit its
declared type stops the run: the reader is wrong or the dump is, and a column of silently coerced
nulls would hide which.

WHAT IS LOST, deliberately: an escaped newline, tab or NUL inside a string becomes a space.

Run:  python -m pharos_gap_kg.sqldump inventory pharos319
      python -m pharos_gap_kg.sqldump extract pharos319 --tables target,protein,disease
      python -m pharos_gap_kg.sqldump --selftest
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterator

import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq

from pharos_gap_kg import paths, registry

# One row: "(" then any run of plain bytes or whole quoted strings, then ")". Possessive, so a
# malformed row fails at once instead of backtracking.
ROW = re.compile(rb"\(((?:[^()'\\]++|'(?:[^'\\]++|\\.)*+')*+)\)", re.S)
ESCAPE = re.compile(rb"\\(.)", re.S)
_TO_SPACE = {b"n", b"r", b"t", b"0", b"Z"}
COLUMN = re.compile(rb"^\s+`([^`]+)`\s+([a-zA-Z]+)")
CHUNK_BYTES = 96 << 20

INT_TYPES = {"int", "bigint", "smallint", "tinyint", "mediumint", "integer", "year"}
FLOAT_TYPES = {"decimal", "double", "float", "numeric", "real"}
BINARY_TYPES = {"blob", "tinyblob", "mediumblob", "longblob", "binary", "varbinary", "bit"}


def _unescape(m: re.Match) -> bytes:
    return b" " if m.group(1) in _TO_SPACE else m.group(0)


def arrow_type(sql_type: str) -> pa.DataType:
    t = sql_type.lower()
    if t in INT_TYPES:
        return pa.int64()
    if t in FLOAT_TYPES:
        return pa.float64()
    return pa.string()


def lines(src: Path | io.BufferedIOBase) -> Iterator[bytes]:
    """Lines of a dump: a .gz path is decompressed by the system gzip, anything else is read as given."""
    if isinstance(src, Path) and src.suffix == ".gz":
        proc = subprocess.Popen(["gzip", "-dc", str(src)], stdout=subprocess.PIPE, bufsize=8 << 20)
        try:
            yield from proc.stdout
        finally:
            proc.stdout.close()
            rc = proc.wait()
        if rc not in (0, -13):
            raise SystemExit(f"gzip -dc {src.name} exited {rc}: the file is truncated or not gzip")
    elif isinstance(src, Path):
        with src.open("rb") as f:
            yield from f
    else:
        yield from src


def events(src) -> Iterator[tuple[str, str, object]]:
    """('create', table, [(column, sql_type), ...]) then ('insert', table, payload_bytes) per INSERT line."""
    name, cols = None, None
    for line in lines(src):
        if cols is not None:
            m = COLUMN.match(line)
            if m:
                cols.append((m.group(1).decode(), m.group(2).decode()))
            elif line.startswith(b")"):
                yield "create", name, cols
                name, cols = None, None
            continue
        if line.startswith(b"CREATE TABLE `"):
            name, cols = line.split(b"`", 2)[1].decode(), []
        elif line.startswith(b"INSERT INTO `"):
            table = line.split(b"`", 2)[1].decode()
            at = line.find(b" VALUES ")
            if at < 0:
                raise SystemExit(f"an INSERT into {table} has no VALUES clause: not a mysqldump this reader knows")
            yield "insert", table, line[at + 8:]


def inventory(src) -> dict:
    """Every table: its columns and declared types, its rows, and the size and SHA-256 of its INSERT payload.

    The digest is over the VALUES text in dump order, so two dumps give the same digest for a table
    only if it holds the same rows, with the same ids, written in the same order. Equal digests mean
    an identical table. Unequal digests mean only "compare the rows"; a re-ordered table differs here.
    """
    out: dict[str, dict] = {}
    hashes: dict[str, "hashlib._Hash"] = {}
    for kind, table, x in events(src):
        if kind == "create":
            out[table] = {"columns": [{"name": c, "type": t} for c, t in x], "rows": 0, "payload_bytes": 0}
        else:
            d = out.setdefault(table, {"columns": [], "rows": 0, "payload_bytes": 0})
            d["rows"] += sum(1 for _ in ROW.finditer(x))
            d["payload_bytes"] += len(x)
            hashes.setdefault(table, hashlib.sha256()).update(x)
    for table, d in out.items():
        d["payload_sha256"] = hashes[table].hexdigest() if table in hashes else None
    return out


def _parse(rows_csv: bytes, cols: list[tuple[str, str]], table: str) -> pa.Table:
    names = [c for c, _ in cols]
    try:
        return pacsv.read_csv(
            pa.BufferReader(rows_csv),
            read_options=pacsv.ReadOptions(column_names=names, block_size=32 << 20),
            parse_options=pacsv.ParseOptions(delimiter=",", quote_char="'", double_quote=False, escape_char="\\"),
            convert_options=pacsv.ConvertOptions(column_types={c: arrow_type(t) for c, t in cols},
                                                 null_values=["NULL"], strings_can_be_null=True,
                                                 quoted_strings_can_be_null=False))
    except pa.ArrowInvalid as e:
        raise SystemExit(f"REFUSED table {table}: a value does not fit its declared type or a row has the "
                         f"wrong number of fields ({e}). Nothing was coerced; fix the reader.") from e


def _rows(payload: bytes) -> tuple[bytes, int]:
    """The rows of one INSERT payload as newline-separated CSV, and how many there are."""
    found = ROW.findall(payload)
    body = b"\n".join(found)
    if b"\\" in body:
        body = ESCAPE.sub(_unescape, body)
    return body.decode("utf-8", "replace").encode("utf-8"), len(found)


def extract(src, wanted: set[str], out_dir: Path) -> dict:
    """Write each wanted table as `<out_dir>/<table>.parquet`; returns {table: {rows, columns}}."""
    out_dir.mkdir(parents=True, exist_ok=True)
    schema: dict[str, list[tuple[str, str]]] = {}
    done: dict[str, dict] = {}
    cur, buf, size, n, writer = None, [], 0, 0, None

    def flush(final: bool) -> None:
        nonlocal buf, size, writer, n, cur
        if cur is None:
            return
        if buf:
            t = _parse(b"\n".join(buf), schema[cur], cur)
            if writer is None:
                writer = pq.ParquetWriter(out_dir / f"{cur}.parquet", t.schema, compression="zstd")
            writer.write_table(t)
            buf, size = [], 0
        if final:
            if writer is not None:
                writer.close()
                got = pq.ParquetFile(out_dir / f"{cur}.parquet").metadata.num_rows
                if got != n:
                    raise SystemExit(f"REFUSED table {cur}: {n:,} rows were read and {got:,} were written")
            done[cur] = {"rows": n, "columns": [{"name": c, "type": t} for c, t in schema[cur]]}
            cur, writer, n = None, None, 0

    for kind, table, x in events(src):
        if kind == "create":
            schema[table] = x
            continue
        if table not in wanted:
            if cur is not None:
                flush(final=True)
            continue
        if table != cur:
            flush(final=True)
            bad = [c for c, t in schema[table] if t.lower() in BINARY_TYPES]
            if bad:
                raise SystemExit(f"REFUSED table {table}: binary columns {bad} are not read by this reader")
            cur = table
        body, k = _rows(x)
        buf.append(body)
        size += len(body)
        n += k
        if size >= CHUNK_BYTES:
            flush(final=False)
    flush(final=True)
    for t in sorted(wanted - set(done)):
        if t in schema:   # declared and empty
            pq.write_table(pa.table({c: pa.array([], arrow_type(ty)) for c, ty in schema[t]}), out_dir / f"{t}.parquet")
            done[t] = {"rows": 0, "columns": [{"name": c, "type": ty} for c, ty in schema[t]]}
    missing = sorted(wanted - set(done))
    if missing:
        raise SystemExit(f"REFUSED: the dump has no table named {missing}")
    return done


def _source(key: str) -> tuple[Path, dict]:
    p = registry.require(key)
    side = p.with_name(p.name + ".lineage.json")
    if not side.is_file():
        raise SystemExit(f"{p.name} has no lineage record beside it; fetch it with:  {registry.SOURCES[key].fetch_command}")
    return p, json.loads(side.read_text())


def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    dump = (
        b"-- MySQL dump\nDROP TABLE IF EXISTS `t`;\nCREATE TABLE `t` (\n"
        b"  `id` int(11) NOT NULL AUTO_INCREMENT,\n  `name` varchar(255) DEFAULT NULL,\n"
        b"  `score` decimal(10,4) DEFAULT NULL,\n  `note` text,\n  PRIMARY KEY (`id`),\n  KEY `k` (`name`)\n"
        b") ENGINE=InnoDB DEFAULT CHARSET=utf8;\nLOCK TABLES `t` WRITE;\n"
        b"INSERT INTO `t` VALUES (1,'plain',0.5000,'a, b'),(2,'paren (x),(y)',NULL,'it\\'s'),(3,NULL,-1.25,'');\n"
        b"INSERT INTO `t` VALUES (4,'back\\\\slash',1e-05,'two\\nlines'),(5,'NULL',2,'caf\xc3\xa9 \"q\"');\n"
        b"UNLOCK TABLES;\nCREATE TABLE `u` (\n  `a` int(11) NOT NULL\n) ENGINE=InnoDB;\n"
        b"INSERT INTO `u` VALUES (7),(8);\nCREATE TABLE `empty` (\n  `a` int(11) NOT NULL\n) ENGINE=InnoDB;\n")
    inv = inventory(io.BytesIO(dump))
    ck("inventory finds three tables", set(inv) == {"t", "u", "empty"})
    ck("columns are the column lines, not the keys", [c["name"] for c in inv["t"]["columns"]] == ["id", "name", "score", "note"])
    ck("a row containing '),(' in a string is one row (5 rows, not 6)", inv["t"]["rows"] == 5)
    ck("a declared table with no INSERT has zero rows", inv["empty"]["rows"] == 0)
    same = inventory(io.BytesIO(dump))
    moved = inventory(io.BytesIO(dump.replace(b"(7),(8)", b"(7),(9)")))
    ck("the same table read twice has one digest; an empty table has none",
       inv["t"]["payload_sha256"] == same["t"]["payload_sha256"] and inv["empty"]["payload_sha256"] is None)
    ck("one changed value changes that table's digest and no other",
       moved["u"]["payload_sha256"] != inv["u"]["payload_sha256"] and moved["t"]["payload_sha256"] == inv["t"]["payload_sha256"])
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        done = extract(io.BytesIO(dump), {"t", "empty"}, Path(d))
        t = pq.read_table(Path(d) / "t.parquet").to_pylist()
        ck("five rows written, and the count is checked against the file", done["t"]["rows"] == 5 and len(t) == 5)
        ck("declared types become column types", pq.read_schema(Path(d) / "t.parquet").types == [pa.int64(), pa.string(), pa.float64(), pa.string()])
        ck("commas and parentheses inside a string survive", t[0]["note"] == "a, b" and t[1]["name"] == "paren (x),(y)")
        ck("an escaped quote is a quote", t[1]["note"] == "it's")
        ck("an unquoted NULL is null; a quoted 'NULL' is the text", t[2]["name"] is None and t[1]["score"] is None and t[4]["name"] == "NULL")
        ck("an empty string is empty, not null", t[2]["note"] == "")
        ck("an escaped backslash is one backslash", t[3]["name"] == "back\\slash")
        ck("an escaped newline becomes a space (declared loss)", t[3]["note"] == "two lines")
        ck("exponent floats and UTF-8 survive", abs(t[3]["score"] - 1e-05) < 1e-12 and t[4]["note"] == 'café "q"')
        ck("an unwanted table is not written", not (Path(d) / "u.parquet").exists())
        ck("a wanted empty table is written with its schema", pq.read_table(Path(d) / "empty.parquet").num_rows == 0)
        bad = dump.replace(b"(1,'plain',0.5000,'a, b')", b"('x','plain',0.5000,'a, b')")
        try:
            extract(io.BytesIO(bad), {"t"}, Path(d) / "bad")
            ck("a value that does not fit its declared type is refused", False)
        except SystemExit:
            ck("a value that does not fit its declared type is refused", True)
        try:
            extract(io.BytesIO(dump), {"nope"}, Path(d) / "none")
            ck("a table the dump does not hold is refused", False)
        except SystemExit:
            ck("a table the dump does not hold is refused", True)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("inventory")
    i.add_argument("key", choices=sorted(registry.SOURCES))
    e = sub.add_parser("extract")
    e.add_argument("key", choices=sorted(registry.SOURCES))
    e.add_argument("--tables", required=True, help="comma-separated table names")
    a = ap.parse_args(argv)
    src, lin = _source(a.key)
    t0 = time.time()
    if a.cmd == "inventory":
        inv = inventory(src)
        out = paths.DOCS_DATA / f"{a.key}_tables.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps({"source": a.key, "sha256": lin["sha256"], "published": lin["published"],
                                   "tables": inv}, indent=1) + "\n")
        for name, d in sorted(inv.items(), key=lambda kv: -kv[1]["payload_bytes"]):
            print(f"  {name:38s} {d['rows']:>12,} rows  {d['payload_bytes'] / 1e6:>9.1f} MB  {len(d['columns']):>3} cols")
        print(f"{len(inv)} tables, {sum(d['rows'] for d in inv.values()):,} rows, {time.time() - t0:.0f} s -> {out}")
    else:
        wanted = {t.strip() for t in a.tables.split(",") if t.strip()}
        out_dir = paths.tables(a.key)
        done = extract(src, wanted, out_dir)
        man_p = out_dir / "_manifest.json"
        man = json.loads(man_p.read_text()) if man_p.is_file() else {"source": a.key, "sha256": lin["sha256"], "tables": {}}
        if man["sha256"] != lin["sha256"]:
            raise SystemExit(f"REFUSED: {out_dir} holds tables from a different file ({man['sha256'][:12]}); remove it first")
        man["tables"].update(done)
        man_p.write_text(json.dumps(man, indent=1) + "\n")
        for name, d in sorted(done.items()):
            print(f"  {name:38s} {d['rows']:>12,} rows")
        print(f"{len(done)} tables in {time.time() - t0:.0f} s -> {out_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
