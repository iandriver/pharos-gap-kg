"""The read-only pages the server shows in a browser: the questions (/), and how tables from papers were added (/add).

The pages hold no data of their own and call no model. / asks, through /api/<question>, the same functions the MCP
tools answer from, with the same arguments, so a browser and an MCP client get the same answer. /add reads
walkthrough_<release>.json, which `walkthrough` builds from the committed records of "Adding a table"
(docs/data/add_*_<release>.json) and from the sheets each run read, by name and size only: no cell of a study's file
is copied into it.

Run:  python -m pharos_gap_kg.web walkthrough pharos319     # write src/pharos_gap_kg/web/walkthrough_pharos319.json
      python -m pharos_gap_kg.web --selftest
"""
from __future__ import annotations

import asyncio
import inspect
import json
import math
import re
import sys
from pathlib import Path
from urllib.parse import parse_qs

HERE = Path(__file__).parent / "web"
STATIC = {"/": ("index.html", "text/html; charset=utf-8"), "/add": ("add.html", "text/html; charset=utf-8"),
          "/style.css": ("style.css", "text/css; charset=utf-8"), "/app.js": ("app.js", "text/javascript; charset=utf-8"),
          "/add.js": ("add.js", "text/javascript; charset=utf-8")}
CSP = ("default-src 'none'; script-src 'self'; style-src 'self'; connect-src 'self'; img-src 'self' data:; base-uri 'none'; "
       "form-action 'self'; frame-ancestors 'none'")
MAX_QUERY = 2048            # a question's arguments are a few dozen characters


# ---- /api/<question>: the MCP tools' own functions, read from their signatures ---------------------

def params(fn) -> dict[str, tuple[str, bool]]:
    """Each argument of a tool function: (its type, whether it must be given). The server's annotations are strings."""
    out = {}
    for name, p in inspect.signature(fn).parameters.items():
        ann = str(p.annotation)
        kind = "bool" if "bool" in ann else "int" if "int" in ann else "str"
        out[name] = (kind, p.default is inspect.Parameter.empty)
    return out


def _clean(x):
    """JSON a browser can parse: NaN and infinities become null."""
    if isinstance(x, float) and not math.isfinite(x):
        return None
    if isinstance(x, dict):
        return {str(k): _clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_clean(v) for v in x]
    return x


def arguments(fn, query: dict[str, list[str]]) -> tuple[dict | None, str | None]:
    """The keyword arguments a query string gives a tool function, or the reason it gives none."""
    spec = params(fn)
    extra = sorted(set(query) - set(spec))
    if extra:
        return None, f"unknown argument(s) {extra}; this question takes {sorted(spec)}"
    kw = {}
    for name, (kind, required) in spec.items():
        if name not in query:
            if required:
                return None, f"the argument {name!r} is required"
            continue
        v = query[name][-1].strip()
        if kind == "int":
            if v == "":
                continue
            try:
                kw[name] = int(v)
            except ValueError:
                return None, f"{name} must be a whole number"
        elif kind == "bool":
            kw[name] = v.casefold() in ("1", "true", "yes", "on")
        else:
            kw[name] = v or None
            if kw[name] is None and required:
                return None, f"the argument {name!r} is required"
    return kw, None


def answer(answers: dict, question: str, query: dict[str, list[str]]) -> tuple[int, dict]:
    """(HTTP status, body) for one question asked by URL."""
    fn = answers.get(question)
    if fn is None:
        return 404, {"error": f"no question {question!r}; the questions are {sorted(answers)}"}
    kw, why = arguments(fn, query)
    if why:
        return 400, {"error": why}
    try:
        return 200, _clean(fn(**kw))
    except (ValueError, TypeError, KeyError) as e:
        return 400, {"error": f"the question could not be answered with these arguments: {type(e).__name__}"}


# ---- the router in front of the MCP app ------------------------------------------------------------

def _walkthrough_file(release: str) -> Path:
    return HERE / f"walkthrough_{release}.json"


class Pages:
    """GET and HEAD for the pages, their two scripts and stylesheet, /walkthrough.json and /api/<question>; every other
    request goes to the MCP app. The Guard in front of it counts these requests as it counts the others."""

    def __init__(self, app, answers: dict, release: str):
        self.app, self.answers = app, answers
        self.files = {p: ((HERE / f).read_bytes(), t) for p, (f, t) in STATIC.items()}
        w = _walkthrough_file(release)
        body = w.read_bytes() if w.is_file() else json.dumps({"release": release, "studies": [],
                                                              "note": f"no walkthrough is built for {release}"}).encode()
        self.files["/walkthrough.json"] = (body, "application/json; charset=utf-8")

    async def _send(self, send, status: int, body: bytes, ctype: str, head: bool, cache: str) -> None:
        h = [(b"content-type", ctype.encode()), (b"x-content-type-options", b"nosniff"), (b"referrer-policy", b"no-referrer"),
             (b"cache-control", cache.encode()), (b"content-length", str(len(body)).encode())]
        if ctype.startswith("text/html"):
            h.append((b"content-security-policy", CSP.encode()))
        await send({"type": "http.response.start", "status": status, "headers": h})
        await send({"type": "http.response.body", "body": b"" if head else body})

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http" or scope["method"] not in ("GET", "HEAD"):
            return await self.app(scope, receive, send)
        path, head = scope["path"], scope["method"] == "HEAD"
        if path in self.files:
            body, ctype = self.files[path]
            return await self._send(send, 200, body, ctype, head, "no-cache")
        if path.startswith("/api/"):
            raw = scope.get("query_string", b"")
            if len(raw) > MAX_QUERY:
                status, out = 414, {"error": "the question's arguments are too long"}
            else:
                q = parse_qs(raw.decode("utf-8", "replace"), keep_blank_values=True)
                status, out = await asyncio.to_thread(answer, self.answers, path[len("/api/"):], q)
            return await self._send(send, status, json.dumps(out, default=str).encode(), "application/json; charset=utf-8", head, "no-store")
        return await self.app(scope, receive, send)


# ---- the walkthrough of "Adding a table", from the committed records --------------------------------

# Which runs the page shows, in order, and the sentence that says what each one is. Every number on the page comes from the records.
RUNS = [("A", False, "The Alzheimer disease study, whose table had been built by hand before: do the stages reproduce it?"),
        ("B", False, "A Parkinson disease study no part of the entry had read: the registered run."),
        ("B", True, "The Parkinson disease study again, after its result, with a way back from the mapping to the sheet: development, not a test."),
        ("C", False, "A type 1 diabetes study no part of the entry had read."),
        ("C", True, "The type 1 diabetes study's logged replies, checked again under the new check of the rows: development, not a test."),
        ("D", False, "A refractive error study no part of the entry had read, chosen by a fixed rule from the GWAS Catalog, under the new check.")]
MAPPING_KEYS = ("map", "reason", "subject", "object", "predicate", "association_category", "knowledge_level", "agent_type", "qualifiers",
                "edge_properties", "keep_if", "split_by", "by_value", "notes")


def _short(model: str) -> str:
    return {"anthropic/claude-haiku-5.5": "Claude Haiku 5.5", "anthropic/claude-sonnet-5.5": "Claude Sonnet 5.5", "hand": "by hand"}.get(model, model)


def _tries(tries: list[dict], stage: str) -> list[dict]:
    out = []
    for t in tries:
        p = t.get("proposal")
        if stage == "mapping" and isinstance(p, dict):
            p = {k: p[k] for k in MAPPING_KEYS if k in p}
        out.append({"model": _short(t["model"]), "attempt": t.get("attempt"), "proposal": p, "accepted": t["accepted"], "problems": t["problems"]})
    return out


def _sheets(release: str, study: str) -> list[dict]:
    from pharos_gap_kg import add
    f = add.run_dir(release, study) / "card1.json"
    if not f.is_file():
        return []
    return [{"name": s["name"], "rows": s["rows"], "columns": s["columns"]} for s in json.loads(f.read_text())["sheets"]]


def _outcome(study: str, rerun: bool, r: dict, k: dict) -> str:
    j = r.get("join") or {}
    if study == "A":
        a = k["A"]
        return (f"Reproduced: the model's mapping makes the same {a['model_pairs']} pairs of a graph gene and Alzheimer disease as the table built "
                f"by hand. It was compared, not joined, since the table built by hand is already beside the graph.")
    if study == "B" and not rerun:
        return ("Not added. The model chose a sheet that names no gene; given it, the mapping stage said the table should not become edges, "
                "and the stages as registered had no way back to the sheet.")
    if study == "C" and not rerun:
        m = re.search(r"(\d+) of (\d+) rows \(([\d.]+%)\)", r["stage2"]["tries"][-2]["problems"][0])
        bad = f"{m.group(1)} of {m.group(2)} rows, {m.group(3)}, resolved, under the 95% the check then required"
        return (f"Nothing joins. The model chose the loci table, and code refused the mapping that read its genes by symbol ({bad}). "
                f"The next mapping made a local node of each name, which the check of the rows counted as resolved; it links no gene of the graph, "
                f"and the join refuses it. That check was changed after this result.")
    k = k.get(f"{study}_rerun" if rerun else study, {})            # this run's checks
    first = (r.get("verdict", {}).get("rules") or [{}])[0].get("first_figure")
    share = f" ({first['ok']} of {first['rows']} rows give an edge of the graph)" if first else ""
    reading = k.get("reading") or ""
    cat = ("The comparison with the GWAS Catalog cannot be made." if reading.startswith("not computable")
           else f"Against the genes the GWAS Catalog maps to the study's variants by position, the reading is {reading}.")
    band = ""
    if k.get("rows_band"):
        lo, hi = k["rows_band"]
        band = (f" Its {k['rows']:,} rows lie within {lo:,} to {hi:,}, half to twice the {k['abstract_loci']:,} the abstract reports, so by the stated check it is a "
                f"table of the study's loci." if k.get("a_table_of_the_studys_loci") else
                f" Its {k['rows']:,} rows lie outside {lo:,} to {hi:,}, half to twice the {k['abstract_loci']:,} the abstract reports, so by the stated check it is "
                f"not a table of the study's loci.")
    grain = ""
    if reading == "disagree" and k.get("table_genes") and k.get("of_them_graph_genes") and k["table_genes"] < k["of_them_graph_genes"]:
        grain = (f" The table gives one gene per locus and the Catalog maps each of its {k['catalog_associations']} associations, so the table could "
                 f"hold at most {k['table_genes'] / k['of_them_graph_genes']:.0%} of the Catalog's {k['of_them_graph_genes']} graph genes; it holds "
                 f"{k['share_of_catalog_genes_the_table_holds'][0]:.0%}.")
    return (f"It would join {j.get('rows', 0)} links from {j.get('targets', 0)} graph genes{share}.{band} {cat}{grain} "
            f"Not joined: it joins only on the person's word.")


def walkthrough(release: str) -> dict:
    """The page's data: each run's description, the sheets it read, every proposal with code's verdict, and what would join."""
    from pharos_gap_kg import add, paths
    k = json.loads((paths.DOCS_DATA / f"add_compare_{release}.json").read_text())
    studies = []
    for study, rerun, what in RUNS:
        p = add.record_path(release, study, rerun)
        if not p.is_file():
            continue
        r = json.loads(p.read_text())
        s = add.STUDIES[study]
        name = f"{study}_rerun" if rerun else study
        back = r.get("back_to_the_sheet")
        checks = k.get(name, {})
        out = {"id": name, "study": study, "rerun": rerun, "what": what, "citation": s["short"], "disease": s["about"]["disease"],
               "licence": s["about"]["licence"], "description": r["description"], "file": r["file"], "sheets": _sheets(release, study),
               "check_8": add.CHECK8[r.get("check_8", "resolved")],
               "stage1": _tries(r["stage1"]["tries"], "sheet"), "stage2": _tries((r.get("stage2") or {}).get("tries", []), "mapping"),
               "back_to_the_sheet": None if not back else {"reason": back.get("reason"), "stage1": _tries(back["stage1"]["tries"], "sheet"),
                                                          "stage2": _tries((back.get("stage2") or {}).get("tries", []), "mapping"),
                                                          "stopped": back.get("stopped")},
               "sheet": r.get("sheet"), "header_row": r.get("header_row"), "rows": r.get("rows"),
               "traits_on_the_card": r.get("traits_on_the_card"), "knowledge_source": r.get("knowledge_source"),
               "mapping": {kk: r["mapping"][kk] for kk in MAPPING_KEYS if kk in (r.get("mapping") or {})} if r.get("mapping") else None,
               "join": r.get("join") and {kk: r["join"][kk] for kk in ("rows", "targets", "traits", "table_rows_with_no_subject", "their_values")},
               "joined": bool(r.get("joined")), "outcome": _outcome(study, rerun, r, k)}
        if release == "pharos319" and study == "D":
            replay = json.loads((paths.DOCS_DATA / f"cheng_reference_replay_{release}.json").read_text())
            if replay["source_sha256"] != r["file_sha256"] or replay["full"]["graph_genes"] != r["join"]["targets"]:
                raise ValueError("The Cheng reference replay does not match the recorded walkthrough")
            differences = ", ".join(f"{row['symbol']} ({', '.join(row['full'])} with full tables)"
                                    for row in replay["resolver_differences"])
            out["reference_note"] = (f"This run used the full Pharos reference tables. Replaying the same mapping "
                                     f"with the downloadable bundle gives {replay['bundle']['graph_genes']} genes. "
                                     f"Symbol lookup differences: {differences}.")
        if study == "A":
            out["check"] = {"against": "the table built by hand", **{kk: checks.get(kk) for kk in ("hand_pairs", "model_pairs", "reproduced")}}
        elif checks:
            out["check"] = {"against": "the genes the GWAS Catalog reports for the study's associations" if "catalog_reported_genes" in checks
                            else "the genes the GWAS Catalog maps by position to the study's associations",
                            **{kk: checks.get(kk) for kk in ("catalog_study", "abstract_loci", "rows_band", "a_table_of_the_studys_loci", "catalog_associations",
                                                            "of_them_graph_genes", "table_genes", "shared", "share_of_table_genes_the_catalog_maps",
                                                            "share_of_catalog_genes_the_table_holds", "reading")}}
        studies.append(out)
    return {"release": release, "built_from": [add.record_path(release, st, rr).name for st, rr, _ in RUNS if add.record_path(release, st, rr).is_file()]
            + [f"add_compare_{release}.json"], "brief_sheet": add.BRIEF1, "brief_mapping": add.brief2(), "studies": studies}


def write(release: str) -> Path:
    import os
    f = _walkthrough_file(release)
    tmp = f.with_suffix(".tmp")
    tmp.write_text(json.dumps(walkthrough(release), indent=1, ensure_ascii=False) + "\n")
    os.replace(tmp, f)
    return f


# ---- self-test ---------------------------------------------------------------------------------------

def _selftest() -> int:
    ok = True

    def ck(name: str, cond: bool) -> None:
        nonlocal ok
        ok &= bool(cond)
        print(f"  {'PASS' if cond else 'FAIL'}  {name}")

    def why(target: str, trait: str, as_of: int | None = None) -> dict:
        return {"target": target, "trait": trait, "as_of": as_of}

    def predict(node: str, as_of: int | None = None, top: int = 10, unlinked: bool = False) -> dict:
        return {"node": node, "top": top, "unlinked": unlinked, "x": float("nan")}
    A = {"why": why, "predict": predict}
    ck("a tool's arguments are read from its signature", params(why) == {"target": ("str", True), "trait": ("str", True), "as_of": ("int", False)})
    ck("a URL's arguments become the function's, typed",
       answer(A, "why", {"target": ["TREM2"], "trait": ["AD"], "as_of": ["2015"]}) == (200, {"target": "TREM2", "trait": "AD", "as_of": 2015}))
    ck("a boolean and an int are read, and NaN becomes null",
       answer(A, "predict", {"node": ["x"], "top": ["5"], "unlinked": ["true"]}) == (200, {"node": "x", "top": 5, "unlinked": True, "x": None}))
    ck("a missing argument, an unknown one, a bad number and an unknown question are refused",
       answer(A, "why", {"target": ["TREM2"]})[0] == 400 and answer(A, "why", {"target": ["a"], "trait": ["b"], "x": ["1"]})[0] == 400
       and answer(A, "why", {"target": ["a"], "trait": ["b"], "as_of": ["soon"]})[0] == 400 and answer(A, "nope", {})[0] == 404)
    ck("every page file is there", all((HERE / f).is_file() for f, _ in STATIC.values()))
    text = "".join((HERE / f).read_text() for f in ("index.html", "add.html"))
    ck("the pages load nothing from another site", "http://" not in text.replace("http://127.0.0.1", "") and "src=\"https://" not in text
       and "href=\"https://fonts" not in text)
    js = (HERE / "app.js").read_text() + (HERE / "add.js").read_text()
    ck("the scripts never write a server's text as HTML", "innerHTML" not in js and "insertAdjacentHTML" not in js and "document.write" not in js)
    return 0 if ok else 1


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if "--selftest" in argv:
        return _selftest()
    if len(argv) == 2 and argv[0] == "walkthrough":
        print(write(argv[1]))
        return 0
    raise SystemExit(__doc__)


if __name__ == "__main__":
    sys.exit(main())
