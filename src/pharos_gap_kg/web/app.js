"use strict";
// The questions page. Every answer is the server's own JSON from /api/<question>, the function an MCP client's tool call reaches.
// Text from the server is only ever set as text (textContent), never parsed as HTML.

const $ = (s, el = document) => el.querySelector(s);

// Put children into an element, leaving out the missing ones (append() would write null as text).
function put(el, ...kids) {
  el.replaceChildren(...kids.flat(Infinity).filter((c) => c !== null && c !== undefined && c !== false && c !== ""));
  return el;
}

function h(tag, attrs, ...kids) {
  const e = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs || {})) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") e.className = v;
    else if (k === "text") e.textContent = v;
    else if (k.startsWith("on")) e.addEventListener(k.slice(2), v);
    else e.setAttribute(k, v === true ? "" : String(v));
  }
  for (const c of kids.flat(Infinity)) {
    if (c === null || c === undefined || c === false || c === "") continue;
    e.append(c instanceof Node ? c : document.createTextNode(String(c)));
  }
  return e;
}

async function api(question, args) {
  const u = new URL("/api/" + question, location.origin);
  for (const [k, v] of Object.entries(args || {})) if (v !== undefined && v !== null && v !== "") u.searchParams.set(k, v);
  let r, body;
  try {
    r = await fetch(u, { headers: { accept: "application/json" } });
    body = await r.json();
  } catch (e) {
    return { ok: false, body: { error: "The server could not be reached. Try again in a moment." }, url: u.pathname + u.search };
  }
  if (r.status === 429) body = { error: "Too many questions in a minute from here. Wait a few seconds and ask again." };
  return { ok: r.ok, body, url: u.pathname + u.search };
}

// ---- reading an answer ---------------------------------------------------------------------------

const human = (k) => k.replace(/_/g, " ").replace(/^./, (c) => c.toUpperCase());
const kindOf = (id, kind) => (kind === "target" || /^UniProtKB:/.test(id) ? "target" : kind === "drug" || /^DrugCentral:/.test(id) ? "drug" : "trait");
const n = (x) => (typeof x === "number" ? x.toLocaleString("en-US") : x);

function say(text, cls) {
  return text ? h("p", { class: "say" + (cls ? " " + cls : "") , text }) : null;
}
function note(label, text, cls) {
  if (!text) return null;
  return h("p", { class: "note" + (cls ? " " + cls : "") }, h("b", { text: label + ". " }), text);
}

function nodeLink(row) {
  const id = row.id || row.drug_id;
  const label = row.symbol || row.name || row.drug || id;
  if (!id || /^pharosgapkg:/.test(id)) return document.createTextNode(label || "");
  const kind = kindOf(id, row.kind);
  const a = h("a", { href: "?" + new URLSearchParams({ id, name: label, kind }), text: label, title: id });
  a.addEventListener("click", (ev) => { ev.preventDefault(); selectNode({ id, name: label, kind }); window.scrollTo({ top: $("#node").offsetTop - 12, behavior: "smooth" }); });
  return a;
}

function tiedBy(t) {
  if (!t || typeof t !== "object") return String(t ?? "");
  const parts = [];
  if (t.gwas_link) parts.push(`GWAS link (${n(t.gwas_link.studies)} studies, first ${t.gwas_link.first_reported})`);
  if (t.sources) parts.push(`${t.sources.length} sources: ${t.sources.join(", ")}`);
  for (const [k, v] of Object.entries(t)) if (k !== "gwas_link" && k !== "sources") parts.push(`${human(k)}: ${cell(v)}`);
  return parts.join("; ");
}

function cell(v, key) {
  if (v === null || v === undefined) return "";
  if (typeof v === "boolean") return v ? "yes" : "no";
  if (typeof v === "number") return PLAIN.has(key) ? String(v) : n(v);
  if (Array.isArray(v)) {
    if (!v.length) return "";
    const s = v.map((x) => (x && typeof x === "object" ? x.symbol || x.name || x.id || JSON.stringify(x) : String(x)));
    return s.length > 8 ? s.slice(0, 8).join(", ") + ` and ${s.length - 8} more` : s.join(", ");
  }
  if (typeof v === "object") {
    if (key === "tied_by") return tiedBy(v);
    if (v.name || v.symbol) return v.name || v.symbol;
    return Object.entries(v).map(([k, x]) => `${human(k)}: ${typeof x === "object" ? JSON.stringify(x) : x}`).join("; ");
  }
  return String(v);
}

const COLS = [["rank", "Rank"], ["symbol", "Gene"], ["genes", "Genes"], ["name", "Name"], ["drug", "Drug"], ["study", "Study"], ["source", "Source"],
  ["condition", "Condition"], ["trait", "Trait"], ["pharos_level", "Level"], ["kind", "Kind"], ["first_reported", "First reported"], ["year", "Year"],
  ["studies_in_release", "Studies"], ["studies", "Studies"], ["pmid", "PubMed"], ["score", "Score"], ["action", "Action"],
  ["is_the_drugs_mechanism", "Its mechanism"], ["mechanism_targets", "Mechanism targets"], ["direction", "Expression"], ["log2_fold_change", "log2 fold change"],
  ["p", "p"], ["sources", "Sources"], ["tied_by", "Tied by"], ["lacks", "Lacks"], ["nearest_route", "Nearest route"], ["criteria_to_clear", "Criteria to clear"],
  ["promoted_in_its_band", "Promoted in its band"], ["distance", "Distance"], ["rank_among_gwas_candidates", "Rank among GWAS candidates"],
  ["explained_by", "Already explained by"], ["same_locus", "Same locus"], ["genes_with_the_link", "Genes with the link"], ["proposal", "What to test"],
  ["locus", "Locus"], ["source_evidence", "Source's evidence"]];
const LINKED = new Set(["symbol", "name", "drug"]);
const PLAIN = new Set(["first_reported", "year", "since", "pmid", "rank"]);      // numbers that are years or ids: no thousands separator

function table(rows, opts = {}) {
  if (!rows || !rows.length) return null;
  const present = COLS.filter(([k]) => rows.some((r) => r && r[k] !== undefined && r[k] !== null && !(Array.isArray(r[k]) && !r[k].length)));
  let cols = present.slice(0, opts.max || 7);
  if (!cols.length) return null;
  if (cols.some(([k]) => k === "symbol") ) cols = cols.filter(([k]) => k !== "name" || !rows.every((r) => r.symbol));
  const head = h("tr", null, cols.map(([, l]) => h("th", { scope: "col", text: l })));
  const body = rows.map((r) => h("tr", null, cols.map(([k]) => {
    const v = r[k];
    if (LINKED.has(k) && (r.id || r.drug_id) && (k === "symbol" || (k === "name" && !r.symbol) || k === "drug")) return h("td", null, nodeLink(r));
    if (k === "pmid" && v) return h("td", { class: "num" }, h("a", { href: `https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(v)}/`, rel: "noopener noreferrer", text: v }));
    if (k === "genes" && Array.isArray(v)) return h("td", null, v.map((g, i) => [i ? ", " : "", g]));
    const num = typeof v === "number";
    return h("td", { class: num ? "num" : null }, k === "pharos_level" && v ? h("span", { class: "lvl", text: v }) : cell(v, k));
  })));
  return h("div", { class: "table" }, h("table", null, h("thead", null, head), h("tbody", null, body)));
}

function list(title, rows, more, extra) {
  if (!rows || !rows.length) return null;
  return h("div", { class: "list" }, h("h3", null, title, more ? h("span", { class: "count", text: `${n(rows.length)} shown, ${n(more)} more` }) : null), extra || null, table(rows));
}

function otherSources(o) {
  if (!o) return null;
  const expr = o.expression_observations;
  const exprRows = Array.isArray(expr) ? expr : expr && expr.largest;
  return h("div", { class: "list" },
    h("h3", null, "Beside the GWAS links, as of the release"),
    say(o.answer),
    list("Asserted with", o.asserted_with, o.asserted_more),
    list("Sources that assert it", o.asserted_by),
    list("Drugs that act on it", o.drugs, o.drugs_more),
    list("Drugs indicated for it", o.drugs_indicated, o.drugs_more),
    list("Through a drug", o.through_a_drug),
    list("Rows a pruning rule withheld", o.withheld),
    exprRows && exprRows.length ? list("Expression observations", exprRows.map((r) => ({ ...r, ...(r.condition && typeof r.condition === "object" ? { condition: r.condition.name } : {}) }))) : null,
    note("Expression", expr && expr.reading));
}

function afterRelease(a) {
  if (!a || !a.answer) return null;
  return h("div", { class: "list" }, h("h3", null, "Added after the release"), note("Added", a.answer, "added"), note("Note", a.note), table(a.top || a.links));
}

function gapsList(g, ctx) {
  const count = g.targets ?? g.loci ?? g.rows;
  const box = h("div", { class: "list" });
  const rows = h("div");
  const show = (top) => put(rows, table(top) || h("p", { class: "small", text: "None." }));
  show(g.top);
  let btn = null;
  if (g.more > 0 && ctx && ctx.disease) {
    btn = h("button", { class: "ghost", type: "button", text: `Show up to 25` });
    btn.addEventListener("click", async () => {
      btn.disabled = true; btn.textContent = "Loading…";
      const r = await api("gaps", { disease: ctx.disease, kind: g.kind, top: 25 });
      const got = r.ok && (r.body.gaps || []).find((x) => x.kind === g.kind);
      if (got) { show(got.top); btn.remove(); } else { btn.textContent = (r.body && r.body.error) || "Could not load"; }
    });
  }
  put(box,
    h("h3", null, human(g.kind), count !== undefined ? h("span", { class: "count", text: `${n(count)}${g.genes ? ` loci, ${n(g.genes)} genes` : ""}` }) : null),
    say(g.what),
    g.closes_with ? h("p", { class: "closes" }, h("b", { text: "Closes with: " }), g.closes_with) : null,
    g.how_often_this_closed ? h("p", { class: "closes" }, h("b", { text: "How often it has closed: " }), g.how_often_this_closed) : null,
    note("Reading", g.reading), rows, btn);
  return box;
}

const HANDLED = new Set(["found", "answer", "reading", "caveat", "method", "note", "calibration", "source", "node", "target", "trait", "condition", "as_of",
  "top", "more", "candidates", "first", "studies", "studies_more", "targets", "targets_more", "indications", "indications_more", "gaps", "counts", "other_sources",
  "after_the_release", "problem", "matches", "claim", "knowledge", "edge", "admitted_by", "read_under", "not_listed"]);

function render(b, ctx) {
  const out = [];
  if (b.error) return [h("p", { class: "err", text: b.error })];
  out.push(say(b.answer, "lead"));
  if (b.found === false && Array.isArray(b.candidates) && b.candidates.length) {
    out.push(h("p", { class: "small", text: "Which is meant?" }), h("ul", { class: "matches" }, b.candidates.slice(0, 12).map((m) => h("li", null, matchButton(m, ctx && ctx.pick)))));
  }
  out.push(note("Reading", b.reading), note("Caveat", b.caveat), note("Note", b.note), note("Calibration", b.calibration));
  if (b.method) out.push(h("details", null, h("summary", { text: "How the candidates are ranked" }), say(b.method)));
  if (b.counts && b.gaps) {
    out.push(h("div", { class: "list" }, h("h3", null, "The six lists"), h("div", { class: "table" }, h("table", null, h("tbody", null,
      Object.entries(b.counts).map(([k, v]) => h("tr", null, h("td", { text: human(k) }), h("td", { class: "num", text: n(v) }))))))));
  }
  if (Array.isArray(b.gaps)) for (const g of b.gaps) out.push(g.top !== undefined ? gapsList(g, ctx) : h("div", { class: "list" }, table([g])));
  out.push(list(ctx && ctx.q === "whats_new" ? "First reported since" : "Best-supported links", ctx && ctx.q === "whats_new" ? b.first : b.top, b.more));
  if (b.by_year) out.push(h("p", { class: "small", text: "By year: " + Object.entries(b.by_year).map(([y, c]) => `${y}: ${n(c)}`).join(", ") }));
  out.push(list("Candidates", b.candidates && b.found !== false ? b.candidates : null));
  out.push(list("The GWAS studies behind the link", b.studies, b.studies_more));
  if (b.node && b.node.kind === "drug") { out.push(list("Targets it acts on", b.targets, b.targets_more), list("Indicated for", b.indications, b.indications_more)); }
  for (const [k, v] of Object.entries(b)) {          // any other list in the answer, as a table
    if (HANDLED.has(k) || !Array.isArray(v) || !v.length || typeof v[0] !== "object") continue;
    out.push(list(human(k), v));
  }
  out.push(otherSources(b.other_sources), afterRelease(b.after_the_release));
  return out.filter(Boolean);
}

function tools(res, claude, kgx) {
  const raw = h("details", null, h("summary", { text: "The full answer, as the server gave it (JSON)" }), h("pre", { text: JSON.stringify(res.body, null, 1) }));
  return [raw, h("div", { class: "tools" },
    h("span", null, "As JSON: ", h("a", { href: res.url, text: res.url })),
    kgx ? h("span", null, h("a", { href: kgx, text: "The graph behind this answer (KGX)" })) : null,
    claude ? h("span", null, "In Claude, with the connector: ", h("q", { text: claude })) : null)];
}

// ---- the node and its questions -----------------------------------------------------------------

const QUESTIONS = {
  trait: [["gaps", "What is missing"], ["known", "What is known"], ["predict", "Candidate genes"], ["whats_new", "New since"]],
  target: [["known", "What is known"], ["gaps", "What it lacks"], ["predict", "Traits it may be linked to next"], ["whats_new", "New since"]],
  drug: [["known", "What is known"]],
};
let current = null;

function matchButton(m, pick) {
  const label = m.symbol || m.name || m.id;
  const count = m.links !== undefined ? `${n(m.links)} GWAS links` : m.rows_in_other_sources !== undefined ? `${n(m.rows_in_other_sources)} rows in other sources` : "";
  return h("button", { type: "button", onclick: () => (pick || selectNode)({ id: m.id, name: label, kind: kindOf(m.id, m.kind), detail: m.kind }) },
    h("span", null, label, " ", h("span", { class: "k", text: m.id })), h("span", { class: "k", text: [m.kind, count].filter(Boolean).join(" · ") }));
}

async function lookUp(text, push = true) {
  const box = $("#matches");
  put(box, h("p", { class: "loading", text: "Looking up…" }));
  if (push) history.pushState({}, "", "?" + new URLSearchParams({ q: text }));
  const r = await api("find", { text, kind: "any", limit: 8 });
  if (!r.ok) return put(box, h("p", { class: "err", text: r.body.error }));
  const b = r.body;
  if (b.exact_match && b.matches && b.matches.length) {
    const rest = b.matches.slice(1);
    put(box, rest.length ? h("details", null, h("summary", { text: `${n(rest.length + (b.more || 0))} other names contain "${text}"` }),
      h("ul", { class: "matches" }, rest.map((m) => h("li", null, matchButton(m))))) : null);
    return selectNode({ id: b.matches[0].id, name: b.matches[0].symbol || b.matches[0].name, kind: kindOf(b.matches[0].id, b.matches[0].kind), detail: b.matches[0].kind });
  }
  put(box, say(b.answer), b.matches && b.matches.length ? h("ul", { class: "matches" }, b.matches.map((m) => h("li", null, matchButton(m)))) : null,
    b.more ? h("p", { class: "small", text: `${n(b.more)} more names match; a longer name narrows the list.` }) : null);
}

function selectNode(node, view, push = true) {
  current = node;
  const qs = QUESTIONS[node.kind] || QUESTIONS.trait;
  view = qs.some(([q]) => q === view) ? view : qs[0][0];
  $("#node").hidden = false;
  $("#node-name").textContent = node.name || node.id;
  $("#node-meta").textContent = [node.detail || node.kind, node.id].filter(Boolean).join(" · ");
  const tabs = $("#tabs");
  put(tabs, ...qs.map(([q, label]) => h("button", { type: "button", role: "tab", "aria-selected": q === view ? "true" : "false", "data-q": q,
    onclick: () => ask(q, true) }, label)));
  ask(view, push);
}

function sinceYear() {
  const v = parseInt(new URLSearchParams(location.search).get("since") || "2019", 10);
  return Number.isFinite(v) ? Math.min(2023, Math.max(2005, v)) : 2019;
}

async function ask(q, push) {
  const node = current;
  for (const b of $("#tabs").children) b.setAttribute("aria-selected", b.dataset.q === q ? "true" : "false");
  const box = $("#answer");
  let since = sinceYear();
  if (push) history.pushState({}, "", "?" + new URLSearchParams({ id: node.id, name: node.name, kind: node.kind, view: q, ...(q === "whats_new" ? { since } : {}) }));
  const head = [];
  if (q === "whats_new") {
    const inp = h("input", { type: "number", min: 2005, max: 2023, value: since, "aria-label": "Year" });
    const go = h("button", { class: "ghost", type: "button", text: "Show" });
    go.addEventListener("click", () => { history.replaceState({}, "", "?" + new URLSearchParams({ id: node.id, name: node.name, kind: node.kind, view: q, since: inp.value })); ask(q, false); });
    head.push(h("p", { class: "since" }, "Links first reported after", inp, go, h("span", { text: "(the data run to 2023)" })));
  }
  put(box, ...head, h("p", { class: "loading", text: "Asking the graph…" }));
  const args = q === "gaps" ? (node.kind === "target" ? { target: node.id } : { disease: node.id })
    : q === "whats_new" ? { node: node.id, since } : q === "predict" ? { node: node.id, top: 10 } : { node: node.id, limit: 10 };
  const r = await api(q, args);
  if (current !== node) return;
  const name = node.name || node.id;
  const claude = { gaps: node.kind === "target" ? `Using pharos-gap-kg, what does ${name} lack?` : `Using pharos-gap-kg, what are the gaps for ${name}?`,
    known: `Using pharos-gap-kg, what is known about ${name}?`,
    predict: node.kind === "target" ? `Using pharos-gap-kg, which traits might ${name} be linked to next?` : `Using pharos-gap-kg, which genes might be linked to ${name} next?`,
    whats_new: `Using pharos-gap-kg, what is new for ${name} since ${since}?` }[q];
  const kgx = "/api/subgraph?" + new URLSearchParams(q === "gaps" ? (node.kind === "target" ? { question: "gaps", a: node.id } : { question: "gaps", disease: node.id })
    : q === "whats_new" ? { question: q, a: node.id, b: since } : { question: q, a: node.id });
  put(box, ...head, ...render(r.body, { q, disease: node.kind === "target" ? null : node.id, pick: (m) => selectNode(m) }), ...(r.ok ? tools(r, claude, kgx) : []));
}

async function why(target, trait, push = true) {
  const box = $("#why-answer");
  if (push) history.pushState({}, "", "?" + new URLSearchParams({ target, trait }));
  put(box, h("p", { class: "loading", text: "Asking the graph…" }));
  const r = await api("why", { target, trait });
  const pick = (m) => {
    const side = m.kind === "target" ? "#why-target" : "#why-trait";
    $(side).value = m.id;
    why($("#why-target").value.trim(), $("#why-trait").value.trim());
  };
  const kgx = "/api/subgraph?" + new URLSearchParams({ question: "why", a: target, b: trait });
  put(box, ...render(r.body, { q: "why", pick }), ...(r.ok ? tools(r, `Using pharos-gap-kg, why is ${target} linked to ${trait}?`, r.body.found === false ? null : kgx) : []));
}

async function status() {
  const r = await api("status", {});
  const el = $("#status");
  if (!r.ok) return put(el, h("span", { class: "err", text: r.body.error }));
  const b = r.body, f = b.built_from || {};
  put(el, say(b.answer),
    f.file ? note("Built from", `${f.file}, published ${f.published}, SHA-256 ${f.sha256}`) : null,
    b.checks ? note("Checks", b.checks.join("; ")) : null,
    b.left_out && b.left_out.length ? h("details", null, h("summary", { text: "What it leaves out" }), h("ul", null, b.left_out.map((x) => h("li", { text: x })))) : null,
    b.reading_a_link ? note("Reading a link", b.reading_a_link) : null);
}

function fromUrl() {
  const p = new URLSearchParams(location.search);
  if (p.get("id")) selectNode({ id: p.get("id"), name: p.get("name") || p.get("id"), kind: kindOf(p.get("id"), p.get("kind")) }, p.get("view"), false);
  else if (p.get("q")) { $("#q").value = p.get("q"); lookUp(p.get("q"), false); }
  if (p.get("target") && p.get("trait")) { $("#why-target").value = p.get("target"); $("#why-trait").value = p.get("trait"); why(p.get("target"), p.get("trait"), false); }
}

document.addEventListener("DOMContentLoaded", () => {
  $("#mcp-url").textContent = location.origin + "/mcp";
  $("#ask").addEventListener("submit", (ev) => { ev.preventDefault(); const t = $("#q").value.trim(); if (t) lookUp(t); });
  for (const c of document.querySelectorAll("#examples .chip")) c.addEventListener("click", () => { $("#q").value = c.dataset.q; lookUp(c.dataset.q); });
  $("#why").addEventListener("submit", (ev) => {
    ev.preventDefault();
    const t = $("#why-target").value.trim(), d = $("#why-trait").value.trim();
    if (t && d) why(t, d);
  });
  window.addEventListener("popstate", () => { $("#node").hidden = true; put($("#matches")); put($("#why-answer")); fromUrl(); });
  fromUrl();
  status();
});
