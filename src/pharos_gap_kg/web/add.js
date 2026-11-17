"use strict";
// The walkthrough of "Adding a table", from /walkthrough.json, which is built from the committed records of the runs.
// Text from the records is only ever set as text (textContent), never parsed as HTML.

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

const n = (x) => (typeof x === "number" ? x.toLocaleString("en-US") : x);
const pct = (x) => `${Math.round(x * 100)}%`;
function problem(p) {      // "4_domain_range [value]: text" -> "Check 4 (domain range, value): text"
  const m = /^(\d)_([a-z_]+)(?: \[([^\]]*)\])?: (.*)$/s.exec(p);
  return m ? h("li", null, h("b", { text: `Check ${m[1]} (${m[2].replace(/_/g, " ")}${m[3] ? ", " + m[3] : ""}): ` }), m[4]) : h("li", { text: p });
}

function verdict(ok) {
  return h("span", { class: "verdict " + (ok ? "yes" : "no"), text: ok ? "accepted" : "refused" });
}

function kv(pairs) {
  const dl = h("dl", { class: "kv" });
  for (const [k, v] of pairs) if (v !== null && v !== undefined && v !== "") dl.append(h("dt", { text: k }), h("dd", null, v));
  return dl;
}

function traitName(run, id) {
  const t = (run.traits_on_the_card || []).find((x) => x.id === id);
  return t ? `${t.name} (${id})` : id;
}

function side(run, s) {
  if (!s) return "";
  if (s.constant) return `every row: ${traitName(run, s.constant)}, as ${s.category}`;
  const cols = (s.columns || []).map((c) => `"${c}"`).join(", then ");
  const how = s.resolver === "gene_symbol" ? "read as a gene symbol" : s.resolver && s.resolver.startsWith("name:") ? `made into a local node with no identifier (${s.resolver})` : `read with ${s.resolver}`;
  return `the column ${cols}, ${how}${s.separator ? `, split on "${s.separator}"` : ""}, as ${s.category}`;
}

function mappingFields(run, p) {
  if (!p) return null;
  if (p.map === false) return kv([["Proposes", "that this table should not become edges"], ["Its reason", h("span", { class: "quote", text: p.reason || "" })]]);
  const keep = p.keep_if ? `${p.keep_if.column} ${p.keep_if.matches ? `matches ${p.keep_if.matches}` : p.keep_if.in ? `is one of ${p.keep_if.in.join(", ")}` : `is at least ${p.keep_if.at_least}`}: ${p.keep_if.reason || ""}` : null;
  const split = p.split_by ? `${p.split_by}: ${Object.entries(p.by_value || {}).map(([v, r]) => (r && r.skip ? `${v} left out (${r.skip})` : `${v} mapped`)).join("; ")}` : null;
  return kv([["Subject", side(run, p.subject)], ["Object", side(run, p.object)], ["Edge", `${p.predicate || ""} in ${p.association_category || ""}`],
    ["Knowledge", [p.knowledge_level, p.agent_type].filter(Boolean).join(", ")], ["Keeps a row if", keep], ["Split by", split],
    ["Carries", (p.edge_properties || []).join(", ")], ["The model's note", p.notes ? h("span", { class: "quote", text: p.notes }) : null]]);
}

function tries(run, list, stage) {
  return list.map((t) => h("div", { class: "try " + (t.accepted ? "yes" : "no") },
    h("div", { class: "who" }, h("b", { text: `${t.model}${t.attempt ? `, attempt ${t.attempt}` : ""}` }), verdict(t.accepted)),
    stage === "sheet" ? (t.proposal ? kv([["Sheet", t.proposal.sheet], ["Column names in row", t.proposal.header_row],
      ["The model's note", t.proposal.notes ? h("span", { class: "quote", text: t.proposal.notes }) : null]]) : null) : mappingFields(run, t.proposal),
    t.problems && t.problems.length ? h("ul", null, t.problems.map(problem)) : null));
}

function sheets(run, chosen) {
  if (!run.sheets || !run.sheets.length) return null;
  return h("details", null, h("summary", { text: `The file's ${run.sheets.length} sheets, by name and size` }),
    h("ul", { class: "sheets" }, run.sheets.map((s) => h("li", { class: chosen.includes(s.name) ? "chosen" : null, text: `${s.name}: ${n(s.rows)} rows, ${n(s.columns)} columns` }))));
}

function checkBlock(run) {
  const c = run.check;
  if (!c) return null;
  if (c.reproduced !== undefined) {
    return h("div", { class: "stage" }, h("h3", { text: "The check fixed before the run" }),
      h("p", null, `Against ${c.against}: ${n(c.model_pairs)} pairs from the model's mapping, ${n(c.hand_pairs)} by hand; `, c.reproduced ? "the same pairs." : "not the same pairs."));
  }
  const rows = [];
  if (c.rows_band) rows.push(["Rows of the sheet", `${n(run.rows)}; the abstract reports ${n(c.abstract_loci)}, so a table of the study's loci has ${n(c.rows_band[0])} to ${n(c.rows_band[1])}: ${c.a_table_of_the_studys_loci ? "it is" : "it is not"}`]);
  if (c.catalog_study) rows.push(["GWAS Catalog study", `${c.catalog_study}, ${n(c.catalog_associations)} associations`]);
  if (c.share_of_table_genes_the_catalog_maps && c.share_of_catalog_genes_the_table_holds) {
    const a = c.share_of_table_genes_the_catalog_maps, b = c.share_of_catalog_genes_the_table_holds;
    rows.push(["Genes shared", `${n(c.shared)}: ${pct(a[0])} of the table's ${n(c.table_genes)} [${pct(a[1])}, ${pct(a[2])}], ${pct(b[0])} of the Catalog's ${n(c.of_them_graph_genes)} graph genes [${pct(b[1])}, ${pct(b[2])}]`]);
  }
  rows.push(["Reading", c.reading]);
  return h("div", { class: "stage" }, h("h3", { text: "The check fixed before the run" }), h("p", { class: "small", text: `Against ${c.against}. Both shares at least 80% read as agree, either below 50% as disagree, and between as partial.` }), kv(rows));
}

function show(run) {
  const back = run.back_to_the_sheet;
  const chosen = [run.stage1.find((t) => t.accepted)?.proposal?.sheet, run.sheet].filter(Boolean);
  const box = $("#run");
  put(box, 
    h("h2", { text: `${run.citation}: ${run.disease}` }),
    h("p", { class: "small", text: run.what }),
    h("p", { class: "outcome", text: run.outcome }),
    h("div", { class: "stage" }, h("h3", { text: "What the person gives" }),
      h("p", { class: "quote", text: run.description }),
      kv([["File", run.file], ["Licence", run.licence], ["Check of the rows", run.check_8]])),
    h("div", { class: "stage" }, h("h3", { text: "1. The sheet" }), sheets(run, chosen), tries(run, run.stage1, "sheet")),
    run.stage2.length ? h("div", { class: "stage" }, h("h3", { text: "2. The mapping" }), tries(run, run.stage2, "mapping")) : null,
    back ? h("div", { class: "stage" }, h("h3", { text: "Back to the sheet" }),
      h("p", null, "The mapping stage said this table should not become edges. Its reason went back to the sheet stage once, and the model chose again."),
      tries(run, back.stage1, "sheet"),
      back.stopped ? h("p", { class: "small", text: `The stages stopped: ${back.stopped}.` }) : null,
      back.stage2.length ? [h("h3", { text: "2. The mapping, on the new sheet" }), tries(run, back.stage2, "mapping")] : null) : null,
    h("div", { class: "stage" }, h("h3", { text: "3. The person" }),
      run.join ? [
        h("p", null, `The accepted mapping reads "${run.sheet}" (${n(run.rows)} rows) and would make ${n(run.join.rows)} links from ${n(run.join.targets)} graph genes to ${(run.join.traits || []).map((t) => traitName(run, t)).join(", ")}.`),
        run.reference_note ? h("p", { text: run.reference_note }) : null,
        run.join.table_rows_with_no_subject ? h("p", { class: "small", text: `${n(run.join.table_rows_with_no_subject)} rows give no graph gene; among their values: ${(run.join.their_values || []).join(", ")}.` }) : null,
        h("p", { class: "small", text: run.joined ? "The person joined it." : run.study === "A" ? "Compared with the table built by hand, not joined: that table is already beside the graph." : "Not joined. It joins only when the person says so." })]
        : h("p", { class: "small", text: "No accepted mapping makes edges, so there is nothing to join." })),
    checkBlock(run));
}

function label(run) {
  return { A: "Alzheimer disease", B: "Parkinson disease", C: "Type 1 diabetes", D: "Refractive error" }[run.study] + (run.rerun ? ", again" : "");
}

document.addEventListener("DOMContentLoaded", async () => {
  let w;
  try { w = await (await fetch("/walkthrough.json")).json(); } catch (e) { put($("#run"), h("p", { class: "err", text: "The records could not be loaded." })); return; }
  $("#brief-sheet").textContent = w.brief_sheet || "";
  $("#brief-mapping").textContent = w.brief_mapping || "";
  $("#built-from").textContent = w.built_from ? `Built from the committed records ${w.built_from.join(", ")}.` : "";
  const runs = w.studies || [];
  if (!runs.length) { put($("#run"), h("p", { class: "small", text: w.note || "No runs are recorded for this release." })); return; }
  const bar = $("#runs");
  const pick = (id, push) => {
    const run = runs.find((r) => r.id === id) || runs[runs.length - 1];
    for (const b of bar.children) b.setAttribute("aria-selected", b.dataset.id === run.id ? "true" : "false");
    if (push) history.replaceState({}, "", "?" + new URLSearchParams({ run: run.id }));
    show(run);
  };
  put(bar, ...runs.map((r) => h("button", { type: "button", role: "tab", "data-id": r.id, onclick: () => pick(r.id, true) },
    label(r), h("small", { text: r.rerun ? "development" : r.study === "A" ? "built by hand before" : "unseen" }))));
  pick(new URLSearchParams(location.search).get("run") || runs[runs.length - 1].id, false);
});
