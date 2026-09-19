# Documentation Map

Last checked: 2026-09-05

Use this file to decide which documents are current guidance and which ones are
kept only for provenance.

## Current Entry Points

| Need | File |
|---|---|
| How to run the repo and where outputs live | [`../README.md`](../README.md) |
| Plain-language team summary | [`team_summary/LTH_results_summary.html`](team_summary/LTH_results_summary.html) |
| Prioritized next-analysis roadmap | [`next_analysis_roadmap.md`](next_analysis_roadmap.md) |
| Analysis and diagnostic inventory | [`analysis_diagnostic_inventory.md`](analysis_diagnostic_inventory.md) |
| Full phenotype results and limitations | [`../RESULTS.docx`](../RESULTS.docx) |
| Current figure catalog | [`../figures/FIGURE_INDEX.md`](../figures/FIGURE_INDEX.md) |
| Data-source status | [`provenance/data_integration_status.md`](provenance/data_integration_status.md) |
| Trinity/Ross Hauru colony context | [`provenance/trinity_conn_hauru_context.md`](provenance/trinity_conn_hauru_context.md) |
| RNA-seq handoff and analysis notes | [`rnaseq/README.md`](rnaseq/README.md) |
| RNA-seq expression-phenotype plan | [`rnaseq/expression_integration_analysis_plan.md`](rnaseq/expression_integration_analysis_plan.md) |

## Provenance, Not Current Guidance

| Folder | Status |
|---|---|
| `provenance/drive_methods/` | Google Drive method exports and photo-source notes. These are useful for tracing decisions, but the root README and analysis scripts define the current pipeline. |
| `../notes/archive/` | Historical plans, meeting notes, and resolved analysis decisions. Keep them for provenance; do not use them as the current design without checking against the root README and `code/_run_all.R`. |

## Current Terminology

- Use **source patch** for A, C, and D in reader-facing docs.
- Use `thicket` or `genet` only when referring to legacy column names or code
  internals.
- Treat A, C, and D as source-patch labels, not confirmed genetic individuals,
  until final SNP/kinship results are available.
- Mark Rachael Bay's 2026-09-02 SNP cluster file as preliminary wherever it is
  used.
- Treat the Trinity/Ross Hauru table as preliminary candidate-match context only. Source
  patch C is closest to `Apul-115` by coordinate, but exact identity requires
  comparing LTH RNA-seq SNPs to Trinity/Ross whole-genome data.
- Use *Acropora cf. pulchra* as the current taxonomy caveat for manuscript-facing
  wording, while leaving historical repo filenames and code labels unchanged.
- Use **healing** for early tissue re-covering of the wound and
  **regeneration** for later rebuilding steps such as tip extension and new
  radial corallite budding.
- Treat the current late temperature-log coverage and microscope healing sample
  size as open coauthor-summary checks until Molly's follow-up is reconciled.
