# Documentation Map

Last checked: 2026-09-03

Use this file to decide which documents are current guidance and which ones are
kept only for provenance.

## Current Entry Points

| Need | File |
|---|---|
| How to run the repo and where outputs live | [`../README.md`](../README.md) |
| Plain-language team summary | [`team_summary/LTH_results_summary.html`](team_summary/LTH_results_summary.html) |
| Full phenotype results and limitations | [`../RESULTS.docx`](../RESULTS.docx) |
| Data-source status | [`provenance/data_integration_status.md`](provenance/data_integration_status.md) |
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
