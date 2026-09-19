# Main Figure Plan

Reviewed 2026-09-16. This is a proposed manuscript structure, not a record of newly fitted models or completed figures. Original coauthor drafts remain unchanged. Verbatim draft extracts are in `SOURCE_EXTRACTS_2026-09-16.md`; they contain outdated claims and are not publication-ready text.

## The Main Result

Heat affected wound covering and later regrowth differently. Wounds regained tissue cover quickly at both temperatures in the separate microscope experiment. In the main experiment, all 12 wounded fragments at 28 C had budded new corallites by Day 15, compared with 4 of 12 at 31 C. Heat also reduced coral condition, with smaller changes in fragments from source patch C.

This supports a paper about different responses of early wound repair, later rebuilding, and whole-fragment condition to sustained heat. It does not yet establish a genetic trade-off, heritability, an energy-allocation mechanism, or permanent failure to regenerate. Tissue covering and later budding were measured in different cohorts, so their contrast is not a demonstrated within-fragment sequence across the full sample.

## Proposed Main Figures

| Figure | Question and panels | Existing material | Changes needed |
|---|---|---|---|
| 1. Early wound covering and later regrowth | A: small, labeled photo sequence defining tissue cover, tip rebuilding, and new-corallite budding. B: percentage with tissue cover by day, separate microscope experiment. C: percentage that had budded new corallites by day, main experiment. | Drive regeneration-stage artwork; `39b_stage_timing_cumulative_option`; `11_microscope_trait_trajectories`; summary photo sequence. | Label the cohorts separately: 8 per temperature for tissue cover, 12 for budding. Show 12/12 versus 4/12 at Day 15. Describe cumulative first observations, not whether a feature remained visible. Do not connect the two cohorts as individual trajectories. |
| 2. Whole-fragment condition under heat and clipping | Four panels: photosynthesis efficiency over time, color over time, symbiont density across sampling days, and skeletal mass gain over the study. | `02_pam_fvfm_trajectory`, `03_color_trajectory`, `06_symbiont_density_by_day`, `05_buoyant_weight_growth`. | Use consistent temperature colors and wound-state symbols/lines. Keep natural units and observed data. Distinguish repeated measurements from destructive sampling. Show tank-aware uncertainty, not fragment standard errors presented as treatment confidence intervals. No chlorophyll panel: it was not measured. |
| 3. Source patches differed in their response to heat | Three panels showing 28 versus 31 C values within A, C, and D for final photosynthesis efficiency, color, and symbiont density. Add growth as a fourth panel only if space allows, explicitly noting weaker evidence. | `13_genet_response_panel`, `19f_source_physiology_heat_penalties`, current fitted contrasts. | Rebuild with source-patch labels, individual observations, and model-based comparisons with uncertainty. Photosynthesis and color are Day 14; symbionts are Day 15. Growth covers the full interval. State whether estimates average wound states; inspect wound-specific contrasts before pooling. Avoid the standardized composite as the main evidence. |

Use three main figures for the current phenotype story. A fourth figure should be chosen after the gene-expression analysis produces a defensible result, not reserved for a claim that has not been tested. If genetic assignments become central, reassess Figure 3 rather than automatically adding another exploratory figure.

The existing `16_manuscript_fig1` mixes temperature validation, condition, morphology, and a multivariate summary. It is an overview, not the strongest opening manuscript figure. The Drive draft's repair-first, physiology-second order provides a better starting point.

## Supporting Figures

- Experimental design and temperature validation: distinguish the 192-fragment main experiment from the separate 16-fragment microscope experiment. Show four tanks per temperature.
- Full regeneration-stage curves, exact observation days, and stage means. Means among fragments that reached a stage cannot represent the whole treatment when many did not reach it; keep reached/total counts visible. The final heated-stage mean uses only four fragments.
- Wound-specific and source-patch-specific results, including growth and model sensitivity checks.
- Trade-off analyses: useful tests of an alternative explanation, but no clear growth cost of regeneration has emerged. The source-level comparisons have only three patches, and growth versus growth lost shares a quantity between axes.
- Preliminary DNA-marker structure and its directly matched symbiont analysis. Current DNA-marker data do not identify the Day 15 fragments with the broader phenotype record.
- Multivariate condition plots, standardized heat-penalty summaries, and diagnostic plots.
- Preliminary Hauru colony comparisons only as clearly labeled context. Coordinate proximity does not confirm colony identity, and a short heat-test threshold is not a safety threshold for prolonged exposure.

## Where the Drafts Stand

| Source | Useful content | What needs reconciliation |
|---|---|---|
| [Drive working draft](https://docs.google.com/document/d/1cswu3uOFNRTIa0vtrWiJqXRWI8eUC1Y-0mlKb8TjSY0/edit), modified September 9 | Coauthor framing, detailed husbandry and measurement methods, repair-first figure order. | Results are mostly placeholders. Title assumes a genetic trade-off. Chlorophyll is described despite not being run. Collected tissues and sequenced libraries need separating. Main morphology cohort is missing from methods. Husbandry values need confirmation, not automatic transfer. |
| [Drive manuscript folder](https://drive.google.com/drive/folders/1r5wpaurHmVS4KEJ9uL4re_xFYI9Lk7sI) | Writing and Figures subfolders. Figure 1 folder contains regeneration-stage PNG and Illustrator artwork. | Figure 2 physiology folder was empty at review. Methods and Results contain competing Figure 1 uses; renumber after agreeing on the structure. |
| `manuscript/Manuscript_LTH.docx` | Most developed local Methods and Results narrative; model and design descriptions. | Several claims are stale or too strong; see corrections below. Abstract, introduction, expression results, and discussion are incomplete. |
| `RESULTS.docx`, internally dated June 30 | Detailed phenotype narrative and statistical reporting. | Not authoritative merely because it says so. Sample counts, terminology, and some tests lag current code and outputs. |
| `docs/team_summary/LTH_results_summary.Rmd` and rendered HTML | Most current plain-language synthesis, cohort distinctions, endpoint timing, limitations, and expanded plots. | Broader than a paper's main figure set. Reconcile claims with model-specific outputs rather than copying every panel. |

## Corrections Before Manuscript Assembly

1. Use **source patch A, C, and D**, not confirmed genets. Remove claims of heritable variation. Preliminary markers indicate that patch labels are not interchangeable with genetic identity.
2. Replace "never regenerated" and "never rebuilt skeleton" with the specific observation: fewer fragments had budded new corallites by Day 15. Many heated fragments formed or extended tips. Later recovery was not measured.
3. Keep direct tissue-cover observations from the separate microscope experiment distinct from the main morphology scores. Wound smoothing and central-polyp formation are not interchangeable with a direct tissue-cover score.
4. Describe 48 Day 15 fragments, including 24 wounded fragments, separately from 144 early destructive samples and the separate 16-fragment microscope experiment. Photosynthesis efficiency and color end on Day 14; growth, final symbionts, and regeneration extend to Day 15.
5. Replace the generic four-response repeated-measures model description with response-specific formulas. Growth has one endpoint per fragment; destructively sampled symbionts are not repeated observations of the same fragment.
6. Do not label the script 13 likelihood-ratio tests as isolated temperature-by-patch tests. They compare all patch interactions with an additive-patch model, including patch interactions with wound state and time. For growth, the current table gives chi-square 12.2876 with 6 degrees of freedom, p = 0.05585; the manuscript's 6.5 and p = 0.37 do not match it. Report heat-specific contrasts separately.
7. Report lower estimated growth with its uncertainty: approximately 34% lower mean mass gain, with the tank-level permutation p = 0.057. Do not turn this into either clear evidence of no effect or an unequivocally established reduction.
8. Do not present the narrow interval around the existing 1.32-fold budding-time estimate without noting that its uncertainty does not account for shared tanks. Resolve this before making the model estimate a headline. The observed 4/12 versus 12/12 comparison is descriptive, not a substitute tank-adjusted test.
9. Replace the old temperature-compliance statement with the verified record: treatment-wide daily means averaged 27.78 and 30.87 C over Days 0-15, a 3.09 C difference. This does not assert that every tank-day was within 0.3 C of target.
10. Remove chlorophyll results and distinguish tissue collected from libraries sequenced. The current 144-library expression design uses Days 1, 3, and 10, not Day 15. Do not predeclare expression as the lead result while results are pending.
11. Explain that little additional whole-fragment response to clipping does not establish the absence of a local physiological or gene-expression response. Likewise, an uncertain association does not prove regeneration is cost-free.

## Assembly Order and Checks

First reconcile the cohort/timing table and trait definitions, then assemble Figure 1. Rebuild Figures 2 and 3 from reproducible scripts with the agreed labels and panel order. Update Methods around the actual sampling units and response-specific models. Write Results in the same order as the figures, then choose the title and discuss expression findings when available.

Before using fitted estimates, check model-specific convergence, residual behavior, separation or censoring, influential fragments/tanks, and the treatment-level replication. Use `docs/analysis_diagnostic_inventory.md` to locate completed checks and gaps; this review did not rerun the models or certify their diagnostics.

Numerical checks for this plan: `output/tables/13_genet_anova.csv`, `05_buoyant_weight_tank_test.csv`, `37_tradeoff_main_results.csv`, and the current team-summary source. Preserve the 115/121 score-review caveat and exclusion sensitivity in supporting material until resolved; do not silently change their scores to make the stage sequence cleaner.
