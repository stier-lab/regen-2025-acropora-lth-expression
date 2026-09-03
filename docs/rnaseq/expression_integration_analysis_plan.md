# Tentative RNA-seq expression-phenotype integration plan

Updated: 2026-09-02

This is a planning document for the LTH heat x wound RNA-seq analysis. It is
written before the expression counts arrive so the analysis does not drift toward
the phenotype story only after we see the transcriptome. The core rule is:

**Run the expression analysis genome-wide first; use physiology and regeneration
results to define contrasts and interpret modules, not to pre-select winners.**

The plan should be revised after the raw-data provenance audit and sequencing QC,
but any revision should record what changed and why.

## 1. Current design

The current RNA-seq design has 144 margin libraries:

| Factor | Levels |
|---|---|
| Temperature | 28 C, 31 C |
| Wound state | wounded, unwounded |
| Day | 1, 3, 10 |
| Source thicket | A, C, D |
| Tank | 4 tanks per temperature |
| Tissue position | margin only |

Useful repo anchors:

| Need | File |
|---|---|
| Phenotype source of truth | `output/tables/20_master_results.csv` |
| Paper-ready phenotype summary | `output/tables/20_master_results_paper_ready.csv` |
| Per-library phenotype covariates | `output/tables/31_rnaseq_phenotype_covariates.csv` |
| Raw library lookup | `output/tables/31_rnaseq_library_lookup_raw.csv` |
| Preliminary SNP cluster covariates | `output/tables/32_prelim_snp_rnaseq_covariates.csv` |
| Preliminary SNP join/design audit | `output/tables/32_prelim_snp_join_audit.csv`, `output/tables/32_prelim_snp_response_joinability.csv`, `output/tables/32_prelim_snp_design_balance.csv` |
| RNA-seq background notes | `docs/rnaseq/README.md` |
| Coral expression literature synthesis | `docs/rnaseq/coral_expression_literature_synthesis.md` |
| Candidate gene evidence table | `docs/rnaseq/candidate_genes_reference.csv` |
| Literature search log | `docs/rnaseq/literature_search_log.md` |

Important limits:

- There is no Day 0 expression baseline in the 144-library design.
- Day 15 tissue may exist, but it is not part of the current 144-library analysis.
- The destructive RNA-seq fragments should not be assigned later individual
  regeneration outcomes unless that exact fragment was followed. Link expression
  to phenotype through treatment, day, source thicket, tank, and pre-defined phenotype
  summaries.
- A, C, and D are field source-thicket labels. They are not yet matched to the
  Cunning CBASS genets.
- The preliminary SNP cluster file received from Rachael Bay on 2026-09-02 is
  explicitly provisional. It can be used for exploratory checks and planning,
  but final expression models should use the final SNP set, genetic PCs, or
  kinship matrix once delivered.

## 2. Bias guardrails

Before differential expression:

1. Freeze the phenotype anchors listed in Section 3.
2. Freeze the sample-inclusion and gene-filtering rules.
3. Document the count matrix basis: host genes only, symbiont reads only, or a
   host/symbiont split.
4. Record genome/annotation version, read-processing workflow, and mapping rules.
5. Verify sample IDs against the plate map, raw library lookup, phenotype
   covariate table, and preliminary/final SNP metadata.
6. Check whether any plate, extraction batch, library batch, lane, or well-position
   variable is confounded with temperature, wound state, day, source thicket, or
   tank.
7. Decide before testing whether candidate genes are used only for interpretation
   or also for formal gene-set scoring.

Reporting rule:

- Label tests in advance as **confirmatory**, **secondary**, or **exploratory**.
- Report effect sizes and uncertainty, not only adjusted p-values.
- If the transcriptome contradicts the phenotype story, treat that as a result,
  not as a problem to tune away.

## 3. Locked phenotype anchors

These are the organismal patterns the expression data can confirm, extend, or
revise. They should not be edited after seeing the counts unless the phenotype
pipeline itself changes.

| Anchor | Phenotype result | Expression question |
|---|---|---|
| Regeneration decoupling | Tissue closure proceeds, but new corallite regeneration is suppressed at 31 C | Does heat suppress the skeletal/regenerative program more than early wound closure? |
| Day 10 transition | New corallite formation diverges most clearly by Day 10 | Does the Day 10 wounded-margin expression state differ between 28 C and 31 C? |
| Chronic heat stress | 31 C is chronic/sublethal rather than an acute heat-shock challenge | Do expression patterns show sustained stress/acclimation rather than only an acute HSP spike? |
| Source/genotype-linked resilience | Source C is the most resilient phenotype; preliminary SNPs show C is coherent while A/D are mixed | Does C show a smaller heat-induced expression shift, or a different protective state, after accounting for final genetic structure? |
| Whole-organism physiology | Heat affects PAM, color, symbiont density, and growth | Do expression modules track photochemistry, pigment, symbiont density, growth, or multivariate condition? |

## 4. Competing hypotheses

### H1. Phase-Decoupling Hypothesis

Heat does not shut down all repair. Instead, early tissue closure and wound
response proceed, while the later skeletal/regeneration program is delayed or
suppressed at 31 C.

Predictions:

- Wounded vs unwounded margins differ strongly at Day 1 and/or Day 3 at both
  temperatures.
- Day 10 wounded margins at 28 C show stronger expression of skeletal matrix,
  biomineralization, growth-factor, or corallite-patterning processes than Day 10
  wounded margins at 31 C.
- Early wound, immune, ECM, or proliferation signals are less temperature-sensitive
  than later skeletal/regeneration signals.

Falsifiers:

- The wound response is equally weak at all days.
- Heat suppresses early wound-response genes just as strongly as later
  regeneration genes.
- Day 10 wounded margins do not differ by temperature after accounting for batch
  and tank.

### H2. General-Stress Hypothesis

Heat creates a broad stress state that reduces many physiological and repair
processes, with regeneration appearing more affected only because it is the most
energy-demanding endpoint.

Predictions:

- 31 C samples show broad expression shifts in proteostasis, oxidative stress,
  metabolism, mitochondrial function, and symbiosis-linked pathways in both
  wounded and unwounded margins.
- Wounded and unwounded heat responses share many genes/modules.
- Physiology-linked modules correlate with PAM, color, symbiont density, and
  growth more strongly than with wound state.

Falsifiers:

- Heat effects are mostly restricted to wounded Day 10 samples.
- Unwounded 31 C margins show little sustained expression shift.

### H3. Source/Genotype-Resilience Hypothesis

Source C resists chronic heat stress better than A or D, either through a
smaller heat-induced transcriptome shift or a different protective expression
state. The preliminary SNP file suggests C may also be genetically coherent,
whereas A and D are mixed source labels; the final genetic analysis should
separate source-label effects from genotype/kinship effects.

Predictions:

- C has a smaller 31 C vs 28 C expression displacement than A/D.
- C differs from A/D in stress-response, redox, proteostasis, symbiosis, or
  calcification-associated modules.
- C's module scores correlate with the pre-defined source-resilience summaries.

Falsifiers:

- C is not transcriptionally distinct after accounting for tank, day, wound state,
  plate, and final genetic structure.
- C has the same or larger heat-induced expression displacement than A/D.

Important language limit:

- Without Day 0 libraries, do not call this a true constitutive baseline test.
  We can test sustained state in unwounded margins at Days 1, 3, and 10, not a
  pre-experiment baseline.

### H4. Post-Transcriptional/Cellular Constraint Hypothesis

The phenotype may not map cleanly onto mRNA abundance. Regeneration could fail
because of cellular architecture, resource allocation, symbiont physiology, or
post-transcriptional regulation rather than a simple on/off expression program.

Predictions:

- Strong organismal regeneration differences occur with weak or diffuse host
  gene-expression contrasts.
- Module-level expression tracks general stress but not new corallite formation.
- Symbiont fraction or symbiont expression explains part of the physiological
  pattern.

Falsifiers:

- A coherent Day 10 regeneration module clearly separates 28 C and 31 C wounded
  margins and tracks the organismal regeneration phenotype.

## 5. Primary confirmatory contrasts

These are the contrasts to specify before looking at differential expression.
The exact software can be chosen by the analyst, but the design must account for
plate/batch and the tank-level temperature treatment.

Preferred model family:

- Use a count-aware RNA-seq framework that can accommodate blocking or random
  effects when possible, such as limma-voom with `duplicateCorrelation` or
  variancePartition/dream.
- A DESeq2 fixed-effect model is acceptable only if the design matrix is
  estimable and the tank-level treatment structure is handled transparently.

Minimum biological design terms:

```text
expression ~ temperature * wound * day + source_thicket + plate + tank/blocking
```

For source- or genotype-specific tests, fit targeted models rather than relying only on one
large omnibus model:

```text
expression ~ temperature * source_thicket + day + wound + plate + tank/blocking
expression ~ temperature * wound * source_thicket + day + plate + tank/blocking
expression ~ temperature * final_genetic_PC_or_kinship + day + wound + plate + tank/blocking
```

Primary contrasts:

| Contrast | Purpose | Confirmatory interpretation |
|---|---|---|
| wounded vs unwounded at Day 1 and Day 3 | Positive control for margin wound response | Confirms the sampled tissue captured a wound signal |
| 31 C wounded Day 10 vs 28 C wounded Day 10 | Test heat effect at the regeneration transition | Supports H1 if regeneration/calcification programs are lower at 31 C |
| temperature x wound at Day 10 | Test whether heat changes the wound response specifically | Separates a wound-specific heat effect from general heat stress |
| temperature x day within wounded margins | Test whether heat changes the temporal wound-to-regeneration trajectory | Supports H1 if divergence increases by Day 10 |
| 31 C vs 28 C within unwounded margins | Test chronic heat state without acute wound response | Supports H2/H3 if sustained heat modules appear in unwounded tissue |
| temperature x source thicket | Test whether C differs from A/D in heat response | Supports H3 if C has smaller or qualitatively different heat response; final SNP PCs/kinship determine whether this is source-label or genotype-linked |

## 6. Secondary module and pathway analyses

Run these after the genome-wide tests.

1. Gene-set enrichment or GO enrichment for each primary contrast.
2. A priori module/gene-set scoring for wound response, ECM remodeling,
   proliferation, oxidative stress, proteostasis, metabolism, symbiosis,
   pigmentation, calcification, biomineralization, and corallite patterning.
3. WGCNA or another co-expression module approach, with module eigengenes tested
   against:
   - temperature
   - wound state
   - day
   - source thicket
   - final SNP PCs/kinship or genotype cluster
   - PAM/FvFm
   - color score
   - symbiont density
   - growth
   - source-resilience score
4. Host/symbiont split or host:symbiont fraction analysis if mapping supports it.

Candidate genes from `docs/rnaseq/candidate_genes_reference.csv` should be used
as interpretation aids and for pre-declared gene-set scoring. The table includes
evidence tiers, source systems, LTH use cases, and caveats. Do not use it as a
filter that hides the genome-wide result.

## 7. Exploratory integration

Label these analyses explicitly as exploratory unless they are written into a
pre-analysis addendum before counts are inspected.

- Expression PCA/MDS axes vs organismal PCA axes.
- Module eigengenes vs individual physiology traits.
- Regularized models asking whether expression predicts phenotype summaries.
- Network modules associated with the closed-but-not-regenerated phenotype.
- Final SNP calling from host RNA-seq reads to model genetic relatedness and, if
  possible, match A/C/D-derived samples to Cunning CBASS genets.
- Preliminary SNP clusters/PCs from `data/raw/rnaseq/PRELIM_LTH_genoclusters.csv`
  are useful for exploratory balance checks and same-fragment symbiont-density
  checks, but should not be treated as final genotype calls.
- Parallel symbiont expression or symbiont-fraction analysis, if read mapping
  supports it.

Exploratory results can motivate the Discussion and future work, but they should
not be written as if they were the original test.

## 8. Provenance and QC checklist

Before biological modeling, write a short conventions file for the RNA-seq data.
Minimum checks:

| Gate | Question |
|---|---|
| Sample identity | Do count-matrix columns match library IDs and fragment IDs exactly? |
| Blanks and controls | What do empty wells, NTCs, library blanks, ERCC/spike-ins, and anchors mean in the delivered data? |
| Duplicates | Are duplicate library IDs technical replicates, repeated sequencing, or errors? |
| Dropped records | Which libraries, genes, contigs, or samples were removed before the count matrix was delivered? |
| Timepoints | Do day labels match collection dates and the 144-library design? |
| Batch | Are plate, well, extraction batch, library batch, lane, and sequencing run balanced against the biological factors? |
| Mapping basis | Are counts host-only, symbiont-only, mixed, or separated by reference? |
| Annotation | Which reference genome/transcriptome and gene annotation version were used? |
| Outliers | Are outlier exclusions based on pre-defined QC metrics rather than biological contrast strength? |

Suggested deliverables:

- `docs/rnaseq/rnaseq_data_conventions.md`
- `output/tables/rnaseq_sample_qc.csv`
- `output/tables/rnaseq_sample_exclusions.csv`
- `output/tables/rnaseq_contrast_registry.csv`
- `output/tables/rnaseq_expression_phenotype_links.csv`

## 9. Minimum figures

The first expression paper or analysis handoff should include:

1. Sample QC PCA/MDS colored by temperature, wound, day, source thicket, tank,
   plate, and final genetic structure.
2. Primary contrast summary for wound-response positive control.
3. Day 10 heat effect in wounded margins.
4. Heat x wound x day expression summary, preferably as module scores or an
   interpretable heatmap.
5. Source C vs A/D heat-response summary, with final SNP PCs/kinship added; use
   the 2026-09-02 preliminary SNP clusters only as a flagged exploratory layer.
6. Expression-module to phenotype map using the pre-defined organismal traits.

## 10. Decision log

Fill this in as choices are made.

| Date | Decision | Rationale | Person |
|---|---|---|---|
| 2026-09-02 | Separate confirmatory contrasts from exploratory phenotype integration | Avoids choosing RNA-seq results after seeing the transcriptome | Adrian/Codex draft |
| 2026-09-02 | Candidate-gene evidence table uses tiers A-D | Keeps direct Acropora evidence separate from reviews, web-only sources, and lab-only candidates | Adrian/Codex draft |
| 2026-09-02 | Add Rachael Bay's preliminary SNP cluster file as an exploratory covariate layer | The file joins cleanly to all 144 RNA-seq libraries, but Rachael marked it preliminary and flagged low-coverage singleton concerns | Adrian/Codex draft |
