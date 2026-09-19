# Next Analysis Roadmap

Last checked: 2026-09-09

This roadmap organizes the analyses that would most improve the project after the
current phenotype summary. It is a planning document, not a results document.
The phenotype pipeline is current; RNA sequencing (RNA-seq) expression counts
are not integrated yet. Rachael Bay's 2026-09-02 single-nucleotide polymorphism
(SNP) cluster file is preliminary and should be replaced by final SNP, genetic
similarity, or genetic relatedness (kinship) results before final genetic
wording.

Molly's 2026-09-03 forwarded email adds Trinity Conn's preliminary Hauru colony table from
the Ross Cunning / Trinity Conn work. That file is now stored as external
context only. The coordinate screen makes `Apul-115` the first candidate to
compare against source patch C, but exact identity requires comparing LTH
RNA-seq SNPs to Trinity/Ross whole-genome genotypes. Molly also noted that
Hauru clones can span different thickets and large distances, and that Trinity
mentioned a chromosomal inversion in four Hauru individuals. That makes
coordinate matching useful for prioritizing candidates, not for assigning
identity.

Molly's 2026-09-04 and 2026-09-08 follow-ups add three near-term gates before
these results move from coauthor-summary language into supplement or manuscript
language: audit temperature-log coverage, keep the cohort
language straight, and decide on any extra Day 15 genetics only after final
genetic identity and current RNA-seq expression quality control are stable. The
cohort language is now: main physiology cohort = 48 fragments total; wounded
stage-timing subset = 12 fragments per temperature; microscope-photo cohort =
16 all-wounded fragments for coenosarc healing only, with 14 of 16 scored on Day
1.

## Decision Order

The next analyses should happen in this order:

| Step | Analysis | Decision it answers |
|---|---|---|
| 0 | Coauthor-summary cleanup | Are the temperature logs complete, and are the stage-timing cohorts labeled correctly? |
| 1 | Final sample identity and genetics | Can A, C, and D be treated as genetic individuals, or only as source patches? |
| 2 | RNA-seq quality control | Which samples, genes, batches, and host/symbiont reads are safe to analyze? |
| 3 | Genome-wide expression contrasts | What changed molecularly with heat, wounding, day, and source patch? |
| 4 | Expression-to-phenotype links | Do gene-expression patterns explain heat-stalled regeneration, source patch C, or symbiont retention? |
| 5 | Supporting checks | Do Day 15 genetics, multistate healing, environment, and external heat-limit comparisons sharpen the story enough to justify the added work? |

This order keeps the project from using the strongest phenotype result to
cherry-pick the transcriptome. It also prevents source-patch labels A, C, and D
from being treated as confirmed genetic individuals before the final genetic
analysis is done.

## Step 0: Close Coauthor-Summary Cleanup

**Status, 2026-09-11:** the recovered Apex exports now cover all eight tanks through Day 15. The summary uses cumulative stage attainment and omits averages for stages not reached by everyone at both temperatures. See [the recovery and staging audit](provenance/molly_followup_2026-09-11.md). Exact clock-hour exposure calculations still need the logger time setting confirmed; later YSI spot checks remain unavailable.

**Why it comes first:** these are not new biological analyses. They are
provenance and wording checks that determine how confidently the current
summary can be shared as a coauthor-facing interpretation.

**Outputs to make:**

- Completed: recovered Apex exports, XML-to-compiled-sheet agreement checks, and daily tank coverage. YSI coverage is tracked separately.
- A final microscope-photo statement for coenosarc cover: 16 all-wounded
  fragments total, with 14 of 16 scored on Day 1.
- A final stage vocabulary: early healing = coenosarc cover; regeneration =
  later rebuilding of the polyp, wound surface, tip, tip extension, and new
  radial corallites.
- Completed: cumulative percentage-reached plots retain all fragments. Mean days and middle-half ranges appear only for stages reached by everyone at both temperatures. The violin remains an optional descriptive view.
- A Day 15 genetics decision point: genotype-only endpoint samples if the
  immediate question is cluster identity; sequence extra endpoint RNA only if
  endpoint gene expression is needed.

**Where it fits:** `code/36_molly_followup_checks.R` now writes the temperature
coverage audit, the healing/regeneration stage-timing tables, and the companion
figures used in the team summary.

## Step 1: Lock Identity Before Mechanism

**Analysis:** final sample identity, SNP, and genetic-relatedness audit.

**Why it comes first:** the source-patch result is currently one of the clearest
phenotype patterns, but the preliminary SNP file says source patch C is coherent
while source patches A and D are mixed. Final wording depends on whether that
pattern holds in the complete genetic analysis.

**Outputs to make:**

- A join audit showing which genetic rows match which fragment IDs.
- A final genetic-similarity or kinship table.
- A candidate-match table that compares LTH RNA-seq SNP calls with Trinity/Ross
  whole-genome genotypes, if those reference data become available.
- A flag for whether inversion structure or long-distance clone structure
  changes the interpretation of the preliminary clusters.
- Balance plots showing how final genetic groups fall across temperature, wound
  state, day, tank, and source patch.
- A short wording decision: "source patch," "genetic individual," or "genetic
  group."

**Where it fits:** `code/32_prelim_snp_phenotype_integration.R` is the current
preliminary layer. Extend it or replace it with a final-genetics script when the
final file arrives.

## Step 2: Make The RNA-seq Data Trustworthy

**Analysis:** RNA-seq quality control and count-basis decision.

**Why it comes next:** expression results are only interpretable after sample
identity, mapping, batch, and host-coral versus algal-symbiont read handling are
documented.

**Outputs to make:**

- Sample quality-control table.
- Mapping summary and host:symbiont read fraction.
- Batch-balance checks for plate, extraction batch, library batch, lane, well
  position, tank, day, wound state, temperature, and source patch.
- A count-basis decision: host-only expression, symbiont-only expression, or a
  paired host/symbiont analysis.

**Where it fits:** future `code/rnaseq/01_*` scripts. Background is in
`docs/rnaseq/README.md`.

## Step 3: Run The Main Expression Tests Genome-Wide

**Analysis:** genome-wide expression contrasts.

**Why it matters:** this is the central missing result. The phenotype analysis
shows what the corals did; RNA-seq can show which biological programs changed.

**Primary contrasts:**

- Heat response across all margin samples.
- Wound response at Days 1, 3, and 10.
- Heat-by-wound-by-day response, with Day 10 wounded margins as the key
  regeneration contrast.
- Source patch C versus source patches A and D under heat, after accounting for
  final genetic structure when possible.

**Outputs to make:**

- Gene-expression contrast tables and diagnostic plots.
- A ranked list of biological processes from the genome-wide results.
- A short "phenotype match" table saying whether each expression result supports,
  revises, or fails to explain the phenotype result.

**Where it fits:** future `code/rnaseq/02_*` scripts. Contrast logic is in
`docs/rnaseq/expression_integration_analysis_plan.md`.

## Step 4: Link Expression Back To Phenotype Carefully

**Analysis:** expression modules and phenotype links. Modules are groups of genes
that rise and fall together.

**Why it matters:** individual genes will be noisy. Coordinated gene groups can
ask whether heat stress, wound response, symbiont loss, source-patch resilience,
or regeneration failure share the same signal.

**Good links right now:**

- Symbiont density by fragment ID for the 144 RNA-seq libraries.
- Pre-defined source-patch heat-penalty scores from
  `output/tables/19_genet_resilience_summary.csv`.
- Treatment/day/wound summaries for healing and regeneration.
- Whole-coral condition summaries by source patch and treatment.

**Links to avoid for now:**

- Assigning Day 15 regeneration outcomes to destructive Day 1, 3, or 10 RNA-seq
  fragments unless that exact fragment was followed.
- Calling A, C, and D "genets" in reader-facing text before final genetic
  identity is resolved.
- Selecting a few candidate genes before the genome-wide expression analysis.

**Where it fits:** future `code/rnaseq/03_*` and `code/rnaseq/04_*` scripts. The
per-library phenotype covariates already exist at
`output/tables/31_rnaseq_phenotype_covariates.csv`.

## Step 5: Add Supporting Analyses Where They Clarify

These analyses are useful, but they should not distract from final genetics and
RNA-seq expression.

### Multistate Healing-To-Regeneration Model

The current result shows closure and regeneration are separate. A multistate
model would estimate movement among stages instead of testing each endpoint on
its own:

```text
open wound -> sealed wound -> rebuilt tip -> new corallites
```

Best output: a transition diagram and treatment/source-patch transition rates.
Likely home: future `code/36_multistate_regeneration.R`.

### Host/Symbiont Split

The phenotype story includes symbiont loss, and source patch C retains symbionts
better. RNA-seq can ask whether that is visible in host:symbiont read fraction,
host expression, symbiont expression, or all three.

Best output: host:symbiont read fractions and symbiont-expression summaries
tested against measured symbiont density.
Likely home: future `code/rnaseq/05_host_symbiont_split.R`.

### Extra Day 15 Genetics

Extra Day 15 genetics should answer a direct linkage question, not just add
another data layer. Molly clarified that the Day 15 corals were tracked through
photosynthesis score, color, skeletal growth, symbiont density, and all
regeneration stages. If those same fragments can be genotyped, genotype-only
data may be enough to test cluster identity against endpoint performance.
Extra Day 15 RNA sequencing is the broader option and is useful only if endpoint
gene expression is worth the added cost and analysis.

Best output: a sample-join table showing exactly which Day 15 samples connect
to endpoint phenotypes, followed by a genotype-only versus RNA-seq
recommendation.
Likely home: extend `code/36_molly_followup_checks.R` for the decision audit,
then add a future RNA-seq script only if sequencing proceeds.

### Tank And Environment Sensitivity

The treatment worked, but reviewer-facing checks need a clean test that tank,
temperature logger, Apex/YSI records, and known husbandry issues do not explain
the biology.

Best output: a sensitivity table showing whether the main effects survive
tank/environment checks.
Likely home: existing `code/sensitivity/`, with targeted additions only if
needed.

### Cunning Acute-Heat Genotype Match

If reference genotypes are available, this could connect acute heat tolerance
from Cunning's coral bleaching automated stress system (CBASS) assay to chronic
heat-plus-wound resilience here. ED50 is the temperature where photosynthesis
efficiency drops halfway.

The current preliminary Trinity/Ross coordinate screen is useful but not enough. Source
patch C is closest to `Apul-115` (10.5 m; ED50 = 36.380 C; clonal group 2;
dominant symbiont A). Source patches A and D are both closest to `Apul-111`, and
several nearby colonies are in clonal group 2. If C eventually matches
`Apul-115`, then its strong chronic LTH resilience would pair with an acute
ED50 that is within the Hauru range but not the highest in the table. That would
make the acute-versus-chronic comparison especially useful.

Best output: a match table and a test of acute CBASS ED50 versus chronic LTH
source-patch or genetic ranking.
Likely home: extend `code/sensitivity/26b_trinity_hauru_context.R` after final
SNP files and Trinity/Ross or Cunning reference genotypes are available.

## Hypotheses To Carry Forward

These hypotheses should be tested against each other rather than treated as
separate stories that can all be true by default.

| Hypothesis | Plain-language mechanism | What would support it | What would weaken it |
|---|---|---|---|
| Phase-Decoupling Hypothesis | Heat allows wounds to seal but blocks the later program that rebuilds a new skeletal tip. | Wounded margins respond at Days 1-3 at both temperatures, but Day 10 wounded margins at 31 C lack skeletal or corallite-building expression seen at 28 C. | Heat suppresses early wound response just as strongly as later regeneration, or Day 10 wounded margins do not differ by temperature. |
| General-Stress Hypothesis | Heat creates a broad stress state, and regeneration fails because it is the most demanding endpoint. | Heat shifts wounded and unwounded margins in similar stress, metabolism, redox, and symbiosis-linked modules. | Heat effects are mostly restricted to wounded Day 10 samples. |
| Source-Patch/Genetic-Resilience Hypothesis | Source patch C is less heat sensitive because it has a different protective state or a smaller response to chronic heat. | C shows a smaller heat-induced expression shift, or protective modules that track the pre-defined C resilience score, after accounting for final genetic structure. | C is not distinct after tank, day, wound state, plate, and final genetic structure are included. |
| Cross-Method Concordance Hypothesis | Acute heat limits and chronic heat-plus-wound resilience rank coral backgrounds similarly. | Final SNP matching places C, or C's final genetic group, with higher acute ED50 and smaller LTH heat penalty. | C matches `Apul-115` or another colony with ordinary ED50 but still has strong chronic resilience, or no reliable match is possible. |
| Symbiont-Mediation Hypothesis | Part of C's advantage comes through holding onto algal symbionts. | Host:symbiont read fractions, symbiont expression, and measured symbiont density all show smaller heat penalties in C or its final genetic group. | C's host expression remains distinct but symbiont read fraction and symbiont density do not differ. |
| Technical-Structure Alternative | Some apparent biology is actually sample structure, batch, tank, or source-patch imbalance. | Main patterns weaken after sample identity, tank, batch, coverage, and environment checks. | Main patterns remain after those checks and the design-balance audits look clean. |

## Practical Implementation Order

1. Finish the coauthor-summary cleanup from Molly's 2026-09-04 and 2026-09-08
   notes: retain the recovered Apex coverage audit, keep the stage
   cohorts straight, and keep the Day 15 genotype-only versus RNA-seq question
   as a decision rather than a default next step.
2. Add final genetic files under `data/raw/rnaseq/` with a README note describing
   where they came from and whether they replace `PRELIM_LTH_genoclusters.csv`.
3. Extend the SNP integration script or create a new final-genetics script that
   writes an audit trail to `output/tables/` and figures to `figures/`.
4. Add RNA-seq quality-control scripts before any biological expression tests.
5. Run the primary expression contrasts from
   `docs/rnaseq/expression_integration_analysis_plan.md`.
6. Run module or pathway summaries only after the genome-wide contrasts are
   documented.
7. Add the multistate healing model if the regeneration section needs a cleaner
   single figure or if reviewers ask how the stages connect.
8. Compare LTH SNPs to Trinity/Ross whole-genome genotypes if possible; use
   coordinates only to prioritize candidates (`Apul-115`, `Apul-111`, clonal
   group 2). Add Cunning CBASS matching only after final genotypes and reference
   genotypes are in hand.
