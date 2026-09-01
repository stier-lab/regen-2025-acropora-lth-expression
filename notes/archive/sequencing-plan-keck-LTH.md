# **Proposal for Long-Term Heating Experiment Analysis of Heat x Wounding**

# **LTH Main Gene Expression Analysis Plan — Margin-Only (Days 1, 3, 10\)**

> **Archive/current-status note:** this file mixes proposal language with the
> final 144-library design. Use the "as built" statements below for analysis:
> margin-only RNA-seq libraries from D1, D3, and D10; 4 tanks per temperature;
> thicket/genet labels A, C, and D; and no D15 libraries in the current design.
> Any QC threshold still marked provisional must be confirmed with the Bay lab.

Updated plan incorporating four tanks per temperature, each with three thicket/genet labels (A, C, D). This represents the main gene expression analysis for margin-only biopsies collected at three timepoints (Days 1, 3, and 10) after wounding. Day 1 captures the wound response after 24 h.

## **1\. Experimental Factors**

- Temperature: 28 °C and 31 °C
- Wound: wounded (about 1 cm clipped from the growing tip) vs no-wound control. Verify the exact control handling before calling it "sham-handled."
- Timepoints: Day 1 (24 h post-wound margin), Day 3, Day 10. Day 15 samples exist but are not in the current 144-library analysis design.
- Position: margin (M) only for all RNA-seq libraries.
- Thicket/genet labels: A, C, D per tank, with 1 fragment per thicket label per day per wound state.
- Tanks: 4 per temperature, 8 total.

![][image1]

## **2\. Sampling Design and Counts**

At each timepoint (Days 1, 3, 10), both wounded (W) and no-wound control (U) margins are sampled from all four tanks per temperature, with three thicket/genet labels (A, C, D) represented in each tank.

* Tank sets by temperature:

  * 28 °C: tanks 3, 6, 9, 12

  * 31 °C: tanks 4, 5, 10, 11

Counts per timepoint:

* Per temperature × day × wound: 4 tanks × 3 genotypes \= 12 samples

* Per temperature × day (W \+ U): 12 W \+ 12 U \= 24 samples

* Across both temperatures per day: 24 (28 °C) \+ 24 (31 °C) \= 48 samples

Total across all days:

* 3 days × 48 samples/day \= 144 libraries

This design gives balanced replication across temperature, wound status, timepoint, and thicket/genet label, yielding a total of 144 libraries.

In summary: 3 thicket/genet labels (A, C, D) x 2 wound states (U, W) x 4 tanks per temperature x 2 temperatures (28, 31) x 3 days x 1 position (margin) = 144 libraries.

## **3\. Balanced Plating and Sequencing**

Samples will be distributed across two 96-well plates. The as-built sample count is 72 primary coral RNA samples per plate, plus 8 fixed control wells; optional anchor wells should be listed only if they were actually used. Plates are balanced so no single factor (temperature, wound, day, tank, thicket/genet label) appears only on one plate or lane.

### **Plate Details**

Each 96-well plate contains 72 primary coral RNA samples, along with 8 fixed control wells and up to 8 optional technical anchor wells reserved for cross-plate normalization if used. Unused wells should remain labeled as unused in the plate files.

Composition per plate:

* Per day: Each plate includes 24 samples, divided evenly across temperatures and wound treatments:

  * 12 samples from the 28 °C treatment (6 wounded, 6 unwounded)

  * 12 samples from the 31 °C treatment (6 wounded, 6 unwounded)

* Across three sampling days (1, 3, and 10): 24 samples × 3 days \= 72 primary libraries per plate.

* Across both plates: 72 × 2 \= 144 total coral libraries.

The remaining wells per plate are:

Fixed controls (8 per plate):

* 2 RNA extraction no-template controls (NTCs) to detect contamination.

* 2 library preparation blanks (no RNA added) to detect cross-sample carry-over.

* 2 ERCC or spike-in RNA standards to monitor library construction efficiency and sequencing consistency.

* 2 technical replicate samples drawn from the coral RNA pool ("anchor" samples) that will be used for cross-plate and cross-run batch correction, if included.

Optional technical anchors:

If actually included, list their sample IDs and wells. If not included, keep these wells marked unused. These replicates estimate plate or run effects during downstream data analysis.

Plate design rationale:

Distributing each day's 28 °C and 31 °C samples evenly across the two plates ensures that temperature, wound status, and day are not confounded with plate or run. Each plate contains all major treatment combinations. Both plates will be pooled and sequenced together on the same NovaSeq run to reduce sequencing-lane effects. Extraction and library-prep batches should follow the same mixed design, with both temperatures, both wound states, multiple tanks, and all thicket/genet labels represented in each batch.

This layout gives balanced coverage across biological and technical sources of variation, simplifying statistical modeling and reducing the need for correction after sequencing.

### **Extraction, Library Prep, and Sequencing details**

### EZNA Total RNA extraction kit starting with DNA/RNA Shield-preserved samples

### NEBNext Ultra Directional RNA Library prep kit. We’ve used ¼ reactions successfully so that’s probably what I would do here.

### Sequence on a single NovaSeq 25B lane, which outputs 3.1 billion read pairs. Expected reads per sample must be recalculated from the final number of libraries and controls; do not use the 192-sample estimate for the 144-library design.

### 

### Anchor/Control references 

### [Construction of RNA reference materials for improving the quantification of transcriptomic data](https://www.nature.com/articles/s41596-024-01111-x?error=cookies_not_supported&code=4996ef96-dbfb-4a94-b751-67ba0549c0ef#:~:text=for%20constructing%20reference%20datasets,applications%2C%20including%20definition%20of%20performance)

[Multi-omics data integration using ratio-based quantitative profiling with Quartet reference materials](https://www.nature.com/articles/s41587-023-01934-1?error=cookies_not_supported&code=7d0e1e8e-0b22-4879-bf22-e57608ca8981#:~:text=relationships%20among%20the%20family%20members,profiling%20with%20common%20reference%20materials)

### **4\. Factor Balancing**

To ensure that no single experimental factor (e.g., temperature, wound status, tank, or thicket/genet label) aligns with a specific plate, batch, or sequencing lane, all samples should be mixed across plates and processing steps. This prevents technical artifacts such as batch or plate effects from mimicking biological differences.

Each Temp x Day x Wound combination includes 12 samples (4 tanks x 3 thicket/genet labels). These 12 samples are divided evenly across the two 96-well plates: six on Plate 1 and six on Plate 2. This guarantees that every treatment condition is represented on both plates.

Both plates therefore contain equal representation of all experimental factors: temperature, wound state, day, tank, and thicket/genet label. Any minor technical variation can then be modeled without being confounded with the biological design.

Extraction and library preparation batches should follow this same mixed structure. Each extraction batch should include RNA from both temperatures, both wound conditions, multiple tanks, and all thicket/genet labels. This prevents any single batch or processing day from being dominated by one treatment.

To further minimize bias from tank-plate associations, thicket/genet-label assignments are rotated across plates and days. For each temperature, genotypes are assigned to plates using a structured but alternating pattern, ensuring that all tanks and thicket/genet labels appear on both plates across the time series:

* 28 °C tanks: 3, 6, 9, 12

* 31 °C tanks: 4, 5, 10, 11

Example daily mapping:

* Tanks 3 & 9 (28 °C) and 4 & 10 (31 °C): genotypes *A* and *C* → Plate 1; genotype *D* → Plate 2

* Tanks 6 & 12 (28 °C) and 5 & 11 (31 °C): genotype *D* → Plate 1; genotypes *A* and *C* → Plate 2

On the next sampling day, the assignment will be flipped (A/C ↔ D), and permuted again on the third day to prevent any consistent tank-to-plate pairing over time.

Under this layout, each plate on every day contains:

* Equal numbers of samples from 28 °C and 31 °C treatments

* Equal numbers of Wounded (W) and Unwounded (U) fragments

* Representation from all four tanks at each temperature

* All three genotypes (A, C, D) represented at both temperatures

This design keeps extraction, library prep, and sequencing plate from being confounded with biological treatments. Any residual technical variation can be modeled statistically without collinearity, improving power to detect biological effects of temperature, wounding, and day.

### 

### **Controls Per Plate (provisional until confirmed)**

- 2 extraction NTCs
- 2 library blanks/NTCs
- 2 ERCC or spike-in controls
- 2 technical replicates (cross-plate anchors; list exact sample IDs and wells if included)

## **4\. RNA and Library Quality Control**

- Provisional target: RNA input 50-100 ng per sample; RIN >= 7 or DV200 > 70%.
- Provisional target: >= 20 million paired-end reads per library.
- Provisional target: >= 50% of reads assigned to host genes; record host:symbiont ratio for each library.
- Metadata must include plate, well, extraction batch, library batch, lane, tank, thicket/genet label, temperature, wound, day, and QC metrics.

Here is the plate-layout spreadsheet: [LTH\_PlateLayout\_with\_IDs](https://docs.google.com/spreadsheets/d/1lHt-ksUymLimWg1dv1NikETNOX0RyZxqF3yclXXt8P4/edit?usp=sharing). The repo-local CSV exports in `data/raw/plate_layout/` are the analysis inputs.

**5\. Hypotheses**

1\. 	**Early Stress and Immune Activation (Day 1):** Immediately after wounding, we expect coral tissue at the margin to activate genes that protect cells from stress and initiate healing. These include genes that neutralize reactive oxygen species (ROS) and genes that strengthen cell membranes and secrete protective proteins into the wound area. No-wound control corals should show little wound-specific change at this stage. Warm-water (31 °C) corals may show a stronger and longer stress response, suggesting that heat amplifies injury sensitivity.

2\. 	**Temperature Modifies the Healing Program:** At 28 °C, we expect corals to begin transitioning from defense to repair within the first few days. Genes involved in tissue growth, calcium deposition, and extracellular matrix (ECM) rebuilding should increase. At 31 °C, these same repair signals may remain lower or delayed, while stress and protein-damage pathways stay elevated.

3\. 	**Progression from Defense to Repair to Recovery (Day 1 to 3 to 10\)**: We expect healing to unfold in distinct stages. Day 1: shock and defense (injury detection, inflammatory response, ROS control, immune activation, cell migration, and proliferation to re-establish the epithelial barrier). Day 3: structural rebuilding (collagen production and extracellular matrix reformation, including matrix metalloproteases). Day 10: remodeling and strengthening (skeletal organic matrix and calcification proteins). By Day 10, corals at 28 °C may look transcriptionally similar to no-wound controls, indicating recovery of the growing-tip state, while 31 °C corals may still express stress-associated genes.

4\. 	**Temperature x Wound Interaction:** We expect heat to prolong the stress phase and delay the shift toward baseline calcification. Warm-exposed corals may show prolonged inflammatory expression, whereas ambient-exposed corals should activate genes tied to cell proliferation and calcification earlier.

5\. 	**Whole-Coral (Systemic) Effects of Heat:** Even in unwounded fragments, we expect high temperature to reduce expression of genes involved in mitochondrial energy production and increase lipid-metabolism genes—reflecting a general metabolic slowdown. This reallocation of energy likely leaves fewer resources available for regeneration when injury occurs.

6\. 	**Integration with Physiology:** We expect gene clusters related to photosynthesis, pigment production, and energy metabolism to correlate with measured physiology: PAM fluorescence, color, and buoyant weight. Stronger correlations at 28 °C would suggest expression and phenotype change together; weaker ones at 31 °C would suggest heat separates molecular repair from visible recovery.

7\. 	**Thicket/genet differences in recovery strategy:** We expect some thicket/genet labels to ramp up repair genes more quickly and suppress stress genes faster, resulting in faster healing and better physiological recovery. Other thicket/genet labels may sustain prolonged oxidative-stress signatures, showing sensitivity to combined wounding and thermal stress.

 

 

**Table 1\.** This table summarizes representative genes drawn from Han et al. (2025) and related coral regeneration studies. Verify that Han et al. (2025) is in the project bibliography before citing it. Each functional category groups genes with similar biological roles relevant to coral wound healing and thermal response. Use these gene sets for gene set enrichment analyses (GSEA) or module interpretation after the genome-wide analysis, not as a filter before differential expression.

Oxidative Stress & Defense: Genes involved in detoxifying reactive oxygen species (ROS), maintaining redox balance, and protecting cellular proteins from damage during the early stress response.

Extracellular Matrix (ECM) & Remodeling: Genes mediating collagen synthesis, matrix turnover, protease regulation, and tissue adhesion during structural rebuilding of the wound.

Biomineralization & Skeletal Organic Matrix (SOM): Genes associated with secretion of organic matrix proteins, calcium handling, and skeletal deposition that restore coral skeletal integrity.

Metabolism & Energy Reallocation: Genes reflecting shifts in mitochondrial respiration, ATP synthesis, and lipid or carbohydrate metabolism that redistribute energetic resources during healing and under heat stress.

Signaling & Cell Communication: Genes coordinating cell proliferation, inflammatory control, and intercellular signaling to synchronize stress response and tissue regeneration.

Together, these categories provide a functional framework for interpreting coral gene-expression trajectories across wound-healing stages and temperature treatments.

| Category | Example Genes (Marker List) | Biological Role / Interpretation |
| ----- | ----- | ----- |
| **Oxidative Stress & Defense** | *CAT*, *AOS-LOX*, *CYBA*, *PXDN*, *PRDX4*, *GSTK1*, *SODC*, *GPX1*, *GSTM3*, *HSP70*, *HSP90AA1*, *DNAJB1*, *TXN*, *TXNRD1* | Detoxify ROS, maintain redox balance, protect proteins from heat and oxidative damage |
| **Extracellular Matrix (ECM) & Remodeling** | *COL1A2*, *COL4A1*, *THBS2*, *LAMC1*, *DAG1*, *SDC1*, *SDC4*, *MMP24*, *MMP25*, *TIMP1*, *CTSL*, *PZP*, *PSAP*, *PLOD1* | Structural rebuilding, cell adhesion, protease/inhibitor balance, wound closure |
| **Biomineralization & Skeletal Organic Matrix (SOM)** | *GALAXIN*, *SAARP1*, *SAARP2*, *MUCIN*, *CDCP1*, *PRG4*, *CARP4*, *CARP5*, *CA2*, *BMP1*, *CARP7*, *EF-hand Ca²+-binding* | Deposition of new skeleton, matrix secretion, CaCO₃ nucleation |
| **Metabolism & Energy Reallocation** | *ATP5A1*, *ATP5B*, *COX1*, *COX3*, *SDHB*, *NDUFA9*, *ACSL4*, *SCP2*, *CPT1A*, *PDHA1*, *PFKFB3*, *LDHA*, *FABP1*, *PPARA* | ATP synthesis, mitochondrial respiration, lipid utilization, glycolytic adjustment |
| **Signaling & Cell Communication** | *TGFBR1*, *TGFB1*, *FGF2*, *EGFR*, *NOTCH1*, *JUN*, *MAPK3*, *NFKBIA*, *IL1R*, *STAT3*, *CREB3L1* | Coordinate stress response, inflammation, and cell proliferation during repair |

###  

[image1]: img/sequencing-plan-keck-LTH_image1.png
