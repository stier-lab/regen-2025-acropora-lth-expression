> **Archive status:** this is an early field plan, not the current analysis
> authority. Use `README.md`, `RESULTS.docx`, `data/DATA_DICTIONARY.docx`, and
> `docs/rnaseq/README.md` for the current phenotype results and the current
> 144-library RNA-seq design. Items below that mention chlorophyll-a, microbiome,
> symbiont genotyping, D15 RNA-seq, or 32 C heat are planning notes unless later
> files confirm they happened.

Outstanding [Project Tasks](https://docs.google.com/spreadsheets/u/0/d/171cWktNYW2Ci-U9L5AbjZgwh_QRAMHJB2N6CzWsbgJ0/edit):

1. Laboratory analysis on slurry from D1, D3, D10, D15
   1. Chlorophyll
   2. AFDW
   3. Zooxanthellae
2. Laboratory analysis on RNA cryotube samples
   1. RNA extractions
   2. RNAseq

### **1\. Objective**

To determine how wounding and elevated temperature (up to 31 C in the current repo) affect physiological and molecular responses in *Acropora pulchra*, including:

* Near-wound vs. far-from-wound gene expression
* Symbiodinium density
* Photosynthetic efficiency (Fv/Fm)
* Skeletal growth (buoyant weight)
* Differences among thickets/genet labels A, C, and D
* Coral microbiome composition
* Symbiont genotyping

### **2\. Experimental Design (15-Day Experiment with Day 1, Day 3, Day 10, and Day 15 Sampling)**

* **Design:** 2 × 2 factorial
* **Factors: 1\. Wounding:** no-wound control vs. wounded (about 1 cm clipped from the growing tip; verify exact no-wound handling before calling controls "sham")
  2\. **Temperature:** Ambient (28 C) vs. Elevated (31 C)
  3\. **Thicket/genet label:** 3 parent thickets, evenly represented across treatments
* **Sampling**: Time points: 24 h, Day 3, Day 10, Day 15
* **Photography**: See “[microscope photograph experimental plan](https://docs.google.com/document/d/1c0i4uat4QWKnpx4hts-h5JiyA4ry10ZulxelFwECvuM/edit?usp=sharing)”
* **Fv/Fm:** PAM measurements starting at baseline (pre-wounding), targeting corals destined for final physiological sampling. Two PAM spots per colony: one at the tip and one at the base.
  1. We will monitor colonies designated for final physiological sampling (buoyant weight, gene expression, symbiont density) using PAM fluorometry (Fv/Fm) starting prior to wounding and continuing weekly throughout the experiment. For each colony, we will take two measurements at each timepoint: one on the side of the colony (intact tissue) and one at the tip or wound site. This approach allows us to track localized versus systemic changes in photosynthetic efficiency. 

#### 

#### 

#### **Treatment Levels and Tank Setup**

* **Treatment levels:**

  * Ambient \+ Control  
  * Ambient \+ Wounded  
  * Heated \+ Control  
  * Heated \+ Wounded

* **Tank Design:**

  * 8 tanks: 4 ambient tanks; 4 heated tanks; 3 thickets per tank.  
  * **Within each tank:**  
    * For each wound treatment in each tank, assign replicates from 3 thickets.

  * Across all tanks, assignments are rotated so that when it comes time for time-series sampling (sacrificed colonies) and for final-only sampling, every treatment combination will have each of the three thickets represented at each timepoint.

### 

### **3\. Total Number of Coral Fragments Needed**

* Each thicket/genet label, wounded and unwounded, is sampled at each time point from a tank assigned to the correct temperature treatment. Verify whether tank choice was random or haphazard before describing it in methods.
* Time series: 8 tanks × 2 wound states × 4 timepoints × 3 thickets = 192 fragments, 24 per tank.
* Current RNA-seq design: 144 selected margin libraries = 8 tanks × 2 wound states × 3 RNA-seq days (D1, D3, D10) × 3 thickets. D15 tissue exists in this plan but is not part of the current 144-library analysis design.

| Factor | Levels |
| ----- | ----- |
| Temperatures | 2 (Ambient, Heated) |
| Treatments | 2 (Control, Wounded) |
| Timepoints | 4 collected (Day 1, 3, 10, 15); 3 in current RNA-seq design (Day 1, 3, 10) |
| Thickets | 3  |

### 

### **4\. Biopsy Sampling:**

* Each colony is sampled in 5 segments. Three segments were planned for molecular analysis (RNA-seq, 16S rRNA microbiome, and host/symbiont genotyping), and two segments were planned for chlorophyll-a and zooxanthellae. Current repo status: symbiont counts were analyzed; chlorophyll-a was planned but not run. Verify microbiome and symbiont-genotyping status before citing them as completed assays.

![][image1]

### 

### **5\. Tissue Biopsy Plan**

| Variable | Method |
| ----- | ----- |
| **Gene expression** | RNA-seq; current 144-library design uses Day 1, Day 3, and Day 10 margin samples |
| **Symbiont density** | Hemocytometer counts |
| **Photosynthetic efficiency** | Weekly PAM fluorometry (Fv/Fm) |
| **Skeletal growth** | Buoyant weight (Day 0, Day 15\) |
| **Host genotyping** | DNA from tissue punches; planned, verify current status |
| **Microbiome** | 16S rRNA amplicon sequencing; planned, verify current status |
| **Symbiont ID** | ITS/qPCR panel for clade or strain identification; planned, verify current status |

### 

### **6\. Timeline**

| Dates | Phase | Details |
| ----- | ----- | ----- |
| **May 18–19** | **Collect & Glue** | Collect coral fragments and glue to plugs. Early prep before mounting. Treat for AEFs with Levamisole |
| **May 19–26** | **Extended Acclimation @ 28°C** | Corals acclimate at ambient temperature (28°C) for 8 days. |
| **May 27–June 2** | **Heat ramp** | Increase temp ~1 C/day in heated tanks. Current repo treats 31 C as the elevated setpoint; verify this row if it was originally written as 32 C. Ambient tanks held at 28 C. |
| **May 30** | **Initial PAM \+ Color card** |  |
| **June 3 (Tue)** | **Pre-Wound PAM (Baseline)** | PAM on **final-only colonies**: one spot on side, one on tip/wound site. |
| **June 4 (Wed)** | **Day 0: Wounding \+ BW \+ Post-Wound PAM \+ color card** | Wound treatments applied and tip collected; buoyant weight measured on final-only colonies. PAM repeated immediately post-wounding (same colonies, side and tip).  Color card photos and tripod photos of final-only (physio corals) |
| **June 5 (Thu)** | **Day 1 Sampling** | 24 hr biopsy: gene expression, microbiome, symbiont ID, host DNA. |
| **June 7 (Sat)** | **Day 3 Sampling \+ PAM \+ color card** | Biopsy (if coenosarc coverage restored). **PAM on final-only colonies**. |
| **June 10** | **PAM \+ color card** |  |
| **June 13** | **PAM \+ color card** |  |
| **June 14** | **Day 10 Sampling** | biopsy: gene expression, microbiome, symbiont ID, host DNA. |
| **June 16** | **PAM \+ color card** |  |
| **June 18** | **PAM \+ color card** |  |
| **June 19**  | **Day 15 Sampling \+ BW** | biopsy: gene expression, microbiome, symbiont ID, host DNA. |

**Personnel timeline:**

**Adrian: May 10- May 27, June 16-28**

* May 9 flight arrives at night, pick up in morning of May 10 (9am ferry)

**Ashley: May 9 \- May 24**

* Flight arrives morning of May 9; Molly pick up morning of May 9.

**Molly: May 1- June 20, July 2 (?) \- July 27**

**Craig: May 10 \- May 27**

* Flight arrives 5:30am, take 9am ferry

**Michelle: May 10 \- July 27**

* Flight arrives 5:30am, take 9am ferry

**Hayden Vega: June 18 \- August 1**

[image1]: img/Experimental_Plan_Gene_Expression_image1.png
