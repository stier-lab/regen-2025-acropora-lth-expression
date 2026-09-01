# **LTH\_expression\_by\_temperature\_2025 \- README**

> **Archive status:** this is the original Drive README. Use the repo root
> `README.md` and `RESULTS.docx` for the current analysis. This file is kept for
> provenance and has been edited only to remove contradictions that would mislead
> a reader.

## **Project Information**

Title: LTH: expression\_by\_temperature\_2025

Principal Investigator: Adrian Stier (UCSB), Ashley Seifert (UKY), Craig Osenberg (UGA), Rachael Bay (UC Davis)

Collaborators: Molly Brzezinski (UCSB), Michelle Diminuco (UGA)

Project Dates: May \- July 2025

Contact: Adrian Stier \- astier@ucsb.edu

## **Research Objectives**

This project investigates how corals recover from a tip-removal wound under heat. Specifically, we compare wounded and no-wound *A. pulchra* fragments at ambient (28 C) and elevated (31 C) temperature over about 2 weeks. The study aims to:

• Assess differences in gene expression during wound healing and skeletal regeneration stages along a single coral fragment

• Determine the effect of wound-by-temperature on physiology measures (symbiont density, buoyant weight, PAM, color; chlorophyll-a was planned but not run in the current repo)

• Characterize visible and microscope-scale wound-recovery stages

## **Data Sources & Collection Methods**

### **Experimental Design**

• Coral Samples: 

  \- 208 total *A. pulchra* fragments (about 6 cm tall with no obvious branches) from 3 source thickets were used across physiology, gene expression, and microscope photography

* 192 *A. pulchra* fragments formed the main destructive sampling set; the current RNA-seq design selects 144 margin libraries from D1, D3, and D10.
  * 48 *A. pulchra* fragments from 3 thickets (field labels A, C, D) were tracked for non-destructive physiology.
* 16 *A. pulchra* fragments from thickets A and C were used for daily microscope photography only.

• Treatment Levels:

  \- Wounded (yes/no) \- About 1 cm was clipped off the growing branch tip of the nubbin. Verify exact no-wound handling before calling controls "sham."

  \- Temperature (28/31 C) \- Temperature was raised 1 C/day until the target treatment temperature was reached. The wound treatment was applied after corals were held at temperature treatments for 7 days.

• Standardization:

  \- A caliper was used in conjunction with a band-saw to standardize the wound depth

 \- Aqualogic/Apex system was used to control and monitor temperature continuously. Verify which system controlled temperature and which logged it before writing the final methods.

### **Response Variables and Protocols** 

| Response Variable | Method | Destructive? | Notes |
| :---- | :---- | :---- | :---- |
| Growth | Buoyant weight | no | 15-day percent skeletal mass change |
| Photochemical efficiency | PAM | no | PAM settings: sat int (6), damp (2), out gain (6), measure int (6) for all data  |
| Pigmentation | Color card | no | Color card from Siebeck et al. 2006: https://coralwatch.org/wp-content/uploads/2023/11/Siebeck-et-al-2006-Monitoring-coral-bleaching.pdf |
| Morphological characterization | Visual assessment (gross and microscopic) | no | Microscope photos taken at 1.25x and 3.2x zoom, BF only Characteristics \- Coenosarc coverage over wound bed  \- Hole in center of wound \- Polyp in center hole \- Wound smoothed over \- Obvious pigment over wound \- Tip exists \- Tip extension \- New corallites on tip |
| Chlorophyll concentration | \-Biopsy by band-saw  \- airbrushing \- wax dipping \- chlorophyll content | yes | Planned but not run in the current repo |
| Symbiont density | \- Biopsy by band-saw \- airbrushing \- wax dipping \- symbiont count | yes | NA |
| Gene expression | \- Biopsy by band-saw \- RNA extraction | yes | Current 144-library design uses Day 1, Day 3, and Day 10 margin samples; verify whether any Day 0 or Day 15 tissues will be sequenced |

### **Sampling and Final Biopsy count**

• Main destructive sampling set: 192 fragments. Current RNA-seq design: 144 selected margin libraries from D1, D3, and D10. Non-destructive physiology subset: 48 fragments. Photo-only set: 16 fragments.

## **Notes on Execution**

• Project Discussion: Adrian, Ashley, Craig, Rachael, Molly, and Michelle were involved in project planning and experimental design

• Experiment Execution: Adrian, Ashley, Craig, Molly, and Michelle collected corals. Molly and Michelle maintained all aquaria and performed all PAM, buoyant weight, color card, biopsies, airbrushing, and microscope photography

• Data management: Molly organized all raw data for physiological and morphological parameters  and performed QA/QC

• Analysis: Molly analyzed all microscope photos for physiological characterization

## **Key Results & Current Status**

• The dataset includes time-series wound recovery metrics, buoyant weight, color card, PAM, and morphological characterization for a subset of corals. See `RESULTS.docx` and `data/DATA_DICTIONARY.docx` for the current analysis counts.

• The dataset includes symbiont density. Chlorophyll-a was planned but not run. Transcriptomic data are pending for the selected RNA-seq libraries.

• Gene expression samples will be processed in Rachael Bay's lab at UC Davis.

## **Contact for Collaboration**

For access to raw data or collaboration, please contact astier@ucsb.edu.

Last Updated: 14 Oct 2025
