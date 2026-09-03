# RNA-seq / Preliminary SNP Inputs

Last updated: 2026-09-03

This folder holds local RNA-seq-adjacent inputs received from the Bay lab.

## `PRELIM_LTH_genoclusters.csv`

Source: email attachment from Rachael A. Bay (`rbay@ucdavis.edu`) received
2026-09-02 14:06 PDT, subject `Re: LTH data`.

Status: **preliminary**. Rachael described these as preliminary clustering
results from an incomplete SNP set. She noted that a few individuals are assigned
to their own cluster and that some of those individuals also had lower values in
the `cov` column, so not all cluster assignments should be treated as final.
When the full SNP set arrives, replace or supersede this file rather than
treating it as confirmed genet identity.

Current use:

- Read by `code/32_prelim_snp_phenotype_integration.R`.
- Joined to `output/tables/31_rnaseq_phenotype_covariates.csv` by parsed coral
  fragment ID.
- Used for exploratory summaries of preliminary SNP cluster structure, design
  balance, and same-fragment symbiont density.

Important limits:

- A/C/D remain field source-patch labels, not verified genetic individuals.
- The current SNP file covers the 144 RNA-seq libraries. Those fragment IDs
  match the symbiont-density biopsy data but do not match the PAM, color-card,
  buoyant-weight, gross morphology, or microscope morphology fragment sets.
- Coverage bands created by `code/32` are descriptive only. They are not a Bay
  lab QC threshold.
