# QA/QC — flagged samples and tanks

Observations carried over from Molly's original exploratory analysis
(`code/archive/molly_original/`), verified against the processed data, plus the
sensitivity analysis that tests whether they change any conclusion.

Verified 2026-06-04. Sensitivity check: `code/sensitivity/22_sensitivity_flagged.R`.

## Flagged individual colonies

| Coral ID | Treatment | Wound | Genet | Tank | Molly's note | Verified |
|----------|-----------|-------|-------|------|--------------|----------|
| **121** | 28 °C | wounded | c | 9 | "looks weird" in morphology trajectories | present; wounded genet-c |
| **116** | 31 °C | wounded | c | 11 | "looks weird" in morphology trajectories | present; wounded genet-c |

Both are wounded genet-c fragments. Molly flagged them as odd in the
per-individual time-to-event morphology plots (irregular trait-onset ordering).
They are retained in all analyses because the notes do not identify a data-entry
error. The sensitivity check below tests the main treatment conclusions after
dropping them.

## Flagged tanks

| Tank | Treatment | Mean growth (%) | Molly's note | Verified |
|------|-----------|-----------------|--------------|----------|
| **3** | 28 °C | 3.93 | "slow grower compared to other ambient tanks" | **confirmed** — other 28 °C tanks 6.74–7.01% |
| **11** | 31 °C | 5.15 | "negative [wound] effect on growth" | tank mean is actually the *highest* of the 31 °C tanks; the negative pattern Molly saw is the *within-tank* wound contrast, not the tank mean |

Tank 3 is a genuine low outlier among ambient tanks (≈40% lower growth than the
other three 28 °C tanks). Because `(1 | tank)` is in the random-effects
structure for every repeated-measures model, this tank-level variation is
already accounted for; for buoyant weight (single endpoint, OLS) tank 3 is one
of four ambient replicates and pulls the ambient mean down slightly. Keeping tank
3 is conservative for the heat-growth contrast. The sensitivity check refits
growth without tank 3 to confirm the treatment effect holds. Tank 11 is not
dropped because the current audit did not identify it as a low/high tank outlier;
it records a within-tank wound contrast instead.

## Bottom line

These are documented, not excluded. The current audit did not find evidence that
any flag is a data-entry error; they are treated as real variation among corals
or tanks. The sensitivity analysis
(`code/sensitivity/22_sensitivity_flagged.R`, `output/tables/22_sensitivity_flagged.csv`)
re-runs the key tests with the flagged colonies and tank 3 removed and reports
whether the direction and significance category of the main treatment results
change. Add the exact rows/p-values from `22_sensitivity_flagged.csv` here when
this note is used in a manuscript-supporting audit.
