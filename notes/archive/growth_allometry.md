# Growth metric: % skeletal mass change (and why areal calcification is not used)

Original decision 2026-06-04 (areal calcification); **reversed in June 2026** after M.
Brzezinski confirmed whole-fragment surface area was never measured. Add the
exact confirmation date and medium before citing the personal communication.
Implemented in `code/05_buoyant_weight.R`; used by the growth models in `12`,
`13`, the PCA in `15`/`16`, and the resilience dashboard in `19`.

> **Status:** growth is reported as **% skeletal mass change** (primary) and
> **specific growth rate, SGR** (robustness). An areal calcification rate
> (mg CaCO₃ cm⁻² d⁻¹) is **not** used — see "Why areal calcification cannot be
> used" below. The heat conclusion points the same direction under both metrics,
> but the primary tank-level growth test remains a trend (p = 0.057).

## The question

How should "growth" of an *Acropora pulchra* fragment be estimated from buoyant
weight? Two options were considered: % mass change (ΔM / M_initial) and areal
calcification rate (mass gain per calcifying surface area).

## Why areal calcification cannot be used

Mass gain per surface area per day = 1000 · ΔM(g) / SA(cm²) / days. The numerator (ΔM) is the
mass gain of the **whole** buoyant-weighed fragment. The only surface area we
have (`wax_clean$sa_curve_cm2`, from wax dipping) is for the **small
sub-fragment** taken for the symbiont slurry — the rest of each fragment was
chopped for transcriptomics, so whole-fragment SA was never measured (M.
Brzezinski, pers. comm.). Dividing whole-fragment mass gain by a sub-fragment's
SA mixes two different objects, so the ratio is not a valid areal rate. Growth is
therefore expressed as % mass change and SGR, which need no surface area.

(The same wax SA *is* valid for symbiont density in `code/06`, where the cell
count and the SA come from the same sub-fragment.)

## Checks for size bias in % mass change

% mass change normalizes by skeletal mass and assumes growth is proportional to
current mass. That is an approximation for a surface-mediated process, but the
bias is weak here because fragments were size-standardized (7–10 g) and mass
gain is independent of size (n = 48, all biopsied terminally at day 15):

| Check | Result | Reading |
|---|---|---|
| slope of `ln(Mf) ~ ln(Mi)` | b = 0.97 | growth approximately proportional over this range |
| ΔM vs initial mass | r = 0.12, p = 0.42 | absolute gain independent of size |
| ΔM vs surface area | r = 0.06 | gain does not track sub-fragment SA across individuals |
| initial mass by treatment | 8.27 vs 8.08 g | balanced — no treatment bias |
| initial mass by genet | c = 8.96 vs a/d 7.7–7.9 g | genet c fragments are larger |

The heat effect is the same under % mass change and SGR (treatment F ≈ 23,
p ≈ 1×10⁻⁵; `output/tables/05b_growth_metric_comparison.csv`).

## Result

Mean growth fell from 6.10% at 28 °C to 4.03% at 31 °C — a 34% reduction; SGR
gives the same (0.39 vs 0.26% d⁻¹, 33% reduction). The inferential test is a
tank-level exact permutation (temperature was randomized at the tank level):
28 °C − 31 °C = 2.07 percentage points, p = 0.057.

Genet × treatment interaction for growth is non-significant under both metrics
(n = 4 per genet × treatment cell — underpowered); genet c lost less growth under
heat (≈−17% vs ≈−39% and ≈−43% for a and d), consistent with its resilience in
PAM/color/symbionts, but not statistically significant for growth.

SGR formula to keep with the method: `100 * (ln(final dry mass) - ln(initial dry mass)) / days`.
