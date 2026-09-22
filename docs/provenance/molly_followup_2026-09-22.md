# Molly follow-up: September 22, 2026

Source: Molly's September 21, 15:04 Mo'orea-time email, "Re: updated LTH summary HTML", message `1a0c6a5011773f84`, thread `1a0639048ef4f103`.

## Requests and Actions

- Check tip-present then absent sequences and missing pigment calls: audited every wounded fragment in the repository dataset. No raw/processed scores changed and no morphology models refitted.
- Clarify the photographs: separate microscope experiment illustrates stages but cannot confirm individual scores in the main experiment. Updated the summary; the earlier whole-fragment color-photo review remains insufficient to assign replacement scores.
- Compare photosynthesis readings 1 cm from tip and base, including temperature and clipping differences: added reproducible script 42, a paired-data table, tank-level tests, repeated-measures models, diagnostic outputs, and a summary figure.
- Assess averaging: retained the existing two-location mean as a fragment summary, but removed the statement that a nonsignificant location test justifies treating positions as equivalent.

## Morphology Audit

The repository's 24 wounded fragments each have 16 daily records. IDs 121 (28 C, tank 9) and 143 (31 C, tank 4) have tip presence only on Day 8, followed by absence on Days 9-15. ID 116 (31 C, tank 11) has a one-day absence on Day 9 after earlier presence, then presence again. These are review flags, not proof of data-entry error or biological loss. ID 115's budding-before-extension record remains a separate, unresolved scoring-rule question.

Pigment is missing for 23/24 wounded fragments on Day 9: 11/12 ambient and 12/12 heated. ID 172 has the sole recorded value (absence). All other wounded-fragment days have pigment calls. Unwounded-fragment pigment blanks are outside this wounded-only audit and must not be counted as missing wound observations. Missing pigment values remain NA. The six-stage summary does not include pigment, though other morphology analyses do.

Tables: `42_morphology_fragment_review.csv`, `42_morphology_missing_calls.csv`. These audit the repository snapshot; they are not an assertion that Molly has completed her field-note review or that the live sheet has been rescored.

## Photosynthesis Design and Checks

There are 672 valid readings: 48 fragments x seven visits x two positions. Each ID/day has one top and one bottom reading; there are no missing pairs or duplicate position keys. Visits are Days -1, 0, 3, 6, 9, 12, 14. Day -1 precedes clipping, not necessarily heat exposure. Analyses of post-clipping trajectories use Days 0-14.

The raw export contains spreadsheet formulas in the ratio column; recover those values from Y/1000, as the established cleaner does. CSV rows omit trailing fields, mostly notes, and include empty rows. Script 42 saves parser issues, checks all retained measurement fields, and verifies that all 336 reconstructed two-position averages match `pam_clean.rds` to numerical precision. No primary physiology inputs changed.

Define paired gap = near-tip minus near-base Fv/Fm (photosynthesis efficiency). Model the gap within each fragment, avoiding independent treatment of the two locations. The mixed model uses temperature x wound x categorical day plus source patch, with tank and fragment random intercepts. The primary heat comparison instead averages the gap across the six post-clipping visits within each tank and compares the four heated with four ambient tanks across all 70 possible 4:4 label allocations. Its validity assumes exchangeable tank assignments. Confidence intervals are Welch t intervals for the eight tank summaries, not exact permutation intervals.

The whole-trajectory test centers each tank's six-visit gap trajectory on its own average, then compares the sum of squared differences between the temperature-group mean curves across the same 70 tank-label allocations. This isolates differences in time pattern rather than the overall offset. Secondary visit-specific tests have Holm-adjusted p-values, including the pre-clipping comparison as descriptive context.

## Results and Interpretation

Mean post-clipping tip-minus-base gap: -0.01244 at 28 C and -0.02135 at 31 C. The heated-minus-ambient gap is -0.00891 (95% Welch interval -0.02335 to 0.00553; exact p=0.20). Temperature differences in gap trajectories are uncertain (tank-trajectory p=0.14286). Day-14 minus Day-0 gap change differs by -0.02733 between temperatures (p=0.11429).

Both positions show the overall heat response: the heated-minus-ambient difference in Day-0-to-Day-14 change is -0.12429 near the tip and -0.09696 near the base (each exact p=0.02857). These correlated supporting contrasts are not independent replications of the heat treatment. The two-location mean gives -0.110625 (p=0.02857). With eight tanks the permutation p-values are necessarily coarse.

The paired-gap LMM finds little evidence for a clipping main effect (p=0.532), temperature x clipping (p=0.895), or temperature x clipping x day (p=0.406). Its temperature x day term is p=0.0295, unlike the tank-trajectory p=0.143. Retain this discrepancy: the more optimistic fragment-level model is not sufficient to claim that heat changes the location pattern. No equivalence margin was prespecified, so absence of a clear interaction does not prove that positions are interchangeable.

## Diagnostics and Sensitivity

- Paired-gap LMM converged without warnings and was not singular.
- Residual/Q-Q plots were visually reviewed: tails are heavier than Gaussian, with unequal spread, especially under heat. Maximum absolute standardized residual is 3.67. Do not label all assumptions passed.
- Adjacent-visit residual correlation is -0.063; this is a diagnostic summary, not an independence proof.
- A continuous-time correlation model with temperature-specific residual variance gives a mean heat-gap estimate -0.00880 (95% interval -0.02394 to 0.00633; p=0.204), similar to the primary tank comparison.
- Omitting one tank at a time keeps the mean heat-gap estimate negative (-0.0133 to -0.0058). This checks direction, not significance under each omission.
- Using each tank's median gap reduces the estimated difference to -0.0055 (exact p=0.457). The weak heat-gap result is not strengthened by reducing extreme-reading influence.

The summary can retain the two-position mean for overall condition because the heat decline occurs at both sites. The spatial comparison must remain available and cannot be described as showing no location difference. Diagnostics are recorded as partial goodness-of-fit coverage rather than an all-pass result.

## Reproduction and Remaining Work

Run `Rscript code/42_pam_location_and_score_audit.R`, then script 38 and render the team summary. This machine used `/usr/local/bin/Rscript` (R 4.5.2); the Homebrew R installation lacks the project packages. No full pipeline rerun was required because raw data and primary derived physiology/morphology datasets were unchanged.

Still needed from Molly: field-note confirmation for 121/143, the Day 9 observation context for 116 and missing pigment calls, and whether budding can precede visible extension under her scoring criteria. Do not fill gaps or enforce monotonic stage scores without that evidence.

After user authorization on September 22, the shared Drive HTML was replaced in place (file `18z9zDpKl690ZjBZrhjjIrdtDvllNllPK`). Its returned MD5 `b1895a53e7ad37df7cdaf39cdbcfbbc9` matches the local HTML (8,840,470 bytes). The follow-up was sent to Molly in the existing thread, Gmail message `1a0c9ffd0abee1fd`, thread `1a0639048ef4f103`. The correspondence file records the sent body. Plain-text and HTML alternatives and RFC reply headers were verified before sending.
