# Photo and score review: fragments 115 and 121

## Finding

No scores changed. The source records agree, but the available photos do not
resolve the disputed stage observations. Fragment 121 is the stronger scoring
inconsistency; fragment 115 may reflect different thresholds for recognizing
budding versus vertical extension.

## Source-score checks

All 96 values for tip existence, tip extension, and new corallites across the
32 visits for these two fragments agree between the raw CSV, archived Excel
workbook, live Google Sheet snapshot, and processed data. The three traits
have no blanks in these rows and there are no duplicate ID/day keys.
The live sheet has a different row order, so comparisons use ID/day, not row
position. The live sheet has no comments or cell notes in the checked Q:U
cells for these visits. Ordinary notes mention checks/treatments for flatworms
on Days 3, 4, and 8 but do not explain the late-stage scores.

Reproduce with `Rscript code/sensitivity/40b_stage_source_audit.R`.
Live source: `1_6CCmYHGpocI7q9XDzaJHya5NaR-SqArgd24aOSpuoM`, `data` tab,
snapshot retrieved 2026-09-15 and saved in
`data/external/molly_stage_review_live_2026-09-15.csv`.

## NAS scope and image identity

Searched recursively for filenames containing 115 or 121 under:

- `/Volumes/Stier_Lab/Moorea 2025/Moorea 2025/LTH 2025/A. pulchra`
- `/Volumes/Stier_Lab/coral-regeneration/Projects/17. LTH_expression_by_temperature_2025`

The first location contains 14 color-card JPGs per fragment: two views on
30 May and 4, 7, 10, 13, 16, and 18 June 2025. The project archive returned no
additional ID-matching files. This is a search of the relevant project
locations, not a claim that no other photos exist anywhere on the NAS.

Images visually inspected under the first location's `Color Card Photos`:

| Fragment | Day | Files inspected |
|---|---|---|
| 121 | 0 | `06042025/06042025_121a.JPG` |
| 121 | 12 | `06162025/06162025_121a.JPG` |
| 121 | 14 | `06182025/06182025_121a.JPG`, `06182025/06182025_121b.JPG` |
| 115 | 0 | `06042025/06042025_115a.JPG` |
| 115 | 12 | `06162025/06162025_115a.JPG` |
| 115 | 14 | `06182025/06182025_115a.JPG`, `06182025/06182025_115b.JPG` |

The handwritten base labels agree with the filenames. These are whole-fragment
side views. Their angle and image detail do not establish whether a small
terminal feature is a newly budded corallite, a pre-existing lateral corallite,
or enough axial growth to satisfy the extension score. No image was enhanced
or used to assign a replacement score. Day 0 views are baseline comparisons,
not proof of whether photography preceded or followed that day's clipping.

The `Ambient (28C)/D15` microscope folder contains IDs 9-16, not 121.
The microscope protocol explicitly describes a separate 16-fragment experiment
photographed in July, whereas these fragments were scored in June. It is not
valid to substitute those microscope images for the main-experiment fragments.
The nearby `Various photos/Pulchra top view photos` folder uses other specimen
labels and offers no documented match to 115 or 121.

## Interpretation by fragment

- **121, 28 C, tank 9, source C:** tip present only on Day 8; tip extension
  never present; new corallites present on Day 15 despite tip existence being
  marked absent that day. The codebook defines those buds as being along the
  axial corallite, making this combination internally inconsistent. The located
  photo series ends on Day 14, before the disputed positive observation. It
  cannot establish which Day 15 score needs correction.
- **115, 31 C, tank 5, source C:** tip present from Day 8; new corallites on
  Day 14; extension first recorded on Day 15. A tip was therefore present when
  buds were scored. The codebook distinguishes tip presence from vertical
  extension, so the one-day ordering is not by itself proof of a scoring error.
  The Day 14 side views are insufficient for a confident reassignment.

## Does the main observation depend on these fragments?

Original budding counts are 12/12 at 28 C versus 4/12 at 31 C. Omitting both
flagged fragments gives 11/11 versus 3/11. The direction of this descriptive
contrast does not depend on these records. This is not a model refit, a test
of significance, or a recommendation to exclude either fragment.

## Remaining decision

Ask Molly to identify the observation/photo supporting 121's Day 15 budding
score and whether the tip-existence entries after Day 8 need correction. For
115, confirm whether budding can be scored before visible vertical extension
under her scoring rule. Preserve original observations until that review.
No email sent; no shared source data or fitted models changed.
