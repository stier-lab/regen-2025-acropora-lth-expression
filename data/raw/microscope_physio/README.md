# Microscope Physio Characterization

Raw export of the LTH microscope-only wound-healing scoring sheet.

Source Google Sheet: `Microscope physio characterization`

Drive URL: https://docs.google.com/spreadsheets/d/1rB5VFJqQov0ZXsPjRKnjwVDOtiPxgbfMMDm9l_9Y-Q4/edit

Drive file ID: `1rB5VFJqQov0ZXsPjRKnjwVDOtiPxgbfMMDm9l_9Y-Q4`

Drive modified time at export: `2025-10-07T18:34:33.762Z`

Exported into repo: `2026-09-01`

This is the separate 16-fragment microscope/photo cohort, not the main
`physio_morphology` scoring table. It contains daily observations for wounded
photo-only corals on days 0-15, with two temperatures, two thicket/genet labels
(`a`, `c`), and four tanks.

Important conventions:

- `sample` is `photo` for every row.
- `wound` is `yes` for every row; there is no unwounded microscope control.
- `hole_in_center` and `polyp_in_hole` are identical in the export and are
  combined downstream as `axial_polyp_formation`.
- `pigment_over_wound` is scored through day 7 and then blank from days 8-15.
  The pipeline treats those blanks as not scored, not as `no`.
- `use_photo` marks image series that are useful for illustrative figures; it is
  not a biological response variable.

The corresponding analysis script is `code/11_microscope_physio.R`.
