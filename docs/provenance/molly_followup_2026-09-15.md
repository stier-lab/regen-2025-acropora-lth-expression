# Molly follow-up: 15 September 2026

Source: Molly Brzezinski, "Re: updated LTH summary HTML", 15 September 2026,
14:21 Mo'orea time. Gmail message `1a0a797724c913c5`.

## Temperature sheet and clock

Molly believes the live app displayed local Mo'orea time throughout. She
asked whether the complete combined sheet could be used on its own and noted
that export record IDs were not measurement timestamps.

`code/08_apex_temperature.R` now checks this explicitly. Her sheet spans
27 May through 20 June 2025 and contains 28,800 unique readings. All 26,496
original-file readings for the eight study tanks from 28 May through 19 June
are present with matching displayed time, probe, and value. Zero are missing.
The previous check also found all 32,256 sheet rows matched the original files;
3,456 rows were repeats. Export record IDs are not used as measurement keys.

The pipeline retains the XML files as its source and the sheet as an independent
check. Either reproduces the assigned study-period record after deduplication.
Displayed times remain unchanged. The exported -7 label versus local UTC-10
remains a metadata question, not evidence of inconsistent temperature readings.
Molly's recollection is recorded as such, not as independently verified clock
configuration. Existing three-hour sensitivity results remain unchanged.

Output: `output/tables/08_apex_compiled_standalone_check.csv`.

## Average-day view

Molly asked to include the average-day plot for group discussion, rather than
remove it entirely. Section 3 now includes an explicit mean-day companion
alongside the percentage-by-day figure. Every point includes reached/total
counts; open symbols distinguish means conditional on reaching the stage.
Lines show the 25th-75th percentile of observed days, not confidence intervals.
No completion days are assigned to non-reachers. The late heated mean uses
only 4/12 fragments and is not a full-group timing estimate. The existing
individual-observation violin view is retained in the expandable section.

Reproduction: `Rscript code/40_molly_stage_mean_review.R`, after script 36.
Added to `code/_run_all.R`. Outputs: `figures/40_stage_mean_review.{png,pdf}`
and `output/tables/40_stage_mean_review.csv`.

## Stage-sequence review

The check uses the same fragment IDs and first-positive observations as the
existing stage summary. It flags corallite budding before, or without,
recorded tip extension; it does not impose a biological stage order.

- Fragment 121, 28 C, tank 9, source C: new corallites on Day 15; no tip
  extension scored on Days 0-15. Tip present on Day 8, absent Days 9-15.
- Fragment 115, 31 C, tank 5, source C: new corallites on Day 14; first tip
  extension on Day 15.

The existing cleaned-data notes do not explain these stage combinations.
Molly offered to check notes for the ambient fragment. Both records need
photo/field-note review before any correction. Original scores, 4/12 versus
12/12 budding counts, and fitted models are unchanged.

Outputs: `output/tables/40_stage_sequence_review.csv` and
`output/tables/40_stage_sequence_original_scores.csv`.

No email was sent and no shared Drive file was replaced for this revision.
