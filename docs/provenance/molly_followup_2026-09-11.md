# Molly follow-up: temperature recovery and stage timing

Updated 2026-09-11. Addresses Molly's 2026-09-10 email, subject
`Re: updated LTH summary HTML`, Gmail message `1a08da4d9cd6fc1d`.

## Source files and conventions

- Recovery folder: https://drive.google.com/drive/folders/1AgT7ktS-6_2p-vaQ77NuxEPvQov1MCrU
- Compiled sheet: https://docs.google.com/spreadsheets/d/1Pxp90_PNhOqG2cCt23RUHTzwNvS9qyuyci51J3TBKc0/edit
- Dated assignments: https://docs.google.com/spreadsheets/d/1xyc3Q-9P1butAde2C4Byz1mPXCfrfSg4UYuWIcwIvv4/edit
- CSV snapshots are in `data/raw/apex/recovered_2026-09-10/`.
  The compiled snapshot preserves columns F:I, the timestamp, probe name,
  probe type, and raw value, for all 32,256 source rows. The assignments
  snapshot preserves populated rows of the Longform tab, columns A:F.
- Four original XML files were copied unchanged from the synced Coral-Regeneration
  `Methods/Apex_temp_download/Moorea Summer 2025 APEX temp data` folder:
  `datalog_05.20.2025_14days.xml`, `datalog_06022025_14days.xml`,
  `datalog_06162025_14days.xml`, and `datalog_06302025_14days.xml`.
  The first is a byte-identical copy of an existing repo export. Source hashes,
  actual date spans, serials, and clock settings are in `08_apex_source_manifest.csv`.
- File names are export labels, not reliable coverage bounds. Read the dates
  inside each record. The new June 2 file spans June 2-16; June 16 spans June
  16-30. The older June 19 files contained only part of June 19.

## Temperature audit

`code/08_apex_temperature.R` now parses raw temperature readings before any
aggregation. The durable key is controller serial + displayed timestamp + probe.
Repeated keys must agree on temperature; a conflicting value stops the run.
The compiled sheet is checked against XML and does not add duplicate readings.

All 32,256 compiled rows match XML values. There are 3,456 repeated compiled
rows, leaving 28,800 unique readings. The original malformed `datalog.xml`
cannot be parsed and is explicitly logged in `08_apex_parse_failures.csv`;
valid exports cover the full study without it. Missing timestamps or numeric
values stop the run. Only named tank-temperature probes are used. Values over
60 are converted from Fahrenheit before aggregation, consistent with the
existing temperature-unit convention. Readings outside 15-40 C are reported
in the excluded-readings table rather than silently removed from daily means.

The dated assignment sheet confirms tanks 3, 6, 9, 12 at 28 C and tanks
4, 5, 10, 11 at 31 C for the main LTH experiment, May 28-June 19.
Other dates include other experiments or periods without main-LTH assignments.
Those dates are not labelled as continued main-experiment treatment exposure.

All eight tanks have 144 unique readings, spaced ten minutes apart, in every
day from Day 0 through Day 15 (June 4-19), and also on Day 16. This is 128
tank-days during Days 0-15. Mean daily tank temperatures are 27.7756 C and
30.8655 C, a 3.0899 C difference. Ranges of daily tank means are 27.3492-28.0988 C
and 30.3796-31.2002 C. Temperature averages exclude Day 16 and later.

Dates retain the logger's displayed clock to match Molly's assignments.
XML records a -7 offset, whereas Mo'orea is -10. UTC in the R timestamp is
only a storage convention for displayed wall time; it is not an inferred
physical time zone. A sensitivity calculation shifting the clock back three
hours changes each treatment mean by less than 0.002 C, with the difference
still 3.09 C. Confirm the controller clock before exact clock-hour exposure
calculations. Later YSI hand-meter records are still unavailable; their absence
does not create a gap in the recovered Apex record.

## Stage-timing audit

`code/36_molly_followup_checks.R` builds one row per wounded fragment and stage,
restricting observations to Days 0-15. `code/39_stage_timing_visual_options.R`
checks duplicate fragment-day keys, counts missing daily scores, and exports
the plotted cumulative counts. A missing score is not recoded as absence.

The first panel uses the separate microscope-photo experiment: 8 wounded
fragments per temperature. The other five panels use the main physiology
experiment's 12 wounded fragments per temperature. They are separate cohorts,
not sequential measurements of the same 16 or 24 individuals. The dated tank
sheet places the microscope-photo experiment later in the season.

Every fragment has final scored follow-up through Day 15. The only missing
daily scores among these six plotted traits are two heated tissue-cover
scores on Day 1. Both were positive on Day 2. Thus the small tissue-cover
first-observation difference cannot establish a biological delay under heat.

| Stage observed at least once by Day 15 | 28 C | 31 C | Timing averages shown? |
|---|---|---|---|
| Living tissue covers wound (separate experiment) | 8/8 | 8/8 | Yes, first observation only; Day 1 missingness noted |
| Main polyp forms | 12/12 | 12/12 | Yes |
| Wound surface remodels | 12/12 | 12/12 | Yes |
| Tip forms | 12/12 | 12/12 | Yes |
| Tip extends | 11/12 | 11/12 | No |
| New skeletal cups bud | 12/12 | 4/12 | No |

The same four of six displayed stages are complete in both temperatures.
One is healing and three are regeneration; this figure is not a display of
six regeneration traits in addition to healing. No stage was dropped to match
the wording of Molly's email.

The denominator stays 8 or 12 as appropriate. Curves show documented first
attainment, not percent currently positive. A fragment stays counted after its
first positive score even if a later score is negative. On Day 15 itself, tip
presence is 11/12 in both treatments, and tip extension is 11/12 versus 9/12.
Those endpoint-presence values previously appeared under a misleading
"reached by" label. The summary now explains the difference explicitly.

Means, medians, and middle-half ranges are suppressed in both treatments if
either treatment is incomplete for that stage. Otherwise an apparently short
average can just reflect the early finishers. Optional violin shapes and
individual dots still describe observed events, with non-events visible.
They are not estimates of the full distribution of completion times.

The existing Weibull timing model in `code/14_morphology_kaplan.R` uses intervals
between negative and positive observations and right-censors unfinished
fragments. It already retains non-events. Its 1.32 time ratio is a model-based
estimate, not the four heated finishers' mean. It adjusts for source patch,
but not tank clustering. The summary foregrounds observed Day 15 attainment
and states that limitation instead of presenting its small P value as decisive.
No survival models were refitted for this presentation/provenance update.

## Reproduction and checks

Run `code/08_apex_temperature.R`, `code/36_molly_followup_checks.R`, and
`code/39_stage_timing_visual_options.R`, then render
`docs/team_summary/LTH_results_summary.Rmd` with `rmarkdown::render()`.
The main-figure temperature panel uses the same corrected assignment window.
The diagnostic inventory and next-analysis roadmap were updated to remove the
resolved late-Apex gap. Previous September 6 searches remain historical records.

The summary HTML was rendered and both primary replacement figures inspected.
The shared Drive HTML was replaced and the companion email sent to Molly on
2026-09-11. The sent-message record is in `docs/correspondence/`.
