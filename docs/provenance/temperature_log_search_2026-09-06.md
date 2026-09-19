> **Superseded for Apex coverage on 2026-09-11.** Molly recovered additional files on September 10. They now fill the late-study gap for all eight tanks. The search below records what was available on September 6; see [the recovery audit](molly_followup_2026-09-11.md) for current status. Later YSI spot checks remain unavailable.

# Temperature Log Search, 2026-09-06

Purpose: check whether the missing late Apex logger or YSI hand-meter records are present outside the current processed repo files.

Experimental Day 0 is 2025-06-04. The current coverage audit covers Days 0-16, or 2025-06-04 through 2025-06-20.

## Checked Locations

- Repo raw files: `data/raw/apex/` and `data/raw/ysi/`
- Synced Google Drive project folder:
  `/Users/adrianstier/Library/CloudStorage/GoogleDrive-astier@ucsb.edu/My Drive/Stier Lab/People/Adrian Stier/Projects/In Progress/Coral-Regeneration/Projects/17. LTH_expression_by_temperature_2025/data/APEX temp data and YSI`
- Mounted NAS project folder:
  `/Volumes/Stier_Lab/coral-regeneration/Projects/17. LTH_expression_by_temperature_2025/data/APEX temp data and YSI`
- Mounted NAS field-photo folder:
  `/Volumes/Stier_Lab/Moorea 2025/Moorea 2025/LTH 2025/A. pulchra`
- Google Drive search for `APEX temp data and YSI`, `LTH Apex YSI temperature`, `datalog`, and `ysi daily log`.
- Google Drive-native sheets named `ysi daily log ` and `APEX XML data`.

## Files Found

The repo, synced Google Drive folder, and mounted NAS project folder contain the same Apex XML files:

- `datalog.xml`
- `datalog_05.20.2025_14days.xml`
- `datalog_06.01.2025_14days.xml`
- `datalog_06.12.2025_14days.xml`
- `datalog_06.19.2025_14days.xml`
- `datalog_06.19.2025_7days.xml`

The repo and Drive project folder contain one YSI hand-meter sheet:

- `ysi daily log `

Google Drive search also found an older `APEX XML data` sheet, but its visible cells did not contain usable 2025 Apex log rows. It also found an older YSI sheet that runs through 2025-05-26, before the LTH experimental window used here.

## What The Files Cover

The processed Apex daily file has complete eight-tank coverage for:

- Day 0 through Day 8: 2025-06-04 through 2025-06-12
- Day 15: 2025-06-19

It has no complete Apex coverage for:

- Day 9 through Day 14: 2025-06-13 through 2025-06-18
- Day 16: 2025-06-20

The June 19 Apex XML files are byte-identical to each other. They do not fill the 2025-06-13 through 2025-06-18 gap.

The repo YSI raw files and the active Drive YSI sheet have records through 2025-06-06. Within the Day 0-16 experimental window, that gives full YSI tank coverage for Days 1-2 only. No later YSI hand-meter file was found in the checked project locations.

## Remaining Uncertainty

This search checked the repo, synced Google Drive project folder, mounted NAS project folder, and the documented LTH field-photo folder. A fully unbounded walk across the entire NAS was attempted but was too slow to complete usefully. If late temperature records exist, they are likely outside the checked project folders or require a fresh export from the Apex controller and the later YSI hand-meter sheet/file.

Molly's 2026-09-08 follow-up says she will connect the Apex controller to wifi when she collects the Moorea samples and check whether it stored the missing logger window. Treat that controller check as the next real source for Days 9-14 and Day 16.

## Action

Do not call the temperature record supplement-ready yet. If late Apex or YSI records exist elsewhere, import them and rerun:

- `code/08_apex_temperature.R`
- `code/09_ysi_water_chem.R`
- `code/36_molly_followup_checks.R`
