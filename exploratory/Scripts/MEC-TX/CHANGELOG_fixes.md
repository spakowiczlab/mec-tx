# MEC-TX Fix Integration — Changelog
# Date: March 12, 2026
# All 5 fixes from mectx_fixes.R integrated into production-ready files.

## Files to Deploy on OSC

### NEW FILES (copy to helpers/)
1. `standardise_status.R`  → `BASE_R/helpers/standardise_status.R`
2. `get_focus_cohort.R`    → `BASE_R/helpers/get_focus_cohort.R`
3. `dominant_exclusive.R`  → `BASE_R/helpers/dominant_exclusive.R`

### UPDATED FILES (replace existing)
4. `tx_normalize.R`        → `BASE_R/tx_normalize.R`
5. `tx_pooled_analysis.R`  → `BASE_R/tx_pooled_analysis.R`

---

## What Changed — File by File

### standardise_status.R (Fix 4 — extracted as standalone helper)
- Standalone function, exported
- Auto-detects: numeric 0/1, "dead"/"alive", "deceased"/"living", etc.
- Converts to integer 0=alive, 1=dead
- Informative message on detection + conversion

### get_focus_cohort.R (Fix 1 — new helper)
- Pure data extraction: returns tibble of sample IDs + focus_share
- Supports modes: only, concurrent, dominant
- No plotting dependency — reusable for ad-hoc cohort building
- Used by: Fig 4 IO resistance analysis, any standalone cohort query

### dominant_exclusive.R (Fix 2 — new helper)
- Assigns ONE mutually exclusive regimen label per patient
- Specificity-first hierarchy: 3-type > 2-type > single-agent
- Eliminates overlap between e.g. Chemo[DOMINANT] and Chemo+IO[DOMINANT]
- Used by: tx_normalize() (Fix 3), Fig 3 regimen distribution

### tx_normalize.R (Fixes 3 + 4 integrated)
- NEW PARAMETER: `dominant_regimen_share = 0.20`
- NEW VALIDATION: check dominant_regimen_share in (0, 1]
- FIX 4: calls `standardise_status()` on NSCLC_metadata before join
  - Applied to metadata's status column, not Cluster_surv
  - Only runs if "status" column exists in metadata
- FIX 3: adds `dominant_regimen` column to output
  - Named `dominant_regimen` (not `treatment_group`) to avoid collision
    with the existing per-row `treatment_group` column
  - Builds temp intervals from df_clean (AgeAtTreatmentStart.mod → start_year,
    AgeAtMedStop → end_year, treatment_group → type)
  - Calls dominant_exclusive() internally
  - Patients with no qualifying types get "Ancillary/Supportive only"
- NEW DEPENDENCY: requires dominant_exclusive.R and standardise_status.R

### tx_pooled_analysis.R (Fix 5 integrated)
- 3-stage n_cohort audit trail printed on every run:
  - Stage 1: n_raw (patients passing mode filter)
  - Stage 2: n_timeline (patients in plot segments)
  - Stage 3: n_km (focus-dominant subset) ← reported as n_cohort
- Explanatory note when n_km < n_raw (patients dropped because their
  dominant type wasn't in focus_types)
- Return list now includes `n_raw` alongside `n_cohort`
- Unicode arrows (→, ✗) replaced with ASCII (->, x) for OSC compatibility

---

## Source Order on OSC
```r
source(file.path(BASE_R, "helpers/standardise_status.R"))   # Fix 4
source(file.path(BASE_R, "helpers/dominant_exclusive.R"))    # Fix 2
source(file.path(BASE_R, "helpers/get_focus_cohort.R"))      # Fix 1
source(file.path(BASE_R, "tx_normalize.R"))                  # Fixes 3+4
# ... other sources ...
source(file.path(BASE_R, "tx_pooled_analysis.R"))            # Fix 5
```

## Post-Deploy Checklist
- [ ] Source all helpers before tx_normalize.R
- [ ] Run test_script.R (62 → expect some new passes if tests reference dominant_regimen)
- [ ] Verify: check overlap between chemo_io_dominant and chemo_dominant IDs → should be 0
- [ ] Update test_script.R with new column expectations (dominant_regimen in tx_normalize output)
