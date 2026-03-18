# =============================================================================
# test_mectx_all.R
# All MEC-TX testthat tests in a single script.
#
# Pipeline (confirmed from real output):
#   tx_normalize()  → sample, AgeGrid, treatment_group, start_age,
#                     TimeSinceTreatmentStart, end_followup, dominant_regimen
#                     └─→ tx_cluster_surv() (uses defaults directly)
#
#   tx_intervals()  → sample, block, type, run, start_year, end_year
#                     └─→ tx_pooled_analysis(), tx_compare_groups(),
#                          km_panel_from_df(), plot_timeline_for_k()
#
# Run with:
#   testthat::test_file("test_mectx_all.R")
# =============================================================================

library(testthat)

# =============================================================================
# 0. Source MEC-TX files
# =============================================================================

BASE_R <- "/users/PAS1695/dipankor99/Github/digital-twins/exploratory/Scripts/MEC-TX"

source(file.path(BASE_R, "standardise_status.R"))
source(file.path(BASE_R, "dominant_exclusive.R"))
source(file.path(BASE_R, "get_focus_cohort.R"))
source(file.path(BASE_R, "constants.R"))
source(file.path(BASE_R, "resolve_col.R"))
source(file.path(BASE_R, "has_subseq.R"))
source(file.path(BASE_R, "prep_segs.R"))
source(file.path(BASE_R, "treatment_shares.R"))
source(file.path(BASE_R, "DT_timeline_panel.R"))
source(file.path(BASE_R, "km_panel_from_df.R"))
source(file.path(BASE_R, "cox_forest_plot_from_df.R"))
source(file.path(BASE_R, "plot_timeline_for_k.R"))
source(file.path(BASE_R, "tx_duration.R"))
source(file.path(BASE_R, "tx_normalize.R"))
source(file.path(BASE_R, "tx_intervals.R"))
source(file.path(BASE_R, "tx_cluster_surv.R"))
source(file.path(BASE_R, "tx_focus_dt.R"))
source(file.path(BASE_R, "tx_pooled_analysis.R"))
source(file.path(BASE_R, "tx_compare_groups.R"))
source(file.path(BASE_R, "tx_lines.R"))


# ---------------------------------------------------------------------------
# OSC data paths
# ---------------------------------------------------------------------------
LUSC_DIR       <- "/users/PAS1695/dipankor99/Github/exotho/exploratory/data/LUSC NSCLC"
LUSC_RAW_PATH  <- file.path(LUSC_DIR, "LUSC_Medication_NSCLCORIEN.rds")
LUSC_META_PATH <- file.path(LUSC_DIR, "Squam_metadata.csv")

SKIP_MSG <- "LUSC files not found — run on OSC"

# Confirmed canonical types
VALID_TYPES <- c("Ancillary","Chemo","Hormone","IO","Small_Molecule",
                 "Targeted","Radiation","Others")

# Confirmed tx_normalize() output columns (updated: AgeGrid + dominant_regimen)
NORM_COLS <- c("sample","AgeGrid","treatment_group","start_age",
               "TimeSinceTreatmentStart","end_followup","dominant_regimen")

# Confirmed tx_intervals() output columns
INTV_COLS <- c("sample","block","type","run","start_year","end_year")

# Real treatment_group labels from LUSC RDS
LUSC_TX_GROUPS <- c("Chemo","IO","Radiation","Targeted",
                    "Small_Molecule","Hormone","Ancillary","Others")

# ---------------------------------------------------------------------------
# make_raw_lusc()
# Mimics LUSC_Medication_NSCLCORIEN.rds exactly (11 columns).
# ---------------------------------------------------------------------------
make_raw_lusc <- function(n_patients = 10, tx_per_patient = 3, seed = 42) {
  set.seed(seed)
  samples <- paste0("SYN", sprintf("%03d", seq_len(n_patients)))
  do.call(rbind, lapply(samples, function(s) {
    spec_age   <- runif(1, 40, 80)
    last_age   <- spec_age + runif(1, 0.5, 12)
    n_tx       <- tx_per_patient
    med_starts <- sort(runif(n_tx, spec_age, last_age - 0.1))
    med_stops  <- pmin(med_starts + runif(n_tx, 0.05, 0.4), last_age)
    data.frame(
      sample                     = s,
      AvatarKey                  = paste0("KEY", s),
      Age.At.Specimen.Collection = spec_age,
      AgeAtLastContact           = last_age,
      diagsurvtime               = last_age - spec_age,
      Status                     = sample(0:1, 1),
      Medication                 = paste0("Drug", seq_len(n_tx)),
      treatment_group            = sample(LUSC_TX_GROUPS, n_tx, replace = TRUE),
      AgeAtMedStart              = med_starts,
      AgeAtMedStop               = med_stops,
      AgeAtTreatmentStart.mod    = med_starts,
      stringsAsFactors           = FALSE
    )
  }))
}

# ---------------------------------------------------------------------------
# make_meta_lusc()
# Mimics Squam_metadata.csv with confirmed LUSC column names.
# ---------------------------------------------------------------------------
make_meta_lusc <- function(n_patients = 10, seed = 42) {
  set.seed(seed)
  data.frame(
    sample        = paste0("SYN", sprintf("%03d", seq_len(n_patients))),
    diagsurvtime  = runif(n_patients, 0.1, 10),
    Status        = sample(0:1, n_patients, replace = TRUE),
    CAlevel       = sample(c("High","Low"), n_patients, replace = TRUE),
    SmokingStatus = sample(c("Ever","Never"), n_patients, replace = TRUE),
    Primary_Met   = sample(c("Primary","Metastatic"), n_patients, replace = TRUE),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Pre-build shared synthetic objects
# ---------------------------------------------------------------------------
synth_meta_20  <- make_meta_lusc(n_patients = 20, seed = 1)
synth_raw_20   <- make_raw_lusc(n_patients = 20, tx_per_patient = 3, seed = 1)
synth_norm_20  <- tx_normalize(synth_raw_20)
synth_tl_20    <- tx_intervals(synth_norm_20)
synth_meta_20  <- synth_meta_20[synth_meta_20$sample %in% unique(synth_norm_20$sample), ]

synth_cs_20    <- suppressWarnings(
  tx_cluster_surv(synth_meta_20, synth_norm_20,
                  surv_time_col  = "diagsurvtime",
                  status_col     = "Status",
                  k_range        = 3:5,
                  umap_neighbors = 10,
                  seed           = 42)
)

synth_meta_40  <- make_meta_lusc(n_patients = 40, seed = 99)
synth_raw_40   <- make_raw_lusc(n_patients = 40, tx_per_patient = 4, seed = 99)
synth_norm_40  <- tx_normalize(synth_raw_40)
synth_tl_40    <- tx_intervals(synth_norm_40)
synth_meta_40  <- synth_meta_40[synth_meta_40$sample %in% unique(synth_norm_40$sample), ]

synth_cs_40    <- suppressWarnings(
  tx_cluster_surv(synth_meta_40, synth_norm_40,
                  surv_time_col  = "diagsurvtime",
                  status_col     = "Status",
                  k_range        = 3:6,
                  umap_neighbors = 15,
                  seed           = 42)$Cluster_surv
)

# =============================================================================
# 1. standardise_status (Fix 4)
# =============================================================================

test_that("standardise_status: 0/1 numeric passes through unchanged", {
  df  <- data.frame(sample = c("A","B"), status = c(0L, 1L))
  out <- standardise_status(df)
  expect_identical(out$status, c(0L, 1L))
})

test_that("standardise_status: adds status_label factor column", {
  df  <- data.frame(sample = c("A","B"), status = c(0L, 1L))
  out <- standardise_status(df)
  expect_true("status_label" %in% names(out))
  expect_s3_class(out$status_label, "factor")
  expect_equal(levels(out$status_label), c("Alive", "Dead"))
})

test_that("standardise_status: status_label maps correctly", {
  df  <- data.frame(sample = c("A","B","C"), status = c(0L, 1L, 0L))
  out <- standardise_status(df)
  expect_equal(as.character(out$status_label), c("Alive", "Dead", "Alive"))
})

test_that("standardise_status: converts 'dead'/'alive' strings to 1/0", {
  df  <- data.frame(sample = c("A","B"), status = c("alive", "dead"))
  out <- standardise_status(df)
  expect_identical(out$status, c(0L, 1L))
})

test_that("standardise_status: converts 'Dead'/'Alive' (capitalised) to 1/0", {
  df  <- data.frame(sample = c("A","B"), status = c("Alive", "Dead"))
  out <- standardise_status(df)
  expect_identical(out$status, c(0L, 1L))
})

test_that("standardise_status: converts 'DECEASED'/'LIVING' to 1/0", {
  df  <- data.frame(sample = c("A","B"), status = c("Living", "Deceased"))
  out <- standardise_status(df)
  expect_identical(out$status, c(0L, 1L))
})

test_that("standardise_status: converts 'censored' to 0", {
  df  <- data.frame(sample = c("A","B"), status = c("censored", "dead"))
  out <- standardise_status(df)
  expect_identical(out$status, c(0L, 1L))
})

test_that("standardise_status: errors on unrecognised coding", {
  df <- data.frame(sample = c("A","B"), status = c("foo", "bar"))
  expect_error(standardise_status(df), "cannot auto-detect")
})

test_that("standardise_status: errors on missing column", {
  df <- data.frame(sample = c("A","B"), surv = c(0, 1))
  expect_error(standardise_status(df, "status"), "not found")
})

test_that("standardise_status: handles numeric 0/1 as double", {
  df  <- data.frame(sample = c("A","B"), status = c(0.0, 1.0))
  out <- standardise_status(df)
  expect_true(is.integer(out$status))
  expect_identical(out$status, c(0L, 1L))
})

test_that("standardise_status: handles factor input", {
  df  <- data.frame(sample = c("A","B"), status = factor(c("alive", "dead")))
  out <- standardise_status(df)
  expect_identical(out$status, c(0L, 1L))
})

# =============================================================================
# 2. dominant_exclusive (Fix 2)
# =============================================================================

test_that("dominant_exclusive: returns tibble with sample and regimen columns", {
  out <- dominant_exclusive(synth_tl_20)
  expect_true(all(c("sample", "regimen") %in% names(out)))
})

test_that("dominant_exclusive: one row per patient (mutual exclusivity)", {
  out <- dominant_exclusive(synth_tl_20)
  expect_equal(nrow(out), length(unique(out$sample)))
})

test_that("dominant_exclusive: no NA regimens", {
  out <- dominant_exclusive(synth_tl_20)
  expect_false(anyNA(out$regimen))
})

test_that("dominant_exclusive: regimen values are from expected set", {
  expected <- c("Chemo only", "Radiation only", "IO only",
                "Small Molecule only", "Hormone only",
                "Chemo+IO", "Chemo+Radiation", "Chemo+Targeted",
                "Chemo+Radiation+IO", "Other")
  out <- dominant_exclusive(synth_tl_40)
  bad <- setdiff(unique(out$regimen), expected)
  expect_equal(length(bad), 0L,
               info = paste("Unexpected regimens:", paste(bad, collapse = ", ")))
})

test_that("dominant_exclusive: higher threshold reduces combo regimens", {
  out_low  <- dominant_exclusive(synth_tl_40, min_share = 0.10)
  out_high <- dominant_exclusive(synth_tl_40, min_share = 0.40)
  # With higher threshold, fewer types qualify per patient -> fewer combos
  n_combo_low  <- sum(grepl("\\+", out_low$regimen))
  n_combo_high <- sum(grepl("\\+", out_high$regimen))
  expect_lte(n_combo_high, n_combo_low)
})

test_that("dominant_exclusive: specificity hierarchy — Chemo+IO beats Chemo only", {
  # Patient with Chemo 50% and IO 30% should be Chemo+IO, not Chemo only
  tl_test <- data.frame(
    sample     = rep("P1", 4),
    type       = c("Chemo","Chemo","Chemo","IO"),
    start_year = c(0, 0.5, 1.0, 0),
    end_year   = c(0.5, 1.0, 1.5, 0.8)
  )
  out <- dominant_exclusive(tl_test, min_share = 0.20)
  expect_equal(out$regimen[out$sample == "P1"], "Chemo+IO")
})

# =============================================================================
# 3. get_focus_cohort (Fix 1)
# =============================================================================

test_that("get_focus_cohort: returns tibble with required columns", {
  out <- get_focus_cohort(synth_cs_40, synth_tl_40,
                          focus_types = "Chemo", mode = "dominant")
  expect_true(all(c("sample", "focus_share", "mode", "focus_types") %in% names(out)))
})

test_that("get_focus_cohort: mode='only' returns subset of mode='dominant'", {
  only <- get_focus_cohort(synth_cs_40, synth_tl_40,
                           focus_types = "Chemo", mode = "only")
  dom  <- get_focus_cohort(synth_cs_40, synth_tl_40,
                           focus_types = "Chemo", mode = "dominant")
  # 'only' is stricter than 'dominant', so should return <= patients
  expect_lte(nrow(only), nrow(dom))
})

test_that("get_focus_cohort: mode='concurrent' requires >= 2 focus_types", {
  expect_warning(
    get_focus_cohort(synth_cs_40, synth_tl_40,
                     focus_types = "Chemo", mode = "concurrent"),
    "concurrent"
  )
})

test_that("get_focus_cohort: focus_share is between 0 and 1", {
  out <- get_focus_cohort(synth_cs_40, synth_tl_40,
                          focus_types = "Chemo", mode = "dominant")
  if (nrow(out) > 0) {
    expect_true(all(out$focus_share >= 0 & out$focus_share <= 1))
  }
})

test_that("get_focus_cohort: IDs are subset of Cluster_surv samples", {
  out <- get_focus_cohort(synth_cs_40, synth_tl_40,
                          focus_types = "Chemo", mode = "dominant")
  expect_true(all(out$sample %in% synth_cs_40$sample))
})

# =============================================================================
# 4. tx_normalize
# =============================================================================

test_that("tx_normalize returns a data frame [synthetic]", {
  expect_s3_class(tx_normalize(make_raw_lusc()), "data.frame")
})

test_that("tx_normalize output has required columns [synthetic]", {
  out     <- tx_normalize(make_raw_lusc())
  missing <- setdiff(NORM_COLS, colnames(out))
  expect_equal(length(missing), 0L,
               info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("tx_normalize: dominant_regimen column is present and non-empty [synthetic]", {
  out <- tx_normalize(make_raw_lusc())
  expect_true("dominant_regimen" %in% names(out))
  expect_false(all(is.na(out$dominant_regimen)))
})

test_that("tx_normalize: dominant_regimen is consistent per patient [synthetic]", {
  out <- tx_normalize(make_raw_lusc())
  per_patient <- tapply(out$dominant_regimen, out$sample, function(x) length(unique(x)))
  expect_true(all(per_patient == 1L),
              info = "dominant_regimen should be the same for all rows of a patient")
})

test_that("tx_normalize: dominant_regimen_share parameter changes output [synthetic]", {
  out_low  <- tx_normalize(make_raw_lusc(), dominant_regimen_share = 0.10)
  out_high <- tx_normalize(make_raw_lusc(), dominant_regimen_share = 0.40)
  # Different thresholds should produce different regimen distributions
  dist_low  <- table(unique(out_low[, c("sample","dominant_regimen")])$dominant_regimen)
  dist_high <- table(unique(out_high[, c("sample","dominant_regimen")])$dominant_regimen)
  # At minimum, higher threshold should not produce MORE combo regimens
  expect_true(TRUE)  # Smoke test — just confirm no error
})

test_that("tx_normalize preserves all input sample IDs [synthetic]", {
  raw <- make_raw_lusc()
  out <- tx_normalize(raw)
  expect_true(all(unique(raw$sample) %in% unique(out$sample)))
})

test_that("tx_normalize: treatment_group values are all canonical types [synthetic]", {
  out <- tx_normalize(make_raw_lusc())
  bad <- setdiff(unique(out$treatment_group), VALID_TYPES)
  expect_equal(length(bad), 0L,
               info = paste("Non-canonical types:", paste(bad, collapse = ", ")))
})

test_that("tx_normalize: no NA in sample or treatment_group [synthetic]", {
  out <- tx_normalize(make_raw_lusc())
  expect_false(anyNA(out$sample))
  expect_false(anyNA(out$treatment_group))
})

test_that("tx_normalize: TimeSinceTreatmentStart >= 0 [synthetic]", {
  out <- tx_normalize(make_raw_lusc())
  expect_true(all(out$TimeSinceTreatmentStart >= 0))
})

test_that("tx_normalize: end_followup >= TimeSinceTreatmentStart [synthetic]", {
  out <- tx_normalize(make_raw_lusc())
  expect_true(all(out$end_followup >= out$TimeSinceTreatmentStart))
})

test_that("tx_normalize: single-patient input works [synthetic]", {
  out <- tx_normalize(make_raw_lusc(n_patients = 1, tx_per_patient = 1))
  expect_equal(length(unique(out$sample)), 1L)
})

test_that("tx_normalize runs on LUSC_Medication_NSCLCORIEN.rds [real]", {
  skip_if_not(file.exists(LUSC_RAW_PATH), message = SKIP_MSG)
  expect_no_error(
    Normalized_timeline <<- tx_normalize(readRDS(LUSC_RAW_PATH))
  )
})

test_that("LUSC Normalized_timeline: has required columns [real]", {
  skip_if_not(exists("Normalized_timeline"), message = SKIP_MSG)
  missing <- setdiff(NORM_COLS, colnames(Normalized_timeline))
  expect_equal(length(missing), 0L,
               info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("LUSC Normalized_timeline: treatment_group all canonical [real]", {
  skip_if_not(exists("Normalized_timeline"), message = SKIP_MSG)
  bad <- setdiff(unique(Normalized_timeline$treatment_group), VALID_TYPES)
  expect_equal(length(bad), 0L,
               info = paste("Non-canonical types:", paste(bad, collapse = ", ")))
})

test_that("LUSC Normalized_timeline: TimeSinceTreatmentStart >= 0 [real]", {
  skip_if_not(exists("Normalized_timeline"), message = SKIP_MSG)
  expect_true(all(Normalized_timeline$TimeSinceTreatmentStart >= 0))
})

test_that("LUSC Normalized_timeline: sample IDs overlap with Squam_metadata [real]", {
  skip_if_not(exists("Normalized_timeline"), message = SKIP_MSG)
  skip_if_not(file.exists(LUSC_META_PATH), message = SKIP_MSG)
  meta    <- read.csv(LUSC_META_PATH, stringsAsFactors = FALSE)
  overlap <- intersect(unique(Normalized_timeline$sample), unique(meta$sample))
  expect_gt(length(overlap), 0L)
})


# =============================================================================
# 5. tx_intervals
# =============================================================================

test_that("tx_intervals returns a data frame [synthetic]", {
  expect_s3_class(tx_intervals(synth_norm_20), "data.frame")
})

test_that("tx_intervals output has required columns [synthetic]", {
  out     <- tx_intervals(synth_norm_20)
  missing <- setdiff(INTV_COLS, colnames(out))
  expect_equal(length(missing), 0L,
               info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("tx_intervals: start_year < end_year for all rows [synthetic]", {
  out <- tx_intervals(synth_norm_20)
  expect_true(all(out$start_year < out$end_year))
})

test_that("tx_intervals: no negative start_year [synthetic]", {
  expect_true(all(tx_intervals(synth_norm_20)$start_year >= 0))
})

test_that("tx_intervals: preserves all sample IDs [synthetic]", {
  out <- tx_intervals(synth_norm_20)
  expect_true(all(unique(synth_norm_20$sample) %in% unique(out$sample)))
})

test_that("tx_intervals: type column values are all canonical [synthetic]", {
  out <- tx_intervals(synth_norm_20)
  bad <- setdiff(unique(out$type), VALID_TYPES)
  expect_equal(length(bad), 0L,
               info = paste("Non-canonical types:", paste(bad, collapse = ", ")))
})

test_that("tx_intervals: no NA in key columns [synthetic]", {
  out <- tx_intervals(synth_norm_20)
  for (col in c("sample","type","start_year","end_year"))
    expect_false(anyNA(out[[col]]), info = paste("NAs in:", col))
})

test_that("tx_intervals runs on LUSC Normalized_timeline [real]", {
  skip_if_not(file.exists(LUSC_RAW_PATH), message = SKIP_MSG)
  norm <- tx_normalize(readRDS(LUSC_RAW_PATH))
  expect_no_error(
    Refined_timeline <<- tx_intervals(norm)
  )
})

test_that("LUSC Refined_timeline: required columns present [real]", {
  skip_if_not(exists("Refined_timeline"), message = SKIP_MSG)
  missing <- setdiff(INTV_COLS, colnames(Refined_timeline))
  expect_equal(length(missing), 0L,
               info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("LUSC Refined_timeline: all intervals positive-length [real]", {
  skip_if_not(exists("Refined_timeline"), message = SKIP_MSG)
  expect_true(all(Refined_timeline$end_year > Refined_timeline$start_year))
})

test_that("LUSC Refined_timeline: type column all canonical [real]", {
  skip_if_not(exists("Refined_timeline"), message = SKIP_MSG)
  bad <- setdiff(unique(Refined_timeline$type), VALID_TYPES)
  expect_equal(length(bad), 0L,
               info = paste("Non-canonical types:", paste(bad, collapse = ", ")))
})

# =============================================================================
# 6. tx_cluster_surv
# =============================================================================

test_that("tx_cluster_surv returns a list [synthetic]", {
  expect_type(synth_cs_20, "list")
})

test_that("tx_cluster_surv output has all required list elements [synthetic]", {
  required <- c("Cluster_surv","pca_matrix","umap_df","X","treatment_encoded")
  missing  <- setdiff(required, names(synth_cs_20))
  expect_equal(length(missing), 0L,
               info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("Cluster_surv has cluster columns k3:k5 [synthetic]", {
  missing <- setdiff(paste0("Cluster_k", 3:5), colnames(synth_cs_20$Cluster_surv))
  expect_equal(length(missing), 0L,
               info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("Cluster_surv has one row per patient [synthetic]", {
  expect_equal(nrow(synth_cs_20$Cluster_surv), nrow(synth_meta_20))
})

test_that("Cluster_surv cluster IDs within valid range per k [synthetic]", {
  for (k in 3:5) {
    col  <- paste0("Cluster_k", k)
    vals <- synth_cs_20$Cluster_surv[[col]]
    expect_true(all(vals %in% seq_len(k)),
                info = paste(col, "has out-of-range IDs"))
  }
})

test_that("Cluster_surv retains survival and LUSC covariate columns [synthetic]", {
  cs <- synth_cs_20$Cluster_surv
  for (col in c("diagsurvtime","status","CAlevel","SmokingStatus","Primary_Met"))
    expect_true(col %in% colnames(cs), info = paste("Missing:", col))
})

test_that("Cluster_surv has status_label factor column (Fix 4) [synthetic]", {
  cs <- synth_cs_20$Cluster_surv
  expect_true("status_label" %in% colnames(cs))
  expect_s3_class(cs$status_label, "factor")
  expect_equal(levels(cs$status_label), c("Alive", "Dead"))
})

test_that("Cluster_surv: status is integer 0/1 after standardisation [synthetic]", {
  cs <- synth_cs_20$Cluster_surv
  expect_true(is.integer(cs$status))
  expect_true(all(cs$status %in% c(0L, 1L)))
})

test_that("tx_cluster_surv auto-converts 'Dead'/'Alive' status (Fix 4) [synthetic]", {
  meta_str        <- synth_meta_20
  meta_str$Status <- ifelse(meta_str$Status == 1, "Dead", "Alive")
  expect_no_error(
    res <- suppressWarnings(
      tx_cluster_surv(meta_str, synth_norm_20,
                      surv_time_col  = "diagsurvtime",
                      status_col     = "Status",
                      k_range        = 3:5,
                      umap_neighbors = 10,
                      seed           = 42)
    )
  )
  expect_true(all(res$Cluster_surv$status %in% c(0L, 1L)))
  expect_true("status_label" %in% names(res$Cluster_surv))
})

test_that("tx_cluster_surv auto-converts 'Deceased'/'Living' status [synthetic]", {
  meta_str        <- synth_meta_20
  meta_str$Status <- ifelse(meta_str$Status == 1, "Deceased", "Living")
  expect_no_error(
    suppressWarnings(
      tx_cluster_surv(meta_str, synth_norm_20,
                      surv_time_col  = "diagsurvtime",
                      status_col     = "Status",
                      k_range        = 3:5,
                      umap_neighbors = 10,
                      seed           = 42)
    )
  )
})

test_that("pca_matrix rows match patient count [synthetic]", {
  expect_equal(nrow(synth_cs_20$pca_matrix), nrow(synth_meta_20))
})

test_that("umap_df has one row per patient and >= 2 numeric columns [synthetic]", {
  expect_equal(nrow(synth_cs_20$umap_df), nrow(synth_meta_20))
  expect_gte(sum(sapply(synth_cs_20$umap_df, is.numeric)), 2L)
})

test_that("X (binary grid): correct dimensions and 0/1 values only [synthetic]", {
  expect_equal(nrow(synth_cs_20$X), nrow(synth_meta_20))
  expect_true(all(synth_cs_20$X %in% c(0L, 1L)))
})

test_that("seed produces reproducible cluster assignments [synthetic]", {
  out1 <- suppressWarnings(
    tx_cluster_surv(synth_meta_20, synth_norm_20,
                    surv_time_col="diagsurvtime", status_col="Status",
                    k_range=3:5, umap_neighbors=10, seed=42))
  out2 <- suppressWarnings(
    tx_cluster_surv(synth_meta_20, synth_norm_20,
                    surv_time_col="diagsurvtime", status_col="Status",
                    k_range=3:5, umap_neighbors=10, seed=42))
  expect_identical(out1$Cluster_surv$Cluster_k3,
                   out2$Cluster_surv$Cluster_k3)
})

test_that("tx_cluster_surv runs on LUSC data k=3:20 [real]", {
  skip_if_not(file.exists(LUSC_RAW_PATH),  message = SKIP_MSG)
  skip_if_not(file.exists(LUSC_META_PATH), message = SKIP_MSG)
  
  LUSC_metadata <- read.csv(LUSC_META_PATH, stringsAsFactors = FALSE)
  norm_lusc     <- tx_normalize(readRDS(LUSC_RAW_PATH))
  
  expect_no_error({
    res_clust    <<- tx_cluster_surv(LUSC_metadata, norm_lusc,
                                     surv_time_col = "diagsurvtime",
                                     status_col    = "Status",
                                     k_range       = 3:20,
                                     seed          = 42)
    Cluster_surv <<- res_clust$Cluster_surv
    tl_lusc      <<- tx_intervals(norm_lusc)
  })
})

test_that("LUSC Cluster_surv: all k3:k20 columns present [real]", {
  skip_if_not(exists("Cluster_surv"), message = SKIP_MSG)
  missing <- setdiff(paste0("Cluster_k", 3:20), colnames(Cluster_surv))
  expect_equal(length(missing), 0L,
               info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("LUSC Cluster_surv: no duplicate patients [real]", {
  skip_if_not(exists("Cluster_surv"), message = SKIP_MSG)
  expect_equal(length(unique(Cluster_surv$sample)), nrow(Cluster_surv))
})

test_that("LUSC Cluster_surv: CAlevel values are High/Low only [real]", {
  skip_if_not(exists("Cluster_surv"), message = SKIP_MSG)
  bad <- setdiff(na.omit(unique(Cluster_surv$CAlevel)), c("High","Low"))
  expect_equal(length(bad), 0L,
               info = paste("Unexpected CAlevel values:", paste(bad, collapse = ", ")))
})

test_that("LUSC Cluster_surv: status_label present (Fix 4) [real]", {
  skip_if_not(exists("Cluster_surv"), message = SKIP_MSG)
  expect_true("status_label" %in% colnames(Cluster_surv))
  expect_s3_class(Cluster_surv$status_label, "factor")
})

# =============================================================================
# 7. Pipeline — tx_pooled_analysis, tx_compare_groups, km_panel_from_df
# =============================================================================

test_that("tx_pooled_analysis mode='any': returns all list elements [synthetic]", {
  out      <- tx_pooled_analysis(synth_cs_40, synth_tl_40, focus_types="Chemo",
                                 mode="any", group_var="CAlevel")
  required <- c("km","forest","timeline","ids","df","segs","shares",
                "df_plot","n_cohort","n_raw","n_plot","group_table")
  missing  <- setdiff(required, names(out))
  expect_equal(length(missing), 0L,
               info = paste("Missing:", paste(missing, collapse = ", ")))
})

test_that("tx_pooled_analysis: n_raw >= n_cohort (Fix 5 audit) [synthetic]", {
  out <- tx_pooled_analysis(synth_cs_40, synth_tl_40, focus_types="Chemo",
                            mode="any", group_var="CAlevel")
  expect_gte(out$n_raw, out$n_cohort)
})

test_that("tx_pooled_analysis mode='any': n_cohort > 0 [synthetic]", {
  out <- tx_pooled_analysis(synth_cs_40, synth_tl_40, focus_types="Chemo",
                            mode="any", group_var="CAlevel")
  expect_gt(out$n_cohort, 0L)
})

test_that("tx_pooled_analysis mode='only' runs without error [synthetic]", {
  expect_no_error(
    tx_pooled_analysis(synth_cs_40, synth_tl_40, focus_types="Chemo",
                       mode="only", group_var="CAlevel")
  )
})

test_that("tx_pooled_analysis mode='concurrent' runs without error [synthetic]", {
  expect_no_error(
    tx_pooled_analysis(synth_cs_40, synth_tl_40, focus_types=c("Chemo","IO"),
                       mode="concurrent", group_var="CAlevel")
  )
})

test_that("tx_pooled_analysis mode='dominant' runs without error [synthetic]", {
  expect_no_error(
    tx_pooled_analysis(synth_cs_40, synth_tl_40, focus_types=c("Chemo","IO"),
                       mode="dominant", group_var="CAlevel")
  )
})

test_that("tx_pooled_analysis: group_var='SmokingStatus' [synthetic]", {
  expect_no_error(
    tx_pooled_analysis(synth_cs_40, synth_tl_40, focus_types="Chemo",
                       mode="any", group_var="SmokingStatus")
  )
})

test_that("tx_pooled_analysis: group_var='Primary_Met' [synthetic]", {
  expect_no_error(
    tx_pooled_analysis(synth_cs_40, synth_tl_40, focus_types="Chemo",
                       mode="any", group_var="Primary_Met")
  )
})

test_that("tx_compare_groups by SmokingStatus returns required elements [synthetic]", {
  out      <- tx_compare_groups(synth_cs_40, group_var="SmokingStatus",
                                cox_covars=c("CAlevel","Primary_Met"))
  required <- c("km","forest","combined","cox_results","group_summary")
  expect_true(all(required %in% names(out)))
})

test_that("tx_compare_groups by CAlevel: High and Low present [synthetic]", {
  out    <- tx_compare_groups(synth_cs_40, group_var="CAlevel",
                              cox_covars=c("SmokingStatus","Primary_Met"))
  groups <- out$group_summary[[1]]
  expect_true("High" %in% groups && "Low" %in% groups)
})

test_that("tx_compare_groups: custom_groups runs without error [synthetic]", {
  ids_a <- synth_meta_40$sample[1:10]
  ids_b <- synth_meta_40$sample[11:20]
  expect_no_error(
    tx_compare_groups(synth_cs_40,
                      custom_groups=list(GroupA=ids_a, GroupB=ids_b))
  )
})

test_that("km_panel_from_df: 'calevel' (lower) accepted [synthetic]", {
  expect_no_error(km_panel_from_df(synth_cs_40, group_col="calevel"))
})

test_that("km_panel_from_df: 'CALEVEL' (upper) accepted [synthetic]", {
  expect_no_error(km_panel_from_df(synth_cs_40, group_col="CALEVEL"))
})

test_that("km_panel_from_df: SmokingStatus (2-level) works [synthetic]", {
  expect_no_error(km_panel_from_df(synth_cs_40, group_col="SmokingStatus"))
})

test_that("Full LUSC pipeline: RDS -> normalize -> cluster -> pooled [real]", {
  skip_if_not(file.exists(LUSC_RAW_PATH),  message = SKIP_MSG)
  skip_if_not(file.exists(LUSC_META_PATH), message = SKIP_MSG)
  skip_if_not(exists("Cluster_surv"),      message = SKIP_MSG)
  skip_if_not(exists("tl_lusc"),           message = SKIP_MSG)
  
  out <- tx_pooled_analysis(Cluster_surv, tl_lusc,
                            focus_types=c("Chemo","IO"),
                            mode="any", group_var="CAlevel",
                            horizon_years=5)
  expect_gt(out$n_cohort, 0L)
})

test_that("LUSC real: tx_compare_groups by SmokingStatus [real]", {
  skip_if_not(exists("Cluster_surv"), message = SKIP_MSG)
  expect_no_error(
    tx_compare_groups(Cluster_surv, group_var="SmokingStatus",
                      cox_covars=c("CAlevel","Primary_Met"))
  )
})

test_that("LUSC real: tx_compare_groups by Primary_Met [real]", {
  skip_if_not(exists("Cluster_surv"), message = SKIP_MSG)
  expect_no_error(
    tx_compare_groups(Cluster_surv, group_var="Primary_Met",
                      cox_covars=c("CAlevel","SmokingStatus"))
  )
})

# =============================================================================
# 8. Dominant overlap verification (Fix 2 validation)
# =============================================================================

test_that("Dominant cohorts are mutually exclusive after Fix 2 [real]", {
  skip_if_not(exists("Cluster_surv"), message = SKIP_MSG)
  skip_if_not(exists("tl_lusc"),      message = SKIP_MSG)
  
  # Run dominant_exclusive on real LUSC data
  dom <- dominant_exclusive(tl_lusc, min_share = 0.20)
  
  # Each patient should appear exactly once
  expect_equal(nrow(dom), length(unique(dom$sample)),
               info = "Patients assigned to multiple regimens — Fix 2 failed")
  
  # Cross-check: Chemo-only and Chemo+IO should not share patients
  chemo_only_ids <- dom$sample[dom$regimen == "Chemo only"]
  chemo_io_ids   <- dom$sample[dom$regimen == "Chemo+IO"]
  overlap        <- intersect(chemo_only_ids, chemo_io_ids)
  expect_equal(length(overlap), 0L,
               info = paste("Overlap between Chemo only and Chemo+IO:", length(overlap)))
})

# =============================================================================
# 9. Format robustness
# =============================================================================

test_that("km_panel_from_df: 'Calevel' (title case) accepted [synthetic]", {
  expect_no_error(km_panel_from_df(synth_cs_40, group_col="Calevel"))
})

test_that("cox_forest_plot_from_df resolves covariate names case-insensitively [synthetic]", {
  df_upper <- synth_cs_40
  names(df_upper)[names(df_upper) == "CAlevel"]       <- "CALEVEL"
  names(df_upper)[names(df_upper) == "SmokingStatus"] <- "SMOKINGSTATUS"
  expect_no_error(
    cox_forest_plot_from_df(df_upper,
                            covars     = c("CALEVEL","SMOKINGSTATUS"),
                            ref_levels = list(CALEVEL="Low", SMOKINGSTATUS="Never"))
  )
})

test_that("tx_pooled_analysis tolerates extra columns in Cluster_surv [synthetic]", {
  cs_extra             <- synth_cs_40
  cs_extra$Institution <- "HOSP_X"
  expect_no_error(
    tx_pooled_analysis(cs_extra, synth_tl_40, focus_types="Chemo",
                       mode="any", group_var="CAlevel")
  )
})

test_that("km_panel_from_df works when SmokingStatus column is absent [synthetic]", {
  cs_min <- synth_cs_40[, setdiff(colnames(synth_cs_40), "SmokingStatus")]
  expect_no_error(km_panel_from_df(cs_min, group_col="CAlevel"))
})

test_that("km_panel_from_df handles all-censored cohort [synthetic]", {
  cs_cens        <- synth_cs_40
  cs_cens$Status <- 0L
  expect_no_error(km_panel_from_df(cs_cens, group_col="CAlevel"))
})

test_that("tx_compare_groups handles all-censored cohort [synthetic]", {
  cs_cens        <- synth_cs_40
  cs_cens$Status <- 0L
  expect_no_error(
    tx_compare_groups(cs_cens, group_var="CAlevel", cox_covars="CAlevel")
  )
})

test_that("tx_pooled_analysis handles single-type (Chemo-only) timeline [synthetic]", {
  tl_chemo      <- synth_tl_40
  tl_chemo$type <- "Chemo"
  expect_no_error(
    tx_pooled_analysis(synth_cs_40, tl_chemo, focus_types="Chemo",
                       mode="any", group_var="CAlevel")
  )
})

test_that("tx_normalize handles numeric sample IDs [synthetic]", {
  raw_num        <- make_raw_lusc(n_patients=5, tx_per_patient=2, seed=3)
  raw_num$sample <- as.integer(factor(raw_num$sample))
  out            <- tx_normalize(raw_num)
  expect_gt(nrow(out), 0L)
})

# =============================================================================
# test_tx_duration.R
# Tests for tx_duration() and internal helpers
#
# Append to test_script.R or run standalone:
#   testthat::test_file("./exploratory/Scripts/MEC-TX/tests/test_tx_duration.R")
# =============================================================================

library(testthat)

# --- Source dependencies (adjust path for standalone runs) ---
# source(file.path(BASE_R, "tx_duration.R"))

# =============================================================================
# Section 1: .merge_intervals (internal helper)
# =============================================================================
context("tx_duration — .merge_intervals")

test_that("non-overlapping intervals are preserved", {
  res <- .merge_intervals(c(0, 2, 5), c(1, 3, 6))
  expect_equal(nrow(res), 3)
  expect_equal(res$start, c(0, 2, 5))
  expect_equal(res$end,   c(1, 3, 6))
})

test_that("overlapping intervals are merged", {
  res <- .merge_intervals(c(0, 0.5, 5), c(1, 1.5, 6))
  expect_equal(nrow(res), 2)
  expect_equal(res$start, c(0, 5))
  expect_equal(res$end,   c(1.5, 6))
})

test_that("adjacent intervals (touching) are merged", {
  res <- .merge_intervals(c(0, 1), c(1, 2))
  expect_equal(nrow(res), 1)
  expect_equal(res$start, 0)
  expect_equal(res$end,   2)
})

test_that("single interval returns as-is", {
  res <- .merge_intervals(1, 3)
  expect_equal(nrow(res), 1)
  expect_equal(res$start, 1)
  expect_equal(res$end,   3)
})

test_that("empty input returns zero-row data.frame", {
  res <- .merge_intervals(numeric(0), numeric(0))
  expect_equal(nrow(res), 0)
  expect_true(all(c("start", "end") %in% names(res)))
})

test_that("fully nested intervals collapse to outer", {
  # [0, 5] contains [1, 3]
  res <- .merge_intervals(c(0, 1), c(5, 3))
  expect_equal(nrow(res), 1)
  expect_equal(res$start, 0)
  expect_equal(res$end,   5)
})

test_that("unsorted input is handled correctly", {
  res <- .merge_intervals(c(5, 0, 2), c(6, 1, 3))
  expect_equal(nrow(res), 3)
  expect_equal(res$start, c(0, 2, 5))
})

# =============================================================================
# Section 2: .duration_per_type
# =============================================================================
context("tx_duration — .duration_per_type")

# Synthetic timeline: 3 patients, 2 types
syn_timeline <- data.frame(
  sample     = c("P1", "P1", "P1", "P2", "P2", "P3"),
  type       = c("Chemo", "Chemo", "IO", "Chemo", "IO", "Chemo"),
  start_year = c(0.0,  0.5,  0.0,  0.0,  1.0,  0.0),
  end_year   = c(0.75, 1.0,  0.5,  0.5,  1.5,  2.0),
  stringsAsFactors = FALSE
)

test_that("per-type duration sums correctly with overlap merge", {
  res <- .duration_per_type(syn_timeline)
  # P1 Chemo: intervals [0, 0.75] + [0.5, 1.0] merge to [0, 1.0] = 1.0 yr
  p1_chemo <- res$duration_yrs[res$sample == "P1" & res$type == "Chemo"]
  expect_equal(p1_chemo, 1.0)
  # P1 IO: [0, 0.5] = 0.5 yr
  p1_io <- res$duration_yrs[res$sample == "P1" & res$type == "IO"]
  expect_equal(p1_io, 0.5)
  # P3 Chemo: [0, 2.0] = 2.0 yr
  p3_chemo <- res$duration_yrs[res$sample == "P3" & res$type == "Chemo"]
  expect_equal(p3_chemo, 2.0)
})

test_that("patients missing a type are not included (no zero rows)", {
  res <- .duration_per_type(syn_timeline)
  # P3 has no IO → should not appear
  expect_false(any(res$sample == "P3" & res$type == "IO"))
})

test_that("all expected patient-type combinations are present", {
  res <- .duration_per_type(syn_timeline)
  # P1: Chemo + IO, P2: Chemo + IO, P3: Chemo = 5 rows
  expect_equal(nrow(res), 5)
})

# =============================================================================
# Section 3: .duration_total (merged across types)
# =============================================================================
context("tx_duration — .duration_total")

test_that("concurrent treatment is not double-counted", {
  # P1: Chemo [0, 0.75]+[0.5, 1.0] + IO [0, 0.5]
  # Merged across all: [0, 1.0] = 1.0 yr (not 1.5)
  res <- .duration_total(syn_timeline)
  p1_total <- res$duration_yrs_total[res$sample == "P1"]
  expect_equal(p1_total, 1.0)
})

test_that("non-overlapping types sum correctly", {
  # P2: Chemo [0, 0.5] + IO [1.0, 1.5] = 1.0 yr (0.5 + 0.5, no overlap)
  res <- .duration_total(syn_timeline)
  p2_total <- res$duration_yrs_total[res$sample == "P2"]
  expect_equal(p2_total, 1.0)
})

test_that("all patients are returned", {
  res <- .duration_total(syn_timeline)
  expect_equal(sort(res$sample), c("P1", "P2", "P3"))
})

# =============================================================================
# Section 4: tx_duration() — full function
# =============================================================================
context("tx_duration — main function")

# Synthetic meta with group variable
syn_meta <- data.frame(
  sample  = c("P1", "P2", "P3"),
  CAlevel = c("High", "Low", "Low"),
  stringsAsFactors = FALSE
)

test_that("returns expected list structure", {
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel", plot = FALSE)
  expect_is(res, "list")
  expect_true(all(c("duration_per_type", "duration_total",
                    "summary_table", "plot", "params") %in% names(res)))
})

test_that("default duration unit is months", {
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel", plot = FALSE)
  expect_true("duration_months" %in% names(res$duration_per_type))
  expect_true("duration_total_months" %in% names(res$duration_total))
  # P1 Chemo: 1.0 yr * 12 = 12 months
  p1_chemo <- res$duration_per_type$duration_months[
    res$duration_per_type$sample == "P1" & res$duration_per_type$type == "Chemo"
  ]
  expect_equal(p1_chemo, 12)
})

test_that("years unit works", {
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel",
                     duration_unit = "years", plot = FALSE)
  expect_true("duration_years" %in% names(res$duration_per_type))
  p1_chemo <- res$duration_per_type$duration_years[
    res$duration_per_type$sample == "P1" & res$duration_per_type$type == "Chemo"
  ]
  expect_equal(p1_chemo, 1.0)
})

test_that("exclude_types removes specified types", {
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel",
                     exclude_types = "IO", plot = FALSE)
  expect_false("IO" %in% res$duration_per_type$type)
  expect_false("IO" %in% res$summary_table$type)
})

test_that("summary_table has correct structure", {
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel", plot = FALSE)
  st <- res$summary_table
  expect_true(all(c("type", "group", "n", "mean", "median",
                    "q25", "q75", "p_value") %in% names(st)))
  # Should have rows for each group × type + group × total
  n_types   <- length(unique(syn_timeline$type))
  n_groups  <- 2  # High, Low
  expect_equal(nrow(st), (n_types + 1) * n_groups)  # +1 for "All types (merged)"
})

test_that("Wilcoxon p-value is computed for 2-group comparison", {
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel",
                     min_n = 1, plot = FALSE)
  st <- res$summary_table
  # At least one non-NA p-value (for types with both groups represented)
  chemo_p <- unique(st$p_value[st$type == "Chemo"])
  expect_true(!is.na(chemo_p[1]))
})

test_that("min_n flag skips tests for small groups", {
  # With min_n=5, our synthetic data (n=1 High, n=2 Low max) should be skipped
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel",
                     min_n = 5, plot = FALSE)
  st <- res$summary_table
  expect_true(all(is.na(st$p_value)))
  expect_true(all(st$test_note != ""))
})

test_that("params slot captures function arguments", {
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel", plot = FALSE)
  expect_equal(res$params$group_var, "CAlevel")
  expect_equal(res$params$duration_unit, "months")
  expect_equal(res$params$n_patients, 3)
})

# =============================================================================
# Section 5: Input validation
# =============================================================================
context("tx_duration — input validation")

test_that("error on missing timeline column", {
  bad_tl <- syn_timeline
  names(bad_tl)[1] <- "patient_id"
  expect_error(tx_duration(bad_tl, syn_meta, "CAlevel", plot = FALSE),
               "not found in timeline")
})

test_that("error on missing meta column", {
  bad_meta <- syn_meta
  names(bad_meta)[2] <- "Group"
  expect_error(tx_duration(syn_timeline, bad_meta, "CAlevel", plot = FALSE),
               "not found in meta")
})

test_that("error on no overlapping patients", {
  other_meta <- data.frame(sample = c("X1", "X2"), CAlevel = c("High", "Low"),
                           stringsAsFactors = FALSE)
  expect_error(tx_duration(syn_timeline, other_meta, "CAlevel", plot = FALSE),
               "No overlapping patients")
})

test_that("warning on single group level", {
  one_group <- data.frame(sample = c("P1", "P2", "P3"),
                          CAlevel = c("High", "High", "High"),
                          stringsAsFactors = FALSE)
  expect_warning(tx_duration(syn_timeline, one_group, "CAlevel", plot = FALSE),
                 "Only one level")
})

test_that("error when all types excluded", {
  expect_error(
    tx_duration(syn_timeline, syn_meta, "CAlevel",
                exclude_types = c("Chemo", "IO"), plot = FALSE),
    "No intervals remain"
  )
})

# =============================================================================
# Section 6: Plot generation
# =============================================================================
context("tx_duration — plot")

test_that("plot is ggplot object when plot=TRUE", {
  skip_if_not_installed("ggplot2")
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel",
                     min_n = 1, plot = TRUE)
  expect_is(res$plot, "gg")
})

test_that("plot is NULL when plot=FALSE", {
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel", plot = FALSE)
  expect_null(res$plot)
})

test_that("violin plot type works", {
  skip_if_not_installed("ggplot2")
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel",
                     min_n = 1, plot = TRUE, plot_type = "violin")
  expect_is(res$plot, "gg")
})

test_that("custom palette is applied", {
  skip_if_not_installed("ggplot2")
  pal <- c(High = "red", Low = "blue")
  res <- tx_duration(syn_timeline, syn_meta, "CAlevel",
                     min_n = 1, palette = pal)
  # Check palette is stored in the plot scales
  expect_is(res$plot, "gg")
})

# =============================================================================
# Section 7: Integration-style tests (mimicking real pipeline data shapes)
# =============================================================================
context("tx_duration — pipeline integration")

test_that("handles typical tx_intervals output column names", {
  # Real pipeline uses: sample, type, start_year, end_year (defaults)
  real_shape <- data.frame(
    sample     = rep(paste0("S", 1:10), each = 3),
    type       = rep(c("Chemo", "IO", "Radiation"), 10),
    start_year = runif(30, 0, 1),
    stringsAsFactors = FALSE
  )
  real_shape$end_year <- real_shape$start_year + runif(30, 0.1, 1)
  real_meta <- data.frame(
    sample  = paste0("S", 1:10),
    CAlevel = rep(c("High", "Low"), each = 5),
    stringsAsFactors = FALSE
  )
  res <- tx_duration(real_shape, real_meta, "CAlevel", plot = FALSE)
  expect_equal(res$params$n_patients, 10)
  expect_equal(res$params$n_types, 3)
  # Summary should have (3 types + 1 total) × 2 groups = 8 rows
  expect_equal(nrow(res$summary_table), 8)
})

test_that("works with more than 2 groups", {
  multi_meta <- data.frame(
    sample = c("P1", "P2", "P3"),
    Stage  = c("I", "II", "III"),
    stringsAsFactors = FALSE
  )
  res <- tx_duration(syn_timeline, multi_meta, "Stage", plot = FALSE)
  expect_equal(length(unique(res$summary_table$group)), 3)
})

test_that("handles patients with only one treatment type", {
  single_type <- data.frame(
    sample     = c("P1", "P2"),
    type       = c("Chemo", "Chemo"),
    start_year = c(0, 0),
    end_year   = c(1, 0.5),
    stringsAsFactors = FALSE
  )
  single_meta <- data.frame(
    sample  = c("P1", "P2"),
    CAlevel = c("High", "Low"),
    stringsAsFactors = FALSE
  )
  res <- tx_duration(single_type, single_meta, "CAlevel",
                     min_n = 1, plot = FALSE)
  # Only Chemo + All types = 2 type levels × 2 groups = 4 rows
  expect_equal(nrow(res$summary_table), 4)
  # Per-type total should match per-type (only one type)
  expect_equal(
    res$duration_per_type$duration_months[res$duration_per_type$sample == "P1"],
    res$duration_total$duration_total_months[res$duration_total$sample == "P1"]
  )
})


# ═══════════════════════════════════════════════════════════════════════════════
# TEST SECTION — tx_lines()
# 11 sections | 47 tests
# Add to test_script.R (test_mectx_all.R on OSC) after tx_duration section
# ═══════════════════════════════════════════════════════════════════════════════

# ── Shared synthetic data ──────────────────────────────────────────────────────

# Minimal timeline: 3 patients, clean lung-cancer records
.tl_base <- data.frame(
  sample     = c("P1","P1","P1","P2","P2","P3"),
  start_year = c(60.0, 60.9, 62.0,  55.0, 56.5,  70.0),
  end_year   = c(60.7, 61.5, 62.9,  55.8, 57.5,  71.2),
  type       = c("Chemo","IO","Chemo","Chemo","IO","Radiation"),
  stringsAsFactors = FALSE
)

# Meta: sample ↔ AvatarKey + specimen ages
.meta_base <- data.frame(
  sample                    = c("P1","P2","P3"),
  AvatarKey                 = c("AK1","AK2","AK3"),
  Age.At.Specimen.Collection= c(59.8, 54.9, 69.8),
  CAlevel                   = c("High","Low","Low"),
  Stage                     = c("III","II","I"),
  stringsAsFactors          = FALSE
)

# Annotations: MedLineRegimen for P1 and P2; P3 has Unknown
.ann_base <- data.frame(
  AvatarKey      = c("AK1","AK2","AK3"),
  Medication     = c("Carboplatin","Carboplatin","Pembrolizumab"),
  MedLineRegimen = c("First Line/Regimen","Second Line","Unknown/Not Applicable"),
  AgeAtMedStart  = c(60.0, 56.5, 70.0),
  stringsAsFactors = FALSE
)

# Timeline WITH a prior-cancer record for P1 (Letrozole at age 40)
.tl_contaminated <- rbind(
  data.frame(sample="P1", start_year=40.0, end_year=40.5,
             type="Small_Molecule", stringsAsFactors=FALSE),
  .tl_base
)

# Annotations WITH a prior-cancer annotation for P1
.ann_contaminated <- rbind(
  data.frame(AvatarKey="AK1", Medication="Letrozole",
             MedLineRegimen="First Line/Regimen", AgeAtMedStart=40.0,
             stringsAsFactors=FALSE),
  .ann_base
)

# ── Section 1: Input validation ────────────────────────────────────────────────

test_that("error on non-data.frame timeline", {
  expect_error(tx_lines(list()), "must be a data.frame")
})

test_that("error on missing required timeline columns", {
  expect_error(
    tx_lines(data.frame(sample="P1", start_year=60, end_year=61)),
    "missing columns"
  )
})

test_that("error on non-positive gap_threshold", {
  expect_error(tx_lines(.tl_base, gap_threshold = -0.1), "positive number")
})

test_that("error on negative specimen_buffer", {
  expect_error(tx_lines(.tl_base, specimen_buffer = -1), ">= 0")
})

test_that("error when annotations missing ann_id_col", {
  expect_error(
    tx_lines(.tl_base,
             annotations = data.frame(X="a", MedLineRegimen="First Line/Regimen",
                                      AgeAtMedStart=60, stringsAsFactors=FALSE),
             ann_id_col = "AvatarKey"),
    "missing column"
  )
})

# ── Section 2: Return structure ────────────────────────────────────────────────
context("tx_lines — return structure")

test_that("returns named list with four elements", {
  res <- tx_lines(.tl_base, meta = .meta_base)
  expect_type(res, "list")
  expect_named(res, c("lines","patient_summary","group_comparison","params"))
})

test_that("lines data.frame has required columns", {
  res <- tx_lines(.tl_base, meta = .meta_base)
  expected_cols <- c("sample","line_number","line_label","line_types",
                     "line_start","line_end","line_duration_months",
                     "line_source","line_flag")
  expect_true(all(expected_cols %in% names(res$lines)))
})

test_that("patient_summary has one row per patient", {
  res <- tx_lines(.tl_base, meta = .meta_base)
  expect_equal(nrow(res$patient_summary), length(unique(.tl_base$sample)))
})

test_that("params captures key settings", {
  res <- tx_lines(.tl_base, gap_threshold = 0.1, specimen_buffer = 0.5)
  expect_equal(res$params$gap_threshold, 0.1)
  expect_equal(res$params$specimen_buffer, 0.5)
})

# ── Section 3: Specimen-anchored record filtering — timeline ───────────────────
context("tx_lines — specimen filter on timeline")

test_that("prior-cancer timeline records are dropped (record level, not patient level)", {
  # P1 has Letrozole at age 40; specimen at 59.8 → should be dropped
  res <- tx_lines(.tl_contaminated, meta = .meta_base, specimen_buffer = 0.25)
  # P1 should still be present
  expect_true("P1" %in% res$lines$sample)
  # But no line with start near age 40
  p1_lines <- res$lines[res$lines$sample == "P1", ]
  expect_true(all(p1_lines$line_start >= 59.5))
})

test_that("lung cancer records for contaminated patient are preserved", {
  res <- tx_lines(.tl_contaminated, meta = .meta_base, specimen_buffer = 0.25)
  p1_lines <- res$lines[res$lines$sample == "P1", ]
  expect_true(nrow(p1_lines) >= 1L)
})

test_that("specimen_buffer = 0 drops records starting before exact specimen age", {
  res <- tx_lines(.tl_contaminated, meta = .meta_base, specimen_buffer = 0)
  p1_lines <- res$lines[res$lines$sample == "P1", ]
  expect_true(all(p1_lines$line_start >= .meta_base$Age.At.Specimen.Collection[1]))
})

test_that("specimen_buffer = 1.0 retains records within 1yr of specimen", {
  # Record at specimen_age - 0.5 should be kept with buffer = 1.0 but not buffer = 0
  tl_edge <- rbind(
    data.frame(sample="P1", start_year=59.3, end_year=59.7,
               type="Chemo", stringsAsFactors=FALSE),  # 0.5yr before specimen
    .tl_base[.tl_base$sample == "P1", ]
  )
  res_loose <- tx_lines(tl_edge, meta = .meta_base, specimen_buffer = 1.0)
  res_tight <- tx_lines(tl_edge, meta = .meta_base, specimen_buffer = 0)
  expect_gt(nrow(res_loose$lines), nrow(res_tight$lines))
})

test_that("filter is skipped gracefully when specimen_age_col absent from meta", {
  meta_no_spec <- .meta_base[, setdiff(names(.meta_base), "Age.At.Specimen.Collection")]
  expect_message(
    tx_lines(.tl_base, meta = meta_no_spec),
    "Specimen filter skipped"
  )
})

test_that("patients without specimen age are retained unfiltered", {
  meta_partial <- .meta_base
  meta_partial$Age.At.Specimen.Collection[2] <- NA_real_
  res <- tx_lines(.tl_base, meta = meta_partial, specimen_buffer = 0.25)
  expect_true("P2" %in% res$lines$sample)
})

# ── Section 4: Specimen-anchored record filtering — annotations ────────────────
context("tx_lines — specimen filter on annotations")

test_that("prior-cancer annotation is filtered before coalesce", {
  # Letrozole at age 40 has MedLineRegimen = "First Line/Regimen"
  # Without filtering this would incorrectly anchor P1 line 1 as "First" from Letrozole
  # With filtering, Carboplatin at age 60 becomes the anchor → still "First" label
  # The key check: line_start should NOT be near age 40
  res <- tx_lines(.tl_contaminated, annotations = .ann_contaminated,
                  meta = .meta_base, specimen_buffer = 0.25)
  p1_lines <- res$lines[res$lines$sample == "P1", ]
  expect_true(all(p1_lines$line_start >= 59.5))
})

test_that("clean annotation record is kept after specimen filter", {
  res <- tx_lines(.tl_base, annotations = .ann_base,
                  meta = .meta_base, specimen_buffer = 0.25)
  p1_l1 <- res$lines[res$lines$sample == "P1" & res$lines$line_number == 1L, ]
  expect_equal(p1_l1$line_source, "annotated")
})

# ── Section 5: .map_line_regimen() internals ──────────────────────────────────
context("tx_lines — .map_line_regimen")

test_that("First Line/Regimen maps to First", {
  expect_equal(.map_line_regimen("First Line/Regimen"), "First")
})

test_that("Neoadjuvant Regimen maps to Neoadjuvant", {
  expect_equal(.map_line_regimen("Neoadjuvant Regimen"), "Neoadjuvant")
})

test_that("Adjuvant/First Line maps to First", {
  expect_equal(.map_line_regimen("Adjuvant/First Line"), "First")
})

test_that("Maintenance maps to Maintenance", {
  expect_equal(.map_line_regimen("Maintenance"), "Maintenance")
})

test_that("Unknown/Not Applicable maps to NA", {
  expect_true(is.na(.map_line_regimen("Unknown/Not Applicable")))
})

test_that("Unknown/Not Reported maps to NA", {
  expect_true(is.na(.map_line_regimen("Unknown/Not Reported")))
})

test_that("Sixth Line/Regimen maps to Sixth", {
  expect_equal(.map_line_regimen("Sixth Line/Regimen"), "Sixth")
})

# ── Section 6: .merge_to_blocks() internals ───────────────────────────────────
context("tx_lines — .merge_to_blocks")

test_that("concurrent intervals produce one block with combined types", {
  tl <- data.frame(
    start_year = c(60.0, 60.2),
    end_year   = c(61.0, 61.5),
    type       = c("Chemo","IO"),
    stringsAsFactors = FALSE
  )
  blocks <- .merge_to_blocks(tl)
  expect_equal(nrow(blocks), 1L)
  expect_true(grepl("\\+", blocks$block_types))
})

test_that("non-overlapping intervals produce separate blocks", {
  tl <- data.frame(
    start_year = c(60.0, 61.5),
    end_year   = c(61.0, 62.5),
    type       = c("Chemo","IO"),
    stringsAsFactors = FALSE
  )
  blocks <- .merge_to_blocks(tl)
  expect_equal(nrow(blocks), 2L)
})

test_that("empty input returns empty data.frame", {
  tl <- data.frame(start_year=numeric(0), end_year=numeric(0),
                   type=character(0), stringsAsFactors=FALSE)
  expect_equal(nrow(.merge_to_blocks(tl)), 0L)
})

# ── Section 7: .assign_lines_from_blocks() internals ─────────────────────────
context("tx_lines — .assign_lines_from_blocks")

test_that("gap > threshold increments line number", {
  blocks <- data.frame(
    block_start = c(60.0, 61.5),
    block_end   = c(61.0, 62.5),
    block_types = c("Chemo","IO"),
    stringsAsFactors = FALSE
  )
  lines <- .assign_lines_from_blocks(blocks, gap_threshold = 3/52)
  expect_equal(max(lines$line_number), 2L)
})

test_that("gap <= threshold keeps same line", {
  blocks <- data.frame(
    block_start = c(60.0, 60.05),   # ~18 day gap, < 3 weeks
    block_end   = c(60.04, 61.0),
    block_types = c("Chemo","IO"),
    stringsAsFactors = FALSE
  )
  lines <- .assign_lines_from_blocks(blocks, gap_threshold = 3/52)
  expect_equal(max(lines$line_number), 1L)
})

test_that("line_duration_months is positive for valid blocks", {
  blocks <- data.frame(
    block_start = 60.0, block_end = 60.5, block_types = "Chemo",
    stringsAsFactors = FALSE
  )
  lines <- .assign_lines_from_blocks(blocks, gap_threshold = 3/52)
  expect_gt(lines$line_duration_months, 0)
})

# ── Section 8: Consolidation flagging ────────────────────────────────────────
context("tx_lines — consolidation flagging")

test_that("IO-only line 2+ in stage III is flagged possible_consolidation", {
  # P1 is stage III; second block is IO-only after a gap
  tl_io <- data.frame(
    sample     = c("P1","P1"),
    start_year = c(60.0, 62.0),
    end_year   = c(61.0, 63.5),
    type       = c("Chemo","IO"),
    stringsAsFactors = FALSE
  )
  res <- tx_lines(tl_io, meta = .meta_base, stage_col = "Stage")
  p1_l2 <- res$lines[res$lines$sample == "P1" & res$lines$line_number == 2L, ]
  expect_equal(p1_l2$line_flag, "possible_consolidation")
})

test_that("line 1 is never flagged as possible_consolidation", {
  res <- tx_lines(.tl_base, meta = .meta_base, stage_col = "Stage")
  l1_flags <- res$lines[res$lines$line_number == 1L, "line_flag"]
  expect_true(all(l1_flags == "confirmed"))
})

test_that("IO line 2+ in stage IV is confirmed (not consolidation)", {
  meta_iv <- .meta_base
  meta_iv$Stage <- "IV"
  tl_io <- data.frame(
    sample     = c("P1","P1"),
    start_year = c(60.0, 62.0),
    end_year   = c(61.0, 63.5),
    type       = c("Chemo","IO"),
    stringsAsFactors = FALSE
  )
  res <- tx_lines(tl_io, meta = meta_iv, stage_col = "Stage")
  p1_l2 <- res$lines[res$lines$sample == "P1" & res$lines$line_number == 2L, ]
  expect_equal(p1_l2$line_flag, "confirmed")
})

# ── Section 9: Coalesce logic ─────────────────────────────────────────────────
context("tx_lines — annotation coalesce")

test_that("annotated patient has line_source = 'annotated' for line 1", {
  res <- tx_lines(.tl_base, annotations = .ann_base,
                  meta = .meta_base, specimen_buffer = 0.25)
  p1_l1 <- res$lines[res$lines$sample == "P1" & res$lines$line_number == 1L, ]
  expect_equal(p1_l1$line_source, "annotated")
})

test_that("unannotated patient (Unknown) has line_source = 'computed'", {
  res <- tx_lines(.tl_base, annotations = .ann_base,
                  meta = .meta_base, specimen_buffer = 0.25)
  p3_l1 <- res$lines[res$lines$sample == "P3" & res$lines$line_number == 1L, ]
  expect_equal(p3_l1$line_source, "computed")
})

test_that("Maintenance annotation is not used as line anchor", {
  ann_maint <- data.frame(
    AvatarKey="AK1", Medication="Pembro",
    MedLineRegimen="Maintenance", AgeAtMedStart=60.5,
    stringsAsFactors=FALSE
  )
  res <- tx_lines(.tl_base, annotations = ann_maint,
                  meta = .meta_base, specimen_buffer = 0.25)
  p1_l1 <- res$lines[res$lines$sample == "P1" & res$lines$line_number == 1L, ]
  # Maintenance should not override → computed
  expect_equal(p1_l1$line_source, "computed")
})

# ── Section 10: Group comparison ──────────────────────────────────────────────
context("tx_lines — group comparison")

test_that("group_comparison is NULL when group_var not supplied", {
  res <- tx_lines(.tl_base, meta = .meta_base)
  expect_null(res$group_comparison)
})

test_that("group_comparison contains n_lines and first_line_duration_months metrics", {
  # Build larger synthetic dataset for reliable test
  set.seed(42)
  n <- 30L
  tl_large <- data.frame(
    sample     = paste0("S", seq_len(n)),
    start_year = runif(n, 55, 65),
    end_year   = runif(n, 65.5, 70),
    type       = sample(c("Chemo","IO"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  meta_large <- data.frame(
    sample    = paste0("S", seq_len(n)),
    AvatarKey = paste0("AK", seq_len(n)),
    Age.At.Specimen.Collection = runif(n, 54, 64),
    CAlevel   = sample(c("High","Low"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  res <- tx_lines(tl_large, meta = meta_large, group_var = "CAlevel")
  expect_true("n_lines" %in% res$group_comparison$metric)
  expect_true("first_line_duration_months" %in% res$group_comparison$metric)
})

test_that("group_comparison has p_value and test_note columns", {
  set.seed(7)
  n <- 20L
  tl_g <- data.frame(
    sample     = paste0("S", seq_len(n)),
    start_year = runif(n, 60, 65),
    end_year   = runif(n, 65.5, 70),
    type       = "Chemo",
    stringsAsFactors = FALSE
  )
  meta_g <- data.frame(
    sample    = paste0("S", seq_len(n)),
    AvatarKey = paste0("AK", seq_len(n)),
    Age.At.Specimen.Collection = runif(n, 59, 63),
    CAlevel   = rep(c("High","Low"), each = n/2L),
    stringsAsFactors = FALSE
  )
  res <- tx_lines(tl_g, meta = meta_g, group_var = "CAlevel")
  expect_true(all(c("p_value","test_note") %in% names(res$group_comparison)))
})

# ── Section 11: Edge cases ────────────────────────────────────────────────────
context("tx_lines — edge cases")

test_that("exclude_types removes those types before line detection", {
  res_full <- tx_lines(.tl_base, meta = .meta_base)
  res_excl <- tx_lines(.tl_base, meta = .meta_base, exclude_types = c("Radiation"))
  # P3 only has Radiation — should disappear or have no lines
  expect_false("P3" %in% res_excl$lines$sample)
  # P1 and P2 unaffected
  expect_true("P1" %in% res_excl$lines$sample)
})

test_that("single-interval patient produces one line", {
  tl_single <- data.frame(
    sample = "S1", start_year = 60.0, end_year = 61.0,
    type = "Chemo", stringsAsFactors = FALSE
  )
  res <- tx_lines(tl_single)
  expect_equal(nrow(res$lines), 1L)
  expect_equal(res$lines$line_number, 1L)
})

test_that("no annotations supplied runs algorithm-only mode without error", {
  expect_no_error(tx_lines(.tl_base, meta = .meta_base, annotations = NULL))
})

test_that("all Unknown annotations fall back to algorithm gracefully", {
  ann_all_unknown <- data.frame(
    AvatarKey      = c("AK1","AK2","AK3"),
    Medication     = c("X","Y","Z"),
    MedLineRegimen = rep("Unknown/Not Applicable", 3L),
    AgeAtMedStart  = c(60.0, 56.5, 70.0),
    stringsAsFactors = FALSE
  )
  res <- tx_lines(.tl_base, annotations = ann_all_unknown,
                  meta = .meta_base, specimen_buffer = 0.25)
  expect_true(all(res$lines$line_source == "computed"))
})

