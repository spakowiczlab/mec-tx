### Function that normalize the timeline so that everything starts from 0.

tx_normalize <- function(
    Modified_medication,
    NSCLC_metadata = NULL
){
  suppressPackageStartupMessages({
    library(dplyr); library(tidyr)
  })
  
  # --- Clean and define mod start (biology: keep intervals, impute missing stop) ---
  df_clean <- Modified_medication %>%
    mutate(
      AgeAtMedStart            = as.numeric(AgeAtMedStart),
      AgeAtMedStop             = as.numeric(AgeAtMedStop),
      AgeAtLastContact         = as.numeric(AgeAtLastContact),
      Age.At.Specimen.Collection = as.numeric(Age.At.Specimen.Collection),
      # If stop missing, treat as ongoing until last follow-up
      AgeAtMedStop = if_else(is.na(AgeAtMedStop),
                             AgeAtLastContact,
                             AgeAtMedStop),
      # Do not count exposure before specimen collection
      AgeAtTreatmentStart.mod = pmax(Age.At.Specimen.Collection,
                                     AgeAtMedStart,
                                     na.rm = TRUE)
    ) %>%
    filter(!is.na(AgeAtMedStart)) %>%
    filter(!is.na(AgeAtTreatmentStart.mod), !is.na(AgeAtMedStop)) %>%
    filter(AgeAtMedStop > AgeAtTreatmentStart.mod)
  
  # --- Expand to biweekly grid (1/26 year ≈ 2 weeks) ---
  # Biweekly resolution chosen over monthly (1/12) because:
  # - 2-week window is the validated concurrency boundary from sensitivity analysis
  # - Approximates standard oncology cycle boundaries (q3w chemo, q3w/q6w IO)
  # - Reduces false concurrent classification from 3-week gaps within monthly bins
  # - NZV filter downstream removes sparse biweekly bins, keeping matrix tractable
  expanded_timeline <- df_clean %>%
    rowwise() %>%
    mutate(biweek_seq = list(seq(AgeAtTreatmentStart.mod, AgeAtMedStop, by = 1/26))) %>%
    ungroup() %>%
    unnest(biweek_seq) %>%
    mutate(biweek_seq = round(biweek_seq, 4)) %>%
    distinct(sample, treatment_group, biweek_seq)
  
  timeline_long_biweekly <- expanded_timeline %>%
    rename(AgeBiweek = biweek_seq) %>%
    select(sample, AgeBiweek, treatment_group)
  
  # --- Normalize so first treatment biweek = 0 ---
  start_age_per_patient <- timeline_long_biweekly %>%
    group_by(sample) %>%
    summarise(start_age = min(AgeBiweek, na.rm = TRUE), .groups = "drop")
  
  timeline_long_norm <- timeline_long_biweekly %>%
    left_join(start_age_per_patient, by = "sample") %>%
    mutate(TimeSinceTreatmentStart = round((AgeBiweek - start_age) * 26) / 26)
  
  # --- Join metadata (case-insensitive + safe) ---
  if (!is.null(NSCLC_metadata)) {
    NSCLC_metadata <- NSCLC_metadata %>% rename_with(tolower)
    keep_cols <- intersect(c("sample", "stage", "status"), names(NSCLC_metadata))
    timeline_long_norm <- timeline_long_norm %>%
      left_join(NSCLC_metadata %>% select(all_of(keep_cols)), by = "sample")
  }
  
  # --- Optional: keep end_followup for later truncation/plotting ---
  timeline_long_norm <- timeline_long_norm %>%
    left_join(
      Modified_medication %>%
        transmute(sample, end_followup = as.numeric(AgeAtLastContact)) %>%
        distinct(),
      by = "sample"
    )
  
  timeline_long_norm
}
