tx_intervals <- function(
    df,
    id_col   = "sample",
    time_col = "TimeSinceTreatmentStart",
    type_col = "treatment_group",
    drop_types = c("None"),
    horizon_years = 5
){
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
  })
  
  # ---- 1) Clean + snap to biweek grid ----
  # Biweekly index (1/26 year ≈ 2 weeks) aligns with validated concurrency
  # window from sensitivity analysis and standard oncology cycle boundaries.
  dat <- df %>%
    transmute(
      sample     = .data[[id_col]],
      t_year_raw = as.numeric(.data[[time_col]]),
      type       = as.character(.data[[type_col]])
    ) %>%
    filter(!is.na(sample), !is.na(t_year_raw), !is.na(type)) %>%
    filter(!(type %in% drop_types)) %>%
    mutate(
      # convert to biweek index (26 biweeks per year)
      m = as.integer(round(t_year_raw * 26)),
      # clamp within horizon
      m = pmax(0L, pmin(m, as.integer(horizon_years * 26)))
    ) %>%
    distinct(sample, m, type)
  
  if (!nrow(dat)) {
    return(tibble::tibble(
      sample     = character(),
      type       = character(),
      block      = integer(),
      start_year = double(),
      end_year   = double()
    ))
  }
  
  # ---- 2) Define active-therapy blocks (no biweek gaps) ----
  biweek_blocks <- dat %>%
    distinct(sample, m) %>%
    arrange(sample, m) %>%
    group_by(sample) %>%
    mutate(
      new_block = m != lag(m, default = first(m)) + 1L,
      block     = cumsum(replace_na(new_block, TRUE))
    ) %>%
    ungroup()
  
  dat2 <- dat %>%
    left_join(biweek_blocks, by = c("sample", "m")) %>%
    arrange(sample, block, m, type)
  
  # ---- 3) Merge consecutive biweeks per type within each block ----
  segs <- dat2 %>%
    group_by(sample, block, type) %>%
    arrange(m, .by_group = TRUE) %>%
    mutate(
      new_run = m != lag(m, default = first(m)) + 1L,
      run     = cumsum(replace_na(new_run, TRUE))
    ) %>%
    group_by(sample, block, type, run) %>%
    summarise(
      start_year = min(m) / 26,
      end_year   = (max(m) + 1L) / 26,  # end-exclusive
      .groups = "drop"
    ) %>%
    arrange(sample, block, start_year, type)
  
  segs
}
