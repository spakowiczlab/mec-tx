# ============================================================
# MEC-TX helper: prep_segs()
# helpers/prep_segs.R
# ============================================================
# ============================================================
# MEC-TX helpers: prep_segs() and treatment_shares()
# helpers/prep_segs.R + helpers/treatment_shares.R
# ============================================================

# -------- segments -> canonical types --------
prep_segs <- function(segs, horizon_years = 5) {

  # valid_types defined locally — no global scope dependency (Bug 4.2 fix)
  local_valid_types <- c(
    "Ancillary", "Chemo", "Hormone", "IO",
    "Small_Molecule", "Targeted", "Radiation", "Others"
  )

  stopifnot(
    "sample"     %in% names(segs),
    "type"       %in% names(segs),
    "start_year" %in% names(segs),
    "end_year"   %in% names(segs)
  )

  segs %>%
    dplyr::mutate(type_raw = as.character(.data$type)) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        stringr::str_detect(.data$type_raw, stringr::regex(
          "radiation|\\brt\\b|xrt|imrt|sbrt|radiother", TRUE)) ~ "Radiation",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "\\bio\\b|immuno|pembro|nivol|atezo|ipi|pd-1|pd-l1|ctla", TRUE)) ~ "IO",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "chemo|chemotherapy|platin|taxel|5fu|gemcitabine|doxo", TRUE)) ~ "Chemo",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "target|tki|inhibitor|egfr|alk|braf|mek|parp|her2|trast", TRUE)) ~ "Targeted",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "hormone|endocrine|androgen|estrogen|tamox|abiraterone|enzalu|aromatase|letro|anastro|fulves", TRUE)) ~ "Hormone",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "small[ _-]?molecule|onco.?drug|oncodrug", TRUE)) ~ "Small_Molecule",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "ancillary|support|palliat|pain control|bisphosphonate", TRUE)) ~ "Ancillary",
        TRUE ~ "Others"
      ),
      t0 = pmax(0, pmin(.data$start_year, horizon_years)),
      t1 = pmax(0, pmin(.data$end_year,   horizon_years))
    ) %>%
    dplyr::filter(
      !is.na(.data$sample),
      !is.na(.data$t0),
      !is.na(.data$t1),
      .data$t1 > .data$t0
    )
}


