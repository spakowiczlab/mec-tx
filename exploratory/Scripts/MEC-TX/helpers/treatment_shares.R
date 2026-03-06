# ============================================================
# MEC-TX helper: treatment_shares()
# helpers/treatment_shares.R
# ============================================================
# -------- per-sample durations & shares (for all valid types) --------
treatment_shares <- function(segs_prepped) {

  # valid_types defined locally — no global scope dependency (Bug 4.2 fix)
  local_valid_types <- c(
    "Ancillary", "Chemo", "Hormone", "IO",
    "Small_Molecule", "Targeted", "Radiation", "Others"
  )

  wide <- segs_prepped %>%
    dplyr::group_by(sample, type) %>%
    dplyr::summarise(dur = sum(t1 - t0), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = type, values_from = dur, values_fill = 0)

  # ensure all treatment columns exist
  for (nm in c(local_valid_types, "None")) {
    if (!nm %in% names(wide)) wide[[nm]] <- 0
  }

  wide %>%
    dplyr::mutate(
      dur_all     = rowSums(dplyr::across(dplyr::all_of(c(local_valid_types, "None")))),
      dur_treated = rowSums(dplyr::across(dplyr::all_of(local_valid_types))),
      dplyr::across(
        dplyr::all_of(local_valid_types),
        ~ ifelse(dur_treated > 0, .x / dur_treated, 0),
        .names = "share_{col}_tx"
      )
    ) %>%
    dplyr::mutate(
      # dominant treated-time type — Bug 4.3 fix: use colnames() directly
      dom_type_tx = {
        shares     <- as.matrix(dplyr::select(., dplyr::starts_with("share_") &
                                                 dplyr::ends_with("_tx")))
        colnames(shares) <- gsub("^share_|_tx$", "", colnames(shares))
        colnames(shares)[max.col(shares, ties.method = "first")]
      }
    )
}
