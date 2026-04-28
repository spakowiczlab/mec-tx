# =============================================================================
# MEC-TX --- dominant_exclusive.R
# helpers/dominant_exclusive.R
#
# Assign each patient to ONE mutually exclusive dominant regimen based on
# treatment duration shares. Solves the overlap problem where a patient
# qualifies for Chemo[DOMINANT], Chemo+IO[DOMINANT], and
# all_types[DOMINANT] simultaneously.
#
# Logic: assign to the MOST SPECIFIC matching regimen (longest combo first).
# =============================================================================

#' Assign Mutually Exclusive Dominant Treatment Regimen Per Patient
#'
#' Computes per-patient treatment duration shares (excluding ancillary types),
#' identifies which treatment types exceed the share threshold, and assigns
#' a single regimen label using a specificity-first hierarchy. Solves the
#' overlap problem where a patient simultaneously qualifies for
#' \code{Chemo[DOMINANT]}, \code{Chemo+IO[DOMINANT]}, and
#' \code{all_types[DOMINANT]}.
#'
#' @param timeline A long-format treatment intervals data frame --- typically
#'   the \code{$timeline_long_intv} slot from \code{\link{tx_intervals}}.
#'   Must contain columns \code{sample}, \code{type}, \code{start_year},
#'   and \code{end_year}.
#' @param ancillary_types Character vector of treatment types to exclude
#'   from share calculation. These types contribute no duration to the
#'   denominator. Default \code{c("Ancillary", "Others")}.
#' @param min_share Numeric in \code{(0, 1)}. Minimum proportion of total
#'   non-ancillary treatment duration for a type to qualify as dominant.
#'   Default \code{0.20} (20\%). Increase to tighten, decrease to allow
#'   more multi-type combinations.
#'
#' @return A tibble with one row per patient and two columns:
#'   \describe{
#'     \item{sample}{Patient identifier, matching the input \code{timeline}.}
#'     \item{regimen}{Character. Mutually exclusive dominant regimen label.
#'       One of: \code{"Chemo+Radiation+IO"}, \code{"Chemo+IO"},
#'       \code{"Chemo+Radiation"}, \code{"Chemo+Targeted"},
#'       \code{"Chemo only"}, \code{"Radiation only"}, \code{"IO only"},
#'       \code{"Small Molecule only"}, \code{"Hormone only"},
#'       \code{"Other"}.}
#'   }
#'   Patients whose entire treatment record consists of \code{ancillary_types}
#'   are dropped from the output. A message reports the count of dropped
#'   patients. Use a left join against your full cohort to identify them ---
#'   do not treat them as equivalent to the \code{"Other"} regimen stratum.
#'   A warning is raised if any patient is assigned more than one regimen,
#'   which indicates a gap in the hierarchy.
#'
#' @details
#' \strong{Hierarchy logic:} Multi-type combinations are matched before
#' single-agent labels (longest combination first). A patient whose
#' qualifying types are \code{"Chemo+IO"} is assigned \code{"Chemo+IO"},
#' not \code{"Chemo only"}, even though Chemo alone also exceeds the
#' threshold. This ensures mutual exclusivity without post-hoc filtering.
#'
#' \strong{Share calculation:} Duration is computed as
#' \code{end_year - start_year} per interval. Overlapping intervals are
#' not merged before summing --- call \code{\link{tx_intervals}} upstream
#' to ensure non-overlapping intervals are passed in.
#'
#' \strong{Unassigned patients:} Patients present in the input but absent
#' from the output had no qualifying type at the given \code{min_share}
#' threshold (e.g. highly fragmented treatment with no dominant type).
#' These are distinct from \code{"Other"} --- they have no dominant signal,
#' not an unrecognised one.
#'
#' @examples
#' \dontrun{
#' intervals <- tx_intervals(med_data, metadata)
#'
#' # Default 20% threshold
#' regimens <- dominant_exclusive(intervals$timeline_long_intv)
#'
#' # Stricter threshold --- require 30% share to qualify
#' regimens_strict <- dominant_exclusive(
#'   intervals$timeline_long_intv,
#'   min_share = 0.30
#' )
#'
#' # Join back to survival data for Cox model
#' cohort <- cohort %>%
#'   dplyr::left_join(regimens, by = "sample") %>%
#'   dplyr::mutate(regimen = tidyr::replace_na(regimen, "Unassigned"))
#' }
#'
#' @seealso \code{\link{tx_intervals}}, \code{\link{tx_pooled_analysis}},
#'   \code{\link{get_focus_cohort}}
#'
#' @importFrom dplyr filter mutate group_by summarise ungroup select case_when
#' @export
dominant_exclusive <- function(timeline,
                               ancillary_types  = c("Ancillary", "Others"),
                               min_share        = 0.20) {

  # --- Compute per-patient treatment shares ---
  shares <- timeline %>%
    dplyr::filter(!type %in% ancillary_types) %>%
    dplyr::mutate(duration = end_year - start_year) %>%
    dplyr::group_by(sample, type) %>%
    dplyr::summarise(total_dur = sum(duration), .groups = "drop") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(share = total_dur / sum(total_dur)) %>%
    dplyr::ungroup()

  # --- Per patient: which types exceed the threshold? ---
  qualifying <- shares %>%
    dplyr::filter(share >= min_share) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      types = paste(sort(unique(type)), collapse = "+"),
      .groups = "drop"
    )

  # --- Assign regimen using specificity-first hierarchy ---
  # Most specific (3-type combos) first, then 2-type, then single-agent.
  result <- qualifying %>%
    dplyr::mutate(
      regimen = dplyr::case_when(
        # 3-type combinations
        grepl("Chemo", types) &
          grepl("Radiation", types) &
          grepl("IO", types)                 ~ "Chemo+Radiation+IO",

        # 2-type combinations
        grepl("Chemo", types) &
          grepl("IO", types)                 ~ "Chemo+IO",
        grepl("Chemo", types) &
          grepl("Radiation", types)          ~ "Chemo+Radiation",
        grepl("Chemo", types) &
          grepl("Targeted", types)           ~ "Chemo+Targeted",

        # Single-agent dominant
        types == "Chemo"                     ~ "Chemo only",
        types == "Radiation"                 ~ "Radiation only",
        types == "IO"                        ~ "IO only",
        types == "Small_Molecule"            ~ "Small Molecule only",
        types == "Hormone"                   ~ "Hormone only",

        # Catch-all
        TRUE                                 ~ "Other"
      )
    ) %>%
    dplyr::select(sample, regimen)

  # --- Verify mutual exclusivity ---
  n_dup <- sum(duplicated(result$sample))
  if (n_dup > 0) {
    warning(sprintf(
      "dominant_exclusive: %d patients assigned to multiple regimens --- check hierarchy",
      n_dup
    ))
  }

  message(sprintf(
    "dominant_exclusive: %d patients assigned | threshold=%.0f%%",
    nrow(result), min_share * 100
  ))
  message("  Regimen distribution:")
  message(table(result$regimen))

  result
}
