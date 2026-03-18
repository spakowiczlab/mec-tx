# ============================================================
# MEC-TX helper: prep_segs()
# helpers/prep_segs.R
# ============================================================

#' @noRd
#'
#' Recode and clip treatment segments to canonical types
#'
#' Validates a raw segments data frame, recodes free-text treatment type
#' labels to the eight canonical MEC-TX types via regex, and clips segment
#' start/end times to \code{[0, horizon_years]}. Zero-duration segments
#' (after clipping) and rows with \code{NA} in \code{sample}, \code{t0},
#' or \code{t1} are dropped.
#'
#' Called upstream by \code{plot_timeline_for_k()} and
#' \code{timeline_panel()} — not intended for direct use.
#'
#' @param segs A data frame with at minimum columns \code{sample},
#'   \code{type} (free-text treatment label), \code{start_year}, and
#'   \code{end_year} (numeric, years since first treatment).
#' @param horizon_years Numeric. Maximum time horizon in years. Segment
#'   endpoints are clipped to \code{[0, horizon_years]} via \code{pmin}/
#'   \code{pmax}. Default \code{5}.
#'
#' @return A data frame with the same rows as \code{segs} (minus dropped
#'   rows) and two new columns: \code{type} (recoded canonical label, one
#'   of \code{"Radiation"}, \code{"IO"}, \code{"Chemo"}, \code{"Targeted"},
#'   \code{"Hormone"}, \code{"Small_Molecule"}, \code{"Ancillary"},
#'   \code{"Others"}) and \code{t0} / \code{t1} (clipped start and end
#'   times in years). The original \code{start_year} / \code{end_year}
#'   columns are retained alongside \code{t0} / \code{t1}.
#'
#' @details
#' \strong{Regex priority:} Rules are applied in \code{case_when()} order —
#' first match wins. Order: Radiation → IO → Chemo → Targeted → Hormone →
#' Small_Molecule → Ancillary → Others. All patterns are
#' case-insensitive. Types already matching a canonical label pass through
#' correctly since the regexes cover the canonical names themselves.
#'
#' \strong{Local valid types:} The canonical type vector is defined
#' internally rather than referencing \code{valid_types} from
#' \code{constants.R}, avoiding a global scope dependency (Bug 4.2 fix).
#'
#' @importFrom dplyr mutate case_when filter
#' @importFrom stringr str_detect regex
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
          "radiation|\\brt\\b|xrt|imrt|sbrt|radiother", TRUE))           ~ "Radiation",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "\\bio\\b|immuno|pembro|nivol|atezo|ipi|pd-1|pd-l1|ctla", TRUE)) ~ "IO",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "chemo|chemotherapy|platin|taxel|5fu|gemcitabine|doxo", TRUE)) ~ "Chemo",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "target|tki|inhibitor|egfr|alk|braf|mek|parp|her2|trast", TRUE)) ~ "Targeted",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "hormone|endocrine|androgen|estrogen|tamox|abiraterone|enzalu|aromatase|letro|anastro|fulves", TRUE)) ~ "Hormone",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "small[ _-]?molecule|onco.?drug|oncodrug", TRUE))              ~ "Small_Molecule",
        stringr::str_detect(.data$type_raw, stringr::regex(
          "ancillary|support|palliat|pain control|bisphosphonate", TRUE)) ~ "Ancillary",
        TRUE                                                              ~ "Others"
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