# ============================================================
# MEC-TX helper: treatment_shares()
# helpers/treatment_shares.R
# ============================================================

#' Compute per-patient treatment type durations and shares
#' Takes a clipped and recoded segments data frame (output of
#' \code{prep_segs}) and computes per-patient treatment duration
#' and share for each of the eight canonical treatment types. Returns a
#' wide-format tibble with one row per patient.
#' Called upstream by \code{\link{timeline_panel}} and
#' \code{\link{plot_timeline_for_k}} --- not intended for direct use.
#' @param segs_prepped A data frame produced by \code{prep_segs}.
#'   Must contain columns \code{sample}, \code{type} (canonical label),
#'   \code{t0}, and \code{t1} (clipped segment start/end in years).
#' @return A wide-format tibble with one row per patient and the following
#'   columns:
#'   \describe{
#'     \item{sample}{Patient identifier.}
#'     \item{Ancillary, Chemo, Hormone, IO, Small_Molecule, Targeted,
#'       Radiation, Others}{Total duration (years) for each treatment type
#'       within the clipped horizon. Zero-filled for types not present.}
#'     \item{dur_all}{Total duration across all types including
#'       \code{"None"}.}
#'     \item{dur_treated}{Total duration across the eight canonical
#'       treatment types (excludes \code{"None"}).}
#'     \item{share_*_tx}{Duration share for each canonical type as a
#'       proportion of \code{dur_treated}. Zero when
#'       \code{dur_treated == 0}. Column names follow the pattern
#'       \code{share_<Type>_tx} (e.g. \code{share_Chemo_tx}).}
#'     \item{dom_type_tx}{Character. The canonical type with the highest
#'       share for this patient. Ties broken by first match in column
#'       order.}
#'   }
#' @details
#' \strong{Input requirement:} \code{segs_prepped} must be the output of
#' \code{prep_segs} --- specifically it must contain \code{t0} and
#' \code{t1} columns (clipped start/end times). Passing raw
#' \code{timeline_long_intv} directly will error because \code{start_year}
#' / \code{end_year} are not the same as \code{t0} / \code{t1}.
#' \strong{Local valid types:} The canonical type vector is defined
#' internally rather than referencing \code{valid_types} from
#' \code{constants.R}, avoiding a global scope dependency (Bug 4.2 fix).
#' \strong{Share calculation:} Shares are computed as a proportion of
#' \code{dur_treated} (non-ancillary, non-None duration), not
#' \code{dur_all}. A patient with only \code{"None"} segments gets
#' \code{share_*_tx = 0} for all types.
#' \strong{Dominant type:} \code{dom_type_tx} uses \code{max.col()} with
#' \code{ties.method = "first"} on the share matrix. Column order follows
#' \code{local_valid_types} (Bug 4.3 fix --- previously used unreliable
#' column name inference).
#' @importFrom dplyr group_by summarise mutate across all_of starts_with
#'   ends_with select
#' @importFrom tidyr pivot_wider
#' @noRd
treatment_shares <- function(segs_prepped) {
  
  # --- Input validation --- must be output of prep_segs() ---  # --- (Thread 8)
  stopifnot(
    "t0" %in% names(segs_prepped),
    "t1" %in% names(segs_prepped)
  )
  
  # valid_types defined locally --- no global scope dependency (Bug 4.2 fix)
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
      # dominant treated-time type --- Bug 4.3 fix: use colnames() directly
      dom_type_tx = {
        shares     <- as.matrix(dplyr::select(., dplyr::starts_with("share_") &
                                                dplyr::ends_with("_tx")))
        colnames(shares) <- gsub("^share_|_tx$", "", colnames(shares))
        colnames(shares)[max.col(shares, ties.method = "first")]
      }
    )
}
