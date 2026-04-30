# =============================================================================
# MEC-TX --- get_focus_cohort.R
# helpers/get_focus_cohort.R
#
# Pure data extraction --- returns sample IDs + metadata for a given
# focus_types + mode combination, without any plotting.
# Reusable helper for ad-hoc cohort building (e.g. IO resistance analysis).
# =============================================================================

#' Extract Focus Cohort Sample IDs
#'
#' Identifies patients matching a specific treatment focus and mode,
#' returning a tibble of sample IDs with associated metadata. Pure data
#' extraction --- no plotting. Designed as a reusable upstream helper for
#' cohort building before survival analysis (e.g. IO-only resistance
#' analysis, concurrent chemoradiation cohort).
#'
#' @param metadata A data frame containing survival and cluster columns
#'   for the analysis cohort. The timeline is scoped to patients present
#'   here before any filtering, ensuring the focus cohort cannot include
#'   samples excluded by upstream filters. Must contain a \code{sample}
#'   column.
#' @param timeline A long-format treatment intervals data frame --- typically
#'   the \code{$timeline_long_intv} slot from \code{\link{tx_intervals}}.
#'   Must contain columns \code{sample}, \code{type}, \code{start_year},
#'   and \code{end_year}.
#' @param focus_types Character vector of treatment types to focus on.
#'   Must match type labels in the \code{type} column of \code{timeline}
#'   (e.g. \code{c("Chemo", "Radiation")}). All specified types must be
#'   present for a patient to qualify.
#' @param mode One of \code{"only"}, \code{"concurrent"}, or
#'   \code{"dominant"}. Controls the inclusion logic --- see Details.
#'   Default \code{"only"}.
#' @param min_share Numeric in \code{(0, 1)}. Minimum treatment duration
#'   share threshold. Used in \code{"only"} mode to exclude patients with
#'   substantial non-focus treatment, and in \code{"dominant"} mode to
#'   require that focus types collectively reach this share. Default
#'   \code{0.20}.
#' @param ancillary_types Character vector of treatment types excluded from
#'   share calculations and mode filtering. Default
#'   \code{c("Ancillary", "Others")}.
#'
#' @return A tibble with one row per qualifying patient and five columns:
#'   \describe{
#'     \item{sample}{Patient identifier.}
#'     \item{focus_share}{Numeric. Combined duration share of all
#'       \code{focus_types} for this patient (sum across types, excluding
#'       ancillary denominator).}
#'     \item{mode}{Character. The \code{mode} argument used, recorded for
#'       traceability.}
#'     \item{focus_types}{Character. \code{focus_types} collapsed with
#'       \code{"+"}, recorded for traceability.}
#'     \item{n_patients}{Integer. Total number of qualifying patients
#'       (same value in every row).}
#'   }
#'   Returns a zero-row tibble (with the same columns) if no patients
#'   qualify. A message is printed reporting the focus, mode, and n.
#'
#' @details
#' \strong{Cohort scoping:} The first operation inside the function is
#' filtering \code{timeline} to patients present in \code{metadata}. This
#' ensures upstream exclusions (stage filters, duplicate removal, minimum
#' follow-up requirements) are respected before any mode logic runs.
#'
#' \strong{Mode definitions:}
#' \describe{
#'   \item{\code{"only"}}{Patient must have ALL \code{focus_types} present
#'     and must have NO non-focus type exceeding \code{min_share}. Use for
#'     pure-modality cohorts (e.g. radiation-only patients).}
#'   \item{\code{"concurrent"}}{Patient must have ALL \code{focus_types}
#'     with at least one pairwise temporal overlap between each type
#'     combination. Requires \code{length(focus_types) >= 2}; returns an
#'     empty result with a warning otherwise. Use for regimens administered
#'     simultaneously (e.g. concurrent chemoradiation --- the standard of
#'     care signal in locally advanced LUSC).}
#'   \item{\code{"dominant"}}{Patient must have ALL \code{focus_types} and
#'     their combined share must meet or exceed \code{min_share}. Less
#'     restrictive than \code{"only"} --- other treatment types may also be
#'     present. Use when focus types need not be exclusive but should
#'     dominate the treatment record.}
#' }
#'
#' \strong{Share calculation:} Duration shares exclude \code{ancillary_types}
#' from the denominator. Overlapping intervals are not merged ---
#' call \code{\link{tx_intervals}} upstream to ensure clean intervals.
#'
#' \strong{Concurrent overlap check:} Overlap is assessed pairwise across
#' all combinations of \code{focus_types}. Two intervals overlap when
#' \code{start1 < end2 & start2 < end1} (strict inequality --- touching
#' endpoints do not count).
#'
#' @examples
#' \dontrun{
#' intervals <- tx_intervals(med_data, metadata)
#'
#' # Concurrent chemoradiation cohort (LUSC standard of care)
#' crad_cohort <- get_focus_cohort(
#'   metadata    = Cluster_surv,
#'   timeline    = intervals$timeline_long_intv,
#'   focus_types = c("Chemo", "Radiation"),
#'   mode        = "concurrent"
#' )
#'
#' # IO-only cohort (no other dominant treatment)
#' io_cohort <- get_focus_cohort(
#'   metadata    = LUAD_metadata,
#'   timeline    = intervals$timeline_long_intv,
#'   focus_types = "IO",
#'   mode        = "only"
#' )
#'
#' # Join back to survival data for downstream analysis
#' io_surv <- Cluster_surv %>%
#'   dplyr::semi_join(io_cohort, by = "sample")
#' }
#'
#' @seealso \code{\link{tx_intervals}}, \code{\link{dominant_exclusive}},
#'   \code{\link{tx_pooled_analysis}}
#'
#' @importFrom dplyr filter mutate group_by summarise ungroup pull left_join
#'   intersect n_distinct n semi_join
#' @importFrom tibble tibble
#' @export
get_focus_cohort <- function(metadata,
                             timeline,
                             focus_types,
                             mode            = c("only", "concurrent", "dominant"),
                             min_share       = 0.20,
                             ancillary_types = c("Ancillary", "Others")) {
  
  mode <- match.arg(mode)
  
  # --- Scope timeline to analysis cohort (Option B fix --- Thread 8) ---
  # Ensures upstream exclusions (stage filters, deduplication, minimum
  # follow-up) are respected before any mode logic runs.
  timeline <- timeline %>%
    dplyr::filter(sample %in% metadata$sample)
  
  # --- Compute per-patient treatment shares (excluding ancillary) ---
  shares <- timeline %>%
    dplyr::filter(!type %in% ancillary_types) %>%
    dplyr::mutate(duration = end_year - start_year) %>%
    dplyr::group_by(sample, type) %>%
    dplyr::summarise(total_dur = sum(duration), .groups = "drop") %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(share = total_dur / sum(total_dur)) %>%
    dplyr::ungroup()
  
  # --- Which patients have ALL focus_types at any point ---
  has_all_focus <- shares %>%
    dplyr::filter(type %in% focus_types) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n_focus_types = dplyr::n_distinct(type),
      .groups       = "drop"
    ) %>%
    dplyr::filter(n_focus_types == length(focus_types))
  
  # --- Mode-specific filtering ---
  if (mode == "only") {
    # Must have ALL focus_types and NO other types above min_share
    non_focus_above_threshold <- shares %>%
      dplyr::filter(!type %in% focus_types, share >= min_share) %>%
      dplyr::pull(sample) %>%
      unique()
    
    ids <- has_all_focus %>%
      dplyr::filter(!sample %in% non_focus_above_threshold) %>%
      dplyr::pull(sample)
    
  } else if (mode == "concurrent") {
    # Focus types must overlap in time
    if (length(focus_types) < 2) {
      warning("get_focus_cohort: concurrent mode requires >= 2 focus_types --- returning empty")
      ids <- character(0)
    } else {
      samples_with_all <- has_all_focus$sample
      concurrent_ids   <- c()
      
      for (s in samples_with_all) {
        s_tl <- timeline %>% dplyr::filter(sample == s, type %in% focus_types)
        # Check pairwise overlap for all focus type combinations
        overlaps <- sapply(focus_types, function(ft1) {
          sapply(focus_types, function(ft2) {
            if (ft1 == ft2) return(TRUE)
            t1 <- s_tl %>% dplyr::filter(type == ft1)
            t2 <- s_tl %>% dplyr::filter(type == ft2)
            if (nrow(t1) == 0 || nrow(t2) == 0) return(FALSE)
            any(outer(
              seq_len(nrow(t1)), seq_len(nrow(t2)),
              FUN = function(i, j) {
                t1$start_year[i] < t2$end_year[j] &
                  t2$start_year[j] < t1$end_year[i]
              }
            ))
          })
        })
        if (all(overlaps)) concurrent_ids <- c(concurrent_ids, s)
      }
      ids <- concurrent_ids
    }
    
  } else if (mode == "dominant") {
    # focus_types must collectively dominate treatment share
    focus_share <- shares %>%
      dplyr::filter(type %in% focus_types) %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(focus_share = sum(share), .groups = "drop") %>%
      dplyr::filter(focus_share >= min_share)
    
    ids <- dplyr::intersect(has_all_focus$sample, focus_share$sample)
  }
  
  # --- Build result tibble with metadata ---
  result <- tibble::tibble(sample = ids) %>%
    dplyr::left_join(
      shares %>%
        dplyr::filter(type %in% focus_types) %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(focus_share = sum(share), .groups = "drop"),
      by = "sample"
    ) %>%
    dplyr::mutate(
      mode        = mode,
      focus_types = paste(focus_types, collapse = "+"),
      n_patients  = dplyr::n()
    )
  
  message(sprintf(
    "get_focus_cohort: focus=%s | mode=%s | n=%d patients",
    paste(focus_types, collapse = "+"), mode, nrow(result)
  ))
  
  result
}
