# =============================================================================
# tx_duration.R
# MEC-TX Package
#
# Compare time-on-treatment by a grouping variable (e.g., CAlevel) within
# each treatment type. Measures total calendar-time exposure per patient
# per type and as a merged total across all types.
#
# Consumes tx_intervals() output directly. No dependency on tx_normalize().
# =============================================================================

#' @noRd
#'
#' Merge overlapping intervals into non-overlapping unions
#'
#' Sorts intervals by start time and greedily merges any that overlap or
#' touch. Used internally by \code{.duration_per_type()} and
#' \code{.duration_total()} to avoid double-counting concurrent treatment.
#'
#' @param starts Numeric vector of interval start times (years).
#' @param ends Numeric vector of interval end times (years). Must be the
#'   same length as \code{starts}.
#'
#' @return A data frame with columns \code{start} and \code{end} giving
#'   the merged non-overlapping intervals. Returns a zero-row data frame
#'   if \code{starts} is empty.
.merge_intervals <- function(starts, ends) {
  if (length(starts) == 0L) {
    return(data.frame(start = numeric(0), end = numeric(0)))
  }
  ord    <- order(starts, ends)
  starts <- starts[ord]
  ends   <- ends[ord]
  
  ms <- starts[1]
  me <- ends[1]
  rs <- numeric()
  re <- numeric()
  
  for (i in seq_along(starts)) {
    if (starts[i] <= me) {
      me <- max(me, ends[i])
    } else {
      rs <- c(rs, ms)
      re <- c(re, me)
      ms <- starts[i]
      me <- ends[i]
    }
  }
  rs <- c(rs, ms)
  re <- c(re, me)
  data.frame(start = rs, end = re)
}


#' @noRd
#'
#' Compute treatment duration per patient per type
#'
#' For each patient × treatment-type combination, merges overlapping
#' intervals within that type via \code{.merge_intervals()}, then sums
#' calendar time. A patient receiving two overlapping chemo courses counts
#' the union once, not twice.
#'
#' @param timeline A data frame from \code{\link{tx_intervals}} containing
#'   per-interval records.
#' @param sample_col Column name for patient ID. Default \code{"sample"}.
#' @param type_col Column name for treatment type. Default \code{"type"}.
#' @param start_col Column name for interval start (years).
#'   Default \code{"start_year"}.
#' @param end_col Column name for interval end (years).
#'   Default \code{"end_year"}.
#'
#' @return A data frame with columns \code{sample}, \code{type}, and
#'   \code{duration_yrs}. One row per patient × type combination that has
#'   at least one interval. Patient × type combinations with no intervals
#'   are omitted (not zero-filled).
.duration_per_type <- function(timeline,
                               sample_col = "sample",
                               type_col   = "type",
                               start_col  = "start_year",
                               end_col    = "end_year") {
  
  patients <- unique(timeline[[sample_col]])
  types    <- unique(timeline[[type_col]])
  
  out <- vector("list", length(patients) * length(types))
  idx <- 0L
  
  for (pid in patients) {
    sub_p <- timeline[timeline[[sample_col]] == pid, , drop = FALSE]
    for (tx in types) {
      sub_pt <- sub_p[sub_p[[type_col]] == tx, , drop = FALSE]
      if (nrow(sub_pt) == 0L) next
      merged <- .merge_intervals(sub_pt[[start_col]], sub_pt[[end_col]])
      dur    <- sum(merged$end - merged$start)
      idx    <- idx + 1L
      out[[idx]] <- data.frame(
        sample           = pid,
        type             = tx,
        duration_yrs     = dur,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out[seq_len(idx)])
}


#' @noRd
#'
#' Compute merged total treatment duration per patient
#'
#' Unions all treatment intervals across all types for each patient,
#' merges overlapping periods via \code{.merge_intervals()}, and sums
#' calendar time. A patient receiving concurrent chemo + IO for 6 months
#' contributes 6 months total, not 12.
#'
#' @param timeline A data frame from \code{\link{tx_intervals}}.
#' @param sample_col Column name for patient ID. Default \code{"sample"}.
#' @param start_col Column name for interval start (years).
#'   Default \code{"start_year"}.
#' @param end_col Column name for interval end (years).
#'   Default \code{"end_year"}.
#'
#' @return A data frame with columns \code{sample} and
#'   \code{duration_yrs_total}. One row per patient.
.duration_total <- function(timeline,
                            sample_col = "sample",
                            start_col  = "start_year",
                            end_col    = "end_year") {
  patients <- unique(timeline[[sample_col]])
  out <- data.frame(
    sample             = patients,
    duration_yrs_total = NA_real_,
    stringsAsFactors   = FALSE
  )
  for (i in seq_along(patients)) {
    sub    <- timeline[timeline[[sample_col]] == patients[i], , drop = FALSE]
    merged <- .merge_intervals(sub[[start_col]], sub[[end_col]])
    out$duration_yrs_total[i] <- sum(merged$end - merged$start)
  }
  out
}


#' Treatment Duration Analysis by Grouping Variable
#'
#' Computes per-patient treatment duration (calendar time) broken down by
#' treatment type and as a merged total across all types. Compares groups
#' using Wilcoxon rank-sum tests (2 groups) or Kruskal-Wallis tests (>2
#' groups) and produces faceted box or violin plots. Concurrent treatment
#' periods are counted once in the merged total — not double-counted.
#'
#' @param timeline A data frame — the direct output of
#'   \code{\link{tx_intervals}}. Must contain columns named in
#'   \code{sample_col}, \code{type_col}, \code{start_col}, and
#'   \code{end_col}.
#' @param meta A data frame of patient-level metadata (e.g.
#'   \code{Cluster_surv} or \code{LUAD_metadata}). Must contain columns
#'   named in \code{sample_col} and \code{group_var}.
#' @param group_var Character string. Name of the grouping column in
#'   \code{meta} (e.g. \code{"CAlevel"}). Must be an exact column name
#'   match.
#' @param sample_col Character string. Patient identifier column name in
#'   both \code{timeline} and \code{meta}. Default \code{"sample"}.
#' @param type_col Character string. Treatment type column name in
#'   \code{timeline}. Default \code{"type"}.
#' @param start_col Character string. Interval start column name in
#'   \code{timeline} (years). Default \code{"start_year"}.
#' @param end_col Character string. Interval end column name in
#'   \code{timeline} (years). Default \code{"end_year"}.
#' @param duration_unit One of \code{"months"} (default) or \code{"years"}.
#'   Controls units in output tables and plot axis labels. Months =
#'   duration in years × 12.
#' @param exclude_types Character vector or \code{NULL}. Treatment types
#'   to exclude from the analysis before computing durations. Default
#'   \code{NULL}.
#' @param min_n Integer. Minimum number of patients per group × type
#'   required to run a statistical test. Types falling below this threshold
#'   are flagged in \code{test_note} but retained in the output.
#'   Default \code{3}.
#' @param plot Logical. If \code{TRUE}, produce a faceted plot.
#'   Default \code{TRUE}.
#' @param plot_type One of \code{"box"} (default) or \code{"violin"}.
#'   Controls the geom used in the plot.
#' @param title Character string or \code{NULL}. Plot title. \code{NULL}
#'   auto-generates \code{"Treatment Duration by <group_var>"}.
#'   Default \code{NULL}.
#' @param palette Named character vector or \code{NULL}. Colours keyed by
#'   group level names. \code{NULL} auto-generates from
#'   \code{grDevices::hcl.colors()}. Default \code{NULL}.
#'
#' @return A named list with five elements:
#'   \describe{
#'     \item{duration_per_type}{Data frame with one row per patient ×
#'       treatment type combination. Columns: \code{sample}, \code{type},
#'       \code{<group_var>}, \code{duration_<duration_unit>}.}
#'     \item{duration_total}{Data frame with one row per patient. Columns:
#'       \code{sample}, \code{<group_var>},
#'       \code{duration_total_<duration_unit>}. Concurrent treatment
#'       periods are merged before summing.}
#'     \item{summary_table}{Data frame with one row per type × group
#'       combination. Columns: \code{type}, \code{group}, \code{n},
#'       \code{mean}, \code{median}, \code{q25}, \code{q75},
#'       \code{p_value} (Wilcoxon or Kruskal-Wallis), \code{test_note}
#'       (reason if test was skipped).}
#'     \item{plot}{A \code{ggplot} object — faceted by treatment type with
#'       p-value annotations. \code{NULL} if \code{plot = FALSE}.}
#'     \item{params}{Named list recording the call parameters:
#'       \code{group_var}, \code{duration_unit}, \code{exclude_types},
#'       \code{min_n}, \code{n_patients}, \code{n_types}.}
#'   }
#'
#' @details
#' \strong{Overlap handling:} Durations are computed on merged
#' non-overlapping intervals. This means a patient receiving two
#' overlapping chemo regimens for 6 months each (3 months concurrent)
#' accumulates 9 months of chemo, not 12. The merged total similarly
#' counts concurrent multi-type treatment once.
#'
#' \strong{Statistical tests:} For exactly 2 groups, Wilcoxon rank-sum
#' with \code{exact = FALSE} (handles ties). For >2 groups,
#' Kruskal-Wallis. Tests are skipped when any group has fewer than
#' \code{min_n} patients for that type — recorded in \code{test_note}.
#'
#' \strong{Output column naming:} Duration columns in the returned data
#' frames are named \code{duration_months} or \code{duration_years}
#' depending on \code{duration_unit}, making it unambiguous when saved
#' to disk or joined to other tables.
#'
#' \strong{No tx_normalize() dependency:} This function consumes
#' \code{tx_intervals()} output directly and does not require the
#' normalised timeline from \code{tx_normalize()}.
#'
#' @examples
#' \dontrun{
#' intervals <- tx_intervals(med_data, metadata)
#'
#' # Default — months, CAlevel comparison, box plot
#' res <- tx_duration(
#'   timeline  = intervals$timeline_long_intv,
#'   meta      = Cluster_surv,
#'   group_var = "CAlevel"
#' )
#' res$summary_table
#'
#' pdf(file.path(out_dir, "tx_duration_lusc.pdf"), width = 14, height = 10)
#' print(res$plot)
#' dev.off()
#'
#' # Violin plot, exclude Others, report in years
#' res2 <- tx_duration(
#'   timeline      = intervals$timeline_long_intv,
#'   meta          = LUAD_metadata,
#'   group_var     = "CAlevel",
#'   duration_unit = "years",
#'   exclude_types = "Others",
#'   plot_type     = "violin"
#' )
#' }
#'
#' @seealso \code{\link{tx_intervals}}, \code{\link{tx_lines}},
#'   \code{\link{tx_pooled_analysis}}
#'
#' @import ggplot2
#' @importFrom stats wilcox.test kruskal.test median quantile
#' @importFrom grDevices hcl.colors
#' @export