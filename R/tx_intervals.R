#' Convert Normalised Treatment Timeline to Treatment Intervals
#'
#' Takes the long-format normalised treatment timeline produced by
#' \code{\link{tx_normalize}} and converts it to a compact interval
#' representation --- one row per contiguous treatment run per patient per
#' type. Intervals are defined on a regular time grid whose resolution
#' must match the \code{grid_weeks} used in \code{\link{tx_normalize}}.
#'
#' @param df A long-format data frame --- the direct output of
#'   \code{\link{tx_normalize}}. Must contain columns named in
#'   \code{id_col}, \code{time_col}, and \code{type_col}.
#' @param id_col Character string. Name of the patient identifier column.
#'   Default \code{"sample"}.
#' @param time_col Character string. Name of the time column (years since
#'   treatment start). Must be numeric. Default
#'   \code{"TimeSinceTreatmentStart"}.
#' @param type_col Character string. Name of the treatment type column.
#'   Default \code{"treatment_group"}.
#' @param drop_types Character vector. Treatment type labels to exclude
#'   before building intervals. Default \code{c("None")} drops untreated
#'   time bins. Set to \code{character(0)} to retain all types.
#' @param horizon_years Numeric. Maximum follow-up horizon in years.
#'   Time points beyond this value are clamped to the horizon before
#'   interval construction. Must match the \code{horizon_years} used in
#'   \code{\link{tx_normalize}}. Default \code{5}.
#' @param grid_weeks Numeric. Time bin width in weeks. Must match the
#'   \code{grid_weeks} used in \code{\link{tx_normalize}} --- mismatches
#'   produce incorrect interval boundaries. Common values: \code{1}
#'   (weekly), \code{2} (biweekly), \code{4} (monthly, default).
#'   Default \code{4}.
#'
#' @return A tibble with one row per contiguous treatment run and six
#'   columns:
#'   \describe{
#'     \item{sample}{Patient identifier.}
#'     \item{type}{Treatment type label.}
#'     \item{block}{Integer. Active-therapy block index per patient.
#'       A new block starts whenever there is a gap of one or more grid
#'       points with no recorded treatment of any type.}
#'     \item{run}{Integer. Contiguous run index within each
#'       \code{sample -- block -- type} combination.}
#'     \item{start_year}{Numeric. Interval start time in years
#'       (grid-aligned, inclusive).}
#'     \item{end_year}{Numeric. Interval end time in years
#'       (grid-aligned, end-exclusive --- i.e. the start of the next bin).}
#'   }
#'   Returns a zero-row tibble with the same columns if no rows remain
#'   after filtering.
#'
#' @details
#' \strong{Grid alignment:} Time values are snapped to the nearest grid
#' point using \code{round(t_year * grid_res)} where
#' \code{grid_res = 52 / grid_weeks}. This ensures interval boundaries
#' align exactly with the bins created by \code{\link{tx_normalize}}.
#'
#' \strong{Block vs run:} A \emph{block} groups all treatment types
#' active within a contiguous set of grid points (no gap). A \emph{run}
#' is a contiguous sequence of the same treatment type within a block.
#' The block structure captures treatment pauses; the run structure
#' separates types within a block.
#'
#' \strong{End-exclusive intervals:} \code{end_year} is set to
#' \code{(max_grid_index + 1) / grid_res}, making intervals
#' end-exclusive. This convention is consistent with standard interval
#' arithmetic and avoids zero-duration segments.
#'
#' \strong{grid_weeks alignment:} If \code{grid_weeks} does not match
#' the value used in \code{\link{tx_normalize}}, interval boundaries will
#' be misaligned and durations will be incorrect. Always pass the same
#' value to both functions.
#'
#' @examples
#' \dontrun{
#' # Standard pipeline
#' norm      <- tx_normalize(med_data, metadata)
#' intervals <- tx_intervals(norm)
#'
#' # Access interval data frame
#' head(intervals)
#'
#' # Biweekly grid --- must match tx_normalize() call
#' norm2      <- tx_normalize(med_data, metadata, grid_weeks = 2)
#' intervals2 <- tx_intervals(norm2, grid_weeks = 2)
#'
#' # Retain None type (untreated bins)
#' intervals_with_none <- tx_intervals(norm, drop_types = character(0))
#' }
#'
#' @seealso \code{\link{tx_normalize}}, \code{tx_duration},
#'   \code{\link{tx_lines}}, \code{\link{tx_pooled_analysis}},
#'   \code{\link{dominant_exclusive}}
#'
#' @importFrom dplyr transmute filter mutate distinct arrange group_by
#'   ungroup left_join summarise lag
#' @importFrom tidyr replace_na
#' @importFrom tibble tibble
#' @export
tx_intervals <- function(
    df,
    id_col        = "sample",
    time_col      = "TimeSinceTreatmentStart",
    type_col      = "treatment_group",
    drop_types    = c("None"),
    horizon_years = 5,
    grid_weeks    = 4
) {
  
  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
  
  # --- 1. df must be a data frame ---
  if (!is.data.frame(df)) {
    stop(
      "tx_intervals(): 'df' must be a data frame.\n",
      "  --- You passed an object of class: ", paste(class(df), collapse = ", "), ".\n",
      "  --- Pass the output of tx_normalize() directly:\n",
      "    tx_intervals(tx_normalize(your_raw_data))"
    )
  }
  
  # --- 2. Required columns must exist ---
  for (col in c(id_col, time_col, type_col)) {
    if (!col %in% names(df)) {
      stop(
        "tx_intervals(): Column '", col, "' not found in 'df'.\n",
        "  --- Columns present: ", paste(names(df), collapse = ", "), "\n",
        "  --- The default expected columns are:\n",
        "      id_col   = 'sample'\n",
        "      time_col = 'TimeSinceTreatmentStart'\n",
        "      type_col = 'treatment_group'\n",
        "  --- If your columns have different names, pass them explicitly:\n",
        "    tx_intervals(df, id_col='your_id', time_col='your_time', type_col='your_type')"
      )
    }
  }
  
  # --- 3. time_col must be numeric ---
  if (!is.numeric(df[[time_col]])) {
    stop(
      "tx_intervals(): Column '", time_col, "' must be numeric.\n",
      "  --- Current class: ", class(df[[time_col]]), "\n",
      "  --- This column should contain time in years since treatment start.\n",
      "  --- If you used tx_normalize(), this column is created automatically."
    )
  }
  
  # --- 4. grid_weeks must be a positive number ---
  if (!is.numeric(grid_weeks) || length(grid_weeks) != 1 || grid_weeks <= 0) {
    stop(
      "tx_intervals(): 'grid_weeks' must be a single positive number.\n",
      "  --- You passed: ", deparse(grid_weeks), "\n",
      "  --- This must match the grid_weeks used in tx_normalize().\n",
      "  --- Default is 4 (monthly). Use 2 for biweekly."
    )
  }
  
  # --- 5. horizon_years must be a positive number ---
  if (!is.numeric(horizon_years) || length(horizon_years) != 1 || horizon_years <= 0) {
    stop(
      "tx_intervals(): 'horizon_years' must be a single positive number.\n",
      "  --- You passed: ", deparse(horizon_years), "\n",
      "  --- Example: horizon_years = 5 (default)"
    )
  }
  
  # --- 6. Warn if drop_types removes all rows ---
  all_types  <- unique(as.character(df[[type_col]]))
  kept_types <- setdiff(all_types, drop_types)
  dropped    <- intersect(all_types, drop_types)
  
  if (length(dropped) > 0) {
    message(
      "tx_intervals(): Dropping treatment type(s): ",
      paste(dropped, collapse = ", ")
    )
  }
  
  if (length(kept_types) == 0) {
    stop(
      "tx_intervals(): All treatment types would be dropped by 'drop_types'.\n",
      "  --- drop_types = ", paste0('"', drop_types, '"', collapse = ", "), "\n",
      "  --- Types in your data: ", paste(all_types, collapse = ", "), "\n",
      "  --- Adjust drop_types so at least one treatment type is retained."
    )
  }
  
  # --- 7. Warn if time_col has all-NA values ---
  if (all(is.na(df[[time_col]]))) {
    stop(
      "tx_intervals(): Column '", time_col, "' is entirely NA.\n",
      "  --- This column must contain numeric time values in years.\n",
      "  --- Check that tx_normalize() completed successfully before calling tx_intervals()."
    )
  }
  
  # ===========================================================================
  # MAIN FUNCTION BODY
  # ===========================================================================
  
  # grid_res: number of grid points per year (e.g. 13 for 4-week grid)
  grid_res <- 52L / grid_weeks
  
  # ---- 1) Clean + snap to grid ----
  dat <- df %>%
    dplyr::transmute(
      sample     = .data[[id_col]],
      t_year_raw = as.numeric(.data[[time_col]]),
      type       = as.character(.data[[type_col]])
    ) %>%
    dplyr::filter(!is.na(sample), !is.na(t_year_raw), !is.na(type)) %>%
    dplyr::filter(!(type %in% drop_types)) %>%
    dplyr::mutate(
      m = as.integer(round(t_year_raw * grid_res)),
      m = pmax(0L, pmin(m, as.integer(horizon_years * grid_res)))
    ) %>%
    dplyr::distinct(sample, m, type)
  
  if (!nrow(dat)) {
    return(tibble::tibble(
      sample     = character(),
      type       = character(),
      block      = integer(),
      run        = integer(),
      start_year = double(),
      end_year   = double()
    ))
  }
  
  # ---- 2) Define active-therapy blocks (no grid gaps) ----
  grid_blocks <- dat %>%
    dplyr::distinct(sample, m) %>%
    dplyr::arrange(sample, m) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(
      new_block = m != dplyr::lag(m, default = dplyr::first(m)) + 1L,
      block     = cumsum(tidyr::replace_na(new_block, TRUE))
    ) %>%
    dplyr::ungroup()
  
  dat2 <- dat %>%
    dplyr::left_join(grid_blocks, by = c("sample", "m")) %>%
    dplyr::arrange(sample, block, m, type)
  
  # ---- 3) Merge consecutive grid points per type within each block ----
  segs <- dat2 %>%
    dplyr::group_by(sample, block, type) %>%
    dplyr::arrange(m, .by_group = TRUE) %>%
    dplyr::mutate(
      new_run = m != dplyr::lag(m, default = dplyr::first(m)) + 1L,
      run     = cumsum(tidyr::replace_na(new_run, TRUE))
    ) %>%
    dplyr::group_by(sample, block, type, run) %>%
    dplyr::summarise(
      start_year = min(m) / grid_res,
      end_year   = (max(m) + 1L) / grid_res,
      .groups    = "drop"
    ) %>%
    dplyr::arrange(sample, block, start_year, type)
  
  segs
}
