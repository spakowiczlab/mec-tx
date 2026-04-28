# ============================================================
# MEC-TX shared helper: case-insensitive column resolver
# helpers/resolve_col.R
# ============================================================

#' Case-insensitive column name resolver
#' Looks up a column in a data frame by name using case-insensitive
#' matching and returns the actual column name as it exists in the data
#' frame. Errors with an informative message if the column is not found
#' or if the name matches multiple columns after case folding.
#' Used throughout the package to handle inconsistent capitalisation in
#' registry data (e.g. \code{"Status"} vs \code{"status"}, \code{"Sample"}
#' vs \code{"sample"}).
#' @param df A data frame to search.
#' @param col_name Character string. The column name to look up
#'   (case-insensitive).
#' @param arg_name Character string. The argument name to report in error
#'   messages, for clearer diagnostics at the call site. Defaults to
#'   \code{col_name}.
#' @return Character scalar. The exact column name as it appears in
#'   \code{names(df)}.
#' @examples
#' df <- data.frame(Sample = 1:3, Status = 0:2)
#' resolve_col(df, "sample")   # returns "Sample"
#' resolve_col(df, "STATUS")   # returns "Status"
#' @noRd
resolve_col <- function(df, col_name, arg_name = col_name) {
  actual <- names(df)[tolower(names(df)) == tolower(col_name)]
  if (length(actual) == 0) {
    stop(sprintf(
      "Column '%s' (arg: %s) not found. Available: %s",
      col_name, arg_name, paste(names(df), collapse = ", ")
    ))
  }
  if (length(actual) > 1) {
    stop(sprintf(
      "Column '%s' (arg: %s) matches multiple columns after case folding: %s",
      col_name, arg_name, paste(actual, collapse = ", ")
    ))
  }
  actual
}
