# ============================================================
# MEC-TX shared helper: case-insensitive column resolver
# helpers/resolve_col.R
# ============================================================

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
