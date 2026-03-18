# ============================================================
# MEC-TX helper: has_subseq()
# helpers/has_subseq.R
# ============================================================

#' @noRd
#'
#' Does \code{pattern} appear as a contiguous subsequence of \code{seq}?
#'
#' Used by \code{tx_focus_dt()} for treatment sequence enforcement.
#' Slides a window of length \code{length(pattern)} along \code{seq} and
#' returns \code{TRUE} on the first exact match.
#'
#' @param seq Character vector. The full sequence to search.
#' @param pattern Character vector. The contiguous subsequence to find.
#'
#' @return Logical scalar. \code{TRUE} if \code{pattern} appears as a
#'   contiguous subsequence of \code{seq}, \code{FALSE} otherwise.
#'   Edge cases: returns \code{TRUE} for zero-length \code{pattern};
#'   returns \code{FALSE} if \code{length(pattern) > length(seq)}.
#'
#' @examples
#' # TRUE  — pattern appears at position 1
#' has_subseq(c("Chemo", "IO", "Chemo"), c("Chemo", "IO"))
#'
#' # FALSE — pattern never appears contiguously
#' has_subseq(c("Chemo", "Chemo"), c("Chemo", "IO"))
has_subseq <- function(seq, pattern) {
  seq     <- as.character(seq)
  pattern <- as.character(pattern)
  n <- length(seq)
  m <- length(pattern)
  if (m == 0L) return(TRUE)
  if (m > n)   return(FALSE)
  for (i in seq_len(n - m + 1L)) {
    if (all(seq[i:(i + m - 1L)] == pattern)) return(TRUE)
  }
  FALSE
}