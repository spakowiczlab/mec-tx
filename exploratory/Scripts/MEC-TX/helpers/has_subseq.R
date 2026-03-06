# ============================================================
# MEC-TX helper: has_subseq()
# helpers/has_subseq.R
# ============================================================

# Does `pattern` appear as a contiguous subsequence of `seq`?
# Used by tx_focus_dt() for sequence enforcement.
#
# Examples:
#   has_subseq(c("Chemo","IO","Chemo"), c("Chemo","IO"))  → TRUE
#   has_subseq(c("Chemo","Chemo"),      c("Chemo","IO"))  → FALSE

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
