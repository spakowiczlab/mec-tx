# =============================================================================
# MEC-TX --- standardise_status.R
# helpers/standardise_status.R
#
# Standardise survival status column to 0/1 integer.
# Detects: numeric 0/1, character "dead"/"alive" (case-insensitive),
#          "deceased"/"living", "censored", etc.
# Convention: 0 = alive/censored, 1 = dead/event.
#
# Output columns:
#   status       --- integer 0/1 (for survival::Surv())
#   status_label --- factor "Alive"/"Dead" (for tables, plots, summaries)
# =============================================================================

#' Standardise Survival Status Column to 0/1 Integer
#'
#' Detects the coding scheme of a survival status column and converts it
#' to the MEC-TX convention: 0 = alive/censored, 1 = dead/event (integer).
#' Adds a \code{status_label} factor column for display use in tables and
#' plots. Issues an informative message reporting what was detected and how
#' it was mapped.
#'
#' @param df A data frame containing a survival status column.
#' @param status_col Character string. Name of the status column to
#'   standardise. Case-sensitive. Default \code{"status"}.
#'
#' @return The input data frame with two modifications:
#'   \describe{
#'     \item{status}{The status column standardised to integer 0/1 and
#'       renamed to \code{"status"} regardless of the input
#'       \code{status_col} name. Safe to pass directly to
#'       \code{survival::Surv()}. Convention: 0 = alive/censored,
#'       1 = dead/event.}
#'     \item{status_label}{A new factor column with levels
#'       \code{c("Alive", "Dead")} added for use in plots, tables, and
#'       summaries. Do not use this column in \code{survival::Surv()} ---
#'       use the integer \code{status} column instead.}
#'   }
#'
#' @details
#' \strong{Detection logic:} Two coding paths are handled:
#' \describe{
#'   \item{Numeric 0/1}{If the column is numeric and all non-NA values are
#'     in \code{\{0, 1\}}, no recoding is performed --- the column is coerced
#'     to integer, renamed to \code{"status"}, and \code{status_label} is
#'     added.}
#'   \item{Character / factor}{Values are lowercased and matched against
#'     known dead patterns (\code{"dead"}, \code{"deceased"}, \code{"died"},
#'     \code{"death"}, \code{"1"}) and alive patterns (\code{"alive"},
#'     \code{"living"}, \code{"censored"}, \code{"censor"}, \code{"0"}).
#'     Any value matching a dead pattern maps to 1; all others map to 0.}
#' }
#' If no known patterns are detected the function stops with an informative
#' error listing the values found, prompting manual recoding upstream.
#'
#' \strong{Column rename:} The output column is always named \code{"status"}
#' (lowercase) regardless of \code{status_col}. This ensures compatibility
#' with all downstream MEC-TX functions which hardcode \code{"status"} in
#' \code{survival::Surv()} calls.
#'
#' \strong{AVATAR-specific note:} The LUSC \code{Cluster_surv} object uses
#' lowercase \code{"status"} while LUAD \code{LUAD_metadata} uses
#' capitalised \code{"Status"}. Pass \code{status_col = "Status"} for LUAD
#' --- the output column will be renamed to lowercase \code{"status"}
#' automatically.
#'
#' \strong{status_label warning:} \code{status_label} is a convenience
#' column for display only. Never pass it to \code{survival::Surv()} ---
#' this will produce incorrect results.
#'
#' @examples
#' df <- data.frame(sample = 1:3, diagsurvtime = c(1.2, 3.4, 2.1),
#'                  status = c(0, 1, 1))
#' standardise_status(df)
#'
#' @seealso \code{\link{tx_normalize}}, \code{resolve_col}
#'
#' @importFrom stats na.omit
#' @export
standardise_status <- function(df, status_col = "status") {
  
  # --- Column must exist ---
  if (!status_col %in% names(df)) {
    stop(sprintf(
      "standardise_status: column '%s' not found. Available columns: %s",
      status_col,
      paste(names(df), collapse = ", ")
    ))
  }
  
  raw      <- df[[status_col]]
  raw_vals <- sort(unique(na.omit(raw)))
  
  # --- Already 0/1 numeric -> just ensure integer type ---
  if (is.numeric(raw) && all(raw_vals %in% c(0, 1))) {
    message(sprintf(
      "standardise_status: '%s' already 0/1 numeric --- no change (0=alive, 1=dead)",
      status_col
    ))
    df[[status_col]] <- as.integer(raw)
    df$status_label  <- factor(
      ifelse(df[[status_col]] == 1, "Dead", "Alive"),
      levels = c("Alive", "Dead")
    )
    names(df)[names(df) == status_col] <- "status"  # --- rename fix (Thread 8)
    return(df)
  }
  
  # --- Character / factor coding -> detect pattern ---
  raw_lower    <- tolower(as.character(raw))
  unique_lower <- sort(unique(na.omit(raw_lower)))
  
  dead_patterns  <- c("dead", "deceased", "died", "death", "1")
  alive_patterns <- c("alive", "living", "censored", "censor", "0")
  
  detected_dead  <- unique_lower[unique_lower %in% dead_patterns]
  detected_alive <- unique_lower[unique_lower %in% alive_patterns]
  
  if (length(detected_dead) == 0 && length(detected_alive) == 0) {
    stop(sprintf(
      paste0(
        "standardise_status: cannot auto-detect coding for '%s'.\n",
        "  Found values: %s\n",
        "  Please recode manually to 0 (alive) / 1 (dead) before calling tx_normalize()."
      ),
      status_col,
      paste(unique_lower, collapse = ", ")
    ))
  }
  
  message(sprintf(
    "standardise_status: '%s' --- detected '%s'=dead, '%s'=alive -> recoding to 1/0",
    status_col,
    paste(detected_dead,  collapse = "/"),
    paste(detected_alive, collapse = "/")
  ))
  
  df[[status_col]] <- as.integer(raw_lower %in% dead_patterns)
  df$status_label  <- factor(
    ifelse(df[[status_col]] == 1, "Dead", "Alive"),
    levels = c("Alive", "Dead")
  )
  names(df)[names(df) == status_col] <- "status"    # --- rename fix (Thread 8)
  df
}
