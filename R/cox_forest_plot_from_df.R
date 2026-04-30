#' Cox Proportional Hazards Forest Plot
#'
#' Fits an adjusted Cox model on a pre-filtered survival dataset and renders
#' a publication-ready forest plot. Handles reference-level encoding,
#' case-insensitive column resolution, numeric covariate rescaling, and
#' events-per-variable (EPV) driven covariate selection automatically.
#' Degrades gracefully --- returns an annotated blank plot rather than
#' erroring when data are insufficient.
#'
#' @param df A data frame containing at minimum \code{diagsurvtime} (numeric,
#'   years), \code{status} (integer, 0/1), and all columns named in
#'   \code{covars}. Typically the \code{$data} slot from
#'   \code{\link{tx_pooled_analysis}} or \code{\link{tx_compare_groups}}.
#'   Column name matching is case-insensitive.
#' @param covars Character vector of covariate column names to include in
#'   the model. The first element is always retained regardless of EPV.
#'   Default: \code{c("CAlevel", "stage_group", "sex", "age", "smokingstatus")}.
#' @param ref_levels Named list specifying the reference level for each
#'   categorical covariate. Keys are matched case-insensitively to column
#'   names. Covariates not listed retain their existing factor ordering.
#'   Default: \code{list(CAlevel = "Low", stage_group = "Local",
#'   smokingstatus = "Never")}.
#' @param title Character string. Plot title. Default \code{"Adjusted Cox"}.
#' @param min_epv Integer. Minimum events-per-variable threshold used to
#'   cap the number of model parameters at \code{floor(events / min_epv)}.
#'   Covariates are dropped in reverse \code{priority} order until the cap
#'   is satisfied. The first covariate is always retained. Default \code{5}.
#' @param priority Character vector or \code{NULL}. Covariates listed first
#'   are retained preferentially when EPV forces dropping. \code{NULL}
#'   uses the \code{covars} order as-is.
#' @param numeric_scale Named list mapping numeric covariate names to a
#'   divisor applied before modelling. The reported HR reflects one unit of
#'   the rescaled variable. Default \code{list(age = 10)} gives HR per
#'   10-year increment.
#' @param numeric_units Named list mapping numeric covariate names to a
#'   unit label used in axis annotations. Default \code{list(age = "years")}.
#' @param numeric_pretty Named list or \code{NULL}. Display name overrides
#'   for numeric covariates on the forest plot axis. \code{NULL} auto-generates
#'   labels from \code{numeric_scale} and \code{numeric_units}.
#' @param base_size Base font size (pt) passed to
#'   \code{ggplot2::theme_minimal()}. Default \code{16}.
#' @param title_size Font size (pt) for the plot title. Default \code{22}.
#' @param axis_title_size Font size (pt) for axis titles. Default \code{14}.
#' @param axis_text_size Font size (pt) for axis tick and covariate labels.
#'   Default \code{12}.
#' @param bold Logical. If \code{TRUE} all text elements use bold weight.
#'   Default \code{TRUE}.
#'
#' @return A \code{ggplot} object. The subtitle reports complete-case n,
#'   event count, and the covariates actually fitted. Returns an annotated
#'   blank plot (still a \code{ggplot}) if columns are missing, fewer than
#'   10 complete cases exist, fewer than 3 events are observed, or the Cox
#'   model fails to converge.
#'
#' @details
#' \strong{Column requirements:} \code{df} must contain \code{diagsurvtime}
#' (time from diagnosis in years) and \code{status} (0 = censored, 1 = event).
#' These names are hardcoded; rename columns upstream if necessary.
#'
#' \strong{EPV selection:} Parameters are budgeted as
#' \code{floor(events / min_epv)}. Categorical covariates consume
#' \code{(n_levels - 1)} parameters each. The first element of \code{covars}
#' (typically \code{CAlevel}) is always included regardless of budget.
#'
#' \strong{Convergence fallback:} If the full selected model fails, covariates
#' are dropped one at a time from the end until the model converges.
#'
#' @examples
#' \dontrun{
#' # Standard usage after tx_pooled_analysis()
#' p <- cox_forest_plot_from_df(
#'   df    = pooled_res$data,
#'   title = "Adjusted Cox --- LUSC Chemoradiation"
#' )
#' pdf(file.path(tempdir(), "cox_forest_lusc.pdf"), width = 10, height = 8)
#' print(p)
#' dev.off()
#'
#' # Override covariate priority and reference levels
#' p2 <- cox_forest_plot_from_df(
#'   df         = pooled_res$data,
#'   priority   = c("CAlevel", "stage_group"),
#'   ref_levels = list(CAlevel = "Low", stage_group = "Local",
#'                     smokingstatus = "Never"),
#'   title      = "Adjusted Cox --- LUAD Chemo + IO"
#' )
#' }
#'
#' @seealso \code{\link{tx_pooled_analysis}}, \code{\link{tx_compare_groups}},
#'   \code{\link{km_panel_from_df}}
#'
#' @import ggplot2
#' @importFrom survival coxph Surv
#' @importFrom broom tidy
#' @importFrom dplyr filter mutate n_distinct
#' @importFrom stats relevel na.omit var as.formula
#' @export
# -------- Cox forest plot --------
cox_forest_plot_from_df <- function(
  df,
  covars          = c("CAlevel", "stage_group", "sex", "age", "smokingstatus"),
  ref_levels      = list(CAlevel = "Low", stage_group = "Local",
                         smokingstatus = "Never"),
  title           = "Adjusted Cox",
  min_epv         = 5,
  priority        = NULL,   # NULL = use covars order as-is
  numeric_scale   = list(age = 10),
  numeric_units   = list(age = "years"),
  numeric_pretty  = NULL,
  base_size       = 16,
  title_size      = 22,
  axis_title_size = 14,
  axis_text_size  = 12,
  bold            = TRUE
) {

  # ---- resolve all covariate column names case-insensitively (Bug 4.4 fix) ----
  covars <- Filter(Negate(is.null), lapply(covars, function(cv) {
    tryCatch(resolve_col(df, cv, cv), error = function(e) {
      warning(sprintf("covar '%s' not found --- skipping.", cv))
      NULL
    })
  }))
  covars <- unlist(covars)

  # resolve ref_levels keys to actual column names
  names(ref_levels) <- vapply(names(ref_levels), function(nm) {
    tryCatch(resolve_col(df, nm, nm), error = function(e) nm)
  }, character(1))

  need <- unique(c("diagsurvtime", "status", covars))
  missing_cols <- setdiff(need, names(df))
  if (length(missing_cols)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = .5, y = .5,
                          label = paste("Missing columns:", paste(missing_cols, collapse = ", "))) +
        ggplot2::labs(title = title) +
        ggplot2::theme_void(base_size = base_size)
    )
  }

  d0 <- df[, need, drop = FALSE]

  # apply ref levels
  for (nm in names(ref_levels)) {
    if (!is.null(ref_levels[[nm]]) && nm %in% names(d0)) {
      d0[[nm]] <- factor(d0[[nm]])
      if (ref_levels[[nm]] %in% levels(d0[[nm]])) {
        d0[[nm]] <- stats::relevel(d0[[nm]], ref = ref_levels[[nm]])
      }
    }
  }

  # coerce numeric columns
  safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  for (nm in names(numeric_scale)) {
    if (nm %in% names(d0) && !is.numeric(d0[[nm]])) d0[[nm]] <- safe_num(d0[[nm]])
    if (nm %in% names(d0) && is.numeric(d0[[nm]]) && numeric_scale[[nm]] > 0) {
      d0[[nm]] <- d0[[nm]] / numeric_scale[[nm]]
    }
  }

  d  <- stats::na.omit(d0)
  n  <- nrow(d)
  ev <- sum(d$status == 1, na.rm = TRUE)

  if (n < 10 || ev < 3) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = .5, y = .6,
                          label = sprintf("n=%d, events=%d", n, ev)) +
        ggplot2::annotate("text", x = .5, y = .4,
                          label = "Too few complete cases") +
        ggplot2::labs(title = title) +
        ggplot2::theme_void(base_size = base_size)
    )
  }

  # EPV-based covariate selection
  varies <- function(v) {
    x <- d[[v]]
    if (is.numeric(x)) var(x, na.rm = TRUE) > 0 else dplyr::n_distinct(x) > 1
  }
  cov_keep   <- covars[sapply(covars, varies)]
  max_params <- max(1L, floor(ev / min_epv))

  param_count <- function(v) {
    if (is.numeric(d[[v]])) 1L else max(1L, dplyr::n_distinct(d[[v]]) - 1L)
  }

  if (is.null(priority)) priority <- cov_keep
  cov_keep <- unique(c(intersect(priority, cov_keep), setdiff(cov_keep, priority)))

  pick <- character(0); used <- 0L
  for (v in cov_keep) {
    pv <- as.integer(param_count(v))
    if (used + pv <= max_params) { pick <- c(pick, v); used <- used + pv }
  }
  # always include first covariate (usually the variable of interest)
  if (!covars[1] %in% pick) pick <- c(covars[1], setdiff(pick, covars[1]))

  # fit Cox model
  fit_ok <- FALSE; used_cov <- pick
  while (!fit_ok && length(used_cov)) {
    form   <- as.formula(paste(
      "survival::Surv(diagsurvtime, status) ~",
      paste(used_cov, collapse = " + ")
    ))
    fit    <- try(survival::coxph(form, data = d, ties = "efron"), silent = TRUE)
    fit_ok <- !inherits(fit, "try-error")
    if (!fit_ok) used_cov <- head(used_cov, -1)
  }
  if (!fit_ok) {
    used_cov <- covars[1]
    fit <- survival::coxph(
      as.formula(paste("survival::Surv(diagsurvtime, status) ~", covars[1])),
      data = d, ties = "efron"
    )
  }

  tt <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  if (!nrow(tt)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = .5, y = .5,
                          label = "No estimable coefficients") +
        ggplot2::labs(title = title) +
        ggplot2::theme_void(base_size = base_size)
    )
  }

  # pretty term labels
  pretty_term <- function(term) {
    for (v in used_cov) {
      if (startsWith(term, v)) {
        suf <- gsub("^[:.=]", "", sub(paste0("^", v), "", term))
        if (suf == "") {
          if (v %in% names(numeric_scale)) {
            nm   <- if (!is.null(numeric_pretty[[v]])) numeric_pretty[[v]] else gsub("_", " ", v)
            unit <- if (!is.null(numeric_units[[v]])) numeric_units[[v]] else "units"
            return(sprintf("%s (per %s %s)", nm, numeric_scale[[v]], unit))
          } else return(v)
        } else {
          ref <- if (v %in% names(ref_levels) && !is.null(ref_levels[[v]])) {
            paste0(" (ref=", ref_levels[[v]], ")")
          } else NULL
          return(paste0(v, ": ", suf, ifelse(is.null(ref), "", ref)))
        }
      }
    }
    term
  }

  tt$label <- vapply(tt$term, pretty_term, character(1))
  tt <- tt %>%
    dplyr::filter(!is.na(estimate)) %>%
    dplyr::mutate(label = factor(label, levels = rev(unique(label))))

  xr   <- range(tt$conf.low, tt$conf.high, na.rm = TRUE)
  if (!all(is.finite(xr))) xr <- c(0.5, 2)
  xpad <- exp(log(xr) + c(-0.25, 0.25))

  ggplot2::ggplot(tt, ggplot2::aes(y = label, x = estimate)) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = conf.low, xmax = conf.high),
                           orientation = "y", width = 0.2, linewidth = 1) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_x_log10(limits = xpad) +
    ggplot2::labs(
      x        = "Hazard Ratio (log scale)",
      y        = NULL,
      title    = title,
      subtitle = sprintf(
        "complete cases n=%d (of %d twins), events=%d; covars: %s",
        n, nrow(df), ev, paste(used_cov, collapse = ", ")
      )
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      text       = ggplot2::element_text(face = if (bold) "bold" else "plain"),
      plot.title = ggplot2::element_text(size = title_size),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text  = ggplot2::element_text(size = axis_text_size)
    )
}
