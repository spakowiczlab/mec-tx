# ============================================================
# MEC-TX visualization: cox_forest_plot_from_df()
# visualization/cox_forest_plot_from_df.R
# ============================================================
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
  resolved_covars <- character(length(covars))
  for (i in seq_along(covars)) {
    m <- names(df)[tolower(names(df)) == tolower(covars[i])]
    if (length(m) == 1) {
      resolved_covars[i] <- m
    } else {
      warning(sprintf("covar '%s' not found — skipping.", covars[i]))
      resolved_covars[i] <- NA_character_
    }
  }
  covars <- resolved_covars[!is.na(resolved_covars)]

  # resolve ref_levels keys to actual column names
  names(ref_levels) <- vapply(names(ref_levels), function(nm) {
    m <- names(df)[tolower(names(df)) == tolower(nm)]
    if (length(m) == 1) m else nm
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
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = conf.low, xmax = conf.high),
                             height = 0.2, linewidth = 1) +
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
