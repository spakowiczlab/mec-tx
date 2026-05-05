# ============================================================
# MEC-TX visualization: km_panel_from_df()
# visualization/km_panel_from_df.R
# ============================================================

#' Kaplan-Meier Survival Panel with Optional Risk Table
#'
#' Fits a Kaplan-Meier survival curve stratified by a grouping variable and
#' renders a publication-ready plot with confidence ribbons, log-rank p-value,
#' and an optional at-risk table below the main panel. Degrades gracefully ---
#' returns an annotated blank plot when fewer than two groups are present.
#'
#' @param df A data frame containing at minimum \code{diagsurvtime} (numeric,
#'   years from diagnosis), \code{status} (integer, 0/1), and the column
#'   named in \code{group_col}. Column name matching for \code{group_col}
#'   is case-insensitive.
#' @param group_col Character string. Name of the grouping variable column.
#'   Default \code{"CAlevel"}. When \code{group_col} is \code{"CAlevel"}
#'   and both \code{"High"} and \code{"Low"} levels are present,
#'   \code{ca_cols} from \code{constants.R} is used automatically for
#'   colouring.
#' @param title Character string. Plot title. Default \code{""}.
#' @param horizon_years Numeric or \code{NULL}. If supplied, the x-axis is
#'   clipped to this value via \code{coord_cartesian}. Also used as the
#'   upper bound for \code{risk_times} when \code{risk_times} is
#'   \code{NULL}. Default \code{NULL} (full follow-up range).
#' @param risk_table Logical. If \code{TRUE}, an at-risk table is appended
#'   below the KM panel. Default \code{TRUE}.
#' @param risk_times Numeric vector or \code{NULL}. Time points at which
#'   at-risk counts are displayed. \code{NULL} auto-generates breakpoints
#'   from \code{0:horizon_years} (if supplied) or
#'   \code{pretty(range(diagsurvtime))}. Default \code{NULL}.
#' @param risk_table_height Numeric. Relative height of the at-risk table
#'   panel as a fraction of the KM panel height, passed to
#'   \code{patchwork::plot_layout()}. Default \code{0.28}.
#' @param risk_text_size Numeric. Font size for at-risk count labels.
#'   Default \code{4}.
#' @param base_size Base font size (pt) for the KM panel. Default \code{16}.
#' @param title_size Font size (pt) for the plot title. Default \code{22}.
#' @param axis_title_size Font size (pt) for axis titles. Default \code{14}.
#' @param axis_text_size Font size (pt) for axis tick labels. Default \code{12}.
#' @param legend_title_size Font size (pt) for the legend title.
#'   Default \code{14}.
#' @param legend_text_size Font size (pt) for legend item labels.
#'   Default \code{12}.
#' @param bold Logical. If \code{TRUE} all text elements use bold weight.
#'   Default \code{TRUE}.
#' @param group_colours Named character vector mapping group levels to hex
#'   colours. Must be named with the levels of \code{group_col}. \code{NULL}
#'   auto-generates: uses \code{ca_cols} for \code{CAlevel}, otherwise
#'   cycles through a built-in 8-colour palette. Default \code{NULL}.
#'
#' @return A \code{patchwork} object in all cases. When
#'   \code{risk_table = TRUE}, the KM panel (top) and at-risk table
#'   (bottom) are combined via \code{patchwork::plot_layout()}. When
#'   \code{risk_table = FALSE}, the KM panel is wrapped in
#'   \code{patchwork::wrap_plots()} for a consistent return type. Returns
#'   an annotated blank \code{ggplot} if fewer than two groups are present
#'   in \code{group_col}.
#'
#' @details
#' \strong{Column requirements:} \code{diagsurvtime} and \code{status} are
#' hardcoded. \code{status} must be coded 0 (censored) / 1 (event).
#' Rename columns upstream if necessary.
#'
#' \strong{Time-zero rows:} The function ensures a survival estimate of 1.0
#' at time 0 exists for every group before plotting, preventing stepped
#' curves that start below 1.
#'
#' \strong{Log-rank p-value:} Computed via \code{survival::survdiff()} with
#' 1 degree of freedom. Displayed in the subtitle as \code{signif(p, 3)}.
#' \code{NA} is shown if the chi-square statistic is non-finite.
#'
#' \strong{Colour precedence:} \code{group_colours} argument >
#' \code{ca_cols} auto-detection > built-in palette.
#'
#' \strong{Consistent return type:} Always returns a \code{patchwork}
#' object so downstream code assembling multi-panel figures does not need
#' to branch on \code{risk_table}.
#'
#' @examples
#' set.seed(42)
#' n <- 6
#' spec_ages <- seq(55, 80, by = 5)
#' tx_types <- list(
#'   c('Chemo','IO','Radiation'),
#'   c('Chemo','Targeted','Others'),
#'   c('IO','Radiation','Chemo'),
#'   c('Targeted','Chemo','IO'),
#'   c('Radiation','Others','Chemo'),
#'   c('IO','Targeted','Chemo')
#' )
#' med_data <- do.call(rbind, lapply(seq_len(n), function(i) {
#'   data.frame(
#'     sample                     = paste0('P', i),
#'     Age.At.Specimen.Collection = spec_ages[i],
#'     AgeAtLastContact           = spec_ages[i] + 3,
#'     diagsurvtime               = 3,
#'     Status                     = i %% 2L,
#'     Medication                 = c('DrugA','DrugB','DrugC'),
#'     treatment_group            = tx_types[[i]],
#'     AgeAtMedStart              = spec_ages[i] + c(0.1, 0.5, 1.0),
#'     AgeAtMedStop               = spec_ages[i] + c(0.4, 0.9, 1.3),
#'     AgeAtTreatmentStart.mod    = spec_ages[i] + c(0.1, 0.5, 1.0),
#'     stringsAsFactors           = FALSE
#'   )
#' }))
#' meta <- data.frame(
#'   sample       = paste0('P', seq_len(n)),
#'   diagsurvtime = rep(3, n),
#'   Status       = seq_len(n) %% 2L,
#'   CAlevel      = rep(c('High','Low'), n/2),
#'   stringsAsFactors = FALSE
#' )
#' norm        <- tx_normalize(med_data)
#' cluster_res <- tx_cluster_surv(meta, norm, k_range = 2,
#'                                umap_neighbors = 5,
#'                                min_feature_variance = 0)
#' p <- km_panel_from_df(
#'   df        = cluster_res$Cluster_surv,
#'   group_col = 'CAlevel',
#'   title     = 'KM by CAlevel'
#' )
#' class(p)
#'
#' @seealso \code{\link{cox_forest_plot_from_df}},
#'   \code{\link{tx_pooled_analysis}}, \code{\link{tx_compare_groups}}
#'
#' @import ggplot2
#' @importFrom survival survfit survdiff Surv
#' @importFrom broom tidy
#' @importFrom dplyr filter arrange bind_rows n_distinct
#' @importFrom tibble tibble
#' @importFrom patchwork plot_layout wrap_plots
#' @importFrom stats pchisq as.formula
#' @export
km_panel_from_df <- function(
    df,
    group_col         = "CAlevel",
    title             = "",
    horizon_years     = NULL,
    risk_table        = TRUE,
    risk_times        = NULL,
    risk_table_height = 0.28,
    risk_text_size    = 4,
    base_size         = 16,
    title_size        = 22,
    axis_title_size   = 14,
    axis_text_size    = 12,
    legend_title_size = 14,
    legend_text_size  = 12,
    bold              = TRUE,
    group_colours     = NULL
) {
  
  # ---- resolve group_col case-insensitively (Bug 4.4 fix) ----
  group_col <- resolve_col(df, group_col, "group_col")
  
  stopifnot(all(c("diagsurvtime", "status", group_col) %in% names(df)))
  
  df <- df %>% dplyr::filter(!is.na(.data[[group_col]]))
  df[[group_col]] <- factor(df[[group_col]])
  lvls <- levels(df[[group_col]])
  
  if (dplyr::n_distinct(df[[group_col]]) < 2) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = .5, y = .5,
                          label = paste("Only one", group_col, "group in subset"), size = 5) +
        ggplot2::labs(title = title, x = "Years", y = "Survival") +
        ggplot2::theme_minimal(base_size = 12)
    )
  }
  
  # ---- auto-generate colours if not supplied ----
  if (is.null(group_colours)) {
    if (tolower(group_col) == "calevel" &&
        all(c("High", "Low") %in% lvls)) {
      group_colours <- ca_cols
    } else {
      pal <- c("#4E79A7","#F28E2B","#E15759","#76B7B2",
               "#59A14F","#EDC948","#B07AA1","#FF9DA7")
      group_colours <- setNames(pal[seq_along(lvls)], lvls)
    }
  }
  
  form <- stats::as.formula(
    paste("survival::Surv(diagsurvtime, status) ~", group_col)
  )
  fit  <- survival::survfit(form, data = df)
  td   <- broom::tidy(fit)
  td$group <- sub(".*=", "", td$strata)
  td$group <- factor(td$group, levels = lvls)
  
  # ensure time-0 rows exist for each group
  for (lv in lvls) {
    if (!any(td$time == 0 & td$group == lv)) {
      td <- dplyr::bind_rows(
        tibble::tibble(
          time = 0, estimate = 1, conf.low = 1, conf.high = 1,
          strata = paste0(group_col, "=", lv),
          n.risk = NA_real_, n.event = NA_real_, n.censor = NA_real_,
          group = factor(lv, levels = lvls)
        ),
        td
      )
    }
  }
  td <- td %>% dplyr::arrange(group, time)
  
  lr   <- survival::survdiff(form, data = df)
  pval <- if (is.finite(lr$chisq)) {
    stats::pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  } else NA_real_
  
  km <- ggplot2::ggplot(td, ggplot2::aes(x = time, y = estimate,
                                         colour = group, fill = group)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = conf.low, ymax = conf.high),
                         alpha = .20, linewidth = 0) +
    ggplot2::geom_step(linewidth = 1.1) +
    ggplot2::scale_colour_manual(values = group_colours, drop = FALSE) +
    ggplot2::scale_fill_manual(values = group_colours, drop = FALSE) +
    ggplot2::labs(
      x        = "Years",
      y        = "Survival",
      title    = title,
      subtitle = sprintf(
        "%s grouping: provided | log-rank p = %s",
        group_col,
        ifelse(is.na(pval), "NA", signif(pval, 3))
      )
    ) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      text            = ggplot2::element_text(face = if (bold) "bold" else "plain"),
      plot.title      = ggplot2::element_text(size = title_size),
      axis.title      = ggplot2::element_text(size = axis_title_size),
      axis.text       = ggplot2::element_text(size = axis_text_size),
      legend.title    = ggplot2::element_text(size = legend_title_size),
      legend.text     = ggplot2::element_text(size = legend_text_size)
    )
  
  if (!is.null(horizon_years)) {
    km <- km + ggplot2::coord_cartesian(xlim = c(0, horizon_years))
  }
  
  # -- consistent return type fix (Thread 8) ---------------------------------
  if (!risk_table) return(patchwork::wrap_plots(km))
  
  if (is.null(risk_times)) {
    risk_times <- if (!is.null(horizon_years)) {
      0:horizon_years
    } else {
      pretty(range(df$diagsurvtime, na.rm = TRUE), n = 6)
    }
  }
  
  sm <- summary(fit, times = risk_times, extend = TRUE)
  rt <- tibble::tibble(
    time   = sm$time,
    group  = factor(sub(".*=", "", as.character(sm$strata)), levels = lvls),
    n_risk = sm$n.risk
  )
  
  risk <- ggplot2::ggplot(rt, ggplot2::aes(x = time, y = group,
                                           label = n_risk, colour = group)) +
    ggplot2::geom_text(size = risk_text_size, show.legend = FALSE) +
    ggplot2::scale_colour_manual(values = group_colours, drop = FALSE) +
    ggplot2::scale_x_continuous(breaks = risk_times, limits = range(risk_times)) +
    ggplot2::labs(x = NULL, y = "At risk") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "none"
    )
  
  (km + ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.text.x  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )) /
    risk +
    patchwork::plot_layout(heights = c(1, risk_table_height))
}
