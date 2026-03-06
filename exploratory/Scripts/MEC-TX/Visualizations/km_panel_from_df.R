# ============================================================
# MEC-TX visualization: km_panel_from_df()
# visualization/km_panel_from_df.R
# ============================================================
# ============================================================
# MEC-TX visualization: km_panel_from_df() and cox_forest_plot_from_df()
# visualization/km_panel_from_df.R
# visualization/cox_forest_plot_from_df.R
# ============================================================

# -------- KM with optional risk table --------
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
  group_colours     = NULL   # named colour vector; auto-generated if NULL
) {

  # ---- resolve group_col case-insensitively (Bug 4.4 fix) ----
  actual_group_col <- names(df)[tolower(names(df)) == tolower(group_col)]
  if (length(actual_group_col) == 1) {
    group_col <- actual_group_col
  } else {
    stop(sprintf(
      "group_col '%s' not found in df. Available: %s",
      group_col, paste(names(df), collapse = ", ")
    ))
  }

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
    # use ca_cols for CAlevel, otherwise generate a palette
    if (tolower(group_col) == "calevel" &&
        all(c("High", "Low") %in% lvls)) {
      group_colours <- ca_cols
    } else {
      pal <- c("#4E79A7","#F28E2B","#E15759","#76B7B2",
               "#59A14F","#EDC948","#B07AA1","#FF9DA7")
      group_colours <- setNames(pal[seq_along(lvls)], lvls)
    }
  }

  form <- as.formula(paste("survival::Surv(diagsurvtime, status) ~", group_col))
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

  if (!is.null(horizon_years)) km <- km + ggplot2::coord_cartesian(xlim = c(0, horizon_years))
  if (!risk_table) return(km)

  if (is.null(risk_times)) {
    risk_times <- if (!is.null(horizon_years)) {
      0:horizon_years
    } else {
      pretty(range(df$diagsurvtime, na.rm = TRUE), n = 6)
    }
  }

  sm <- summary(fit, times = risk_times, extend = TRUE)
  rt <- tibble::tibble(
    time    = sm$time,
    group   = factor(sub(".*=", "", as.character(sm$strata)), levels = lvls),
    n_risk  = sm$n.risk
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


