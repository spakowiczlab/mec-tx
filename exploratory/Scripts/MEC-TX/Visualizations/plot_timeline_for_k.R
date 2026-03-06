# ============================================================
# MEC-TX visualization: plot_timeline_for_k()
# visualization/plot_timeline_for_k.R
# ============================================================

plot_timeline_for_k <- function(
  kc,
  Cluster_surv,
  segs,
  out_file          = NULL,   # if NULL, plot is returned but not saved
  ncols             = 3,
  horizon_years     = 6,
  base_size         = 14,
  title_size        = 20,
  axis_title_size   = 13,
  axis_text_size    = 11,
  legend_title_size = 12,
  legend_text_size  = 11,
  bold              = TRUE,
  # ordering within each cluster facet
  # "first_start"    = order by earliest treatment start
  # "dominant_share" = order by dominant treatment type share
  order_by          = c("first_start", "dominant_share")
) {

  order_by <- match.arg(order_by)

  # local colour palette — consistent with rest of MEC-TX package
  local_tx_cols <- c(
    Ancillary      = "#E1BE6A",
    Chemo          = "#FDB863",
    Hormone        = "#DC267F",
    IO             = "#2CA02C",
    Small_Molecule = "#76B7B2",
    Targeted       = "#4E79A7",
    Radiation      = "#6A51A3",
    Others         = "#8C8C8C"
  )

  # ---- recode & clip using prep_segs() — no duplicated logic ----
  segs_clip <- prep_segs(segs, horizon_years = horizon_years)

  # ---- attach cluster labels ----
  cl_df <- Cluster_surv %>%
    dplyr::transmute(
      sample,
      cluster = as.character(.data[[kc]])
    ) %>%
    dplyr::filter(!is.na(cluster))

  df <- segs_clip %>%
    dplyr::inner_join(cl_df, by = "sample")

  if (nrow(df) == 0) {
    warning(paste0("No segments to plot for ", kc))
    return(invisible(NULL))
  }

  # ---- cluster counts for facet labels ----
  cl_counts <- df %>%
    dplyr::distinct(cluster, sample) %>%
    dplyr::count(cluster, name = "n_patients")

  # ---- patient ordering within each cluster ----
  if (order_by == "first_start") {
    order_tbl <- df %>%
      dplyr::group_by(cluster, sample) %>%
      dplyr::summarise(first_start = min(t0), .groups = "drop") %>%
      dplyr::arrange(cluster, first_start, sample) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(y = dplyr::row_number()) %>%
      dplyr::ungroup()

  } else {
    # dominant treated-time share within horizon
    dur <- df %>%
      dplyr::group_by(cluster, sample, type) %>%
      dplyr::summarise(dur = sum(t1 - t0), .groups = "drop")

    tot <- dur %>%
      dplyr::group_by(cluster, sample) %>%
      dplyr::summarise(total = sum(dur), .groups = "drop")

    dom <- dur %>%
      dplyr::group_by(cluster, sample) %>%
      dplyr::slice_max(dur, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(tot, by = c("cluster", "sample")) %>%
      dplyr::mutate(
        dom_share = dplyr::if_else(total > 0, dur / total, 0)
      )

    earliest <- df %>%
      dplyr::group_by(cluster, sample) %>%
      dplyr::summarise(earliest = min(t0), .groups = "drop")

    order_tbl <- dom %>%
      dplyr::left_join(earliest, by = c("cluster", "sample")) %>%
      dplyr::arrange(cluster, dplyr::desc(dom_share), earliest, sample) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(y = dplyr::row_number()) %>%
      dplyr::ungroup()
  }

  # ---- merge ordering and build facet labels ----
  df <- df %>%
    dplyr::left_join(order_tbl %>% dplyr::select(cluster, sample, y),
                     by = c("cluster", "sample")) %>%
    dplyr::left_join(cl_counts, by = "cluster") %>%
    dplyr::mutate(
      cluster_num = suppressWarnings(as.integer(cluster)),
      cluster_lab = paste0(
        "c", ifelse(is.na(cluster_num), cluster, cluster_num),
        " (n=", n_patients, ")"
      ),
      cluster_lab = factor(cluster_lab, levels = unique(cluster_lab)),
      type        = factor(type, levels = names(local_tx_cols))
    )

  # ---- palette subset to types present in data ----
  present_types <- intersect(names(local_tx_cols), unique(as.character(df$type)))
  pal_use       <- local_tx_cols[present_types]

  # ---- build plot ----
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = t0, xend = t1, y = y, yend = y, color = type)
  ) +
    ggplot2::geom_segment(linewidth = 0.55, lineend = "butt") +
    ggplot2::facet_wrap(~ cluster_lab, scales = "free_y", ncol = ncols) +
    ggplot2::scale_color_manual(
      values = pal_use,
      drop   = FALSE,
      guide  = ggplot2::guide_legend(ncol = 1, title = "Treatment Type")
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, horizon_years),
      breaks = 0:horizon_years,
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      title = paste0("Treatment timelines — ", kc),
      x     = "Time since first treatment (years)",
      y     = "Patient"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background   = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      panel.spacing      = grid::unit(0.7, "lines"),
      axis.text.y        = ggplot2::element_blank(),
      axis.ticks.y       = ggplot2::element_blank(),
      legend.position    = "right",
      text               = ggplot2::element_text(
        face = if (bold) "bold" else "plain"),
      plot.title         = ggplot2::element_text(size = title_size),
      axis.title         = ggplot2::element_text(size = axis_title_size),
      axis.text          = ggplot2::element_text(size = axis_text_size),
      legend.title       = ggplot2::element_text(size = legend_title_size),
      legend.text        = ggplot2::element_text(size = legend_text_size)
    )

  # ---- save if out_file provided ----
  if (!is.null(out_file)) {
    ggplot2::ggsave(
      out_file, p,
      width     = 20,
      height    = 12,
      dpi       = 150,
      limitsize = FALSE
    )
    message(sprintf("Saved: %s", out_file))
  }

  invisible(p)
}
