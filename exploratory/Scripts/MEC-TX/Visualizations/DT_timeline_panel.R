# ============================================================
# MEC-TX visualization: timeline_panel()
# visualization/timeline_panel.R
# ============================================================

timeline_panel <- function(
  segs_prepped,
  share_df,
  twin_ids,
  title             = "Top twins",
  horizon_years     = 5,
  focus_types       = NULL,   # used for ordering; NULL = order by first share col
  base_size         = 16,
  title_size        = 22,
  axis_title_size   = 14,
  axis_text_size    = 12,
  legend_title_size = 14,
  legend_text_size  = 12,
  bold              = TRUE
) {

  # local constants — no global scope dependency
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

  # ---- restrict to selected twins and attach shares ----
  segs_plot <- segs_prepped %>%
    dplyr::semi_join(tibble::tibble(sample = twin_ids), by = "sample") %>%
    dplyr::left_join(
      share_df %>%
        dplyr::select(sample, dplyr::starts_with("share_"), dur_treated),
      by = "sample"
    )

  # ---- dynamic ordering by focus_types share (Bug fix: was hardcoded Radiation) ----
  # if focus_types supplied, order by combined focus share
  # otherwise fall back to first share column found
  share_cols <- names(segs_plot)[
    startsWith(names(segs_plot), "share_") & endsWith(names(segs_plot), "_tx")
  ]

  if (!is.null(focus_types) && length(focus_types) >= 1) {
    focus_share_cols <- paste0("share_", focus_types, "_tx")
    focus_share_cols <- intersect(focus_share_cols, names(segs_plot))
  } else {
    focus_share_cols <- share_cols[1]
  }

  # compute ordering score per sample
  order_scores <- segs_plot %>%
    dplyr::distinct(sample, dplyr::across(dplyr::all_of(focus_share_cols))) %>%
    dplyr::mutate(
      order_score = rowSums(dplyr::across(dplyr::all_of(focus_share_cols)),
                            na.rm = TRUE)
    ) %>%
    dplyr::select(sample, order_score)

  segs_plot <- segs_plot %>%
    dplyr::left_join(order_scores, by = "sample") %>%
    dplyr::mutate(
      sample     = forcats::fct_reorder(sample, order_score, .desc = TRUE),
      sample_idx = as.numeric(sample),
      type       = factor(type, levels = names(local_tx_cols)),
      type_idx   = as.numeric(type)
    )

  used_types <- intersect(
    names(local_tx_cols),
    unique(as.character(segs_plot$type))
  )

  # ---- vertical offset per treatment type ----
  offset_scale <- 0.10
  segs_plot <- segs_plot %>%
    dplyr::mutate(
      y_pos = sample_idx + (type_idx - mean(type_idx, na.rm = TRUE)) * offset_scale
    )

  y_breaks <- sort(unique(segs_plot$sample_idx))
  y_labels <- levels(segs_plot$sample)

  ggplot2::ggplot(segs_plot) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x      = t0,
        xend   = t1,
        y      = y_pos,
        yend   = y_pos,
        colour = type
      ),
      linewidth = 1
    ) +
    ggplot2::scale_colour_manual(
      values = local_tx_cols[used_types],
      breaks = used_types,
      drop   = FALSE,
      name   = "Treatment Type"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(0, horizon_years),
      breaks = 0:horizon_years,
      name   = "Years since first treatment"
    ) +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      name   = "Twin sample"
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position = "right",
      text            = ggplot2::element_text(
        face = if (bold) "bold" else "plain"
      ),
      plot.title   = ggplot2::element_text(size = title_size),
      axis.title   = ggplot2::element_text(size = axis_title_size),
      axis.text    = ggplot2::element_text(size = axis_text_size),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text  = ggplot2::element_text(size = legend_text_size)
    )
}
