
## Plotting the timeline
# 2) Plot timelines for all k
plot_timeline_for_k <- function(
    kc, Cluster_normal_Surv, segs,
    out_file,
    ncols = 3,
    x_years = 13,
    # styling
    palette = c(
      Ancillary       = "#E1BE6A",
      Chemo           = "#FDB863",
      Hormone         = "#DC267F",
      IO              = "#2CA02C",
      Small_Molecule  = "#76B7B2",
      Targeted        = "#4E79A7",
      Radiation       = "#6A51A3",
      Other           = "#3182BD",
      None            = "#C0C0C0"
    ),
    base_size = 14,
    title_size = 20,
    axis_title_size = 13,
    axis_text_size = 11,
    legend_title_size = 12,
    legend_text_size = 11,
    bold = TRUE,
    # ordering: "first_start" (default) or "dominant_share"
    order_by = c("first_start", "dominant_share")
) {
  order_by <- match.arg(order_by)
  
  # --- map raw types -> canonical + clip to horizon ---
  segs_clip <- segs %>%
    dplyr::mutate(type_raw = as.character(type)) %>%
    dplyr::mutate(
      type = dplyr::case_when(
        stringr::str_detect(type_raw, stringr::regex("radiation|\\brt\\b|xrt|imrt|sbrt|radiother", TRUE)) ~ "Radiation",
        stringr::str_detect(type_raw, stringr::regex("\\bio\\b|immuno|pembro|nivol|atezo|ipi|pd-1|pd-l1|ctla", TRUE)) ~ "IO",
        stringr::str_detect(type_raw, stringr::regex("chemo|chemotherapy|platin|taxel|5fu|gemcitabine|doxo", TRUE)) ~ "Chemo",
        stringr::str_detect(type_raw, stringr::regex("hormone|endocrine|androgen|estrogen|tamox|abiraterone|enzalu|aromatase|letro|anastro|fulves", TRUE)) ~ "Hormone",
        stringr::str_detect(type_raw, stringr::regex("target|tki|inhibitor|egfr|alk|braf|mek|parp|her2|trast", TRUE)) ~ "Targeted",
        stringr::str_detect(type_raw, stringr::regex("small[ _-]?molecule|onco.?drug|oncodrug", TRUE)) ~ "Small_Molecule",
        stringr::str_detect(type_raw, stringr::regex("ancillary|support|palliat", TRUE)) ~ "Ancillary",
        TRUE ~ "Other"
      ),
      t0 = pmax(0, pmin(start_year, x_years)),
      t1 = pmax(0, pmin(end_year,   x_years))
    ) %>%
    dplyr::filter(!is.na(sample), t1 > t0)
  
  # --- attach cluster labels ---
  cl_df <- Cluster_normal_Surv %>%
    dplyr::transmute(sample, cluster = as.character(.data[[kc]])) %>%
    dplyr::filter(!is.na(cluster))
  
  df <- segs_clip %>% dplyr::inner_join(cl_df, by = "sample")
  if (nrow(df) == 0) {
    warning(paste0("No segments to plot for ", kc))
    return(invisible(NULL))
  }
  
  # cluster counts for facet labels
  cl_counts <- df %>%
    dplyr::distinct(cluster, sample) %>%
    dplyr::count(cluster, name = "n_patients")
  
  # choose ordering within cluster
  if (order_by == "first_start") {
    order_tbl <- df %>%
      dplyr::group_by(cluster, sample) %>%
      dplyr::summarise(first_start = min(t0), .groups = "drop") %>%
      dplyr::arrange(cluster, first_start, sample) %>%
      dplyr::group_by(cluster) %>% dplyr::mutate(y = dplyr::row_number()) %>% dplyr::ungroup()
  } else {
    # dominant treated-time share within horizon
    dur <- df %>%
      dplyr::group_by(cluster, sample, type) %>%
      dplyr::summarise(dur = sum(t1 - t0), .groups = "drop")
    tot <- dur %>% dplyr::group_by(cluster, sample) %>% dplyr::summarise(total = sum(dur), .groups = "drop")
    dom <- dur %>%
      dplyr::group_by(cluster, sample) %>%
      dplyr::slice_max(dur, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::left_join(tot, by = c("cluster","sample")) %>%
      dplyr::mutate(dom_share = dplyr::if_else(total > 0, dur/total, 0))
    earliest <- df %>%
      dplyr::group_by(cluster, sample) %>%
      dplyr::summarise(earliest = min(t0), .groups = "drop")
    order_tbl <- dom %>%
      dplyr::left_join(earliest, by = c("cluster","sample")) %>%
      dplyr::arrange(cluster, dplyr::desc(dom_share), earliest, sample) %>%
      dplyr::group_by(cluster) %>% dplyr::mutate(y = dplyr::row_number()) %>% dplyr::ungroup()
  }
  
  # merge order + tidy factors
  df <- df %>%
    dplyr::left_join(order_tbl, by = c("cluster","sample")) %>%
    dplyr::left_join(cl_counts, by = "cluster") %>%
    dplyr::mutate(
      cluster_num = suppressWarnings(as.integer(cluster)),
      cluster_lab = paste0("c", ifelse(is.na(cluster_num), cluster, cluster_num),
                           " (n=", n_patients, ")"),
      cluster_lab = factor(cluster_lab, levels = unique(cluster_lab))
    )
  
  # x-limit
  xmax <- if (is.null(x_years)) ceiling(max(df$t1, na.rm = TRUE)) else x_years
  
  # palette subset to present types
  present_types <- sort(unique(df$type))
  pal_use <- palette[names(palette) %in% present_types]
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = t0, xend = t1, y = y, yend = y, color = type)) +
    ggplot2::geom_segment(linewidth = 0.55, lineend = "butt") +
    ggplot2::facet_wrap(~ cluster_lab, scales = "free_y", ncol = ncols) +
    ggplot2::scale_color_manual(values = pal_use, drop = FALSE,
                                guide = ggplot2::guide_legend(ncol = 1, title = "Treatment Type")) +
    ggplot2::scale_x_continuous(limits = c(0, xmax), breaks = 0:xmax,
                                expand = ggplot2::expansion(mult = c(0, 0.02))) +
    ggplot2::labs(
      title = paste0("Treatment timelines — ", kc),
      x = "Time Since First Treatment (Years)",
      y = "Patient"
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey92", color = "grey70"),
      panel.spacing    = grid::unit(0.7, "lines"),
      axis.text.y      = ggplot2::element_blank(),
      axis.ticks.y     = ggplot2::element_blank(),
      legend.position  = "right",
      # bold + sizes
      text          = ggplot2::element_text(face = if (bold) "bold" else "plain"),
      plot.title    = ggplot2::element_text(size = title_size, face = if (bold) "bold" else "plain"),
      axis.title    = ggplot2::element_text(size = axis_title_size, face = if (bold) "bold" else "plain"),
      axis.text     = ggplot2::element_text(size = axis_text_size, face = if (bold) "bold" else "plain"),
      legend.title  = ggplot2::element_text(size = legend_title_size, face = if (bold) "bold" else "plain"),
      legend.text   = ggplot2::element_text(size = legend_text_size,  face = if (bold) "bold" else "plain")
    )
  
  ggplot2::ggsave(out_file, p, width = 20, height = 12, dpi = 150, limitsize = FALSE)
  p
}


# 3) Build segments from your long file and save PNGs for k=13,14,15


k_cols <- paste0("Cluster_k", 3:20)


# keep only those that actually exist in your data:
k_cols <- intersect(k_cols, names(Cluster_surv_2))


dir.create("//users/PAS1695/dipankor99/Github/digital-twins/exploratory/Figures/Function_figs/dt_timeline", showWarnings = FALSE)
for (kc in k_cols) {
  plot_timeline_for_k(
    kc,
    Cluster_normal_Surv = Cluster_surv_2,  # or Cluster_normal_Surv_2 if that’s your current object
    segs = Refined_timeline,
    out_file = file.path("/users/PAS1695/dipankor99/Github/digital-twins/exploratory/Figures/Function_figs/dt_timeline", paste0("timelines_", kc, ".png")),
    ncols   = 3,
    x_years = 6
  )
}
