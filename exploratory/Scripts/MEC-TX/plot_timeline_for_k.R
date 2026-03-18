# ============================================================
# MEC-TX visualization: plot_timeline_for_k()
# visualization/plot_timeline_for_k.R
# ============================================================

#' Treatment Timeline Facet Plot by Cluster Assignment
#'
#' Renders a faceted swimlane plot of patient treatment timelines, with one
#' facet per cluster. Each horizontal segment represents a treatment interval
#' coloured by type. Patients within each facet are ordered by earliest
#' treatment start or dominant treatment share. Optionally saves to disk.
#'
#' @param kc Character string. Name of the cluster assignment column in
#'   \code{metadata} (e.g. \code{"Cluster_k3"}). Used both to facet the
#'   plot and to label each panel as \code{c1 (n=X)}, \code{c2 (n=X)},
#'   etc.
#' @param metadata A data frame containing at minimum \code{sample} and
#'   the column named in \code{kc}. Patients with \code{NA} in \code{kc}
#'   are silently dropped. Typically \code{Cluster_surv} (LUSC) or
#'   \code{LUAD_metadata} (LUAD).
#' @param segs A segments data frame passed to \code{\link{prep_segs}} for
#'   recoding and horizon clipping. Must contain columns \code{sample},
#'   \code{type}, \code{t0}, and \code{t1}.
#' @param out_file Character string or \code{NULL}. If supplied, the plot
#'   is saved to this path via \code{pdf()} at 20 × 12 inches. Always use
#'   a \code{.pdf} extension — \code{png()} requires X11/Cairo which is
#'   unavailable on OSC. If \code{NULL}, the plot is returned invisibly
#'   without saving. Default \code{NULL}.
#' @param ncols Integer. Number of columns in the \code{facet_wrap} layout.
#'   Default \code{3}.
#' @param horizon_years Numeric. Maximum x-axis extent in years, passed to
#'   \code{\link{prep_segs}}. Segments extending beyond this value are
#'   clipped. Default \code{6}.
#' @param base_size Base font size (pt). Default \code{14}.
#' @param title_size Font size (pt) for the plot title. Default \code{20}.
#' @param axis_title_size Font size (pt) for axis titles. Default \code{13}.
#' @param axis_text_size Font size (pt) for axis tick labels. Default \code{11}.
#' @param legend_title_size Font size (pt) for the legend title.
#'   Default \code{12}.
#' @param legend_text_size Font size (pt) for legend item labels.
#'   Default \code{11}.
#' @param bold Logical. If \code{TRUE} all text elements use bold weight.
#'   Default \code{TRUE}.
#' @param order_by One of \code{"first_start"} or \code{"dominant_share"}.
#'   Controls row ordering within each cluster facet.
#'   \code{"first_start"} orders by earliest treatment start time
#'   (ascending). \code{"dominant_share"} orders by the proportion of
#'   total treated time accounted for by the single most common treatment
#'   type (descending), with earliest start as a tiebreaker.
#'   Default \code{"first_start"}.
#'
#' @return The \code{ggplot} object, returned \code{invisible()}. If
#'   \code{out_file} is supplied the plot is also written to disk via
#'   \code{pdf()} and a message is printed. Returns \code{invisible(NULL)}
#'   with a warning if no segments remain after joining to \code{metadata}.
#'
#' @details
#' \strong{Colour palette:} Treatment type colours are defined locally
#' inside the function (matching \code{tx_cols} in \code{constants.R}) and
#' do not depend on any global object. Only types present in the data are
#' included in the legend.
#'
#' \strong{Facet labels:} Each panel is labelled \code{cN (n=X)} where
#' \code{N} is the numeric cluster index extracted from the \code{kc}
#' column value and \code{X} is the patient count. Non-numeric cluster
#' values are used as-is.
#'
#' \strong{Y-axis:} Patient identity labels are suppressed on the y-axis
#' for readability. Row position within each facet reflects \code{order_by}
#' only.
#'
#' @examples
#' \dontrun{
#' segs <- prep_segs(intervals$timeline_long_intv)
#'
#' # Plot k=3 clustering, order by dominant share, save as PDF
#' plot_timeline_for_k(
#'   kc       = "Cluster_k3",
#'   metadata = Cluster_surv,
#'   segs     = segs,
#'   out_file = file.path(out_dir, "timelines_k3.pdf"),
#'   order_by = "dominant_share"
#' )
#'
#' # Return plot object only, no saving
#' p <- plot_timeline_for_k(
#'   kc       = "Cluster_k3",
#'   metadata = LUAD_metadata,
#'   segs     = segs
#' )
#' pdf(file.path(out_dir, "timelines_k3_luad.pdf"), width = 20, height = 12)
#' print(p)
#' dev.off()
#' }
#'
#' @seealso \code{\link{prep_segs}}, \code{\link{timeline_panel}},
#'   \code{\link{tx_cluster_surv}}
#'
#' @import ggplot2
#' @importFrom dplyr transmute filter inner_join distinct count group_by
#'   summarise arrange mutate ungroup left_join select slice_max if_else
#'   row_number
#' @importFrom grid unit
#' @export
plot_timeline_for_k <- function(
    kc,
    metadata,                           # ← renamed from Cluster_surv (Thread 8)
    segs,
    out_file          = NULL,
    ncols             = 3,
    horizon_years     = 6,
    base_size         = 14,
    title_size        = 20,
    axis_title_size   = 13,
    axis_text_size    = 11,
    legend_title_size = 12,
    legend_text_size  = 11,
    bold              = TRUE,
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
  
  # ---- attach cluster labels ----                    # ← metadata (Thread 8)
  cl_df <- metadata %>%
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
  
  # ---- save if out_file provided ----          # ← pdf() fix (Thread 8)
  if (!is.null(out_file)) {
    out_file <- sub("\\.(png|jpg|jpeg)$", ".pdf", out_file)
    pdf(out_file, width = 20, height = 12)
    print(p)
    dev.off()
    message(sprintf("Saved: %s", out_file))
  }
  
  invisible(p)
}