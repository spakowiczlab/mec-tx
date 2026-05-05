# ============================================================
# MEC-TX visualization: timeline_panel()
# visualization/timeline_panel.R
# ============================================================

#' Treatment Timeline Panel for Selected Patients
#'
#' Renders a horizontal swimlane plot showing treatment type over time for a
#' specified set of patients (e.g. cluster representatives or "top twins").
#' Each patient occupies one row; treatment types are drawn as coloured
#' horizontal segments with a small vertical offset per type to reduce
#' overplotting. Row order is driven by a focus-type share score.
#'
#' @param segs_prepped A pre-processed segments data frame produced by
#'   \code{prep_segs}. Must contain columns \code{sample},
#'   \code{type}, \code{t0}, and \code{t1} (treatment start and end in
#'   years since first treatment).
#' @param share_df A data frame with one row per patient containing
#'   \code{sample} and any number of \code{share_*_tx} columns (treatment
#'   type duration shares) plus \code{dur_treated}. Typically the
#'   \code{$shares} slot from \code{treatment_shares}.
#' @param twin_ids Character vector of patient \code{sample} identifiers to
#'   display. Rows in \code{segs_prepped} not in \code{twin_ids} are
#'   silently dropped.
#' @param title Character string. Plot title. Default \code{"Top twins"}.
#' @param horizon_years Numeric. Maximum x-axis extent in years.
#'   Default \code{5}.
#' @param focus_types Character vector or \code{NULL}. Treatment type(s) used
#'   to compute the row ordering score. Patients with higher combined share
#'   of these types appear at the top. Must match type labels used in
#'   \code{share_*_tx} column names (e.g. \code{"Radiation"} maps to
#'   \code{share_Radiation_tx}). \code{NULL} falls back to the first
#'   \code{share_*_tx} column found.
#' @param base_size Base font size (pt). Default \code{16}.
#' @param title_size Font size (pt) for the plot title. Default \code{22}.
#' @param axis_title_size Font size (pt) for axis titles. Default \code{14}.
#' @param axis_text_size Font size (pt) for axis tick labels. Default \code{12}.
#' @param legend_title_size Font size (pt) for the legend title.
#'   Default \code{14}.
#' @param legend_text_size Font size (pt) for legend item labels.
#'   Default \code{12}.
#' @param bold Logical. If \code{TRUE} all text elements use bold weight.
#'   Default \code{TRUE}.
#'
#' @return A \code{ggplot} object. Print to display or save with
#'   \code{pdf()} / \code{ggsave()}.
#'
#' @details
#' \strong{Colour palette:} Treatment type colours are defined locally
#' inside the function (matching \code{tx_cols} in \code{constants.R}) and
#' do not depend on any global object. Only types present in the data are
#' included in the legend.
#'
#' \strong{Row ordering:} The ordering score is the row-wise sum of all
#' \code{focus_share_cols} shares per patient. Ties are broken by the
#' factor level order of \code{sample}.
#'
#' \strong{Vertical offset:} Each treatment type within a patient row is
#' offset by \code{0.10 * (type_idx - mean(type_idx))} on the y-axis to
#' visually separate overlapping segments without displacing the row label.
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
#' norm      <- tx_normalize(med_data)
#' intervals <- tx_intervals(norm)
#' segs      <- prep_segs(intervals)
#' sh        <- treatment_shares(segs)
#' p <- timeline_panel(segs, sh, twin_ids = unique(intervals$sample))
#' class(p)
#'
#' @seealso \code{prep_segs}, \code{treatment_shares},
#'   \code{\link{plot_timeline_for_k}}
#'
#' @import ggplot2
#' @importFrom dplyr semi_join left_join select starts_with distinct mutate
#'   across all_of
#' @importFrom tibble tibble
#' @importFrom forcats fct_reorder
#' @export

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

  # local constants --- no global scope dependency
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
