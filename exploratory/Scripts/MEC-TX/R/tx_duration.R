# =============================================================================
# tx_duration.R
# MEC-TX Package
#
# Compare time-on-treatment by a grouping variable (e.g., CAlevel) within
# each treatment type. Measures total calendar-time exposure per patient
# per type and as a merged total across all types.
#
# Consumes tx_intervals() output directly. No dependency on tx_normalize().
# =============================================================================

#' Merge overlapping intervals into non-overlapping unions
#'
#' @param starts Numeric vector of interval start times
#' @param ends   Numeric vector of interval end times
#' @return data.frame with columns `start`, `end` (merged intervals)
#' @keywords internal
.merge_intervals <- function(starts, ends) {
  if (length(starts) == 0L) {
    return(data.frame(start = numeric(0), end = numeric(0)))
  }
  ord   <- order(starts, ends)
  starts <- starts[ord]
  ends   <- ends[ord]
  
  ms <- starts[1]
  me <- ends[1]
  rs <- numeric()
  re <- numeric()
  
  for (i in seq_along(starts)) {
    if (starts[i] <= me) {
      me <- max(me, ends[i])
    } else {
      rs <- c(rs, ms)
      re <- c(re, me)
      ms <- starts[i]
      me <- ends[i]
    }
  }
  rs <- c(rs, ms)
  re <- c(re, me)
  data.frame(start = rs, end = re)
}


#' Compute treatment duration per patient per type
#'
#' For each patient -- treatment-type combination, merges overlapping intervals
#' within that type, then sums calendar time.
#'
#' @param timeline  data.frame from tx_intervals() with columns for sample,
#'                  type, start_year, end_year
#' @param sample_col Column name for patient ID (default "sample")
#' @param type_col   Column name for treatment type (default "type")
#' @param start_col  Column name for interval start (default "start_year")
#' @param end_col    Column name for interval end (default "end_year")
#' @return data.frame: sample, type, duration_yrs
#' @keywords internal
.duration_per_type <- function(timeline,
                               sample_col = "sample",
                               type_col   = "type",
                               start_col  = "start_year",
                               end_col    = "end_year") {
  
  patients <- unique(timeline[[sample_col]])
  types    <- unique(timeline[[type_col]])
  
  out <- vector("list", length(patients) * length(types))
  idx <- 0L
  
  for (pid in patients) {
    sub_p <- timeline[timeline[[sample_col]] == pid, , drop = FALSE]
    for (tx in types) {
      sub_pt <- sub_p[sub_p[[type_col]] == tx, , drop = FALSE]
      if (nrow(sub_pt) == 0L) next
      merged <- .merge_intervals(sub_pt[[start_col]], sub_pt[[end_col]])
      dur    <- sum(merged$end - merged$start)
      idx <- idx + 1L
      out[[idx]] <- data.frame(
        sample       = pid,
        type         = tx,
        duration_yrs = dur,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out[seq_len(idx)])
}


#' Compute merged total treatment duration per patient
#'
#' Unions all treatment intervals (across types) for each patient, merges
#' overlapping periods, and sums calendar time. A patient receiving concurrent
#' chemo + IO for 6 months counts as 6 months total (not 12).
#'
#' @param timeline  data.frame from tx_intervals()
#' @param sample_col Column name for patient ID (default "sample")
#' @param start_col  Column name for interval start (default "start_year")
#' @param end_col    Column name for interval end (default "end_year")
#' @return data.frame: sample, duration_yrs_total
#' @keywords internal
.duration_total <- function(timeline,
                            sample_col = "sample",
                            start_col  = "start_year",
                            end_col    = "end_year") {
  patients <- unique(timeline[[sample_col]])
  out <- data.frame(
    sample             = patients,
    duration_yrs_total = NA_real_,
    stringsAsFactors   = FALSE
  )
  for (i in seq_along(patients)) {
    sub <- timeline[timeline[[sample_col]] == patients[i], , drop = FALSE]
    merged <- .merge_intervals(sub[[start_col]], sub[[end_col]])
    out$duration_yrs_total[i] <- sum(merged$end - merged$start)
  }
  out
}


#' Treatment Duration Analysis by Grouping Variable
#'
#' Computes per-patient treatment duration (calendar time) broken down by
#' treatment type and as a merged total. Compares groups using Wilcoxon
#' rank-sum tests and produces box/violin plots.
#'
#' @param timeline   data.frame --- output of \code{tx_intervals()}.
#'                   Must contain columns for sample, type, start_year, end_year.
#' @param meta       data.frame --- patient-level metadata (e.g., Cluster_surv).
#'                   Must contain \code{sample_col} and \code{group_var}.
#' @param group_var  Character --- column name in \code{meta} for the grouping
#'                   variable (e.g., "CAlevel").
#' @param sample_col Character --- patient ID column (default "sample").
#' @param type_col   Character --- treatment type column (default "type").
#' @param start_col  Character --- interval start column (default "start_year").
#' @param end_col    Character --- interval end column (default "end_year").
#' @param duration_unit Character --- "months" (default) or "years". Controls
#'                   the unit in output tables and plot axis labels.
#' @param exclude_types Character vector --- treatment types to exclude from
#'                   analysis (e.g., \code{c("Others")}). Default NULL.
#' @param min_n      Integer --- minimum patients per group -- type to run a
#'                   test (default 3). Types with fewer are flagged but kept.
#' @param plot       Logical --- produce a plot? (default TRUE)
#' @param plot_type  Character --- "box" (default) or "violin".
#' @param title      Character --- plot title. NULL for auto-generated.
#' @param palette    Named character vector --- colours keyed by group levels.
#'                   NULL for automatic palette.
#'
#' @return A named list:
#' \describe{
#'   \item{duration_per_type}{data.frame: sample, type, group_var, duration
#'         (in \code{duration_unit})}
#'   \item{duration_total}{data.frame: sample, group_var, duration_total}
#'   \item{summary_table}{data.frame: type, group, n, mean, median, q25, q75,
#'         p_value (Wilcoxon)}
#'   \item{plot}{ggplot object (or NULL if \code{plot = FALSE})}
#' }
#'
#' @examples
#' \dontrun{
#' res <- tx_duration(
#'   timeline  = timeline_long_intv,
#'   meta      = Cluster_surv,
#'   group_var = "CAlevel"
#' )
#' res$summary_table
#' res$plot
#' }
#'
#' @export
tx_duration <- function(timeline,
                        meta,
                        group_var,
                        sample_col    = "sample",
                        type_col      = "type",
                        start_col     = "start_year",
                        end_col       = "end_year",
                        duration_unit = "months",
                        exclude_types = NULL,
                        min_n         = 3L,
                        plot          = TRUE,
                        plot_type     = "box",
                        title         = NULL,
                        palette       = NULL) {
  
  # ------------------------------------------------------------------
  # 0. Validate inputs
  # ------------------------------------------------------------------
  stopifnot(is.data.frame(timeline), is.data.frame(meta))
  for (col in c(sample_col, type_col, start_col, end_col)) {
    if (!col %in% names(timeline))
      stop("Column '", col, "' not found in timeline.")
  }
  if (!sample_col %in% names(meta))
    stop("Column '", sample_col, "' not found in meta.")
  if (!group_var %in% names(meta))
    stop("Column '", group_var, "' not found in meta.")
  
  duration_unit <- match.arg(duration_unit, c("months", "years"))
  plot_type     <- match.arg(plot_type, c("box", "violin"))
  multiplier    <- if (duration_unit == "months") 12 else 1
  unit_label    <- if (duration_unit == "months") "Duration (months)" else "Duration (years)"
  
  # ------------------------------------------------------------------
  # 1. Filter types
  # ------------------------------------------------------------------
  if (!is.null(exclude_types)) {
    timeline <- timeline[!timeline[[type_col]] %in% exclude_types, , drop = FALSE]
  }
  if (nrow(timeline) == 0L) stop("No intervals remain after excluding types.")
  
  # Restrict to patients present in both timeline and meta
  shared_ids <- intersect(unique(timeline[[sample_col]]),
                          unique(meta[[sample_col]]))
  if (length(shared_ids) == 0L) stop("No overlapping patients between timeline and meta.")
  timeline <- timeline[timeline[[sample_col]] %in% shared_ids, , drop = FALSE]
  
  # ------------------------------------------------------------------
  # 2. Compute durations
  # ------------------------------------------------------------------
  dur_type <- .duration_per_type(timeline, sample_col, type_col, start_col, end_col)
  dur_total <- .duration_total(timeline, sample_col, start_col, end_col)
  
  # Apply unit conversion
  dur_type$duration  <- dur_type$duration_yrs * multiplier
  dur_total$duration_total <- dur_total$duration_yrs_total * multiplier
  
  # ------------------------------------------------------------------
  # 3. Join group variable from meta
  # ------------------------------------------------------------------
  group_lookup <- meta[, c(sample_col, group_var), drop = FALSE]
  group_lookup <- group_lookup[!duplicated(group_lookup[[sample_col]]), ]
  
  dur_type  <- merge(dur_type,  group_lookup, by = sample_col, all.x = TRUE)
  dur_total <- merge(dur_total, group_lookup, by = sample_col, all.x = TRUE)
  
  # Drop patients with missing group
  dur_type  <- dur_type[!is.na(dur_type[[group_var]]), , drop = FALSE]
  dur_total <- dur_total[!is.na(dur_total[[group_var]]), , drop = FALSE]
  
  groups <- sort(unique(dur_type[[group_var]]))
  if (length(groups) < 2L) {
    warning("Only one level of '", group_var, "' present. No comparison possible.")
  }
  
  # ------------------------------------------------------------------
  # 4. Summary statistics + Wilcoxon tests
  # ------------------------------------------------------------------
  types <- sort(unique(dur_type$type))
  # Include "__TOTAL__" as a pseudo-type for the merged total
  summary_rows <- vector("list", length(types) + 1L)
  
  .summarise_group <- function(vals, grp_label, type_label, n_grp_other) {
    n <- length(vals)
    data.frame(
      type    = type_label,
      group   = grp_label,
      n       = n,
      mean    = if (n > 0) round(mean(vals), 2) else NA_real_,
      median  = if (n > 0) round(median(vals), 2) else NA_real_,
      q25     = if (n > 0) round(quantile(vals, 0.25), 2) else NA_real_,
      q75     = if (n > 0) round(quantile(vals, 0.75), 2) else NA_real_,
      p_value = NA_real_,
      test_note = "",
      stringsAsFactors = FALSE
    )
  }
  
  for (i in seq_along(types)) {
    tx    <- types[i]
    sub   <- dur_type[dur_type$type == tx, , drop = FALSE]
    rows  <- vector("list", length(groups))
    
    # Split by group
    vals_by_group <- lapply(groups, function(g) {
      sub$duration[sub[[group_var]] == g]
    })
    names(vals_by_group) <- groups
    
    # Wilcoxon test (pairwise for 2 groups; Kruskal-Wallis for >2)
    pval      <- NA_real_
    test_note <- ""
    ns <- sapply(vals_by_group, length)
    
    if (length(groups) == 2L && all(ns >= min_n)) {
      wt   <- wilcox.test(vals_by_group[[1]], vals_by_group[[2]])
      pval <- round(wt$p.value, 4)
    } else if (length(groups) > 2L && sum(ns >= min_n) >= 2L) {
      kt   <- kruskal.test(duration ~ get(group_var), data = sub)
      pval <- round(kt$p.value, 4)
    } else {
      test_note <- paste0("skipped (min_n=", min_n, " not met)")
    }
    
    for (j in seq_along(groups)) {
      rows[[j]] <- .summarise_group(vals_by_group[[j]], groups[j], tx, ns)
      rows[[j]]$p_value   <- pval
      rows[[j]]$test_note <- test_note
    }
    summary_rows[[i]] <- do.call(rbind, rows)
  }
  
  # Merged total row
  total_rows <- vector("list", length(groups))
  vals_total_by_group <- lapply(groups, function(g) {
    dur_total$duration_total[dur_total[[group_var]] == g]
  })
  names(vals_total_by_group) <- groups
  ns_total <- sapply(vals_total_by_group, length)
  
  pval_total <- NA_real_
  tnote_total <- ""
  if (length(groups) == 2L && all(ns_total >= min_n)) {
    wt <- wilcox.test(vals_total_by_group[[1]], vals_total_by_group[[2]])
    pval_total <- round(wt$p.value, 4)
  } else if (length(groups) > 2L && sum(ns_total >= min_n) >= 2L) {
    kt <- kruskal.test(duration_total ~ get(group_var), data = dur_total)
    pval_total <- round(kt$p.value, 4)
  } else {
    tnote_total <- paste0("skipped (min_n=", min_n, " not met)")
  }
  
  for (j in seq_along(groups)) {
    total_rows[[j]] <- .summarise_group(
      vals_total_by_group[[j]], groups[j], "All types (merged)", ns_total
    )
    total_rows[[j]]$p_value   <- pval_total
    total_rows[[j]]$test_note <- tnote_total
  }
  summary_rows[[length(types) + 1L]] <- do.call(rbind, total_rows)
  summary_df <- do.call(rbind, summary_rows)
  rownames(summary_df) <- NULL
  
  # ------------------------------------------------------------------
  # 5. Plot
  # ------------------------------------------------------------------
  p <- NULL
  if (plot && requireNamespace("ggplot2", quietly = TRUE)) {
    
    # Build combined data for faceted plot
    plot_df <- dur_type[, c(sample_col, "type", "duration", group_var), drop = FALSE]
    
    total_plot <- dur_total[, c(sample_col, "duration_total", group_var), drop = FALSE]
    names(total_plot)[names(total_plot) == "duration_total"] <- "duration"
    total_plot$type <- "All types (merged)"
    total_plot <- total_plot[, names(plot_df), drop = FALSE]
    
    plot_df <- rbind(plot_df, total_plot)
    
    # Order facets: specific types alphabetically, then "All types (merged)" last
    type_levels <- c(types, "All types (merged)")
    plot_df$type <- factor(plot_df$type, levels = type_levels)
    
    # Palette
    if (is.null(palette)) {
      if (length(groups) <= 8) {
        pal <- grDevices::hcl.colors(max(length(groups), 3), "Dark 3")[seq_along(groups)]
      } else {
        pal <- grDevices::hcl.colors(length(groups), "Dark 3")
      }
      names(pal) <- groups
    } else {
      pal <- palette
    }
    
    # Build p-value labels for facets
    pval_labels <- summary_df[!duplicated(summary_df$type), c("type", "p_value", "test_note")]
    pval_labels$label <- ifelse(
      pval_labels$test_note != "",
      pval_labels$test_note,
      ifelse(
        is.na(pval_labels$p_value),
        "",
        ifelse(
          pval_labels$p_value < 0.001,
          "p < 0.001",
          paste0("p = ", format(pval_labels$p_value, nsmall = 3))
        )
      )
    )
    
    # Compute y-position for annotation (top of each facet)
    max_per_type <- tapply(plot_df$duration, plot_df$type, max, na.rm = TRUE)
    pval_labels$y <- as.numeric(max_per_type[as.character(pval_labels$type)]) * 1.08
    pval_labels$type <- factor(pval_labels$type, levels = type_levels)
    
    if (is.null(title)) {
      title <- paste0("Treatment Duration by ", group_var)
    }
    
    p <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = .data[[group_var]], y = duration, fill = .data[[group_var]])
    )
    
    if (plot_type == "violin") {
      p <- p +
        ggplot2::geom_violin(alpha = 0.6, trim = FALSE) +
        ggplot2::geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA)
    } else {
      p <- p +
        ggplot2::geom_boxplot(alpha = 0.6, outlier.size = 1)
    }
    
    p <- p +
      ggplot2::geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
      ggplot2::facet_wrap(~ type, scales = "free_y") +
      ggplot2::geom_text(
        data    = pval_labels,
        mapping = ggplot2::aes(x = Inf, y = y, label = label),
        inherit.aes = FALSE,
        hjust   = 1.05,
        vjust   = 0,
        size    = 3.2,
        colour  = "grey30"
      ) +
      ggplot2::scale_fill_manual(values = pal) +
      ggplot2::labs(
        title = title,
        x     = group_var,
        y     = unit_label,
        fill  = group_var
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        strip.text       = ggplot2::element_text(face = "bold"),
        legend.position  = "bottom",
        plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold")
      )
  }
  
  # ------------------------------------------------------------------
  # 6. Return
  # ------------------------------------------------------------------
  # Clean output columns
  out_type <- dur_type[, c(sample_col, "type", group_var, "duration"), drop = FALSE]
  names(out_type)[names(out_type) == "duration"] <- paste0("duration_", duration_unit)
  
  out_total <- dur_total[, c(sample_col, group_var, "duration_total"), drop = FALSE]
  names(out_total)[names(out_total) == "duration_total"] <- paste0("duration_total_", duration_unit)
  
  list(
    duration_per_type = out_type,
    duration_total    = out_total,
    summary_table     = summary_df,
    plot              = p,
    params            = list(
      group_var     = group_var,
      duration_unit = duration_unit,
      exclude_types = exclude_types,
      min_n         = min_n,
      n_patients    = length(shared_ids),
      n_types       = length(types)
    )
  )
}
