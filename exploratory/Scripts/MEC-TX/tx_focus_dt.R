# ============================================================
# MEC-TX analysis: tx_focus_dt()
# analysis/tx_focus_dt.R
# ============================================================

#' Digital Twin Focus Analysis for a Specific Treatment Type
#'
#' Identifies patients ("twins") whose treatment history is dominated by
#' one or more specified treatment types, then produces a three-panel
#' composite figure: a treatment timeline swimlane, a Kaplan-Meier
#' survival panel, and an adjusted Cox forest plot. Designed for deep
#' inspection of a treatment-type-specific subgroup, optionally restricted
#' to a single k-means cluster.
#'
#' @param Cluster_surv A data frame — the \code{$Cluster_surv} slot from
#'   \code{\link{tx_cluster_surv}}. Must contain \code{sample}, the column
#'   named in \code{kc}, and all columns in \code{cox_covars}.
#' @param segs A data frame — the interval output of
#'   \code{\link{tx_intervals}}. Must contain \code{sample}, \code{type},
#'   \code{start_year}, and \code{end_year}.
#' @param kc Character string. Name of the cluster assignment column in
#'   \code{Cluster_surv} to use for subsetting (e.g. \code{"Cluster_k14"}).
#'   Default \code{"Cluster_k14"}.
#' @param cl Integer or \code{NULL}. Cluster value to restrict twin
#'   selection to. \code{NULL} auto-selects the cluster with the highest
#'   median focus-type share across all patients. Default \code{NULL}.
#' @param focus_types Character vector. Treatment type(s) defining the
#'   focus cohort. Must be canonical MEC-TX type labels. Patients are
#'   ranked by combined share of these types. Default \code{"Radiation"}.
#' @param group_col Character string. Grouping variable for the KM and
#'   Cox panels. Matched case-insensitively. Default \code{"CAlevel"}.
#' @param horizon_years Numeric. Timeline and KM x-axis extent in years.
#'   Default \code{5}.
#' @param n_twins Integer. Maximum number of twins to display in the
#'   timeline and include in KM/Cox panels. Set to \code{999} to keep all.
#'   Default \code{20}.
#' @param min_share_tx Numeric in \code{[0, 1]}. Minimum combined
#'   \code{focus_types} share threshold for twin selection. Automatically
#'   relaxed in 0.05 decrements if fewer than \code{n_twins} qualify at
#'   the initial threshold. Default \code{0.33}.
#' @param enforce_sequence Logical. If \code{TRUE}, only patients whose
#'   treatment sequence contains \code{seq_pattern} as a contiguous
#'   subsequence are included. Default \code{FALSE}.
#' @param seq_pattern Character vector or \code{NULL}. Subsequence to
#'   enforce when \code{enforce_sequence = TRUE}. \code{NULL} defaults to
#'   \code{focus_types}. Default \code{NULL}.
#' @param sequence_strict Logical. If \code{TRUE}, sequence matching
#'   requires strict temporal ordering of first occurrences (no
#'   re-ordering). If \code{FALSE}, uses \code{\link{has_subseq}} for
#'   contiguous subsequence matching. Default \code{FALSE}.
#' @param start_filter One of \code{"all"}, \code{"single_only"}, or
#'   \code{"combo_only"}. Filters patients by what treatment type(s) they
#'   started with. \code{"single_only"} excludes patients who started two
#'   or more types simultaneously. \code{"combo_only"} requires at least
#'   two focus types to start together. Default \code{"all"}.
#' @param pure_focus_only Logical. If \code{TRUE}, only patients whose
#'   entire treatment sequence consists exclusively of \code{focus_types}
#'   are included. Default \code{FALSE}.
#' @param add_forest Logical. If \code{TRUE}, an adjusted Cox forest plot
#'   is added as the third panel. Default \code{TRUE}.
#' @param cox_covars Character vector. Covariates for the Cox model.
#'   Default \code{c("CAlevel", "stage_group", "sex", "age", "smokingstatus")}.
#' @param cox_ref_levels Named list. Reference levels for Cox covariates.
#'   Default \code{list(CAlevel = "Low", stage_group = "Local", smokingstatus = "Never")}.
#' @param forest_min_epv Integer. Minimum EPV for Cox covariate selection
#'   in the forest panel. Default \code{1} (relaxed — twin cohorts are
#'   small).
#' @param forest_priority Character vector. Covariate priority order for
#'   EPV-based selection. Default matches \code{cox_covars}.
#' @param forest_numeric_scale Named list. Divisors for numeric covariates.
#'   Default \code{list(age = 5)}.
#' @param forest_numeric_units Named list. Unit labels for numeric
#'   covariates. Default \code{list(age = "years")}.
#' @param forest_drop_stage_unknown Logical. If \code{TRUE}, patients with
#'   stage values in \code{forest_stage_unknown_levels} are excluded from
#'   the Cox forest only (not from KM). Default \code{FALSE}.
#' @param forest_stage_unknown_levels Character vector. Stage values
#'   treated as unknown when \code{forest_drop_stage_unknown = TRUE}.
#'   Default \code{c("Unknown", "Unknown/Not Applicable", "Not Applicable")}.
#' @param km_risk_table Logical. Show at-risk table below KM panel.
#'   Default \code{TRUE}.
#' @param km_risk_table_height Numeric. Relative height of risk table
#'   panel. Default \code{0.26}.
#' @param show_km_legend Logical. If \code{FALSE}, the KM legend is
#'   suppressed to save space in the composite layout. Default \code{FALSE}.
#' @param base_size Base font size (pt) applied to all three panels.
#'   Default \code{14}.
#' @param title_size Font size (pt) for panel titles. Default \code{16}.
#' @param widths Numeric vector of length 3. Relative widths of the
#'   timeline, KM, and forest panels in the composite layout. Default
#'   \code{c(1.4, 1, 1)}.
#'
#' @return A \code{patchwork} object combining the three panels side by
#'   side. The selected twin sample IDs are attached as
#'   \code{attr(result, "twin_ids")} for downstream access without
#'   re-running the function.
#'
#' @details
#' \strong{Twin selection pipeline:}
#' \enumerate{
#'   \item Compute treatment type shares via \code{\link{treatment_shares}}.
#'   \item Restrict to the specified cluster (\code{cl}) or auto-select
#'     the cluster with highest median focus share.
#'   \item Apply optional filters: \code{start_filter},
#'     \code{pure_focus_only}, \code{enforce_sequence}.
#'   \item Rank candidates by combined focus share, then by total treated
#'     time.
#'   \item Apply \code{min_share_tx} dominance threshold, relaxing in 0.05
#'     steps if fewer than \code{n_twins} qualify.
#'   \item Take the top \code{n_twins} by focus share.
#' }
#'
#' \strong{Accessing twin IDs:}
#' \code{attr(result, "twin_ids")} returns the sample IDs used in the
#' figure without re-running the full pipeline.
#'
#' \strong{Forest EPV:} Twin cohorts are small by design. \code{forest_min_epv = 1}
#' is the default to prevent all covariates from being dropped. Increase
#' to \code{5} for larger cohorts.
#'
#' \strong{Parameter naming note:} \code{forest_numeric_scale} defaults
#' to \code{list(age = 5)} here (per 5-year increment) rather than the
#' package-wide default of \code{list(age = 10)}, reflecting that twin
#' cohorts are smaller and the 5-year increment gives more estimable HRs.
#'
#' @examples
#' \dontrun{
#' intervals <- tx_intervals(med_data, Cluster_surv)
#' segs      <- intervals$timeline_long_intv
#'
#' # Radiation-focused twins, auto-select cluster
#' p <- tx_focus_dt(
#'   Cluster_surv = res$Cluster_surv,
#'   segs         = segs,
#'   kc           = "Cluster_k14",
#'   focus_types  = "Radiation"
#' )
#' pdf(file.path(out_dir, "radiation_twins_k14.pdf"), width = 18, height = 8)
#' print(p)
#' dev.off()
#'
#' # Retrieve twin IDs for downstream analysis
#' twin_ids <- attr(p, "twin_ids")
#'
#' # Chemo+IO twins in a specific cluster, with sequence enforcement
#' p2 <- tx_focus_dt(
#'   Cluster_surv     = res$Cluster_surv,
#'   segs             = segs,
#'   kc               = "Cluster_k8",
#'   cl               = 3,
#'   focus_types      = c("Chemo", "IO"),
#'   enforce_sequence = TRUE,
#'   seq_pattern      = c("Chemo", "IO")
#' )
#' }
#'
#' @seealso \code{\link{tx_cluster_surv}}, \code{\link{tx_intervals}},
#'   \code{\link{timeline_panel}}, \code{\link{km_panel_from_df}},
#'   \code{\link{cox_forest_plot_from_df}}, \code{\link{has_subseq}}
#'
#' @import ggplot2
#' @importFrom dplyr select filter arrange group_by summarise mutate
#'   left_join semi_join pull slice_head across all_of desc
#' @importFrom patchwork plot_layout
#' @importFrom tibble tibble
#' @importFrom stats relevel median
#' @export

tx_focus_dt <- function(
    Cluster_surv, segs,
    kc            = "Cluster_k14",
    cl            = NULL,
    focus_types   = c("Radiation"),
    group_col     = "CAlevel",
    horizon_years = 5,
    n_twins       = 20,
    min_share_tx  = 0.33,
    
    # sequence options
    enforce_sequence = FALSE,
    seq_pattern      = NULL,
    sequence_strict  = FALSE,
    
    # start-of-treatment filter
    start_filter     = c("all", "single_only", "combo_only"),
    pure_focus_only  = FALSE,
    
    # forest controls
    add_forest                  = TRUE,
    cox_covars                  = c("CAlevel", "stage_group", "sex", "age",
                                    "smokingstatus"),
    cox_ref_levels              = list(CAlevel = "Low", stage_group = "Local",
                                       smokingstatus = "Never"),
    forest_min_epv              = 1,
    forest_priority             = c("CAlevel", "stage_group", "sex", "age",
                                    "smokingstatus"),
    forest_numeric_scale        = list(age = 5),
    forest_numeric_units        = list(age = "years"),
    forest_drop_stage_unknown   = FALSE,
    forest_stage_unknown_levels = c("Unknown", "Unknown/Not Applicable",
                                    "Not Applicable"),
    
    # KM controls
    km_risk_table        = TRUE,
    km_risk_table_height = 0.26,
    show_km_legend       = FALSE,
    
    # layout
    base_size  = 14,
    title_size = 16,
    widths     = c(1.4, 1, 1)
) {
  
  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
  
  # --- 1. Cluster_surv must be a data frame ---
  if (!is.data.frame(Cluster_surv)) {
    stop(
      "tx_focus_dt(): 'Cluster_surv' must be a data frame.\n",
      "  → You passed an object of class: ", paste(class(Cluster_surv), collapse = ", "), ".\n",
      "  → Pass the '$Cluster_surv' element from tx_cluster_surv() output:\n",
      "    res <- tx_cluster_surv(...)\n",
      "    tx_focus_dt(res$Cluster_surv, ...)"
    )
  }
  
  # --- 2. segs must be a data frame with required columns ---
  if (!is.data.frame(segs)) {
    stop(
      "tx_focus_dt(): 'segs' must be a data frame.\n",
      "  → You passed an object of class: ", paste(class(segs), collapse = ", "), ".\n",
      "  → Pass the output of tx_intervals():\n",
      "    tx_focus_dt(Cluster_surv, tx_intervals(norm), ...)"
    )
  }
  
  required_segs <- c("sample", "type", "start_year", "end_year")
  missing_segs  <- setdiff(required_segs, names(segs))
  if (length(missing_segs) > 0) {
    stop(
      "tx_focus_dt(): Required column(s) missing from 'segs':\n",
      paste0("  ✗ ", missing_segs, collapse = "\n"), "\n",
      "  → These columns are produced by tx_intervals().\n",
      "  → Make sure you pass tx_intervals() output, not tx_normalize() output."
    )
  }
  
  # --- 3. kc column must exist in Cluster_surv ---
  if (!kc %in% names(Cluster_surv)) {
    available_k <- grep("^Cluster_k\\d+$", names(Cluster_surv), value = TRUE)
    stop(
      "tx_focus_dt(): Cluster column '", kc, "' not found in 'Cluster_surv'.\n",
      "  → Available cluster columns: ",
      if (length(available_k) > 0) paste(available_k, collapse = ", ")
      else "none found", "\n",
      "  → Example: kc = '", if (length(available_k) > 0) available_k[1] else "Cluster_k5", "'"
    )
  }
  
  # --- 4. focus_types must have at least one valid type ---
  local_valid_types <- c(
    "Ancillary", "Chemo", "Hormone", "IO",
    "Small_Molecule", "Targeted", "Radiation", "Others"
  )
  valid_focus <- intersect(focus_types, local_valid_types)
  if (length(valid_focus) == 0) {
    stop(
      "tx_focus_dt(): None of the supplied 'focus_types' are valid.\n",
      "  → You supplied: ", paste(focus_types, collapse = ", "), "\n",
      "  → Valid types: ", paste(local_valid_types, collapse = ", ")
    )
  }
  if (length(valid_focus) < length(focus_types)) {
    warning(
      "tx_focus_dt(): Ignoring unrecognised focus_types: ",
      paste(setdiff(focus_types, local_valid_types), collapse = ", "), "\n",
      "  → Valid types: ", paste(local_valid_types, collapse = ", ")
    )
  }
  
  # --- 5. group_col must exist in Cluster_surv (case-insensitive) ---
  actual_gc <- names(Cluster_surv)[tolower(names(Cluster_surv)) == tolower(group_col)]
  if (length(actual_gc) != 1) {
    stop(
      "tx_focus_dt(): group_col '", group_col, "' not found in 'Cluster_surv'.\n",
      "  → Columns available: ", paste(names(Cluster_surv), collapse = ", "), "\n",
      "  → group_col is matched case-insensitively."
    )
  }
  
  # --- 6. cl must be a valid cluster level if supplied ---
  if (!is.null(cl)) {
    valid_cls <- unique(Cluster_surv[[kc]])
    if (!cl %in% valid_cls) {
      stop(
        "tx_focus_dt(): cl = ", cl, " is not a valid cluster for '", kc, "'.\n",
        "  → Valid cluster values: ", paste(sort(valid_cls), collapse = ", ")
      )
    }
  }
  
  # --- 7. Sample overlap between Cluster_surv and segs ---
  cs_samples <- unique(as.character(Cluster_surv$sample))
  sg_samples <- unique(as.character(segs$sample))
  overlap    <- intersect(cs_samples, sg_samples)
  if (length(overlap) == 0) {
    stop(
      "tx_focus_dt(): No sample IDs match between 'Cluster_surv' and 'segs'.\n",
      "  → Cluster_surv samples (first 5): ", paste(head(cs_samples, 5), collapse = ", "), "\n",
      "  → segs samples (first 5):         ", paste(head(sg_samples, 5), collapse = ", "), "\n",
      "  → Make sure both come from the same dataset."
    )
  }
  
  # --- 8. n_twins must be a positive integer ---
  if (!is.numeric(n_twins) || length(n_twins) != 1 || n_twins < 1) {
    stop(
      "tx_focus_dt(): 'n_twins' must be a single positive integer.\n",
      "  → You passed: ", deparse(n_twins), "\n",
      "  → Example: n_twins = 20 (default) or n_twins = 999 to keep all."
    )
  }
  
  # --- 9. min_share_tx must be between 0 and 1 ---
  if (!is.numeric(min_share_tx) || length(min_share_tx) != 1 ||
      min_share_tx < 0 || min_share_tx > 1) {
    stop(
      "tx_focus_dt(): 'min_share_tx' must be a number between 0 and 1.\n",
      "  → You passed: ", deparse(min_share_tx), "\n",
      "  → Example: min_share_tx = 0.33 (default — focus types must make up\n",
      "    at least 33% of treatment time)"
    )
  }
  
  # --- 10. enforce_sequence + seq_pattern consistency ---
  if (enforce_sequence && !is.null(seq_pattern)) {
    bad_seq <- setdiff(seq_pattern, local_valid_types)
    if (length(bad_seq) > 0) {
      warning(
        "tx_focus_dt(): seq_pattern contains unrecognised treatment types: ",
        paste(bad_seq, collapse = ", "), "\n",
        "  → These will never match and no patients will pass the sequence filter."
      )
    }
  }
  
  # --- 11. widths must have exactly 3 positive values ---
  if (!is.numeric(widths) || length(widths) != 3 || any(widths <= 0)) {
    stop(
      "tx_focus_dt(): 'widths' must be a numeric vector of length 3 with all positive values.\n",
      "  → You passed: ", deparse(widths), "\n",
      "  → Example: widths = c(1.4, 1, 1)  (timeline | KM | forest)"
    )
  }
  
  # ===========================================================================
  # MAIN FUNCTION BODY
  # ===========================================================================
  
  # ---- arg matching ----
  start_filter <- match.arg(start_filter)
  
  focus_types <- valid_focus
  
  # ---- prep & shares ----
  segs_prep <- prep_segs(segs, horizon_years = horizon_years)
  share_df  <- treatment_shares(segs_prep)
  
  # ---- per-sample sequences & first-time info ----
  seq_tbl <- segs_prep %>%
    dplyr::arrange(sample, t0) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      tx_seq = list({
        v <- as.character(type); rle(v)$values
      }),
      first_time = list({
        sp <- split(t0, type); vapply(sp, min, numeric(1))
      }),
      first_t0    = min(t0),
      types_first = list(unique(type[t0 == first_t0])),
      .groups = "drop"
    )
  
  if (enforce_sequence && is.null(seq_pattern)) seq_pattern <- focus_types
  
  # ---- capture whether cl was user-supplied BEFORE auto-assignment ----
  search_all <- is.null(cl)
  
  # ---- auto-select best cluster if cl not given ----
  if (is.null(cl)) {
    cl <- Cluster_surv %>%
      dplyr::select(sample, cluster = dplyr::all_of(kc)) %>%
      dplyr::left_join(share_df, by = "sample") %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(
        med_focus = median(
          rowSums(dplyr::across(
            dplyr::all_of(paste0("share_", focus_types, "_tx"))
          )),
          na.rm = TRUE
        ),
        .groups = "drop"
      ) %>%
      dplyr::arrange(dplyr::desc(med_focus)) %>%
      dplyr::slice_head(n = 1) %>%
      dplyr::pull(cluster)
  }
  
  # ---- restrict to cluster (or all patients if search_all) ----
  if (search_all) {
    df_cl <- Cluster_surv
  } else {
    df_cl <- Cluster_surv %>%
      dplyr::filter(.data[[kc]] == cl)
  }
  
  # ---- build candidate twin set ----
  cand <- df_cl %>%
    dplyr::select(sample) %>%
    dplyr::left_join(share_df, by = "sample") %>%
    dplyr::left_join(seq_tbl,  by = "sample") %>%
    dplyr::mutate(
      focus_share_tx = rowSums(
        dplyr::across(dplyr::all_of(paste0("share_", focus_types, "_tx"))),
        na.rm = TRUE
      ),
      combo_start_any = vapply(
        types_first,
        function(tt) length(unique(tt)) >= 2L, logical(1)
      ),
      combo_focus_start = vapply(
        types_first,
        function(tt) { tt <- as.character(tt); sum(tt %in% focus_types) >= 2L },
        logical(1)
      ),
      only_focus_types = vapply(
        tx_seq,
        function(v) all(as.character(v) %in% focus_types), logical(1)
      ),
      seq_ok = if (!enforce_sequence || is.null(seq_pattern)) {
        TRUE
      } else if (!sequence_strict) {
        vapply(tx_seq, has_subseq, logical(1), pattern = seq_pattern)
      } else {
        vapply(first_time, function(ft) {
          times <- ft[seq_pattern]
          if (any(is.na(times)) || anyDuplicated(times)) return(FALSE)
          all(diff(times) > 0)
        }, logical(1))
      }
    ) %>%
    dplyr::filter(
      dur_treated > 0,
      seq_ok,
      dplyr::case_when(
        start_filter == "all"         ~ TRUE,
        start_filter == "single_only" ~ !combo_start_any,
        start_filter == "combo_only"  ~ combo_focus_start,
        TRUE                          ~ TRUE
      ),
      if (pure_focus_only) only_focus_types else TRUE
    ) %>%
    dplyr::arrange(dplyr::desc(focus_share_tx), dplyr::desc(dur_treated))
  
  # ---- dominance filter ----
  sel <- cand %>%
    dplyr::filter(
      dom_type_tx %in% focus_types,
      focus_share_tx >= min_share_tx
    ) %>%
    dplyr::pull(sample) %>%
    unique()
  
  # ---- relax threshold if not enough twins ----
  thr <- min_share_tx
  while (length(sel) < min(n_twins, nrow(cand)) && thr > 0) {
    thr <- thr - 0.05
    add <- cand %>%
      dplyr::filter(
        dom_type_tx %in% focus_types,
        focus_share_tx >= thr
      ) %>%
      dplyr::pull(sample)
    sel <- unique(c(sel, add))
  }
  
  # ---- final fallback ----
  if (length(sel) < min(n_twins, nrow(cand))) {
    sel <- unique(c(
      sel,
      cand %>%
        dplyr::filter(dom_type_tx %in% focus_types) %>%
        dplyr::pull(sample)
    ))
  }
  
  twin_ids    <- head(sel, min(n_twins, nrow(cand)))
  df_twins    <- df_cl %>%
    dplyr::semi_join(tibble::tibble(sample = twin_ids), by = "sample")
  focus_label <- paste(focus_types, collapse = "+")
  
  # helper to apply consistent sizes to any ggplot
  apply_sizes <- function(p) {
    p + ggplot2::theme(
      text       = ggplot2::element_text(size = base_size),
      plot.title = ggplot2::element_text(size = title_size, face = "bold")
    )
  }
  
  # ---- timeline panel ----
  p_timeline <- timeline_panel(
    segs_prep, share_df, twin_ids,
    title = sprintf(
      "Top-%d twins (%s focus) in %s cl%s",
      length(twin_ids), focus_label, k_label(kc), cl
    ),
    horizon_years = horizon_years
  ) |> apply_sizes()
  
  # ---- KM panel ----
  p_km <- km_panel_from_df(
    df_twins,
    group_col         = group_col,
    title             = sprintf(
      "KM in %s (cluster %s) — twins only (n=%d, focus=%s)",
      k_label(kc), cl, length(twin_ids), focus_label
    ),
    horizon_years     = horizon_years,
    risk_table        = km_risk_table,
    risk_times        = 0:horizon_years,
    risk_table_height = km_risk_table_height
  ) |> apply_sizes()
  
  if (!show_km_legend) {
    p_km <- p_km + ggplot2::theme(legend.position = "none")
  }
  
  # ---- forest panel ----
  if (add_forest) {
    df_forest <- df_twins
    if (forest_drop_stage_unknown) {
      stage_actual <- names(df_forest)[tolower(names(df_forest)) == "stage"]
      if (length(stage_actual) == 1) {
        df_forest <- df_forest %>%
          dplyr::mutate(!!stage_actual := as.character(.data[[stage_actual]])) %>%
          dplyr::filter(
            !is.na(.data[[stage_actual]]),
            !(.data[[stage_actual]] %in% forest_stage_unknown_levels)
          )
        if (nrow(df_forest)) {
          df_forest[[stage_actual]] <- factor(df_forest[[stage_actual]])
          ref_s <- if (!is.null(cox_ref_levels[["stage"]])) {
            cox_ref_levels[["stage"]]
          } else {
            cox_ref_levels[["Stage"]]
          }
          if (!is.null(ref_s) && ref_s %in% levels(df_forest[[stage_actual]])) {
            df_forest[[stage_actual]] <- stats::relevel(
              df_forest[[stage_actual]], ref = ref_s
            )
          }
        }
      }
    }
    
    p_forest <- cox_forest_plot_from_df(
      df_forest,
      covars        = cox_covars,
      ref_levels    = cox_ref_levels,
      title         = sprintf(
        "Adjusted Cox — twins only (focus=%s%s)",
        focus_label,
        if (forest_drop_stage_unknown) ", Stage unknown excluded" else ""
      ),
      min_epv       = forest_min_epv,
      priority      = forest_priority,
      numeric_scale = forest_numeric_scale,
      numeric_units = forest_numeric_units
    ) |> apply_sizes()
    
  } else {
    p_forest <- ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::labs(title = "Forest (omitted)") |>
      apply_sizes()
  }
  
  # ---- compose row ----
  p_row <- (p_timeline | p_km | p_forest) +
    patchwork::plot_layout(widths = widths, guides = "keep")
  
  attr(p_row, "twin_ids") <- twin_ids
  p_row
}