# ============================================================
# MEC-TX analysis: tx_focus_dt()
# analysis/tx_focus_dt.R
# ============================================================

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
    # "all"         = no restriction
    # "single_only" = first_t0 must have exactly 1 type (no combo start)
    # "combo_only"  = first_t0 must have >=2 focus types simultaneously
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
  
  # ---- arg matching ----
  start_filter <- match.arg(start_filter)
  
  # valid_types defined locally — no global scope dependency (Bug 4.2 fix)
  local_valid_types <- c(
    "Ancillary", "Chemo", "Hormone", "IO",
    "Small_Molecule", "Targeted", "Radiation", "Others"
  )
  
  focus_types <- intersect(focus_types, local_valid_types)
  stopifnot(length(focus_types) >= 1)
  
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
      dplyr::select(sample, cluster = .data[[kc]]) %>%
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
  # search_all = TRUE: search across all k solutions (cl = NULL was passed)
  # search_all = FALSE: restrict to the single specified cluster
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
  
  # ---- relax threshold if not enough twins — dom_type_tx filter preserved ----
  thr <- min_share_tx
  while (length(sel) < min(n_twins, nrow(cand)) && thr > 0) {
    thr <- thr - 0.05
    add <- cand %>%
      dplyr::filter(
        dom_type_tx %in% focus_types,   # Bug fix: keep dominance filter in loop
        focus_share_tx >= thr
      ) %>%
      dplyr::pull(sample)
    sel <- unique(c(sel, add))
  }
  
  # ---- final fallback — dom_type_tx filter preserved ----
  if (length(sel) < min(n_twins, nrow(cand))) {
    sel <- unique(c(
      sel,
      cand %>%
        dplyr::filter(dom_type_tx %in% focus_types) %>%  # Bug fix: keep dominance filter
        dplyr::pull(sample)
    ))
  }
  
  twin_ids  <- head(sel, min(n_twins, nrow(cand)))
  df_twins  <- df_cl %>%
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
          ref_s <- cox_ref_levels[["stage"]] %||% cox_ref_levels[["Stage"]]
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
