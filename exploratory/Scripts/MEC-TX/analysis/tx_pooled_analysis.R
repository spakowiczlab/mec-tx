# ============================================================
# MEC-TX analysis: tx_pooled_analysis()
# analysis/tx_pooled_analysis.R
#
# Wraps the full cohort-building and visualization workflow
# for any treatment combination and any grouping variable.
#
# Modes:
#   "any"        — patients who ever received ALL focus_types
#   "only"       — patients who received focus_types and nothing else
#   "concurrent" — patients where focus_types overlapped in time
#   "dominant"   — cluster-aware; focus_types dominate treatment time
# ============================================================

tx_pooled_analysis <- function(
  Cluster_surv,
  timeline,

  # cohort definition
  focus_types   = c("Chemo", "IO"),
  mode          = c("any", "only", "concurrent", "dominant"),

  # grouping variable (any categorical column in Cluster_surv)
  group_var     = "CAlevel",

  # shared options
  horizon_years = 5,
  min_share_tx  = 0.33,      # dominant mode only

  # concurrent mode: minimum overlap in years to count as concurrent
  min_overlap   = 0,

  # dominant mode: passed to tx_focus_dt
  n_twins       = 999,
  enforce_sequence = FALSE,
  sequence_strict  = FALSE,
  start_filter     = "all",
  pure_focus_only  = FALSE,

  # survival analysis
  cox_covars    = c("stage_group", "sex", "age", "smokingstatus"),
  ref_levels    = NULL,       # NULL = auto-detect from data
  numeric_scale = list(age = 5),
  numeric_units = list(age = "years"),
  min_epv       = 5,

  # visualization
  show_timeline = TRUE,       # include treatment timeline strip
  group_colours = NULL,       # NULL = auto-generated
  horizon_plot  = NULL,       # NULL = uses horizon_years

  # layout
  base_size  = 14,
  title_size = 16,
  widths     = c(1.4, 1, 1)
) {

  # ---- arg matching ----
  mode <- match.arg(mode)
  if (is.null(horizon_plot)) horizon_plot <- horizon_years

  # local valid types — no global dependency
  local_valid_types <- c(
    "Ancillary", "Chemo", "Hormone", "IO",
    "Small_Molecule", "Targeted", "Radiation", "Others"
  )

  focus_types <- intersect(focus_types, local_valid_types)
  stopifnot(
    "focus_types must have at least one valid type" = length(focus_types) >= 1,
    "group_var must be a single string"             = length(group_var) == 1
  )

  # resolve group_var case-insensitively
  group_var_actual <- names(Cluster_surv)[
    tolower(names(Cluster_surv)) == tolower(group_var)
  ]
  if (length(group_var_actual) != 1) {
    stop(sprintf(
      "group_var '%s' not found in Cluster_surv. Available: %s",
      group_var, paste(names(Cluster_surv), collapse = ", ")
    ))
  }
  group_var <- group_var_actual

  focus_label <- paste(focus_types, collapse = "+")

  # =========================================================
  # STEP 1 — Find cohort IDs based on mode
  # =========================================================
  segs_prep <- prep_segs(timeline, horizon_years = horizon_years)

  cohort_ids <- switch(mode,

    # ---- any: received ALL focus_types at any point ----
    "any" = {
      presence <- segs_prep %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(
          has_all = all(focus_types %in% unique(as.character(type))),
          .groups = "drop"
        ) %>%
        dplyr::filter(has_all)
      presence$sample
    },

    # ---- only: received focus_types and nothing else (Ancillary allowed) ----
    "only" = {
      other_types <- setdiff(local_valid_types, c(focus_types, "Ancillary"))
      presence <- segs_prep %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(
          has_all_focus = all(focus_types %in% unique(as.character(type))),
          has_other     = any(as.character(type) %in% other_types),
          .groups = "drop"
        ) %>%
        dplyr::filter(has_all_focus & !has_other)
      presence$sample
    },

    # ---- concurrent: focus_types overlap in time ----
    "concurrent" = {
      # for each pair of focus types, check if their segments overlap
      if (length(focus_types) < 2) {
        stop("concurrent mode requires at least 2 focus_types")
      }

      # get segments for each focus type
      seg_list <- lapply(focus_types, function(ft) {
        segs_prep %>%
          dplyr::filter(as.character(type) == ft) %>%
          dplyr::select(sample, t0, t1) %>%
          dplyr::rename(t0_a = t0, t1_a = t1)
      })

      # find patients with overlapping intervals across all focus_type pairs
      concurrent_ids <- segs_prep %>%
        dplyr::distinct(sample) %>%
        dplyr::pull(sample)

      # check all pairs
      pairs <- combn(length(focus_types), 2, simplify = FALSE)
      for (pair in pairs) {
        fa <- focus_types[pair[1]]
        fb <- focus_types[pair[2]]

        segs_a <- segs_prep %>%
          dplyr::filter(as.character(type) == fa) %>%
          dplyr::select(sample, t0_a = t0, t1_a = t1)

        segs_b <- segs_prep %>%
          dplyr::filter(as.character(type) == fb) %>%
          dplyr::select(sample, t0_b = t0, t1_b = t1)

        # cross join within sample and check overlap
        overlap_ids <- segs_a %>%
          dplyr::inner_join(segs_b, by = "sample") %>%
          dplyr::mutate(
            overlap = pmin(t1_a, t1_b) - pmax(t0_a, t0_b)
          ) %>%
          dplyr::filter(overlap > min_overlap) %>%
          dplyr::distinct(sample) %>%
          dplyr::pull(sample)

        concurrent_ids <- intersect(concurrent_ids, overlap_ids)
      }
      concurrent_ids
    },

    # ---- dominant: cluster-aware via tx_focus_dt across all k ----
    "dominant" = {
      ids <- lapply(paste0("Cluster_k", 3:20), function(kc) {
        # only run if kc exists in Cluster_surv
        if (!kc %in% names(Cluster_surv)) return(NULL)
        tryCatch({
          demo <- tx_focus_dt(
            Cluster_surv  = Cluster_surv,
            segs          = timeline,
            kc            = kc,
            cl            = NULL,
            focus_types   = focus_types,
            group_col     = group_var,
            horizon_years = horizon_years,
            n_twins       = n_twins,
            min_share_tx  = min_share_tx,
            enforce_sequence = enforce_sequence,
            sequence_strict  = sequence_strict,
            start_filter     = start_filter,
            pure_focus_only  = pure_focus_only,
            add_forest       = FALSE
          )
          attr(demo, "twin_ids")
        }, error = function(e) NULL)
      })
      unique(unlist(ids))
    }
  )

  if (length(cohort_ids) == 0) {
    stop(sprintf(
      "No patients found for focus_types = c(%s) with mode = '%s'.",
      paste0('"', focus_types, '"', collapse = ", "), mode
    ))
  }

  message(sprintf(
    "[tx_pooled_analysis] Mode='%s' | focus=%s | n=%d patients found",
    mode, focus_label, length(cohort_ids)
  ))

  # =========================================================
  # STEP 2 — Build cohort data objects
  # =========================================================
  df_cohort <- Cluster_surv %>%
    dplyr::filter(sample %in% cohort_ids)

  segs_cohort  <- segs_prep %>%
    dplyr::filter(sample %in% cohort_ids)

  share_cohort <- treatment_shares(segs_cohort)

  # =========================================================
  # STEP 3 — Summary stats
  # =========================================================
  message(sprintf("  n in Cluster_surv: %d", nrow(df_cohort)))
  message(sprintf("  %s distribution:", group_var))
  print(table(df_cohort[[group_var]], useNA = "ifany"))

  if ("stage_group" %in% names(df_cohort)) {
    message("  Stage distribution:")
    print(table(df_cohort$stage_group, useNA = "ifany"))
  }

  # =========================================================
  # STEP 4 — Build timeline plot data
  # =========================================================
  # order_tbl: patients with focus_types as dominant treatment
  order_tbl <- share_cohort %>%
    dplyr::select(sample, dom_type_tx) %>%
    dplyr::filter(dom_type_tx %in% focus_types) %>%
    dplyr::left_join(
      segs_cohort %>%
        dplyr::group_by(sample) %>%
        dplyr::summarise(first_start = min(t0), .groups = "drop"),
      by = "sample"
    ) %>%
    dplyr::arrange(dom_type_tx, first_start) %>%
    dplyr::mutate(y = dplyr::row_number())

  # df_plot: segments for focus_types only
  df_plot <- segs_cohort %>%
    dplyr::filter(as.character(type) %in% focus_types) %>%
    dplyr::left_join(order_tbl, by = "sample") %>%
    dplyr::filter(!is.na(y)) %>%
    dplyr::left_join(
      df_cohort %>% dplyr::select(sample, dplyr::all_of(group_var)),
      by = "sample"
    ) %>%
    dplyr::mutate(type = factor(type, levels = focus_types))

  # group_strip: sidebar coloured by group_var
  group_strip <- order_tbl %>%
    dplyr::left_join(
      df_cohort %>% dplyr::select(sample, dplyr::all_of(group_var)),
      by = "sample"
    )

  n_timeline_patients <- dplyr::n_distinct(df_plot$sample)

  # =========================================================
  # STEP 5 — Colours for group_var
  # =========================================================
  grp_levels <- sort(unique(as.character(df_cohort[[group_var]])))
  grp_levels <- grp_levels[!is.na(grp_levels)]

  if (is.null(group_colours)) {
    # use ca_cols for CAlevel, otherwise generate palette
    if (tolower(group_var) == "calevel" &&
        all(c("High", "Low") %in% grp_levels)) {
      group_colours <- c(High = "#F28E2B", Low = "#56B4E9")
    } else {
      base_pal <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2",
                    "#59A14F", "#EDC948", "#B07AA1", "#FF9DA7",
                    "#9C755F", "#BAB0AC")
      if (length(grp_levels) <= length(base_pal)) {
        pal <- base_pal[seq_along(grp_levels)]
      } else {
        # more than 10 levels — generate using hcl palette
        message(sprintf(
          "[tx_pooled_analysis] %d levels detected for '%s' — using hcl.colors() palette.",
          length(grp_levels), group_var
        ))
        pal <- grDevices::hcl.colors(length(grp_levels), palette = "Dark 3")
      }
      group_colours <- stats::setNames(pal, grp_levels)
    }
  }

  # auto ref_levels if not supplied
  if (is.null(ref_levels)) {
    ref_levels <- stats::setNames(
      list(grp_levels[1]),
      group_var
    )
  }

  # =========================================================
  # STEP 6 — KM plot
  # =========================================================
  p_km <- km_panel_from_df(
    df_cohort %>% dplyr::filter(sample %in% df_plot$sample),
    group_col     = group_var,
    title         = sprintf("All %s — all clusters pooled (n=%d)",
                            focus_label, n_timeline_patients),
    horizon_years = horizon_plot,
    risk_table    = TRUE,
    group_colours = group_colours
  )

  # =========================================================
  # STEP 7 — Cox forest plot
  # =========================================================
  all_covars <- unique(c(group_var, cox_covars))
  all_covars <- intersect(all_covars, names(df_cohort))

  p_forest <- cox_forest_plot_from_df(
    df_cohort %>% dplyr::filter(sample %in% df_plot$sample),
    covars        = all_covars,
    ref_levels    = ref_levels,
    title         = sprintf("Adjusted Cox — %s pooled (all clusters)", focus_label),
    min_epv       = min_epv,
    priority      = all_covars,
    numeric_scale = numeric_scale,
    numeric_units = numeric_units,
    base_size     = base_size,
    title_size    = title_size
  )

  # =========================================================
  # STEP 8 — Timeline strip plot
  # =========================================================
  if (show_timeline) {

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

    strip_x <- horizon_plot + 0.4

    p_timeline <- ggplot2::ggplot() +
      # group_var sidebar
      ggplot2::geom_tile(
        data = group_strip %>%
          dplyr::mutate(grp = factor(.data[[group_var]], levels = grp_levels)),
        ggplot2::aes(x = strip_x, y = y, fill = grp),
        width = 0.4, height = 0.9, alpha = 0.8
      ) +
      ggplot2::scale_fill_manual(
        values = group_colours,
        name   = group_var,
        na.value = "grey80"
      ) +
      ggnewscale::new_scale_colour() +
      # treatment segments
      ggplot2::geom_segment(
        data = df_plot,
        ggplot2::aes(x = t0, xend = t1, y = y, yend = y, colour = type),
        linewidth = 0.5
      ) +
      ggplot2::scale_colour_manual(
        values = local_tx_cols[focus_types],
        name   = "Treatment Type",
        drop   = FALSE
      ) +
      # dominant type label on left
      ggplot2::geom_text(
        data = order_tbl %>%
          dplyr::group_by(dom_type_tx) %>%
          dplyr::summarise(y_mid = mean(y), .groups = "drop"),
        ggplot2::aes(x = -0.3, y = y_mid, label = dom_type_tx),
        hjust = 1, size = 3.5, fontface = "bold"
      ) +
      ggplot2::scale_x_continuous(
        limits = c(-0.5, strip_x + 0.3),
        breaks = 0:horizon_plot,
        name   = "Years since first treatment"
      ) +
      ggplot2::scale_y_continuous(
        name   = sprintf("Patients (n=%d)", n_timeline_patients),
        breaks = NULL
      ) +
      ggplot2::labs(
        title    = sprintf(
          "Treatment Timelines — %s Dominant Patients (n=%d)",
          focus_label, n_timeline_patients
        ),
        subtitle = sprintf(
          "Ordered by dominant treatment type | Right strip = %s", group_var
        )
      ) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        plot.title         = ggplot2::element_text(size = title_size, face = "bold"),
        plot.subtitle      = ggplot2::element_text(size = base_size - 3),
        axis.text.y        = ggplot2::element_blank(),
        axis.ticks.y       = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        legend.position    = "right"
      )
  } else {
    p_timeline <- NULL
  }

  # =========================================================
  # STEP 9 — Return
  # =========================================================
  list(
    # plots
    km            = p_km,
    forest        = p_forest,
    timeline      = p_timeline,

    # data
    ids           = cohort_ids,
    df            = df_cohort,
    segs          = segs_cohort,
    shares        = share_cohort,
    df_plot       = df_plot,

    # metadata
    mode          = mode,
    focus_types   = focus_types,
    group_var     = group_var,
    n_cohort      = length(cohort_ids),
    n_plot        = n_timeline_patients,
    group_table   = table(df_cohort[[group_var]], useNA = "ifany")
  )
}
