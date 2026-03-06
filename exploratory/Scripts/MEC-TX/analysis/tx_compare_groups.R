# ============================================================
# MEC-TX analysis: tx_compare_groups()
# analysis/tx_compare_groups.R
#
# Flexible KM + adjusted Cox forest for any grouping variable
# or custom sample ID lists. Works with any categorical column
# in Cluster_surv — not limited to CAlevel.
#
# Usage examples:
#   # Compare by smoking status in Chemo+IO cohort
#   tx_compare_groups(df_cohort, group_var = "smokingstatus")
#
#   # Compare two custom user-defined groups
#   tx_compare_groups(
#     Cluster_surv,
#     group_var  = NULL,
#     custom_groups = list(
#       GroupA = c("sample1","sample2","sample3"),
#       GroupB = c("sample4","sample5","sample6")
#     )
#   )
# ============================================================

tx_compare_groups <- function(
  df,                        # data frame with survival data (diagsurvtime, status)
  group_var     = NULL,      # column name for grouping; NULL if using custom_groups
  custom_groups = NULL,      # named list of sample ID vectors for custom groups
  ref_level     = NULL,      # reference level for group_var; NULL = first level

  # survival analysis
  horizon_years = 5,
  cox_covars    = c("stage_group", "sex", "age", "smokingstatus"),
  ref_levels    = NULL,      # ref levels for ALL covariates (list)
  numeric_scale = list(age = 5),
  numeric_units = list(age = "years"),
  min_epv       = 5,

  # visualization
  title         = NULL,      # NULL = auto-generated
  group_colours = NULL,      # NULL = auto-generated
  risk_table    = TRUE,
  show_forest   = TRUE,

  # layout
  base_size  = 14,
  title_size = 16,
  widths     = c(1, 1)
) {

  # =========================================================
  # STEP 1 — Validate inputs and build group column
  # =========================================================

  if (is.null(group_var) && is.null(custom_groups)) {
    stop("Must provide either group_var or custom_groups.")
  }

  if (!is.null(custom_groups) && !is.null(group_var)) {
    message("Both group_var and custom_groups provided — using custom_groups.")
    group_var <- NULL
  }

  # ---- custom groups: build group column from sample ID lists ----
  if (!is.null(custom_groups)) {
    if (!is.list(custom_groups) || is.null(names(custom_groups))) {
      stop("custom_groups must be a named list, e.g. list(GroupA = c('s1','s2'), GroupB = c('s3'))")
    }

    group_col_name <- "custom_group"
    group_map <- purrr::imap_dfr(custom_groups, function(ids, grp_name) {
      tibble::tibble(sample = ids, !!group_col_name := grp_name)
    })

    df <- df %>%
      dplyr::inner_join(group_map, by = "sample")

    group_var <- group_col_name
  }

  # ---- resolve group_var case-insensitively ----
  actual_gv <- names(df)[tolower(names(df)) == tolower(group_var)]
  if (length(actual_gv) != 1) {
    stop(sprintf(
      "group_var '%s' not found in df. Available: %s",
      group_var, paste(names(df), collapse = ", ")
    ))
  }
  group_var <- actual_gv

  # ---- check required survival columns ----
  missing_surv <- setdiff(c("diagsurvtime", "status"), tolower(names(df)))
  if (length(missing_surv) > 0) {
    # try case-insensitive resolve
    for (col in c("diagsurvtime", "status")) {
      actual <- names(df)[tolower(names(df)) == col]
      if (length(actual) == 1 && actual != col) {
        names(df)[names(df) == actual] <- col
      }
    }
  }
  stopifnot(
    "diagsurvtime not found in df" = "diagsurvtime" %in% tolower(names(df)),
    "status not found in df"       = "status"       %in% tolower(names(df))
  )

  # =========================================================
  # STEP 2 — Prepare grouping factor
  # =========================================================
  df <- df %>% dplyr::filter(!is.na(.data[[group_var]]))
  df[[group_var]] <- factor(df[[group_var]])

  # apply ref_level if supplied
  if (!is.null(ref_level) && ref_level %in% levels(df[[group_var]])) {
    df[[group_var]] <- stats::relevel(df[[group_var]], ref = ref_level)
  }

  grp_levels <- levels(df[[group_var]])
  n_groups   <- length(grp_levels)

  if (n_groups < 2) {
    stop(sprintf(
      "group_var '%s' has only 1 level after filtering NAs: %s",
      group_var, paste(grp_levels, collapse = ", ")
    ))
  }

  message(sprintf(
    "[tx_compare_groups] group_var='%s' | levels: %s | n=%d patients",
    group_var, paste(grp_levels, collapse = ", "), nrow(df)
  ))

  # =========================================================
  # STEP 3 — Auto-generate colours
  # =========================================================
  if (is.null(group_colours)) {
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
          "[tx_compare_groups] %d levels detected for '%s' — using hcl.colors() palette.",
          length(grp_levels), group_var
        ))
        pal <- grDevices::hcl.colors(length(grp_levels), palette = "Dark 3")
      }
      group_colours <- stats::setNames(pal, grp_levels)
    }
  }

  # =========================================================
  # STEP 4 — Auto-generate ref_levels for Cox if not supplied
  # =========================================================
  if (is.null(ref_levels)) {
    ref_levels <- stats::setNames(list(grp_levels[1]), group_var)
  }

  # ensure group_var is in ref_levels
  if (!group_var %in% names(ref_levels)) {
    ref_levels[[group_var]] <- grp_levels[1]
  }

  # =========================================================
  # STEP 5 — Auto-generate title
  # =========================================================
  if (is.null(title)) {
    title <- sprintf(
      "Survival by %s (n=%d)", group_var, nrow(df)
    )
  }

  # =========================================================
  # STEP 6 — KM plot
  # =========================================================
  p_km <- km_panel_from_df(
    df,
    group_col     = group_var,
    title         = title,
    horizon_years = horizon_years,
    risk_table    = risk_table,
    risk_times    = 0:horizon_years,
    group_colours = group_colours,
    base_size     = base_size,
    title_size    = title_size
  )

  # =========================================================
  # STEP 7 — Cox forest plot
  # =========================================================
  if (show_forest) {
    all_covars <- unique(c(group_var, cox_covars))
    all_covars <- intersect(all_covars, names(df))

    p_forest <- cox_forest_plot_from_df(
      df,
      covars        = all_covars,
      ref_levels    = ref_levels,
      title         = paste("Adjusted Cox —", title),
      min_epv       = min_epv,
      priority      = all_covars,
      numeric_scale = numeric_scale,
      numeric_units = numeric_units,
      base_size     = base_size,
      title_size    = title_size
    )
  } else {
    p_forest <- NULL
  }

  # =========================================================
  # STEP 8 — Compose panels
  # =========================================================
  if (show_forest && !is.null(p_forest)) {
    p_combined <- (p_km | p_forest) +
      patchwork::plot_layout(widths = widths, guides = "keep")
  } else {
    p_combined <- p_km
  }

  # =========================================================
  # STEP 9 — Tidy Cox results as data frame
  # =========================================================
  cox_results <- tryCatch({
    all_covars <- unique(c(group_var, cox_covars))
    all_covars <- intersect(all_covars, names(df))
    d <- stats::na.omit(df[, unique(c("diagsurvtime", "status", all_covars))])

    # apply ref levels
    for (nm in names(ref_levels)) {
      if (nm %in% names(d)) {
        d[[nm]] <- factor(d[[nm]])
        if (ref_levels[[nm]] %in% levels(d[[nm]])) {
          d[[nm]] <- stats::relevel(d[[nm]], ref = ref_levels[[nm]])
        }
      }
    }

    # scale numerics
    for (nm in names(numeric_scale)) {
      if (nm %in% names(d) && is.numeric(d[[nm]])) {
        d[[nm]] <- d[[nm]] / numeric_scale[[nm]]
      }
    }

    form <- as.formula(paste(
      "survival::Surv(diagsurvtime, status) ~",
      paste(all_covars, collapse = " + ")
    ))
    fit <- survival::coxph(form, data = d, ties = "efron")
    broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE)
  }, error = function(e) {
    warning(sprintf("Cox model failed: %s", e$message))
    NULL
  })

  # =========================================================
  # STEP 10 — Group summary table
  # =========================================================
  group_summary <- df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      n        = dplyr::n(),
      n_events = sum(status == 1, na.rm = TRUE),
      median_surv_years = median(diagsurvtime, na.rm = TRUE),
      .groups  = "drop"
    )

  # =========================================================
  # Return
  # =========================================================
  list(
    # plots
    km          = p_km,
    forest      = p_forest,
    combined    = p_combined,

    # results
    cox_results   = cox_results,
    group_summary = group_summary,

    # data
    df          = df,
    group_var   = group_var,
    grp_levels  = grp_levels,
    n           = nrow(df)
  )
}
