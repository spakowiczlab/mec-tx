# ============================================================
# MEC-TX analysis: tx_compare_groups()
# analysis/tx_compare_groups.R
#
# Flexible KM + adjusted Cox forest for any grouping variable
# or custom sample ID lists. Works with any categorical column
# in Cluster_surv — not limited to CAlevel.
# ============================================================

#' Flexible Survival Comparison by Grouping Variable or Custom Groups
#'
#' Produces a paired Kaplan-Meier panel and adjusted Cox forest plot for
#' any categorical grouping variable or user-defined sample ID lists.
#' Designed to work with any column in the \code{Cluster_surv} output —
#' not limited to \code{CAlevel}. Returns plots, tidy Cox results, and a
#' group summary table.
#'
#' @param df A data frame with survival data. Must contain
#'   \code{diagsurvtime} (numeric, years) and \code{status} (integer 0/1)
#'   columns (case-insensitive). Typically the \code{$Cluster_surv} slot
#'   from \code{\link{tx_cluster_surv}} or the \code{$df} slot from
#'   \code{\link{tx_pooled_analysis}}.
#' @param group_var Character string or \code{NULL}. Name of the grouping
#'   column in \code{df}. Matched case-insensitively. Ignored if
#'   \code{custom_groups} is also supplied (with a message). \code{NULL}
#'   requires \code{custom_groups} to be provided. Default \code{NULL}.
#' @param custom_groups Named list or \code{NULL}. Each element is a
#'   character vector of sample IDs defining a group. Must have at least
#'   two named elements. Sample IDs not found in \code{df} are dropped
#'   with a warning. If both \code{group_var} and \code{custom_groups} are
#'   supplied, \code{custom_groups} takes precedence. Default \code{NULL}.
#' @param ref_level Character string or \code{NULL}. Reference level for
#'   \code{group_var} in the Cox model and KM plot. \code{NULL} uses the
#'   first factor level. Default \code{NULL}.
#' @param horizon_years Numeric. X-axis clip for KM plot and risk table,
#'   in years. Default \code{5}.
#' @param cox_covars Character vector. Adjustment covariates for the Cox
#'   model, in addition to \code{group_var}. Columns not present in
#'   \code{df} are silently dropped. Default
#'   \code{c("stage_group", "sex", "age", "smokingstatus")}.
#' @param ref_levels Named list. Reference levels for all Cox covariates.
#'   \code{NULL} auto-sets \code{group_var} reference to its first factor
#'   level. Default \code{NULL}.
#' @param numeric_scale Named list. Divisors for numeric covariates in
#'   the Cox model. Default \code{list(age = 5)} (HR per 5-year increment)
#'   — intentionally smaller than the package-wide default of 10, since
#'   group comparison cohorts may be subsets with reduced power.
#' @param numeric_units Named list. Unit labels for numeric covariates in
#'   the forest plot. Default \code{list(age = "years")}.
#' @param min_epv Integer. Minimum events-per-variable for Cox covariate
#'   selection. Passed to \code{\link{cox_forest_plot_from_df}}.
#'   Default \code{5}.
#' @param title Character string or \code{NULL}. Plot title. \code{NULL}
#'   auto-generates \code{"Survival by <group_var> (n=<N>)"}. Default
#'   \code{NULL}.
#' @param group_colours Named character vector or \code{NULL}. Hex colours
#'   for each group level. \code{NULL} auto-generates: uses
#'   \code{ca_cols} for \code{CAlevel}, cycles through a 10-colour
#'   palette for up to 10 levels, falls back to
#'   \code{grDevices::hcl.colors()} for larger groupings. Default
#'   \code{NULL}.
#' @param risk_table Logical. If \code{TRUE}, an at-risk table is appended
#'   below the KM panel. Passed to \code{\link{km_panel_from_df}}.
#'   Default \code{TRUE}.
#' @param show_forest Logical. If \code{TRUE}, the Cox forest plot is
#'   computed and included in \code{combined}. Default \code{TRUE}.
#' @param base_size Base font size (pt). Default \code{14}.
#' @param title_size Font size (pt) for plot titles. Default \code{16}.
#' @param widths Numeric vector of length 2. Relative widths of the KM
#'   and forest panels in the combined layout, passed to
#'   \code{patchwork::plot_layout()}. Default \code{c(1, 1)}.
#'
#' @return A named list with nine elements:
#'   \describe{
#'     \item{km}{A \code{patchwork} object — KM panel with optional risk
#'       table.}
#'     \item{forest}{A \code{ggplot} object — Cox forest plot. \code{NULL}
#'       if \code{show_forest = FALSE}.}
#'     \item{combined}{A \code{patchwork} object — KM and forest side by
#'       side. Equal to \code{km} alone if \code{show_forest = FALSE}.}
#'     \item{cox_results}{Tibble from \code{broom::tidy()} with columns
#'       \code{term}, \code{estimate} (HR), \code{conf.low}, \code{conf.high},
#'       \code{p.value}. \code{NULL} if Cox model fails.}
#'     \item{group_summary}{Tibble with one row per group level: \code{n},
#'       \code{n_events}, \code{median_surv_years}.}
#'     \item{df}{The analysis data frame after group column resolution,
#'       NA filtering, and factor relevel — useful for downstream
#'       inspection.}
#'     \item{group_var}{Character. The resolved (actual) group column name
#'       used.}
#'     \item{grp_levels}{Character vector of group factor levels in plot
#'       order.}
#'     \item{n}{Integer. Number of patients in the analysis.}
#'   }
#'   Returns a list of \code{NULL}s (with a warning) if \code{group_var}
#'   has fewer than two non-NA levels after filtering.
#'
#' @details
#' \strong{Group resolution:} When \code{group_var} is supplied, it is
#' matched case-insensitively to column names in \code{df}. When
#' \code{custom_groups} is supplied, a temporary \code{"custom_group"}
#' column is built via an inner join on \code{sample}, restricting
#' \code{df} to only the listed patients.
#'
#' \strong{Cox model:} \code{group_var} is always included as the first
#' covariate regardless of EPV. Additional covariates from \code{cox_covars}
#' are selected by EPV budget. The Cox model is refitted cleanly for the
#' \code{cox_results} slot — independently of the forest plot.
#'
#' \strong{Colour fallback:} For \code{group_var = "CAlevel"} with High/Low
#' levels, \code{ca_cols} from \code{constants.R} is used automatically.
#' For more than 10 group levels, \code{grDevices::hcl.colors()} is used
#' with a message.
#'
#' @examples
#' \dontrun{
#' # Compare by CAlevel in full LUSC cohort
#' res <- tx_compare_groups(
#'   df        = Cluster_surv,
#'   group_var = "CAlevel",
#'   title     = "Overall Survival by CAlevel — LUSC"
#' )
#' pdf(file.path(out_dir, "calevel_compare_lusc.pdf"), width = 14, height = 8)
#' print(res$combined)
#' dev.off()
#'
#' # Compare by smoking status in chemo-only cohort
#' res2 <- tx_compare_groups(
#'   df        = chemo_cohort,
#'   group_var = "smokingstatus",
#'   ref_level = "Never"
#' )
#'
#' # Compare two custom user-defined groups
#' res3 <- tx_compare_groups(
#'   df            = Cluster_surv,
#'   custom_groups = list(
#'     GroupA = c("sample1", "sample2", "sample3"),
#'     GroupB = c("sample4", "sample5", "sample6")
#'   )
#' )
#'
#' # Access tidy Cox results
#' res$cox_results
#' res$group_summary
#' }
#'
#' @seealso \code{\link{tx_pooled_analysis}}, \code{\link{km_panel_from_df}},
#'   \code{\link{cox_forest_plot_from_df}}, \code{\link{tx_cluster_surv}}
#'
#' @import ggplot2
#' @importFrom dplyr filter group_by summarise inner_join n
#' @importFrom patchwork plot_layout
#' @importFrom purrr imap_dfr
#' @importFrom tibble tibble
#' @importFrom survival coxph Surv
#' @importFrom broom tidy
#' @importFrom stats relevel setNames na.omit as.formula
#' @importFrom grDevices hcl.colors
#' @export

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
  
  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
  
  # --- 1. df must be a data frame ---
  if (!is.data.frame(df)) {
    stop(
      "tx_compare_groups(): 'df' must be a data frame.\n",
      "  → You passed an object of class: ", paste(class(df), collapse = ", "), ".\n",
      "  → Pass the '$Cluster_surv' element from tx_cluster_surv() output,\n",
      "    or a cohort data frame from tx_pooled_analysis()$df."
    )
  }
  
  # --- 2. Must provide group_var OR custom_groups, not neither ---
  if (is.null(group_var) && is.null(custom_groups)) {
    stop(
      "tx_compare_groups(): Must provide either 'group_var' or 'custom_groups'.\n",
      "  → group_var example:    group_var = 'CAlevel'\n",
      "  → custom_groups example: custom_groups = list(GroupA = c('s1','s2'), GroupB = c('s3','s4'))"
    )
  }
  
  # --- 3. custom_groups must be a named list if supplied ---
  if (!is.null(custom_groups)) {
    if (!is.list(custom_groups) || is.null(names(custom_groups))) {
      stop(
        "tx_compare_groups(): 'custom_groups' must be a named list.\n",
        "  → You passed: ", class(custom_groups), "\n",
        "  → Example: custom_groups = list(GroupA = c('s1','s2'), GroupB = c('s3','s4'))"
      )
    }
    if (length(custom_groups) < 2) {
      stop(
        "tx_compare_groups(): 'custom_groups' must have at least 2 groups.\n",
        "  → You supplied ", length(custom_groups), " group(s).\n",
        "  → Example: custom_groups = list(GroupA = c('s1','s2'), GroupB = c('s3','s4'))"
      )
    }
    # check all sample IDs in custom_groups exist in df
    if ("sample" %in% names(df)) {
      all_custom_ids <- unlist(custom_groups)
      missing_ids    <- setdiff(all_custom_ids, df$sample)
      if (length(missing_ids) > 0) {
        warning(
          "tx_compare_groups(): ", length(missing_ids),
          " sample ID(s) in 'custom_groups' not found in 'df' — these will be dropped.\n",
          "  → Missing (first 5): ", paste(head(missing_ids, 5), collapse = ", ")
        )
      }
    }
  }
  
  # --- 4. group_var must be a single string if supplied ---
  if (!is.null(group_var) && (!is.character(group_var) || length(group_var) != 1)) {
    stop(
      "tx_compare_groups(): 'group_var' must be a single character string.\n",
      "  → You passed: ", deparse(group_var), "\n",
      "  → Example: group_var = 'CAlevel'"
    )
  }
  
  # --- 5. df must have required survival columns ---
  df_lower <- tolower(names(df))
  for (col in c("diagsurvtime", "status")) {
    if (!col %in% df_lower) {
      stop(
        "tx_compare_groups(): Required survival column '", col, "' not found in 'df'.\n",
        "  → Columns present: ", paste(names(df), collapse = ", "), "\n",
        "  → Required columns (case-insensitive): diagsurvtime, status\n",
        "  → Make sure 'df' comes from tx_cluster_surv()$Cluster_surv or tx_pooled_analysis()$df."
      )
    }
  }
  
  # --- 6. status column must be 0/1 ---
  actual_status <- names(df)[df_lower == "status"][1]
  status_vals   <- unique(na.omit(df[[actual_status]]))
  if (!all(status_vals %in% c(0, 1))) {
    stop(
      "tx_compare_groups(): '", actual_status, "' column must contain only 0 and 1.\n",
      "  → Values found: ", paste(status_vals, collapse = ", "), "\n",
      "  → 0 = censored, 1 = event (e.g. death)."
    )
  }
  
  # --- 7. horizon_years must be positive ---
  if (!is.numeric(horizon_years) || length(horizon_years) != 1 || horizon_years <= 0) {
    stop(
      "tx_compare_groups(): 'horizon_years' must be a single positive number.\n",
      "  → You passed: ", deparse(horizon_years), "\n",
      "  → Example: horizon_years = 5 (default)"
    )
  }
  
  # --- 8. widths must be length 2 with positive values ---
  if (!is.numeric(widths) || length(widths) != 2 || any(widths <= 0)) {
    stop(
      "tx_compare_groups(): 'widths' must be a numeric vector of length 2 with positive values.\n",
      "  → You passed: ", deparse(widths), "\n",
      "  → Example: widths = c(1, 1)  (KM | forest)"
    )
  }
  
  # ===========================================================================
  # STEP 1 — Build group column
  # ===========================================================================
  
  if (!is.null(custom_groups) && !is.null(group_var)) {
    message("tx_compare_groups(): Both group_var and custom_groups provided — using custom_groups.")
    group_var <- NULL
  }
  
  # ---- custom groups: build group column from sample ID lists ----
  if (!is.null(custom_groups)) {
    group_col_name <- "custom_group"
    group_map <- purrr::imap_dfr(custom_groups, function(ids, grp_name) {
      tibble::tibble(sample = ids, !!group_col_name := grp_name)
    })
    df        <- df %>% dplyr::inner_join(group_map, by = "sample")
    group_var <- group_col_name
  }
  
  # ---- resolve group_var case-insensitively ----
  actual_gv <- names(df)[tolower(names(df)) == tolower(group_var)]
  if (length(actual_gv) != 1) {
    stop(
      "tx_compare_groups(): group_var '", group_var, "' not found in 'df'.\n",
      "  → Columns available: ", paste(names(df), collapse = ", "), "\n",
      "  → group_var is matched case-insensitively, so 'calevel' matches 'CAlevel'."
    )
  }
  group_var <- actual_gv
  
  # ---- case-insensitive resolve for survival columns ----
  for (col in c("diagsurvtime", "status")) {
    actual <- names(df)[tolower(names(df)) == col]
    if (length(actual) == 1 && actual != col) {
      names(df)[names(df) == actual] <- col
    }
  }
  
  # ===========================================================================
  # STEP 2 — Prepare grouping factor
  # ===========================================================================
  df <- df %>% dplyr::filter(!is.na(.data[[group_var]]))
  df[[group_var]] <- factor(df[[group_var]])
  
  if (!is.null(ref_level) && ref_level %in% levels(df[[group_var]])) {
    df[[group_var]] <- stats::relevel(df[[group_var]], ref = ref_level)
  }
  
  grp_levels <- levels(df[[group_var]])
  n_groups   <- length(grp_levels)
  
  if (n_groups < 2) {
    warning(sprintf(
      "tx_compare_groups(): group_var '%s' has only 1 level after filtering NAs: %s — skipping.",
      group_var, paste(grp_levels, collapse = ", ")
    ))
    return(list(km = NULL, forest = NULL, combined = NULL,
                cox_results = NULL, group_summary = NULL))
  }
  
  message(sprintf(
    "[tx_compare_groups] group_var='%s' | levels: %s | n=%d patients",
    group_var, paste(grp_levels, collapse = ", "), nrow(df)
  ))
  
  # ===========================================================================
  # STEP 3 — Auto-generate colours
  # ===========================================================================
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
        message(sprintf(
          "[tx_compare_groups] %d levels detected for '%s' — using hcl.colors() palette.",
          length(grp_levels), group_var
        ))
        pal <- grDevices::hcl.colors(length(grp_levels), palette = "Dark 3")
      }
      group_colours <- stats::setNames(pal, grp_levels)
    }
  }
  
  # ===========================================================================
  # STEP 4 — Auto-generate ref_levels for Cox if not supplied
  # ===========================================================================
  if (is.null(ref_levels)) {
    ref_levels <- stats::setNames(list(grp_levels[1]), group_var)
  }
  if (!group_var %in% names(ref_levels)) {
    ref_levels[[group_var]] <- grp_levels[1]
  }
  
  # ===========================================================================
  # STEP 5 — Auto-generate title
  # ===========================================================================
  if (is.null(title)) {
    title <- sprintf("Survival by %s (n=%d)", group_var, nrow(df))
  }
  
  # ===========================================================================
  # STEP 6 — KM plot
  # ===========================================================================
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
  
  # ===========================================================================
  # STEP 7 — Cox forest plot
  # ===========================================================================
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
  
  # ===========================================================================
  # STEP 8 — Compose panels
  # ===========================================================================
  if (show_forest && !is.null(p_forest)) {
    p_combined <- (p_km | p_forest) +
      patchwork::plot_layout(widths = widths, guides = "keep")
  } else {
    p_combined <- p_km
  }
  
  # ===========================================================================
  # STEP 9 — Tidy Cox results as data frame
  # ===========================================================================
  cox_results <- tryCatch({
    all_covars <- unique(c(group_var, cox_covars))
    all_covars <- intersect(all_covars, names(df))
    d <- stats::na.omit(df[, unique(c("diagsurvtime", "status", all_covars))])
    
    for (nm in names(ref_levels)) {
      if (nm %in% names(d)) {
        d[[nm]] <- factor(d[[nm]])
        if (ref_levels[[nm]] %in% levels(d[[nm]])) {
          d[[nm]] <- stats::relevel(d[[nm]], ref = ref_levels[[nm]])
        }
      }
    }
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
  
  # ===========================================================================
  # STEP 10 — Group summary table
  # ===========================================================================
  group_summary <- df %>%
    dplyr::group_by(.data[[group_var]]) %>%
    dplyr::summarise(
      n                 = dplyr::n(),
      n_events          = sum(status == 1, na.rm = TRUE),
      median_surv_years = median(diagsurvtime, na.rm = TRUE),
      .groups           = "drop"
    )
  
  # ===========================================================================
  # Return
  # ===========================================================================
  list(
    km            = p_km,
    forest        = p_forest,
    combined      = p_combined,
    cox_results   = cox_results,
    group_summary = group_summary,
    df            = df,
    group_var     = group_var,
    grp_levels    = grp_levels,
    n             = nrow(df)
  )
}