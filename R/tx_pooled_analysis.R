# ============================================================
# MEC-TX analysis: tx_pooled_analysis()
# analysis/tx_pooled_analysis.R
#
# Wraps the full cohort-building and visualization workflow
# for any treatment combination and any grouping variable.
#
# Modes:
#   "any"        --- patients who ever received ALL focus_types
#   "only"       --- patients who received focus_types and nothing else
#   "concurrent" --- patients where focus_types overlapped in time
#   "dominant"   --- cluster-aware; focus_types dominate treatment time
# ============================================================

#' Pooled Treatment Cohort Analysis
#'
#' Builds a treatment-focused patient cohort using one of four selection
#' modes, then produces a three-panel composite figure: a treatment
#' timeline strip, a Kaplan-Meier survival panel, and an adjusted Cox
#' forest plot. Returns all intermediate data objects for downstream
#' inspection. Implements three-stage \code{n_cohort} transparency
#' reporting.
#'
#' @param Cluster_surv A data frame --- the \code{$Cluster_surv} slot from
#'   \code{\link{tx_cluster_surv}}. Must contain \code{sample},
#'   \code{diagsurvtime}, \code{status}, at least one \code{Cluster_kN}
#'   column, and the column named in \code{group_var}.
#' @param timeline A data frame --- the direct output of
#'   \code{\link{tx_intervals}}. Must contain \code{sample}, \code{type},
#'   \code{start_year}, and \code{end_year}.
#' @param focus_types Character vector. Treatment type(s) defining the
#'   cohort. Must be canonical MEC-TX type labels. Default
#'   \code{c("Chemo", "IO")}.
#' @param mode One of \code{"any"}, \code{"only"}, \code{"concurrent"},
#'   or \code{"dominant"}. Controls patient inclusion logic --- see Details.
#'   Default \code{"any"}.
#' @param group_var Character string. Grouping variable column in
#'   \code{Cluster_surv}. Matched case-insensitively. Default
#'   \code{"CAlevel"}.
#' @param horizon_years Numeric. Analysis horizon in years used for
#'   timeline clipping and segment preparation. Default \code{5}.
#' @param min_share_tx Numeric in \code{[0, 1]}. Minimum combined
#'   \code{focus_types} share for dominant mode. Passed to
#'   \code{\link{tx_focus_dt}}. Default \code{0.33}.
#' @param concurrent_window Numeric. Maximum gap in years between two
#'   treatment intervals to still be classified as concurrent. Default
#'   \code{4/52} (4 weeks). A value of \code{0} means strict overlap
#'   only. Implements the PI-specified definition: a new treatment
#'   starting within 4 weeks of a prior treatment ending is classified
#'   as concurrent. Applied symmetrically --- direction of the gap does
#'   not matter.
#' @param n_twins Integer. Maximum twins per cluster for dominant mode.
#'   \code{999} retains all. Default \code{999}.
#' @param enforce_sequence Logical. Passed to \code{\link{tx_focus_dt}}
#'   for dominant mode. Default \code{FALSE}.
#' @param sequence_strict Logical. Passed to \code{\link{tx_focus_dt}}
#'   for dominant mode. Default \code{FALSE}.
#' @param start_filter One of \code{"all"}, \code{"single_only"},
#'   \code{"combo_only"}. Passed to \code{\link{tx_focus_dt}} for
#'   dominant mode. Default \code{"all"}.
#' @param pure_focus_only Logical. Passed to \code{\link{tx_focus_dt}}
#'   for dominant mode. Default \code{FALSE}.
#' @param cox_covars Character vector. Adjustment covariates for the Cox
#'   model. Columns absent from \code{Cluster_surv} are silently dropped.
#'   Default \code{c("stage_group", "sex", "age", "smokingstatus")}.
#' @param ref_levels Named list or \code{NULL}. Reference levels for Cox
#'   covariates. \code{NULL} auto-sets \code{group_var} to its first
#'   factor level. Default \code{NULL}.
#' @param numeric_scale Named list. Divisors for numeric covariates in
#'   the Cox model. Default \code{list(age = 5)}.
#' @param numeric_units Named list. Unit labels for numeric covariates in
#'   forest plot. Default \code{list(age = "years")}.
#' @param min_epv Integer. Minimum events-per-variable for Cox covariate
#'   selection. Default \code{5}.
#' @param show_timeline Logical. If \code{TRUE}, include the treatment
#'   timeline strip as the first panel. Default \code{TRUE}.
#' @param group_colours Named character vector or \code{NULL}.
#'   Default \code{NULL}.
#' @param horizon_plot Numeric or \code{NULL}. Default \code{NULL}.
#' @param base_size Base font size (pt). Default \code{14}.
#' @param title_size Font size (pt) for panel titles. Default \code{16}.
#' @param widths Numeric vector of length 3. Default \code{c(1.4, 1, 1)}.
#'
#' @return A named list with fifteen elements: \code{km}, \code{forest},
#'   \code{timeline}, \code{ids}, \code{df}, \code{segs}, \code{shares},
#'   \code{df_plot}, \code{mode}, \code{focus_types}, \code{group_var},
#'   \code{n_cohort}, \code{n_raw}, \code{n_plot}, \code{group_table}.
#'
#' @details
#' \strong{Mode definitions:}
#' \describe{
#'   \item{\code{"any"}}{Patients who ever received ALL \code{focus_types}.
#'     No dominance filter --- all n_raw patients used in KM and Cox.}
#'   \item{\code{"only"}}{Patients who received ALL \code{focus_types}
#'     and no other treatment type (Ancillary allowed).
#'     No dominance filter --- all n_raw patients used in KM and Cox.}
#'   \item{\code{"concurrent"}}{Patients where all pairs of
#'     \code{focus_types} were administered within \code{concurrent_window}
#'     years of each other (default 4 weeks).
#'     No dominance filter --- all n_raw patients used in KM and Cox.}
#'   \item{\code{"dominant"}}{Patients whose dominant treatment type is
#'     one of \code{focus_types}. Dominance filter applied --- only
#'     focus-type dominant patients used in KM and Cox.}
#' }
#'
#' \strong{Dominance filter:} Applied ONLY for \code{mode = "dominant"}.
#' For \code{"any"}, \code{"only"}, and \code{"concurrent"}, all patients
#' passing the mode filter are included in KM and Cox, so
#' \code{n_cohort == n_raw} for these three modes.
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
#' meta <- data.frame(
#'   sample       = paste0('P', seq_len(n)),
#'   diagsurvtime = rep(3, n),
#'   Status       = seq_len(n) %% 2L,
#'   CAlevel      = rep(c('High','Low'), n/2),
#'   stringsAsFactors = FALSE
#' )
#' norm         <- tx_normalize(med_data)
#' intervals    <- tx_intervals(norm)
#' cluster_res  <- tx_cluster_surv(meta, norm, k_range = 2,
#'                                 umap_neighbors = 5,
#'                                 min_feature_variance = 0)
#' res <- tx_pooled_analysis(
#'   Cluster_surv = cluster_res$Cluster_surv,
#'   timeline     = intervals,
#'   focus_types  = c('Chemo', 'Radiation'),
#'   group_var    = 'CAlevel'
#' )
#' res$n_cohort
#'
#' @seealso \code{\link{tx_cluster_surv}}, \code{\link{tx_intervals}},
#'   \code{\link{km_panel_from_df}}, \code{\link{cox_forest_plot_from_df}}
#'
#' @import ggplot2
#' @importFrom dplyr filter select mutate group_by summarise left_join
#'   distinct pull inner_join n_distinct all_of arrange row_number
#' @importFrom patchwork plot_layout
#' @importFrom ggnewscale new_scale_colour
#' @importFrom stats setNames
#' @importFrom grDevices hcl.colors
#' @export

tx_pooled_analysis <- function(
    Cluster_surv,
    timeline,
    focus_types       = c("Chemo", "IO"),
    mode              = c("any", "only", "concurrent", "dominant"),
    group_var         = "CAlevel",
    horizon_years     = 5,
    min_share_tx      = 0.33,
    concurrent_window = 4/52,
    n_twins           = 999,
    enforce_sequence  = FALSE,
    sequence_strict   = FALSE,
    start_filter      = "all",
    pure_focus_only   = FALSE,
    cox_covars        = c("stage_group", "sex", "age", "smokingstatus"),
    ref_levels        = NULL,
    numeric_scale     = list(age = 5),
    numeric_units     = list(age = "years"),
    min_epv           = 5,
    show_timeline     = TRUE,
    group_colours     = NULL,
    horizon_plot      = NULL,
    base_size         = 14,
    title_size        = 16,
    widths            = c(1.4, 1, 1)
) {
  
  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
  
  if (!is.data.frame(Cluster_surv))
    stop("tx_pooled_analysis(): 'Cluster_surv' must be a data frame.")
  
  if (!is.data.frame(timeline))
    stop("tx_pooled_analysis(): 'timeline' must be a data frame.")
  
  required_tl <- c("sample", "type", "start_year", "end_year")
  missing_tl  <- setdiff(required_tl, names(timeline))
  if (length(missing_tl) > 0)
    stop("tx_pooled_analysis(): Required column(s) missing from 'timeline': ",
         paste(missing_tl, collapse = ", "))
  
  if (length(grep("^Cluster_k\\d+$", names(Cluster_surv))) == 0)
    stop("tx_pooled_analysis(): No Cluster_kN columns found in 'Cluster_surv'.")
  
  cs_lower <- tolower(names(Cluster_surv))
  for (col in c("sample", "diagsurvtime", "status"))
    if (!col %in% cs_lower)
      stop("tx_pooled_analysis(): Required column '", col,
           "' not found in 'Cluster_surv'.")
  
  local_valid_types <- c("Ancillary", "Chemo", "Hormone", "IO",
                         "Small_Molecule", "Targeted", "Radiation", "Others")
  valid_focus <- intersect(focus_types, local_valid_types)
  if (length(valid_focus) == 0)
    stop("tx_pooled_analysis(): No valid focus_types supplied.")
  if (length(valid_focus) < length(focus_types))
    warning("tx_pooled_analysis(): Ignoring unrecognised focus_types: ",
            paste(setdiff(focus_types, local_valid_types), collapse = ", "))
  
  mode_arg <- match.arg(mode)
  if (mode_arg == "concurrent" && length(valid_focus) < 2)
    stop("tx_pooled_analysis(): mode = 'concurrent' requires >= 2 focus_types.")
  
  if (!is.character(group_var) || length(group_var) != 1)
    stop("tx_pooled_analysis(): 'group_var' must be a single character string.")
  
  group_var_actual <- names(Cluster_surv)[
    tolower(names(Cluster_surv)) == tolower(group_var)]
  if (length(group_var_actual) != 1)
    stop("tx_pooled_analysis(): group_var '", group_var,
         "' not found in 'Cluster_surv'.")
  
  if (length(intersect(as.character(Cluster_surv$sample),
                       as.character(timeline$sample))) == 0)
    stop("tx_pooled_analysis(): No sample IDs match between 'Cluster_surv' and 'timeline'.")
  
  if (!is.numeric(horizon_years) || horizon_years <= 0)
    stop("tx_pooled_analysis(): 'horizon_years' must be a positive number.")
  
  if (!is.numeric(concurrent_window) || concurrent_window < 0)
    stop("tx_pooled_analysis(): 'concurrent_window' must be non-negative.")
  
  # ===========================================================================
  # SETUP
  # ===========================================================================
  
  mode        <- match.arg(mode)
  focus_types <- valid_focus
  group_var   <- group_var_actual
  focus_label <- paste(focus_types, collapse = "+")
  mode_label  <- toupper(mode)
  if (is.null(horizon_plot)) horizon_plot <- horizon_years
  
  # ===========================================================================
  # STEP 1 --- Find cohort IDs
  # ===========================================================================
  
  segs_prep <- prep_segs(timeline, horizon_years = horizon_years)
  
  cohort_ids <- switch(mode,
                       
                       "any" = {
                         segs_prep %>%
                           dplyr::group_by(sample) %>%
                           dplyr::summarise(
                             has_all = all(focus_types %in% unique(as.character(type))),
                             .groups = "drop") %>%
                           dplyr::filter(has_all) %>%
                           dplyr::pull(sample)
                       },
                       
                       "only" = {
                         other_types <- setdiff(local_valid_types, c(focus_types, "Ancillary"))
                         segs_prep %>%
                           dplyr::group_by(sample) %>%
                           dplyr::summarise(
                             has_all_focus = all(focus_types %in% unique(as.character(type))),
                             has_other     = any(as.character(type) %in% other_types),
                             .groups = "drop") %>%
                           dplyr::filter(has_all_focus & !has_other) %>%
                           dplyr::pull(sample)
                       },
                       
                       # gap = max(t0_a, t0_b) - min(t1_a, t1_b)
                       # negative = overlap, positive = sequential gap
                       # qualifies if gap < concurrent_window for every pair
                       "concurrent" = {
                         concurrent_ids <- segs_prep %>%
                           dplyr::distinct(sample) %>%
                           dplyr::pull(sample)
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
                           concurrent_ids <- intersect(concurrent_ids,
                                                       segs_a %>%
                                                         dplyr::inner_join(segs_b, by = "sample",
                                                                           relationship = "many-to-many") %>%
                                                         dplyr::mutate(gap = pmax(t0_a, t0_b) - pmin(t1_a, t1_b)) %>%
                                                         dplyr::filter(gap < concurrent_window) %>%
                                                         dplyr::distinct(sample) %>%
                                                         dplyr::pull(sample))
                         }
                         concurrent_ids
                       },
                       
                       "dominant" = {
                         k_vals <- as.integer(sub("Cluster_k", "",
                                                  grep("^Cluster_k\\d+$", names(Cluster_surv), value = TRUE)))
                         unique(unlist(lapply(paste0("Cluster_k", k_vals), function(kc) {
                           tryCatch({
                             demo <- tx_focus_dt(
                               Cluster_surv     = Cluster_surv,
                               segs             = timeline,
                               kc               = kc,
                               cl               = NULL,
                               focus_types      = focus_types,
                               group_col        = group_var,
                               horizon_years    = horizon_years,
                               n_twins          = n_twins,
                               min_share_tx     = min_share_tx,
                               enforce_sequence = enforce_sequence,
                               sequence_strict  = sequence_strict,
                               start_filter     = start_filter,
                               pure_focus_only  = pure_focus_only,
                               add_forest       = FALSE)
                             attr(demo, "twin_ids")
                           }, error = function(e) NULL)
                         })))
                       }
  )
  
  if (length(cohort_ids) == 0)
    stop(sprintf(
      "tx_pooled_analysis(): No patients found for focus_types = c(%s) with mode = '%s'.",
      paste0('"', focus_types, '"', collapse = ", "), mode))
  
  n_raw <- length(cohort_ids)
  message(sprintf(
    "[tx_pooled_analysis] Mode='%s' | focus=%s | window=%.4f yrs (%.1f wks) | n_raw=%d",
    mode, focus_label, concurrent_window, concurrent_window * 52, n_raw))
  
  # ===========================================================================
  # STEP 2 --- Build cohort objects
  # ===========================================================================
  
  df_cohort    <- Cluster_surv %>% dplyr::filter(sample %in% cohort_ids)
  segs_cohort  <- segs_prep    %>% dplyr::filter(sample %in% cohort_ids)
  share_cohort <- treatment_shares(segs_cohort)
  
  message(sprintf("  n in Cluster_surv: %d", nrow(df_cohort)))
  message(sprintf("  %s distribution:", group_var))
  message(table(df_cohort[[group_var]], useNA = "ifany"))
  if ("stage_group" %in% names(df_cohort)) {
    message("  Stage distribution:")
    message(table(df_cohort$stage_group, useNA = "ifany"))
  }
  
  # ===========================================================================
  # STEP 3 --- Timeline plot data
  #
  # DOMINANCE FILTER applies ONLY for mode = "dominant".
  # For "any", "only", "concurrent" --- all n_raw patients go into
  # KM and Cox. n_cohort == n_raw for these three modes.
  # ===========================================================================
  
  if (mode == "dominant") {
    # Only include patients whose dominant tx is in focus_types
    order_tbl <- share_cohort %>%
      dplyr::select(sample, dom_type_tx) %>%
      dplyr::filter(dom_type_tx %in% focus_types) %>%
      dplyr::left_join(
        segs_cohort %>%
          dplyr::group_by(sample) %>%
          dplyr::summarise(first_start = min(t0), .groups = "drop"),
        by = "sample") %>%
      dplyr::arrange(dom_type_tx, first_start) %>%
      dplyr::mutate(y = dplyr::row_number())
  } else {
    # Include ALL cohort patients --- order by dom_type_tx for visual clarity
    # but do NOT drop patients whose dominant type is outside focus_types
    order_tbl <- share_cohort %>%
      dplyr::select(sample, dom_type_tx) %>%
      dplyr::left_join(
        segs_cohort %>%
          dplyr::group_by(sample) %>%
          dplyr::summarise(first_start = min(t0), .groups = "drop"),
        by = "sample") %>%
      dplyr::arrange(dom_type_tx, first_start) %>%
      dplyr::mutate(y = dplyr::row_number())
  }
  
  # For dominant mode: show only focus_types segments (clean, focused view)
  # For any/only/concurrent: show ALL treatment types so the full treatment
  # context is visible --- patients whose dominant type is outside focus_types
  # are shown with their complete treatment history
  df_plot <- segs_cohort %>%
    {if (mode == "dominant")
      dplyr::filter(., as.character(type) %in% focus_types)
      else .} %>%
    dplyr::left_join(order_tbl, by = "sample") %>%
    dplyr::filter(!is.na(y)) %>%
    dplyr::left_join(
      df_cohort %>% dplyr::select(sample, dplyr::all_of(group_var)),
      by = "sample") %>%
    dplyr::mutate(
      type = if (mode == "dominant") {
        factor(as.character(type), levels = focus_types)
      } else {
        factor(as.character(type), levels = c(
          "Chemo", "IO", "Radiation", "Targeted",
          "Small_Molecule", "Hormone", "Ancillary", "Others"))
      }
    )
  
  group_strip <- order_tbl %>%
    dplyr::left_join(
      df_cohort %>% dplyr::select(sample, dplyr::all_of(group_var)),
      by = "sample")
  
  # KM cohort: dominant mode uses focus-dominant subset only
  # all other modes use full n_raw cohort
  km_ids <- if (mode == "dominant") unique(df_plot$sample) else cohort_ids
  
  n_timeline_patients <- length(km_ids)
  df_km               <- df_cohort %>% dplyr::filter(sample %in% km_ids)
  n_km                <- dplyr::n_distinct(df_km$sample)
  
  # Audit message
  if (mode == "dominant") {
    message(sprintf(
      paste0(
        "  [n_cohort audit] %s [%s] | %s\n",
        "    Stage 1 --- n_raw (pass mode filter)    : %d\n",
        "    Stage 2 --- n_km (focus-dominant for KM): %d  <- n_cohort"),
      focus_label, mode_label, group_var, n_raw, n_km))
    if (n_km < n_raw)
      message(sprintf(
        "    Note: %d patients dropped --- dominant tx not in focus_types ('%s')",
        n_raw - n_km, focus_label))
  } else {
    message(sprintf(
      paste0(
        "  [n_cohort audit] %s [%s] | %s\n",
        "    n_raw = n_cohort = %d (no dominance filter for mode = '%s')"),
      focus_label, mode_label, group_var, n_km, mode))
  }
  
  # ===========================================================================
  # STEP 4 --- Colours
  # ===========================================================================
  
  grp_levels <- sort(unique(as.character(df_cohort[[group_var]])))
  grp_levels <- grp_levels[!is.na(grp_levels)]
  
  if (is.null(group_colours)) {
    if (tolower(group_var) == "calevel" && all(c("High", "Low") %in% grp_levels)) {
      group_colours <- c(High = "#F28E2B", Low = "#56B4E9")
    } else {
      base_pal <- c("#4E79A7","#F28E2B","#E15759","#76B7B2",
                    "#59A14F","#EDC948","#B07AA1","#FF9DA7",
                    "#9C755F","#BAB0AC")
      pal <- if (length(grp_levels) <= length(base_pal)) {
        base_pal[seq_along(grp_levels)]
      } else {
        grDevices::hcl.colors(length(grp_levels), palette = "Dark 3")
      }
      group_colours <- stats::setNames(pal, grp_levels)
    }
  }
  
  if (is.null(ref_levels))
    ref_levels <- stats::setNames(list(grp_levels[1]), group_var)
  
  # ===========================================================================
  # STEP 5 --- KM plot
  # ===========================================================================
  
  p_km <- km_panel_from_df(
    df_km,
    group_col     = group_var,
    title         = sprintf("All %s [%s] --- all clusters pooled (n=%d)",
                            focus_label, mode_label, n_km),
    horizon_years = horizon_plot,
    risk_table    = TRUE,
    group_colours = group_colours)
  
  # ===========================================================================
  # STEP 6 --- Cox forest plot
  # ===========================================================================
  
  all_covars <- intersect(unique(c(group_var, cox_covars)), names(df_cohort))
  
  p_forest <- cox_forest_plot_from_df(
    df_km,
    covars        = all_covars,
    ref_levels    = ref_levels,
    title         = sprintf("Adjusted Cox --- %s [%s] pooled (all clusters)",
                            focus_label, mode_label),
    min_epv       = min_epv,
    priority      = all_covars,
    numeric_scale = numeric_scale,
    numeric_units = numeric_units,
    base_size     = base_size,
    title_size    = title_size)
  
  # ===========================================================================
  # STEP 7 --- Timeline strip
  # ===========================================================================
  
  if (show_timeline) {
    
    local_tx_cols <- c(
      Ancillary      = "#E1BE6A", Chemo     = "#FDB863",
      Hormone        = "#DC267F", IO        = "#2CA02C",
      Small_Molecule = "#76B7B2", Targeted  = "#4E79A7",
      Radiation      = "#6A51A3", Others    = "#8C8C8C")
    
    strip_x   <- horizon_plot + 0.4
    label_tbl <- order_tbl %>%
      dplyr::group_by(dom_type_tx) %>%
      dplyr::summarise(y_mid = mean(y), .groups = "drop")
    
    p_timeline <- ggplot2::ggplot() +
      ggplot2::geom_tile(
        data = group_strip %>%
          dplyr::mutate(grp = factor(.data[[group_var]], levels = grp_levels)),
        ggplot2::aes(x = strip_x, y = y, fill = grp),
        width = 0.4, height = 0.9, alpha = 0.8) +
      ggplot2::scale_fill_manual(
        values = group_colours, name = group_var, na.value = "grey80") +
      ggnewscale::new_scale_colour() +
      ggplot2::geom_segment(
        data = df_plot,
        ggplot2::aes(x = t0, xend = t1, y = y, yend = y, colour = type),
        linewidth = 0.5) +
      ggplot2::scale_colour_manual(
        values = if (mode == "dominant") local_tx_cols[focus_types] else local_tx_cols,
        name   = "Treatment Type", drop = FALSE) +
      ggplot2::geom_text(
        data = label_tbl,
        ggplot2::aes(x = -0.3, y = y_mid, label = dom_type_tx),
        hjust = 1, size = 3.5, fontface = "bold") +
      ggplot2::scale_x_continuous(
        limits = c(-0.5, strip_x + 0.3),
        breaks = 0:horizon_plot,
        name   = "Years since first treatment") +
      ggplot2::scale_y_continuous(
        name = sprintf("Patients (n=%d)", n_km), breaks = NULL) +
      ggplot2::labs(
        title    = sprintf("Treatment Timelines --- %s [%s] Patients (n=%d)",
                           focus_label, mode_label, n_km),
        subtitle = sprintf("Ordered by dominant treatment type | Right strip = %s",
                           group_var)) +
      ggplot2::theme_minimal(base_size = base_size) +
      ggplot2::theme(
        plot.title         = ggplot2::element_text(size = title_size, face = "bold"),
        plot.subtitle      = ggplot2::element_text(size = base_size - 3),
        axis.text.y        = ggplot2::element_blank(),
        axis.ticks.y       = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor   = ggplot2::element_blank(),
        legend.position    = "right")
  } else {
    p_timeline <- NULL
  }
  
  # ===========================================================================
  # STEP 8 --- Return
  # ===========================================================================
  
  list(
    km          = p_km,
    forest      = p_forest,
    timeline    = p_timeline,
    ids         = cohort_ids,
    df          = df_cohort,
    segs        = segs_cohort,
    shares      = share_cohort,
    df_plot     = df_plot,
    mode        = mode,
    focus_types = focus_types,
    group_var   = group_var,
    n_cohort    = n_km,
    n_raw       = n_raw,
    n_plot      = n_timeline_patients,
    group_table = table(df_km[[group_var]], useNA = "ifany")
  )
}
