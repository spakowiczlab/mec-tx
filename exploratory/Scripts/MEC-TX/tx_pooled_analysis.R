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

#' Pooled Treatment Cohort Analysis
#'
#' Builds a treatment-focused patient cohort using one of four selection
#' modes, then produces a three-panel composite figure: a treatment
#' timeline strip, a Kaplan-Meier survival panel, and an adjusted Cox
#' forest plot. Returns all intermediate data objects for downstream
#' inspection. Implements three-stage \code{n_cohort} transparency
#' reporting.
#'
#' @param Cluster_surv A data frame — the \code{$Cluster_surv} slot from
#'   \code{\link{tx_cluster_surv}}. Must contain \code{sample},
#'   \code{diagsurvtime}, \code{status}, at least one \code{Cluster_kN}
#'   column, and the column named in \code{group_var}.
#' @param timeline A data frame — the direct output of
#'   \code{\link{tx_intervals}}. Must contain \code{sample}, \code{type},
#'   \code{start_year}, and \code{end_year}.
#' @param focus_types Character vector. Treatment type(s) defining the
#'   cohort. Must be canonical MEC-TX type labels. Default
#'   \code{c("Chemo", "IO")}.
#' @param mode One of \code{"any"}, \code{"only"}, \code{"concurrent"},
#'   or \code{"dominant"}. Controls patient inclusion logic — see Details.
#'   Default \code{"any"}.
#' @param group_var Character string. Grouping variable column in
#'   \code{Cluster_surv}. Matched case-insensitively. Default
#'   \code{"CAlevel"}.
#' @param horizon_years Numeric. Analysis horizon in years used for
#'   timeline clipping and segment preparation. Default \code{5}.
#' @param min_share_tx Numeric in \code{[0, 1]}. Minimum combined
#'   \code{focus_types} share for dominant mode. Passed to
#'   \code{\link{tx_focus_dt}}. Default \code{0.33}.
#' @param min_overlap Numeric. Minimum overlap duration in years required
#'   for concurrent mode. Default \code{0} (any overlap).
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
#'   the Cox model. Default \code{list(age = 5)} (HR per 5-year increment)
#'   — intentionally smaller than the package-wide default of 10, since
#'   pooled treatment cohorts are subsets of the full cohort and may have
#'   reduced power.
#' @param numeric_units Named list. Unit labels for numeric covariates in
#'   forest plot. Default \code{list(age = "years")}.
#' @param min_epv Integer. Minimum events-per-variable for Cox covariate
#'   selection. Default \code{5}.
#' @param show_timeline Logical. If \code{TRUE}, include the treatment
#'   timeline strip as the first panel. Default \code{TRUE}.
#' @param group_colours Named character vector or \code{NULL}.
#'   \code{NULL} auto-generates: uses \code{ca_cols} for \code{CAlevel},
#'   cycles through a 10-colour palette, or falls back to
#'   \code{grDevices::hcl.colors()} for large groupings.
#'   Default \code{NULL}.
#' @param horizon_plot Numeric or \code{NULL}. X-axis extent for the
#'   timeline and KM plots. \code{NULL} uses \code{horizon_years}.
#'   Default \code{NULL}.
#' @param base_size Base font size (pt). Default \code{14}.
#' @param title_size Font size (pt) for panel titles. Default \code{16}.
#' @param widths Numeric vector of length 3. Relative widths of the
#'   timeline, KM, and forest panels. Default \code{c(1.4, 1, 1)}.
#'
#' @return A named list with fifteen elements:
#'   \describe{
#'     \item{km}{A \code{patchwork} object — KM panel with risk table.}
#'     \item{forest}{A \code{ggplot} object — adjusted Cox forest plot.}
#'     \item{timeline}{A \code{ggplot} object — treatment timeline strip
#'       with group colour strip. \code{NULL} if
#'       \code{show_timeline = FALSE}.}
#'     \item{ids}{Character vector. All patient IDs passing the mode
#'       filter (\code{n_raw}).}
#'     \item{df}{Tibble. Full cohort metadata for all patients in
#'       \code{ids}.}
#'     \item{segs}{Tibble. Clipped treatment segments for all cohort
#'       patients.}
#'     \item{shares}{Tibble. Treatment type duration shares for all
#'       cohort patients.}
#'     \item{df_plot}{Tibble. Subset of segments used in the timeline
#'       plot — focus-type segments for focus-dominant patients only.}
#'     \item{mode}{Character. The \code{mode} argument used.}
#'     \item{focus_types}{Character vector. The validated focus types
#'       used.}
#'     \item{group_var}{Character. The resolved group column name.}
#'     \item{n_cohort}{Integer. Focus-dominant patient count used for
#'       KM, forest, and timeline panels (Stage 3 of audit trail).}
#'     \item{n_raw}{Integer. Patients passing the mode filter before
#'       dominance restriction (Stage 1 of audit trail).}
#'     \item{n_plot}{Integer. Patients with segments in the timeline
#'       plot (Stage 2 of audit trail).}
#'     \item{group_table}{Table. \code{group_var} distribution in the
#'       focus-dominant KM cohort.}
#'   }
#'
#' @details
#' \strong{Mode definitions:}
#' \describe{
#'   \item{\code{"any"}}{Patients who ever received ALL \code{focus_types}
#'     at any point in their treatment history.}
#'   \item{\code{"only"}}{Patients who received ALL \code{focus_types}
#'     and no other treatment type above threshold (Ancillary always
#'     allowed).}
#'   \item{\code{"concurrent"}}{Patients where all pairs of
#'     \code{focus_types} overlapped in time by more than
#'     \code{min_overlap} years. Requires at least 2 focus types.}
#'   \item{\code{"dominant"}}{Patients whose dominant treatment type is
#'     one of \code{focus_types}, aggregated across all cluster
#'     assignments via \code{\link{tx_focus_dt}}.}
#' }
#'
#' \strong{Three-stage n_cohort audit:} A message block reports patient
#' counts at three stages — \code{n_raw} (pass mode filter),
#' \code{n_plot} (have timeline segments), and \code{n_cohort} (focus-
#' dominant, used in KM and forest). The difference between \code{n_raw}
#' and \code{n_cohort} reflects patients whose dominant treatment type was
#' not in \code{focus_types} — this is expected and correct behaviour.
#'
#' \strong{Concurrent mode:} Uses pairwise interval overlap checking
#' across all combinations of \code{focus_types}. A patient must show
#' overlap for every pair to qualify.
#'
#' \strong{Dominant mode:} Calls \code{\link{tx_focus_dt}} for every
#' \code{Cluster_kN} column in \code{Cluster_surv} and pools the twin
#' IDs. This is the most conservative and biologically specific mode.
#'
#' @examples
#' \dontrun{
#' intervals <- tx_intervals(tx_normalize(med_data, metadata))
#'
#' # Concurrent chemoradiation — the LUSC standard-of-care signal
#' res <- tx_pooled_analysis(
#'   Cluster_surv = res_cluster$Cluster_surv,
#'   timeline     = intervals,
#'   focus_types  = c("Chemo", "Radiation"),
#'   mode         = "concurrent",
#'   group_var    = "CAlevel"
#' )
#'
#' # Save composite figure
#' pdf(file.path(out_dir, "chemorad_concurrent_lusc.pdf"),
#'     width = 18, height = 8)
#' print(res$km | res$forest)
#' dev.off()
#'
#' # Access Cox results
#' res$forest
#' res$group_table
#' res$n_cohort
#'
#' # IO-only cohort, LUAD
#' res2 <- tx_pooled_analysis(
#'   Cluster_surv = res_luad$Cluster_surv,
#'   timeline     = intervals_luad,
#'   focus_types  = "IO",
#'   mode         = "only",
#'   group_var    = "CAlevel"
#' )
#' }
#'
#' @seealso \code{\link{tx_cluster_surv}}, \code{\link{tx_intervals}},
#'   \code{\link{tx_compare_groups}}, \code{\link{tx_focus_dt}},
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
  
  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
  
  # --- 1. Cluster_surv must be a data frame ---
  if (!is.data.frame(Cluster_surv)) {
    stop(
      "tx_pooled_analysis(): 'Cluster_surv' must be a data frame.\n",
      "  -> You passed an object of class: ", paste(class(Cluster_surv), collapse = ", "), ".\n",
      "  -> Pass the '$Cluster_surv' element from tx_cluster_surv() output:\n",
      "    res <- tx_cluster_surv(...)\n",
      "    tx_pooled_analysis(res$Cluster_surv, ...)"
    )
  }
  
  # --- 2. timeline must be a data frame ---
  if (!is.data.frame(timeline)) {
    stop(
      "tx_pooled_analysis(): 'timeline' must be a data frame.\n",
      "  -> You passed an object of class: ", paste(class(timeline), collapse = ", "), ".\n",
      "  -> Pass the output of tx_intervals():\n",
      "    tx_pooled_analysis(Cluster_surv, tx_intervals(norm), ...)"
    )
  }
  
  # --- 3. timeline must have required columns ---
  required_tl <- c("sample", "type", "start_year", "end_year")
  missing_tl  <- setdiff(required_tl, names(timeline))
  if (length(missing_tl) > 0) {
    stop(
      "tx_pooled_analysis(): Required column(s) missing from 'timeline':\n",
      paste0("  x ", missing_tl, collapse = "\n"), "\n",
      "  -> These columns are produced by tx_intervals().\n",
      "  -> Make sure you pass tx_intervals() output, not tx_normalize() output."
    )
  }
  
  # --- 4. Cluster_surv must have at least one Cluster_k column ---
  cluster_cols <- grep("^Cluster_k\\d+$", names(Cluster_surv), value = TRUE)
  if (length(cluster_cols) == 0) {
    stop(
      "tx_pooled_analysis(): No cluster assignment columns found in 'Cluster_surv'.\n",
      "  -> Expected columns named 'Cluster_k3', 'Cluster_k4', ...\n",
      "  -> Columns present: ", paste(names(Cluster_surv), collapse = ", "), "\n",
      "  -> Make sure you pass '$Cluster_surv' from tx_cluster_surv() output."
    )
  }
  
  # --- 5. Cluster_surv must have survival columns ---
  cs_lower <- tolower(names(Cluster_surv))
  for (col in c("sample", "diagsurvtime", "status")) {
    if (!col %in% cs_lower) {
      stop(
        "tx_pooled_analysis(): Required column '", col, "' not found in 'Cluster_surv'.\n",
        "  -> Columns present: ", paste(names(Cluster_surv), collapse = ", "), "\n",
        "  -> 'Cluster_surv' must contain: sample, diagsurvtime, status."
      )
    }
  }
  
  # --- 6. focus_types must have at least one valid type ---
  local_valid_types <- c(
    "Ancillary", "Chemo", "Hormone", "IO",
    "Small_Molecule", "Targeted", "Radiation", "Others"
  )
  valid_focus <- intersect(focus_types, local_valid_types)
  if (length(valid_focus) == 0) {
    stop(
      "tx_pooled_analysis(): None of the supplied 'focus_types' are valid.\n",
      "  -> You supplied: ", paste(focus_types, collapse = ", "), "\n",
      "  -> Valid types: ", paste(local_valid_types, collapse = ", ")
    )
  }
  if (length(valid_focus) < length(focus_types)) {
    invalid <- setdiff(focus_types, local_valid_types)
    warning(
      "tx_pooled_analysis(): Ignoring unrecognised focus_types: ",
      paste(invalid, collapse = ", "), "\n",
      "  -> Valid types: ", paste(local_valid_types, collapse = ", ")
    )
  }
  
  # --- 7. concurrent mode requires >= 2 focus_types ---
  mode_arg <- match.arg(mode)
  if (mode_arg == "concurrent" && length(valid_focus) < 2) {
    stop(
      "tx_pooled_analysis(): mode = 'concurrent' requires at least 2 focus_types.\n",
      "  -> You supplied: ", paste(valid_focus, collapse = ", "), "\n",
      "  -> Example: focus_types = c('Chemo', 'IO')"
    )
  }
  
  # --- 8. group_var must be a single string ---
  if (!is.character(group_var) || length(group_var) != 1) {
    stop(
      "tx_pooled_analysis(): 'group_var' must be a single character string.\n",
      "  -> You passed: ", deparse(group_var), "\n",
      "  -> Example: group_var = 'CAlevel'"
    )
  }
  
  # --- 9. group_var must exist in Cluster_surv (case-insensitive) ---
  group_var_actual <- names(Cluster_surv)[
    tolower(names(Cluster_surv)) == tolower(group_var)
  ]
  if (length(group_var_actual) != 1) {
    stop(
      "tx_pooled_analysis(): group_var '", group_var, "' not found in 'Cluster_surv'.\n",
      "  -> Columns available: ", paste(names(Cluster_surv), collapse = ", "), "\n",
      "  -> group_var is matched case-insensitively, so 'calevel' matches 'CAlevel'."
    )
  }
  
  # --- 10. Sample overlap between Cluster_surv and timeline ---
  cs_samples <- unique(as.character(Cluster_surv$sample))
  tl_samples <- unique(as.character(timeline$sample))
  overlap    <- intersect(cs_samples, tl_samples)
  if (length(overlap) == 0) {
    stop(
      "tx_pooled_analysis(): No sample IDs match between 'Cluster_surv' and 'timeline'.\n",
      "  -> Cluster_surv samples (first 5): ", paste(head(cs_samples, 5), collapse = ", "), "\n",
      "  -> timeline samples (first 5):     ", paste(head(tl_samples, 5), collapse = ", "), "\n",
      "  -> Make sure both come from the same dataset and use the same sample IDs."
    )
  }
  
  # --- 11. horizon_years must be positive ---
  if (!is.numeric(horizon_years) || length(horizon_years) != 1 || horizon_years <= 0) {
    stop(
      "tx_pooled_analysis(): 'horizon_years' must be a single positive number.\n",
      "  -> You passed: ", deparse(horizon_years), "\n",
      "  -> Example: horizon_years = 5 (default)"
    )
  }
  
  # ===========================================================================
  # MAIN FUNCTION BODY
  # ===========================================================================
  
  # ---- arg matching ----
  mode <- match.arg(mode)
  if (is.null(horizon_plot)) horizon_plot <- horizon_years
  
  focus_types <- valid_focus
  group_var   <- group_var_actual
  focus_label <- paste(focus_types, collapse = "+")
  mode_label  <- toupper(mode)   # "CONCURRENT", "ONLY", "ANY", "DOMINANT"
  
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
                           
                           overlap_ids <- segs_a %>%
                             dplyr::inner_join(segs_b, by = "sample", relationship = "many-to-many") %>%
                             dplyr::mutate(overlap = pmin(t1_a, t1_b) - pmax(t0_a, t0_b)) %>%
                             dplyr::filter(overlap > min_overlap) %>%
                             dplyr::distinct(sample) %>%
                             dplyr::pull(sample)
                           
                           concurrent_ids <- intersect(concurrent_ids, overlap_ids)
                         }
                         concurrent_ids
                       },
                       
                       # ---- dominant: cluster-aware via tx_focus_dt across all k ----
                       "dominant" = {
                         k_vals <- as.integer(sub("Cluster_k", "", grep("^Cluster_k\\d+$", names(Cluster_surv), value = TRUE)))
                         ids <- lapply(paste0("Cluster_k", k_vals), function(kc) {
                           if (!kc %in% names(Cluster_surv)) return(NULL)
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
      "tx_pooled_analysis(): No patients found for focus_types = c(%s) with mode = '%s'.\n  -> Try mode = 'any' to broaden the cohort, or check that these treatment types exist in your timeline.",
      paste0('"', focus_types, '"', collapse = ", "), mode
    ))
  }
  
  # =========================================================
  # FIX 5 — Stage 1: n_raw (patients passing mode filter)
  # =========================================================
  n_raw <- length(cohort_ids)
  
  message(sprintf(
    "[tx_pooled_analysis] Mode='%s' | focus=%s | n_raw=%d patients pass mode filter",
    mode, focus_label, n_raw
  ))
  
  # =========================================================
  # STEP 2 — Build cohort data objects
  # =========================================================
  df_cohort    <- Cluster_surv %>% dplyr::filter(sample %in% cohort_ids)
  segs_cohort  <- segs_prep    %>% dplyr::filter(sample %in% cohort_ids)
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
  
  df_plot <- segs_cohort %>%
    dplyr::filter(as.character(type) %in% focus_types) %>%
    dplyr::left_join(order_tbl, by = "sample") %>%
    dplyr::filter(!is.na(y)) %>%
    dplyr::left_join(
      df_cohort %>% dplyr::select(sample, dplyr::all_of(group_var)),
      by = "sample"
    ) %>%
    dplyr::mutate(type = factor(type, levels = focus_types))
  
  group_strip <- order_tbl %>%
    dplyr::left_join(
      df_cohort %>% dplyr::select(sample, dplyr::all_of(group_var)),
      by = "sample"
    )
  
  # Focus-dominant subset: patients whose dominant tx is one of the focus types.
  # This is used consistently for timeline, KM, and forest — biologically correct
  # because patients dominated by other types (e.g. IO in a Chemo+Radiation query)
  # will appear in their own appropriate analysis.
  dominant_ids        <- unique(df_plot$sample)
  n_timeline_patients <- length(dominant_ids)
  df_km               <- df_cohort %>% dplyr::filter(sample %in% dominant_ids)
  
  # =========================================================
  # FIX 5 — Stage 2 & 3: n_timeline and n_km with audit trail
  # =========================================================
  n_timeline <- dplyr::n_distinct(df_plot$sample)
  n_km       <- dplyr::n_distinct(df_km$sample)
  
  message(sprintf(
    paste0(
      "  [n_cohort audit] %s [%s] | %s\n",
      "    Stage 1 — n_raw (pass mode filter)     : %d\n",
      "    Stage 2 — n_timeline (in plot segments) : %d\n",
      "    Stage 3 — n_km (focus-dominant for KM)  : %d  <- n_cohort reported"
    ),
    focus_label, mode_label, group_var,
    n_raw, n_timeline, n_km
  ))
  
  if (n_km < n_raw) {
    dropped <- n_raw - n_km
    message(sprintf(
      paste0(
        "    Note: %d patients dropped from KM — their dominant treatment\n",
        "    type was not in focus_types ('%s').\n",
        "    This is correct behaviour: they appear in their own regimen analysis."
      ),
      dropped, focus_label
    ))
  }
  
  # =========================================================
  # STEP 5 — Colours for group_var
  # =========================================================
  grp_levels <- sort(unique(as.character(df_cohort[[group_var]])))
  grp_levels <- grp_levels[!is.na(grp_levels)]
  
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
          "[tx_pooled_analysis] %d levels detected for '%s' — using hcl.colors() palette.",
          length(grp_levels), group_var
        ))
        pal <- grDevices::hcl.colors(length(grp_levels), palette = "Dark 3")
      }
      group_colours <- stats::setNames(pal, grp_levels)
    }
  }
  
  if (is.null(ref_levels)) {
    ref_levels <- stats::setNames(list(grp_levels[1]), group_var)
  }
  
  # =========================================================
  # STEP 6 — KM plot (focus-dominant subset)
  # =========================================================
  p_km <- km_panel_from_df(
    df_km,
    group_col     = group_var,
    title         = sprintf(
      "All %s [%s] — all clusters pooled (n=%d)",
      focus_label, mode_label, n_timeline_patients
    ),
    horizon_years = horizon_plot,
    risk_table    = TRUE,
    group_colours = group_colours
  )
  
  # =========================================================
  # STEP 7 — Cox forest plot (focus-dominant subset)
  # =========================================================
  all_covars <- unique(c(group_var, cox_covars))
  all_covars <- intersect(all_covars, names(df_cohort))
  
  p_forest <- cox_forest_plot_from_df(
    df_km,
    covars        = all_covars,
    ref_levels    = ref_levels,
    title         = sprintf(
      "Adjusted Cox — %s [%s] pooled (all clusters)",
      focus_label, mode_label
    ),
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
      ggplot2::geom_tile(
        data = group_strip %>%
          dplyr::mutate(grp = factor(.data[[group_var]], levels = grp_levels)),
        ggplot2::aes(x = strip_x, y = y, fill = grp),
        width = 0.4, height = 0.9, alpha = 0.8
      ) +
      ggplot2::scale_fill_manual(
        values   = group_colours,
        name     = group_var,
        na.value = "grey80"
      ) +
      ggnewscale::new_scale_colour() +
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
          "Treatment Timelines — %s [%s] Patients (n=%d)",
          focus_label, mode_label, n_timeline_patients
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
    km          = p_km,
    forest      = p_forest,
    timeline    = p_timeline,
    ids         = cohort_ids,          # all patients passing the mode filter
    df          = df_cohort,           # full cohort data (all ids)
    segs        = segs_cohort,
    shares      = share_cohort,
    df_plot     = df_plot,
    mode        = mode,
    focus_types = focus_types,
    group_var   = group_var,
    n_cohort    = n_km,                # Fix 5: focus-dominant count (KM/forest/timeline)
    n_raw       = n_raw,               # Fix 5: total passing mode filter
    n_plot      = n_timeline_patients,
    group_table = table(df_km[[group_var]], useNA = "ifany")
  )
}