#' Cluster Treatment Timelines and Attach Survival Data
#'
#' Takes normalised treatment segments and patient metadata, builds a
#' binary treatment feature matrix, performs PCA dimensionality reduction,
#' clusters patients using k-means (k = 3 to 20 by default), and attaches
#' survival and covariate data to the cluster assignments. UMAP is computed
#' for visualisation only --- clustering runs in PCA space.
#'
#' @param metadata A data frame containing patient-level metadata.
#'   Must contain columns named in \code{sample_col}, \code{surv_time_col},
#'   and \code{status_col} (case-insensitive). Typically \code{Cluster_surv}
#'   (LUSC) or \code{LUAD_metadata} (LUAD).
#' @param timeline_long_norm A long-format normalised treatment timeline ---
#'   the direct output of \code{\link{tx_normalize}}. Must contain columns
#'   named in \code{sample_col}, \code{time_col}, and \code{tx_col}
#'   (case-insensitive).
#' @param sample_col Character string. Name of the patient identifier column
#'   in both \code{metadata} and \code{timeline_long_norm}.
#'   Case-insensitive. Default \code{"sample"}.
#' @param time_col Character string. Name of the treatment time column in
#'   \code{timeline_long_norm}. Case-insensitive. Default
#'   \code{"TimeSinceTreatmentStart"}.
#' @param tx_col Character string. Name of the treatment group column in
#'   \code{timeline_long_norm}. Case-insensitive. Default
#'   \code{"treatment_group"}.
#' @param surv_time_col Character string. Name of the survival time column
#'   (years from diagnosis) in \code{metadata}. Case-insensitive.
#'   Default \code{"diagsurvtime"}.
#' @param status_col Character string. Name of the survival status column
#'   in \code{metadata}. Auto-detected and standardised to integer
#'   0/1 via \code{\link{standardise_status}}. Case-insensitive. Default
#'   \code{"status"}.
#' @param meta_keep Character vector or \code{NULL}. Additional metadata
#'   columns to carry through to the output \code{Cluster_surv} data frame.
#'   \code{NULL} retains all columns. Default \code{NULL}.
#' @param horizon_years Numeric. Maximum follow-up horizon in years. Time
#'   bins beyond this value are excluded. Must match the \code{horizon_years}
#'   used in \code{\link{tx_normalize}}. Default \code{5}.
#' @param grid_weeks Numeric. Time bin width in weeks. Must match the
#'   \code{grid_weeks} used in \code{\link{tx_normalize}}. Common values:
#'   \code{1} (weekly), \code{2} (biweekly), \code{4} (monthly, default).
#'   Default \code{4}.
#' @param include_none Logical. If \code{TRUE}, time bins with no recorded
#'   treatment are filled with \code{"None"} before encoding. Default
#'   \code{TRUE}.
#' @param drop_none_cols Logical. If \code{TRUE}, binary feature columns
#'   corresponding to \code{"None"} treatments are removed before PCA,
#'   reducing noise from untreated time bins. Default \code{TRUE}.
#' @param min_feature_variance Numeric. Near-zero-variance threshold.
#'   Feature columns with variance at or below this value are dropped before
#'   PCA. Default \code{0.01}.
#' @param seed Integer. Random seed for reproducibility of PCA, UMAP, and
#'   k-means. Default \code{42}.
#' @param n_pcs Integer. Maximum number of principal components to retain
#'   for clustering and UMAP. Capped at the actual number of PCs available.
#'   Default \code{50}.
#' @param umap_neighbors Integer. Number of neighbours for UMAP. Must be
#'   strictly less than the number of patients --- use
#'   \code{min(30, n_patients - 1)} as a rule of thumb. Default \code{30}.
#' @param umap_min_dist Numeric. Minimum distance parameter for UMAP
#'   layout. Smaller values produce tighter clusters visually. Default
#'   \code{0.3}.
#' @param k_range Integer vector. Range of k values for k-means clustering.
#'   All values must be \code{>= 2} and \code{< n_patients}. Default
#'   \code{3:20}.
#' @param kmeans_nstart Integer. Number of random starts for k-means.
#'   Higher values improve stability at the cost of runtime. Default
#'   \code{50}.
#'
#' @return A named list with eight elements:
#'   \describe{
#'     \item{Cluster_surv}{Tibble. One row per patient with cluster
#'       assignments for every k in \code{k_range} (\code{Cluster_k3},
#'       \code{Cluster_k4}, ...), plus \code{diagsurvtime}, \code{status}
#'       (integer 0/1), \code{status_label} (factor "Alive"/"Dead"), and
#'       any columns in \code{meta_keep}. This is the primary output for
#'       downstream survival analyses.}
#'     \item{cluster_results}{Tibble. Cluster assignments only --- one row
#'       per patient, all k columns, no survival data.}
#'     \item{umap_df}{Data frame with columns \code{sample}, \code{UMAP1},
#'       \code{UMAP2}. Visualisation only --- not used in clustering.}
#'     \item{pca_matrix}{Numeric matrix. PCA scores for the first
#'       \code{n_pcs} components. Rownames are sample IDs. This is the
#'       space in which clustering is performed.}
#'     \item{pca_var_explained}{Named numeric vector. Proportion of variance
#'       explained by each retained PC.}
#'     \item{X}{Numeric matrix. NZV-filtered binary feature matrix used
#'       as PCA input. Rownames are sample IDs, column names are
#'       \code{<TimeBin>_<TreatmentType>}.}
#'     \item{treatment_encoded}{Tibble. Wide-format multi-hot encoded
#'       treatment matrix before NZV filtering.}
#'     \item{treatment_matrix_ordered}{Tibble. Wide-format treatment matrix
#'       with time bins as columns, ordered chronologically.}
#'   }
#'
#' @details
#' \strong{Pipeline:}
#' \enumerate{
#'   \item Resolve column names case-insensitively via
#'     \code{resolve_col}.
#'   \item Standardise status column via \code{\link{standardise_status}}.
#'   \item Build time-bin matrix: bin \code{TimeSinceTreatmentStart} into
#'     \code{grid_weeks}-wide bins up to \code{horizon_years}.
#'   \item Multi-hot encode treatment types per bin (handles
#'     \code{"+"}-separated combination regimens).
#'   \item Remove near-zero-variance features (\code{var <= min_feature_variance}).
#'   \item PCA on the NZV-filtered binary matrix (\code{scale. = FALSE}).
#'   \item UMAP on PCA scores --- visualisation only.
#'   \item K-means on PCA scores for each k in \code{k_range}.
#'   \item Merge cluster assignments with survival metadata.
#' }
#'
#' \strong{Clustering vs UMAP:} K-means runs in PCA space, not UMAP space.
#' UMAP is computed solely for visualisation. Do not use UMAP coordinates
#' as clustering input.
#'
#' \strong{Parameter alignment:} \code{grid_weeks} and \code{horizon_years}
#' must match the values used in the upstream \code{\link{tx_normalize}}
#' call. Mismatches produce incorrect time bin alignment.
#'
#' \strong{Sample audit:} A message block is printed reporting sample counts
#' at each pipeline stage --- metadata input, timeline encoding, clustering,
#' and survival merge. Samples lost at each stage are identified by ID.
#'
#' \strong{Duplicate handling:} Duplicate sample IDs in \code{metadata}
#' produce a warning; the first occurrence is retained.
#'
#' @examples
#' \dontrun{
#' # Standard LUSC pipeline
#' norm <- tx_normalize(LUSC_med_data, Cluster_surv)
#' res  <- tx_cluster_surv(Cluster_surv, norm)
#'
#' # Access primary output
#' Cluster_surv_k <- res$Cluster_surv
#'
#' # LUAD --- capitalised Status column, extra metadata columns retained
#' norm_luad <- tx_normalize(LUAD_med_data, LUAD_metadata)
#' res_luad  <- tx_cluster_surv(
#'   metadata           = LUAD_metadata,
#'   timeline_long_norm = norm_luad,
#'   status_col         = "Status",
#'   meta_keep          = c("CAlevel", "stage_group", "smokingstatus")
#' )
#'
#' # Adjust umap_neighbors for small cohorts
#' res_small <- tx_cluster_surv(
#'   metadata           = metadata,
#'   timeline_long_norm = norm,
#'   umap_neighbors     = min(30, nrow(metadata) - 1)
#' )
#' }
#'
#' @seealso \code{\link{tx_normalize}}, \code{\link{standardise_status}},
#'   \code{\link{tx_pooled_analysis}}, \code{\link{km_panel_from_df}},
#'   \code{\link{plot_timeline_for_k}}
#'
#' @importFrom dplyr filter select mutate group_by summarise left_join
#'   distinct any_of all_of transmute
#' @importFrom tidyr pivot_wider pivot_longer separate_rows replace_na
#' @importFrom tibble column_to_rownames
#' @importFrom stringr str_extract str_trim
#' @importFrom stats prcomp kmeans var
#' @importFrom umap umap
#' @export

tx_cluster_surv <- function(
    metadata,                              # --- renamed from NSCLC_metadata (Thread 8)
    timeline_long_norm,
    
    # column names (case-insensitive)
    sample_col    = "sample",
    time_col      = "TimeSinceTreatmentStart",
    tx_col        = "treatment_group",
    surv_time_col = "diagsurvtime",
    status_col    = "status",
    meta_keep     = NULL,          # NULL = keep all metadata columns
    
    # horizon / binning
    horizon_years        = 5,
    grid_weeks           = 4,
    include_none         = TRUE,
    drop_none_cols       = TRUE,
    
    # feature filtering
    min_feature_variance = 0.01,
    
    # PCA
    seed   = 42,
    n_pcs  = 50,
    
    # UMAP (visualization only)
    umap_neighbors = 30,
    umap_min_dist  = 0.3,
    
    # clustering (in PCA space)
    k_range        = 3:20,
    kmeans_nstart  = 50
) {
  
  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
  
  # --- 1. Both inputs must be data frames ---
  if (!is.data.frame(metadata)) {
    stop(
      "tx_cluster_surv(): 'metadata' must be a data frame.\n",
      "  -> You passed an object of class: ", paste(class(metadata), collapse = ", "), ".\n",
      "  -> Load your metadata with read.csv() or similar."
    )
  }
  
  if (!is.data.frame(timeline_long_norm)) {
    stop(
      "tx_cluster_surv(): 'timeline_long_norm' must be a data frame.\n",
      "  -> You passed an object of class: ", paste(class(timeline_long_norm), collapse = ", "), ".\n",
      "  -> Pass the output of tx_normalize() directly:\n",
      "    tx_cluster_surv(metadata, tx_normalize(your_raw_data))"
    )
  }
  
  # --- 2. Required columns in timeline_long_norm ---
  for (col in c(time_col, tx_col)) {
    actual <- names(timeline_long_norm)[tolower(names(timeline_long_norm)) == tolower(col)]
    if (length(actual) == 0) {
      stop(
        "tx_cluster_surv(): Column '", col, "' not found in 'timeline_long_norm'.\n",
        "  -> Columns present: ", paste(names(timeline_long_norm), collapse = ", "), "\n",
        "  -> Expected columns from tx_normalize() output:\n",
        "      time_col = 'TimeSinceTreatmentStart'\n",
        "      tx_col   = 'treatment_group'\n",
        "  -> If your columns differ, pass them explicitly:\n",
        "    tx_cluster_surv(..., time_col = 'your_col', tx_col = 'your_col')"
      )
    }
  }
  
  # --- 3. Required columns in metadata ---
  meta_lower <- tolower(names(metadata))
  for (col in c(sample_col, surv_time_col, status_col)) {
    if (!tolower(col) %in% meta_lower) {
      stop(
        "tx_cluster_surv(): Column '", col, "' not found in 'metadata'.\n",
        "  -> Columns present: ", paste(names(metadata), collapse = ", "), "\n",
        "  -> Required metadata columns (case-insensitive):\n",
        "      sample_col    = 'sample'\n",
        "      surv_time_col = 'diagsurvtime'\n",
        "      status_col    = 'status'\n",
        "  -> If your columns differ, pass them explicitly:\n",
        "    tx_cluster_surv(..., surv_time_col = 'your_col', status_col = 'your_col')"
      )
    }
  }
  
  # --- 4. Standardise status column ---
  actual_status <- names(metadata)[meta_lower == tolower(status_col)][1]
  metadata      <- standardise_status(metadata, status_col = actual_status)
  
  # Post-standardisation safety check
  status_vals <- unique(na.omit(metadata[[actual_status]]))
  if (!all(status_vals %in% c(0L, 1L))) {
    stop(
      "tx_cluster_surv(): standardise_status() did not produce clean 0/1 for '",
      actual_status, "'.\n",
      "  -> Values found after standardisation: ", paste(status_vals, collapse = ", "), "\n",
      "  -> Please check your data and report this as a bug."
    )
  }
  
  # --- 5. Sample overlap between metadata and timeline ---
  actual_sample_tl <- names(timeline_long_norm)[
    tolower(names(timeline_long_norm)) == tolower(sample_col)
  ][1]
  actual_sample_md <- names(metadata)[
    tolower(names(metadata)) == tolower(sample_col)
  ][1]
  
  overlap <- intersect(
    unique(as.character(timeline_long_norm[[actual_sample_tl]])),
    unique(as.character(metadata[[actual_sample_md]]))
  )
  if (length(overlap) == 0) {
    stop(
      "tx_cluster_surv(): No sample IDs match between 'timeline_long_norm' and 'metadata'.\n",
      "  -> Timeline samples (first 5): ",
      paste(head(unique(timeline_long_norm[[actual_sample_tl]]), 5), collapse = ", "), "\n",
      "  -> Metadata samples (first 5): ",
      paste(head(unique(metadata[[actual_sample_md]]), 5), collapse = ", "), "\n",
      "  -> Check that sample IDs use the same format in both data frames."
    )
  }
  if (length(overlap) < 10) {
    warning(
      "tx_cluster_surv(): Only ", length(overlap), " sample(s) overlap between ",
      "timeline and metadata.\n",
      "  -> This may produce unreliable clusters.\n",
      "  -> Check that sample IDs match between your two data frames."
    )
  }
  
  # --- 6. umap_neighbors must be < number of overlapping patients ---
  n_patients <- length(overlap)
  if (umap_neighbors >= n_patients) {
    stop(
      "tx_cluster_surv(): 'umap_neighbors' (", umap_neighbors, ") must be less than ",
      "the number of patients (", n_patients, ").\n",
      "  -> Set umap_neighbors to a value less than ", n_patients, ".\n",
      "  -> Rule of thumb: umap_neighbors = min(30, n_patients - 1)\n",
      "  -> Example: tx_cluster_surv(..., umap_neighbors = ",
      min(30L, n_patients - 1L), ")"
    )
  }
  
  # --- 7. k_range must be valid ---
  if (!is.numeric(k_range) || length(k_range) < 1 || any(k_range < 2)) {
    stop(
      "tx_cluster_surv(): 'k_range' must be a numeric vector with all values >= 2.\n",
      "  -> You passed: ", deparse(k_range), "\n",
      "  -> Example: k_range = 3:20"
    )
  }
  if (max(k_range) >= n_patients) {
    stop(
      "tx_cluster_surv(): max(k_range) (", max(k_range), ") must be less than ",
      "the number of patients (", n_patients, ").\n",
      "  -> Reduce k_range, e.g. k_range = 3:", min(n_patients - 1L, 20L)
    )
  }
  
  # --- 8. grid_weeks must be a positive number ---
  if (!is.numeric(grid_weeks) || length(grid_weeks) != 1 || grid_weeks <= 0) {
    stop(
      "tx_cluster_surv(): 'grid_weeks' must be a single positive number.\n",
      "  -> You passed: ", deparse(grid_weeks), "\n",
      "  -> This must match the grid_weeks used in tx_normalize() (default = 4)."
    )
  }
  
  # --- 9. horizon_years must be positive ---
  if (!is.numeric(horizon_years) || length(horizon_years) != 1 || horizon_years <= 0) {
    stop(
      "tx_cluster_surv(): 'horizon_years' must be a single positive number.\n",
      "  -> You passed: ", deparse(horizon_years), "\n",
      "  -> Example: horizon_years = 5 (default)"
    )
  }
  
  # ===========================================================================
  # MAIN FUNCTION BODY
  # ===========================================================================
  
  # ---------------------------------------------------------------------------
  # 0) Resolve column names case-insensitively
  # ---------------------------------------------------------------------------
  sample_col_tl <- resolve_col(timeline_long_norm, sample_col,    "sample_col")
  time_col_tl   <- resolve_col(timeline_long_norm, time_col,      "time_col")
  tx_col_tl     <- resolve_col(timeline_long_norm, tx_col,        "tx_col")
  
  sample_col_md <- resolve_col(metadata, sample_col,    "sample_col")
  surv_col_md   <- resolve_col(metadata, surv_time_col, "surv_time_col")
  status_col_md <- resolve_col(metadata, status_col,    "status_col")
  
  # meta_keep: NULL = keep all columns; otherwise soft resolve
  if (is.null(meta_keep)) {
    keep_meta <- setdiff(
      names(metadata),
      c(sample_col_md, surv_col_md, status_col_md)
    )
  } else {
    keep_meta <- character(0)
    for (m in meta_keep) {
      found <- names(metadata)[tolower(names(metadata)) == tolower(m)]
      if (length(found) == 1) {
        keep_meta <- c(keep_meta, found)
      } else {
        warning(sprintf(
          "meta_keep column '%s' not found in metadata --- skipping.", m
        ))
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # 1) Build time-bin timeline
  # ---------------------------------------------------------------------------
  stopifnot(
    "grid_weeks must be a single positive integer" =
      length(grid_weeks) == 1 && grid_weeks > 0
  )
  res        <- 52 / grid_weeks
  bin_prefix <- sprintf("Week%d", grid_weeks)
  
  dat <- timeline_long_norm %>%
    transmute(
      sample                  = .data[[ sample_col_tl ]],
      TimeSinceTreatmentStart = as.numeric(.data[[ time_col_tl ]]),
      treatment_group         = as.character(.data[[ tx_col_tl ]])
    ) %>%
    filter(!is.na(sample), !is.na(TimeSinceTreatmentStart), !is.na(treatment_group)) %>%
    mutate(
      bin_idx = round(TimeSinceTreatmentStart * res),
      TimeBin = sprintf("%s_%03d", bin_prefix, bin_idx)
    ) %>%
    filter(bin_idx >= 0, bin_idx <= horizon_years * res)
  
  # ---------------------------------------------------------------------------
  # 2) Wide bin matrix
  # ---------------------------------------------------------------------------
  treatment_matrix <- dat %>%
    group_by(sample, TimeBin) %>%
    summarise(
      Treatment = paste(sort(unique(treatment_group)), collapse = "+"),
      .groups   = "drop"
    ) %>%
    pivot_wider(names_from = TimeBin, values_from = Treatment, values_fill = NA)
  
  bin_cols <- setdiff(names(treatment_matrix), "sample")
  treatment_matrix_ordered <- treatment_matrix[,
                                               c("sample", bin_cols[order(as.numeric(stringr::str_extract(bin_cols, "\\d+")))]),
                                               drop = FALSE
  ]
  
  # ---------------------------------------------------------------------------
  # 3) Restrict to horizon, fill gaps, long format
  # ---------------------------------------------------------------------------
  keep_bins <- sprintf("%s_%03d", bin_prefix, 0:round(horizon_years * res))
  keep_cols <- intersect(keep_bins, names(treatment_matrix_ordered))
  
  treatment_short <- treatment_matrix_ordered %>%
    select(sample, all_of(keep_cols)) %>%
    pivot_longer(cols = -sample, names_to = "TimeBin", values_to = "Treatment")
  
  if (include_none) {
    treatment_short <- treatment_short %>%
      mutate(Treatment = replace_na(Treatment, "None"))
  }
  
  # ---------------------------------------------------------------------------
  # 4) Multi-hot encode
  # ---------------------------------------------------------------------------
  treatment_encoded <- treatment_short %>%
    filter(!is.na(Treatment)) %>%
    mutate(Treatment = as.character(Treatment)) %>%
    tidyr::separate_rows(Treatment, sep = "\\+") %>%
    mutate(
      Treatment = stringr::str_trim(Treatment),
      Feature   = paste(TimeBin, Treatment, sep = "_"),
      Value     = 1L
    ) %>%
    select(sample, Feature, Value) %>%
    distinct() %>%
    pivot_wider(
      names_from  = Feature,
      values_from = Value,
      values_fill = list(Value = 0L)
    )
  
  if (drop_none_cols) {
    none_cols <- grep("_None$", names(treatment_encoded), value = TRUE)
    if (length(none_cols)) {
      treatment_encoded <- treatment_encoded %>% select(-all_of(none_cols))
    }
  }
  
  # ---------------------------------------------------------------------------
  # 5) Feature matrix --- remove near-zero-variance columns before PCA
  # ---------------------------------------------------------------------------
  X      <- treatment_encoded %>% tibble::column_to_rownames("sample") %>% as.matrix()
  sl_ids <- rownames(X)
  if (length(sl_ids) < 3) stop("Too few samples after encoding to cluster.")
  
  col_var       <- apply(X, 2, var)
  n_dropped_nzv <- sum(col_var <= min_feature_variance)
  if (n_dropped_nzv > 0) {
    message(sprintf(
      "[tx_cluster_surv] Dropping %d near-zero-variance features (var <= %.3f) before PCA.",
      n_dropped_nzv, min_feature_variance
    ))
    X <- X[, col_var > min_feature_variance, drop = FALSE]
  }
  if (ncol(X) < 2) stop("Too few features remain after NZV filtering.")
  
  # ---------------------------------------------------------------------------
  # 6) PCA
  # ---------------------------------------------------------------------------
  pca_result <- prcomp(X, scale. = FALSE)
  n_take     <- min(n_pcs, ncol(pca_result$x))
  pca_matrix <- pca_result$x[, 1:n_take, drop = FALSE]
  rownames(pca_matrix) <- sl_ids
  
  pca_var_explained <- summary(pca_result)$importance[2, 1:n_take]
  
  # ---------------------------------------------------------------------------
  # 7) UMAP --- visualization ONLY
  # ---------------------------------------------------------------------------
  if (!requireNamespace("umap", quietly = TRUE)) {
    stop("Package 'umap' not installed. Install with install.packages('umap').")
  }
  set.seed(seed)
  umap_result <- umap::umap(
    pca_matrix,
    n_neighbors = umap_neighbors,
    min_dist    = umap_min_dist
  )
  umap_df <- as.data.frame(umap_result$layout)
  colnames(umap_df) <- c("UMAP1", "UMAP2")
  umap_df$sample <- sl_ids
  
  # ---------------------------------------------------------------------------
  # 8) k-means clustering in PCA space
  # ---------------------------------------------------------------------------
  set.seed(seed)
  cluster_results <- data.frame(sample = sl_ids, stringsAsFactors = FALSE)
  for (k in k_range) {
    km <- kmeans(pca_matrix, centers = k, nstart = kmeans_nstart)
    cluster_results[[paste0("Cluster_k", k)]] <- km$cluster
  }
  
  # ---------------------------------------------------------------------------
  # 9) Merge with survival metadata
  #    status_label is included automatically (created by standardise_status)
  # ---------------------------------------------------------------------------
  Surv_data <- metadata %>%
    filter(.data[[ sample_col_md ]] %in% cluster_results$sample) %>%
    select(
      sample       = all_of(sample_col_md),
      diagsurvtime = all_of(surv_col_md),
      status       = all_of(status_col_md),
      any_of("status_label"),
      all_of(keep_meta)
    )
  
  dup_samples <- Surv_data$sample[duplicated(Surv_data$sample)]
  if (length(dup_samples) > 0) {
    warning(sprintf(
      "[tx_cluster_surv] %d duplicate sample(s) in metadata: %s. Keeping first.",
      length(dup_samples), paste(unique(dup_samples), collapse = ", ")
    ))
    Surv_data <- Surv_data %>% filter(!duplicated(sample))
  }
  
  Cluster_surv <- cluster_results %>%
    left_join(Surv_data, by = "sample")
  
  n_missing_surv <- sum(is.na(Cluster_surv$diagsurvtime) | is.na(Cluster_surv$status))
  if (n_missing_surv > 0) {
    message(sprintf(
      "[tx_cluster_surv] %d sample(s) dropped due to missing survival data.",
      n_missing_surv
    ))
  }
  
  Cluster_surv <- Cluster_surv %>%
    filter(!is.na(diagsurvtime), !is.na(status))
  
  # ---------------------------------------------------------------------------
  # Sample audit
  # ---------------------------------------------------------------------------
  n_metadata   <- length(unique(metadata[[ sample_col_md ]]))
  n_timeline   <- length(unique(dat$sample))
  n_encoded    <- nrow(X)
  n_clustered  <- nrow(cluster_results)
  n_surv_merge <- nrow(Cluster_surv)
  
  message(sprintf("
--- Sample audit ---
  In metadata:             %d
  In timeline (encoded):   %d   (lost: %d --- not in timeline)
  After PCA/clustering:    %d
  After survival merge:    %d   (lost: %d --- missing survival data)
  Final:                   %d
",
                  n_metadata,
                  n_timeline,   n_metadata  - n_timeline,
                  n_encoded,
                  n_surv_merge, n_clustered - n_surv_merge,
                  n_surv_merge
  ))
  
  samples_no_timeline <- setdiff(
    unique(metadata[[ sample_col_md ]]),
    unique(dat$sample)
  )
  if (length(samples_no_timeline) > 0) {
    message(sprintf(
      "  Samples in metadata but not in timeline (%d): %s",
      length(samples_no_timeline),
      paste(head(samples_no_timeline, 10), collapse = ", ")
    ))
  }
  
  samples_no_surv <- cluster_results$sample[
    !cluster_results$sample %in% Cluster_surv$sample
  ]
  if (length(samples_no_surv) > 0) {
    message(sprintf(
      "  Samples dropped for missing survival data (%d): %s",
      length(samples_no_surv),
      paste(head(samples_no_surv, 10), collapse = ", ")
    ))
  }
  
  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  list(
    Cluster_surv             = Cluster_surv,
    cluster_results          = cluster_results,
    umap_df                  = umap_df,
    pca_matrix               = pca_matrix,
    pca_var_explained        = pca_var_explained,
    X                        = X,
    treatment_encoded        = treatment_encoded,
    treatment_matrix_ordered = treatment_matrix_ordered
  )
}
