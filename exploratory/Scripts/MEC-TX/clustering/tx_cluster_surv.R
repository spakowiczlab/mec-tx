#' Cluster treatment timelines and attach survival data
#'
#' Takes normalized treatment segments and patient metadata, performs
#' PCA + UMAP dimensionality reduction, clusters patients using
#' hierarchical clustering (k = 3 to 20), and attaches survival
#' and covariate data to the cluster assignments.


tx_cluster_surv <- function(
    NSCLC_metadata,
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
  # ---------------------------------------------------------------------------
  # 0) Resolve column names case-insensitively
  #    resolve_col() is a shared package helper from helpers/resolve_col.R
  # ---------------------------------------------------------------------------
  sample_col_tl <- resolve_col(timeline_long_norm, sample_col,    "sample_col")
  time_col_tl   <- resolve_col(timeline_long_norm, time_col,      "time_col")
  tx_col_tl     <- resolve_col(timeline_long_norm, tx_col,        "tx_col")
  
  sample_col_md <- resolve_col(NSCLC_metadata, sample_col,    "sample_col")
  surv_col_md   <- resolve_col(NSCLC_metadata, surv_time_col, "surv_time_col")
  status_col_md <- resolve_col(NSCLC_metadata, status_col,    "status_col")
  
  # meta_keep: NULL = keep all columns; otherwise soft resolve
  if (is.null(meta_keep)) {
    keep_meta <- setdiff(
      names(NSCLC_metadata),
      c(sample_col_md, surv_col_md, status_col_md)
    )
  } else {
    keep_meta <- character(0)
    for (m in meta_keep) {
      found <- names(NSCLC_metadata)[tolower(names(NSCLC_metadata)) == tolower(m)]
      if (length(found) == 1) {
        keep_meta <- c(keep_meta, found)
      } else {
        warning(sprintf(
          "meta_keep column '%s' not found in NSCLC_metadata — skipping.", m
        ))
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # 1) Build biweek-bin timeline
  #    Biweekly resolution (1/26 year ≈ 2 weeks) aligns with validated
  #    concurrency window from sensitivity analysis and standard oncology
  #    cycle boundaries (q3w chemo, q3w/q6w IO, radiation course transitions).
  # ---------------------------------------------------------------------------
  dat <- timeline_long_norm %>%
    transmute(
      sample                  = .data[[ sample_col_tl ]],
      TimeSinceTreatmentStart = as.numeric(.data[[ time_col_tl ]]),
      treatment_group         = as.character(.data[[ tx_col_tl ]])
    ) %>%
    filter(!is.na(sample), !is.na(TimeSinceTreatmentStart), !is.na(treatment_group)) %>%
    mutate(
      biweek_idx = round(TimeSinceTreatmentStart * 26),
      TimeBin    = sprintf("Biweek_%03d", biweek_idx)
    ) %>%
    filter(biweek_idx >= 0, biweek_idx <= horizon_years * 26)
  
  # ---------------------------------------------------------------------------
  # 2) Wide biweek matrix (for inspection / debugging)
  # ---------------------------------------------------------------------------
  treatment_matrix <- dat %>%
    group_by(sample, TimeBin) %>%
    summarise(
      Treatment = paste(sort(unique(treatment_group)), collapse = "+"),
      .groups   = "drop"
    ) %>%
    pivot_wider(names_from = TimeBin, values_from = Treatment, values_fill = NA)
  
  biweek_cols <- setdiff(names(treatment_matrix), "sample")
  treatment_matrix_ordered <- treatment_matrix[,
                                               c("sample", biweek_cols[order(as.numeric(stringr::str_extract(biweek_cols, "\\d+")))]),
                                               drop = FALSE
  ]
  
  # ---------------------------------------------------------------------------
  # 3) Restrict to horizon, fill gaps, long format
  # ---------------------------------------------------------------------------
  keep_biweeks <- sprintf("Biweek_%03d", 0:(horizon_years * 26))
  keep_cols    <- intersect(keep_biweeks, names(treatment_matrix_ordered))
  
  treatment_short <- treatment_matrix_ordered %>%
    select(sample, all_of(keep_cols)) %>%
    pivot_longer(cols = -sample, names_to = "Biweek", values_to = "Treatment")
  
  if (include_none) {
    treatment_short <- treatment_short %>%
      mutate(Treatment = replace_na(Treatment, "None"))
  }
  
  # ---------------------------------------------------------------------------
  # 4) Multi-hot encode: one binary column per (Biweek x TreatmentType)
  # ---------------------------------------------------------------------------
  treatment_encoded <- treatment_short %>%
    filter(!is.na(Treatment)) %>%
    mutate(Treatment = as.character(Treatment)) %>%
    tidyr::separate_rows(Treatment, sep = "\\+") %>%
    mutate(
      Treatment = stringr::str_trim(Treatment),
      Feature   = paste(Biweek, Treatment, sep = "_"),
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
  # 5) Feature matrix — remove near-zero-variance columns before PCA
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
  # 6) PCA (no scaling — NZV handled above; prcomp is deterministic)
  # ---------------------------------------------------------------------------
  pca_result <- prcomp(X, scale. = FALSE)
  n_take     <- min(n_pcs, ncol(pca_result$x))
  pca_matrix <- pca_result$x[, 1:n_take, drop = FALSE]
  rownames(pca_matrix) <- sl_ids
  
  pca_var_explained <- summary(pca_result)$importance[2, 1:n_take]
  
  # ---------------------------------------------------------------------------
  # 7) UMAP — visualization ONLY, not used for clustering
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
  # ---------------------------------------------------------------------------
  Surv_data <- NSCLC_metadata %>%
    filter(.data[[ sample_col_md ]] %in% cluster_results$sample) %>%
    select(
      sample       = all_of(sample_col_md),
      diagsurvtime = all_of(surv_col_md),
      status       = all_of(status_col_md),
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
  # Sample audit: track where samples are lost
  # ---------------------------------------------------------------------------
  n_metadata   <- length(unique(NSCLC_metadata[[ sample_col_md ]]))
  n_timeline   <- length(unique(dat$sample))
  n_encoded    <- nrow(X)
  n_clustered  <- nrow(cluster_results)
  n_surv_merge <- nrow(Cluster_surv)
  
  message(sprintf("
--- Sample audit ---
  In metadata:             %d
  In timeline (encoded):   %d   (lost: %d — not in timeline)
  After PCA/clustering:    %d
  After survival merge:    %d   (lost: %d — missing survival data)
  Final:                   %d
",
                  n_metadata,
                  n_timeline,   n_metadata  - n_timeline,
                  n_encoded,
                  n_surv_merge, n_clustered - n_surv_merge,
                  n_surv_merge
  ))
  
  samples_no_timeline <- setdiff(
    unique(NSCLC_metadata[[ sample_col_md ]]),
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