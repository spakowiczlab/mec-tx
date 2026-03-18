# tx_lines.R
# Detect lines of therapy from treatment timelines.
#
# Uses a coalesce approach:
#   1. Clinician-annotated MedLineRegimen used where available (non-Unknown)
#   2. Gap-based algorithm for Unknown/Not Applicable records
#
# KEY DESIGN: Specimen-anchored record filtering
#   AVATAR captures all medications across a patient's lifetime, including
#   drugs from prior cancers (breast, lymphoma, transplant etc.). These are
#   filtered at the RECORD level (not patient level) by requiring:
#     AgeAtMedStart >= Age.At.Specimen.Collection - specimen_buffer
#   This retains all lung cancer records while dropping prior-cancer
#   contamination, without excluding any patients from the analysis.
#
# Architecture:
#   .map_line_regimen()        -- MedLineRegimen strings → standardised labels
#   .line_label_to_rank()      -- label → numeric rank (Maintenance/Palliative → NA)
#   .line_number_to_label()    -- numeric → label (1 → "First", 7 → "Line 7")
#   .merge_to_blocks()         -- overlapping/concurrent intervals → blocks
#   .assign_lines_from_blocks()-- walk blocks chronologically, increment on gap
#   .flag_consolidation()      -- mark IO-only lines in stage I–III
#   tx_lines()                 -- exported function

# ── helpers ───────────────────────────────────────────────────────────────────

#' @noRd
#'
#' Map MedLineRegimen strings to standardised line labels
#'
#' Converts free-text \code{MedLineRegimen} values from the AVATAR
#' Medication file to canonical line labels. Returns \code{NA_character_}
#' for Unknown, Not Applicable, and Not Reported entries — these are
#' handled by the gap-based algorithm in \code{.assign_lines_from_blocks()}.
#'
#' @param x Character vector of raw \code{MedLineRegimen} values.
#'
#' @return Character vector. One of \code{"Neoadjuvant"}, \code{"First"},
#'   \code{"Second"}, \code{"Third"}, \code{"Fourth"}, \code{"Fifth"},
#'   \code{"Sixth"}, \code{"Maintenance"}, \code{"Palliative"}, or
#'   \code{NA_character_}. Maintenance and Palliative are modifiers —
#'   not counted as new lines.
#'
#' @importFrom dplyr case_when
.map_line_regimen <- function(x) {
  dplyr::case_when(
    grepl("Neoadjuvant", x, ignore.case = TRUE)                          ~ "Neoadjuvant",
    grepl("Adjuvant/First|First Line/Regimen|^First Line$", x,
          ignore.case = TRUE)                                             ~ "First",
    grepl("Second Line/Regimen|^Second Line$", x, ignore.case = TRUE)    ~ "Second",
    grepl("Third Line/Regimen|^Third Line$",   x, ignore.case = TRUE)    ~ "Third",
    grepl("Fourth Line/Regimen|^Fourth Line$", x, ignore.case = TRUE)    ~ "Fourth",
    grepl("Fifth Line/Regimen|^Fifth Line$",   x, ignore.case = TRUE)    ~ "Fifth",
    grepl("Sixth Line/Regimen|^Sixth Line$",   x, ignore.case = TRUE)    ~ "Sixth",
    grepl("Maintenance",                       x, ignore.case = TRUE)    ~ "Maintenance",
    grepl("Palliative",                        x, ignore.case = TRUE)    ~ "Palliative",
    grepl("Unknown|Not Applicable|Not Reported", x, ignore.case = TRUE)  ~ NA_character_,
    TRUE                                                                  ~ NA_character_
  )
}


#' @noRd
#'
#' Convert line label to integer rank
#'
#' Maps canonical line labels to integer ranks for chronological ordering.
#' Maintenance and Palliative return \code{NA_integer_} — they are
#' continuations of an existing line, not new lines.
#'
#' @param label Character vector of canonical line labels (output of
#'   \code{.map_line_regimen()}).
#'
#' @return Integer vector. \code{0} = Neoadjuvant, \code{1} = First ...
#'   \code{6} = Sixth. \code{NA_integer_} for Maintenance, Palliative,
#'   and unrecognised labels.
.line_label_to_rank <- function(label) {
  ranks <- c(
    Neoadjuvant = 0L,
    First       = 1L, Second = 2L, Third  = 3L,
    Fourth      = 4L, Fifth  = 5L, Sixth  = 6L
  )
  vapply(label, function(l) {
    v <- ranks[l]
    if (is.na(v)) NA_integer_ else as.integer(v)
  }, integer(1L), USE.NAMES = FALSE)
}


#' @noRd
#'
#' Convert numeric line number to descriptive label
#'
#' Maps integer line numbers to human-readable labels for display in
#' output tables and plots.
#'
#' @param n Integer vector of line numbers. Must be \code{>= 0}.
#'
#' @return Character vector. \code{0} → \code{"Neoadjuvant"},
#'   \code{1} → \code{"First"} ... \code{6} → \code{"Sixth"},
#'   \code{7+} → \code{"Line 7"}, \code{"Line 8"}, etc.
.line_number_to_label <- function(n) {
  named <- c("Neoadjuvant", "First", "Second", "Third", "Fourth", "Fifth", "Sixth")
  ifelse(
    n >= 0L & n <= 6L,
    named[n + 1L],
    paste0("Line ", n)
  )
}


#' @noRd
#'
#' Merge overlapping or concurrent treatment intervals into blocks
#'
#' For a single patient, sorts intervals by start time and greedily merges
#' any that overlap or are concurrent. Concurrent treatment types are
#' concatenated with \code{"+"} (e.g. \code{"Chemo+IO"}).
#'
#' @param timeline_pt A data frame for a single patient with columns
#'   \code{start_year}, \code{end_year}, and \code{type}.
#'
#' @return A data frame with columns \code{block_start}, \code{block_end},
#'   and \code{block_types}. Returns a zero-row data frame if
#'   \code{timeline_pt} is empty.
.merge_to_blocks <- function(timeline_pt) {
  if (nrow(timeline_pt) == 0L) {
    return(data.frame(
      block_start = numeric(0), block_end = numeric(0),
      block_types = character(0), stringsAsFactors = FALSE
    ))
  }
  
  df <- timeline_pt[order(timeline_pt$start_year, timeline_pt$end_year), ]
  
  cur_start <- df$start_year[1L]
  cur_end   <- df$end_year[1L]
  cur_types <- df$type[1L]
  blocks    <- list()
  
  for (i in seq_len(nrow(df))[-1L]) {
    if (df$start_year[i] <= cur_end) {
      cur_end   <- max(cur_end, df$end_year[i])
      cur_types <- unique(c(cur_types, df$type[i]))
    } else {
      blocks[[length(blocks) + 1L]] <- data.frame(
        block_start = cur_start,
        block_end   = cur_end,
        block_types = paste(sort(cur_types), collapse = "+"),
        stringsAsFactors = FALSE
      )
      cur_start <- df$start_year[i]
      cur_end   <- df$end_year[i]
      cur_types <- df$type[i]
    }
  }
  
  blocks[[length(blocks) + 1L]] <- data.frame(
    block_start = cur_start,
    block_end   = cur_end,
    block_types = paste(sort(cur_types), collapse = "+"),
    stringsAsFactors = FALSE
  )
  
  do.call(rbind, blocks)
}


#' @noRd
#'
#' Assign line numbers by walking treatment blocks chronologically
#'
#' Walks treatment blocks in chronological order for a single patient and
#' increments the line counter whenever the gap between consecutive blocks
#' exceeds \code{gap_threshold}. Blocks within the same line are
#' consolidated into a single record.
#'
#' @param blocks A data frame — output of \code{.merge_to_blocks()} for
#'   one patient. Must contain \code{block_start}, \code{block_end}, and
#'   \code{block_types}.
#' @param gap_threshold Numeric. Gap in years above which a new line is
#'   triggered.
#'
#' @return A data frame with columns \code{line_number}, \code{line_start},
#'   \code{line_end}, \code{line_types}, \code{line_duration_months},
#'   \code{line_label}, and \code{line_source} (always \code{"computed"}
#'   at this stage — may be overridden by annotation coalesce in
#'   \code{tx_lines()}). Returns a zero-row data frame if \code{blocks}
#'   is empty.
.assign_lines_from_blocks <- function(blocks, gap_threshold) {
  if (nrow(blocks) == 0L) {
    return(data.frame(
      line_number = integer(0), line_start = numeric(0),
      line_end    = numeric(0), line_types = character(0),
      line_duration_months = numeric(0), line_label = character(0),
      line_source = character(0), stringsAsFactors = FALSE
    ))
  }
  
  blocks       <- blocks[order(blocks$block_start), ]
  line_ids     <- integer(nrow(blocks))
  current_line <- 1L
  line_ids[1L] <- current_line
  
  for (i in seq_len(nrow(blocks))[-1L]) {
    gap <- blocks$block_start[i] - blocks$block_end[i - 1L]
    if (gap > gap_threshold) current_line <- current_line + 1L
    line_ids[i] <- current_line
  }
  
  blocks$line_id <- line_ids
  
  out <- do.call(rbind, lapply(unique(line_ids), function(lid) {
    b         <- blocks[blocks$line_id == lid, ]
    all_types <- unique(unlist(strsplit(b$block_types, "\\+")))
    data.frame(
      line_number          = lid,
      line_start           = min(b$block_start),
      line_end             = max(b$block_end),
      line_types           = paste(sort(all_types), collapse = "+"),
      line_duration_months = (max(b$block_end) - min(b$block_start)) * 12,
      line_label           = .line_number_to_label(lid),
      line_source          = "computed",
      stringsAsFactors     = FALSE
    )
  }))
  
  out[order(out$line_number), ]
}


#' @noRd
#'
#' Flag possible consolidation lines
#'
#' Marks IO-only lines after line 1 in early-stage patients as
#' \code{"possible_consolidation"}. These may represent durvalumab or
#' adjuvant pembrolizumab rather than a true new therapy line and require
#' PI review before counting. Line 1 is always \code{"confirmed"}.
#'
#' @param lines_pt A data frame of lines for a single patient — output of
#'   \code{.assign_lines_from_blocks()}.
#' @param continuation_types Character vector. Treatment types that may
#'   represent consolidation (e.g. \code{"IO"}).
#' @param continuation_stages Character vector. Stage values where
#'   consolidation therapy is expected (e.g. \code{c("I","II","III")}).
#' @param stage Character scalar or \code{NA}. Stage value for this
#'   patient.
#'
#' @return The input \code{lines_pt} data frame with a \code{line_flag}
#'   column added. Values: \code{"confirmed"} or
#'   \code{"possible_consolidation"}.
.flag_consolidation <- function(lines_pt, continuation_types,
                                continuation_stages, stage) {
  lines_pt$line_flag <- "confirmed"
  
  if (!is.null(stage) && !is.na(stage) &&
      stage %in% continuation_stages &&
      nrow(lines_pt) > 1L) {
    
    for (i in seq_len(nrow(lines_pt))) {
      if (lines_pt$line_number[i] == 1L) next
      types_in_line <- trimws(unlist(strsplit(lines_pt$line_types[i], "\\+")))
      if (all(types_in_line %in% continuation_types)) {
        lines_pt$line_flag[i] <- "possible_consolidation"
      }
    }
  }
  
  lines_pt
}


# ── exported function ─────────────────────────────────────────────────────────

#' Detect Lines of Therapy from Treatment Timelines
#'
#' Identifies lines of therapy for each patient by combining
#' clinician-annotated \code{MedLineRegimen} labels (where available and
#' non-Unknown) with a gap-based algorithm for unannotated records.
#' Applies specimen-anchored record filtering to remove prior-cancer drug
#' contamination from the AVATAR registry before any line detection.
#'
#' @param timeline A data frame — the direct output of
#'   \code{\link{tx_intervals}}. Must contain columns \code{sample},
#'   \code{start_year}, \code{end_year}, and \code{type}.
#' @param annotations A data frame or \code{NULL}. The AVATAR Medication
#'   annotation file containing \code{MedLineRegimen} labels. Must contain
#'   columns named in \code{ann_id_col}, \code{ann_line_col}, and
#'   \code{ann_start_col}. \code{NULL} runs algorithm-only mode with no
#'   annotation coalescing. Default \code{NULL}.
#' @param meta A data frame or \code{NULL}. Patient-level metadata
#'   (e.g. \code{Cluster_surv} or \code{LUAD_metadata}). Must contain
#'   \code{sample}. Used for: grouping variable join (\code{group_var}),
#'   stage extraction (\code{stage_col}), \code{AvatarKey} crosswalk for
#'   annotation joining, and specimen age for record-level filtering.
#'   Default \code{NULL}.
#' @param group_var Character string or \code{NULL}. Column in \code{meta}
#'   for group comparison (e.g. \code{"CAlevel"}). If supplied, a
#'   Wilcoxon / Kruskal-Wallis comparison of \code{n_lines} and
#'   \code{first_line_duration_months} is computed between groups.
#'   Default \code{NULL}.
#' @param ann_id_col Character string. Patient identifier column in
#'   \code{annotations}. Used for the \code{AvatarKey} → \code{sample}
#'   crosswalk via \code{meta}. Default \code{"AvatarKey"}.
#' @param ann_line_col Character string. Line annotation column in
#'   \code{annotations}. Default \code{"MedLineRegimen"}.
#' @param ann_start_col Character string. Treatment start age column in
#'   \code{annotations} (years). Used for specimen-anchored filtering of
#'   annotation records. Default \code{"AgeAtMedStart"}.
#' @param gap_threshold Numeric. Gap in years above which consecutive
#'   treatment blocks are assigned to a new line. Default \code{3/52}
#'   (~3 weeks). PI-confirmed threshold for AVATAR data.
#' @param continuation_types Character vector. Treatment types that may
#'   represent consolidation / adjuvant therapy rather than a true new
#'   line. IO-only lines containing only these types are flagged
#'   \code{"possible_consolidation"} in stage I–III patients. Default
#'   \code{"IO"}.
#' @param continuation_stages Character vector. Stage values where
#'   consolidation flagging is applied. Default \code{c("I","II","III")}.
#' @param stage_col Character string or \code{NULL}. Column in \code{meta}
#'   containing cancer stage, used for consolidation flagging. \code{NULL}
#'   disables flagging. Default \code{NULL}.
#' @param exclude_types Character vector or \code{NULL}. Treatment types
#'   to remove from \code{timeline} before line detection. Recommended:
#'   \code{c("Ancillary", "Others")} to avoid short ancillary records
#'   triggering spurious line breaks. Default \code{NULL}.
#' @param specimen_age_col Character string. Column in \code{meta}
#'   containing specimen collection age in years. Used for specimen-
#'   anchored record filtering in both \code{timeline} and
#'   \code{annotations}. Default \code{"Age.At.Specimen.Collection"}.
#' @param specimen_buffer Numeric. Years before specimen collection date
#'   to allow medication records. Records with start age earlier than
#'   \code{specimen_age - specimen_buffer} are dropped. Default \code{0.25}
#'   (3 months) — covers the window where lung cancer treatment begins
#'   just before biopsy processing. Increase cautiously; large values
#'   risk re-introducing prior-cancer contamination.
#' @param filter_timeline Logical. If \code{TRUE}, apply specimen-anchored
#'   filtering to \code{timeline} intervals. Set to \code{FALSE} when
#'   \code{timeline} uses relative time from \code{TimeSinceTreatmentStart}
#'   (cannot be compared to absolute specimen age). Annotation filtering
#'   always runs regardless of this flag. Default \code{TRUE}.
#'
#' @return A named list with four elements:
#'   \describe{
#'     \item{lines}{Data frame with one row per patient × line.
#'       Columns: \code{sample}, \code{line_number}, \code{line_label},
#'       \code{line_types}, \code{line_start}, \code{line_end},
#'       \code{line_duration_months},
#'       \code{line_source} (\code{"annotated"} or \code{"computed"}),
#'       \code{line_flag} (\code{"confirmed"} or
#'       \code{"possible_consolidation"}).}
#'     \item{patient_summary}{Data frame with one row per patient.
#'       Columns: \code{sample}, \code{n_lines}, \code{max_line},
#'       \code{first_line_label}, \code{first_line_types},
#'       \code{first_line_duration_months},
#'       \code{n_possible_consolidation}.}
#'     \item{group_comparison}{Data frame of Wilcoxon / Kruskal-Wallis
#'       test results by \code{group_var} for \code{n_lines} and
#'       \code{first_line_duration_months}. \code{NULL} if \code{group_var}
#'       is not supplied.}
#'     \item{params}{Named list recording all function settings for
#'       reproducibility.}
#'   }
#'   Returns a list with an empty \code{lines} data frame and a warning
#'   if no lines are detected.
#'
#' @details
#' \strong{Specimen-anchored record filtering:} AVATAR captures all
#' medications across a patient's lifetime, including drugs for prior
#' cancers. Records are filtered at the individual record level — any
#' record with start age earlier than
#' \code{specimen_age - specimen_buffer} is dropped. This removes
#' prior-cancer contamination (e.g. Letrozole for breast cancer) without
#' excluding any patients from the cohort.
#'
#' \strong{Coalesce logic:} For each patient, if a non-Unknown
#' \code{MedLineRegimen} annotation exists after filtering, it overrides
#' the algorithm-computed label for line 1. Subsequent lines rely on the
#' gap algorithm. Line 1 is the only line that can be annotation-coalesced
#' — higher lines are always algorithm-derived.
#'
#' \strong{Consolidation flagging:} IO-only lines after line 1 in stage
#' I–III patients are flagged \code{"possible_consolidation"} and must be
#' reviewed by the PI before counting as true new therapy lines. These
#' likely represent durvalumab consolidation or adjuvant pembrolizumab.
#'
#' \strong{filter_timeline = FALSE:} Use this when \code{timeline} was
#' produced from \code{tx_intervals()} with relative
#' \code{TimeSinceTreatmentStart} values. These cannot be compared to
#' absolute \code{AgeAtMedStart} values and the filter would drop all
#' records. Annotation filtering still runs.
#'
#' @examples
#' \dontrun{
#' intervals <- tx_intervals(tx_normalize(med_data, metadata))
#'
#' # Algorithm-only mode (no annotations)
#' res <- tx_lines(
#'   timeline  = intervals,
#'   meta      = Cluster_surv,
#'   group_var = "CAlevel",
#'   exclude_types = c("Ancillary", "Others")
#' )
#' res$patient_summary
#' res$group_comparison
#'
#' # With MedLineRegimen annotations (LUAD)
#' res2 <- tx_lines(
#'   timeline     = intervals,
#'   annotations  = LUAD_medication,
#'   meta         = LUAD_metadata,
#'   group_var    = "CAlevel",
#'   stage_col    = "stage_group",
#'   exclude_types = c("Ancillary", "Others"),
#'   filter_timeline = FALSE
#' )
#'
#' # Access lines data frame
#' head(res2$lines)
#'
#' # Consolidation-flagged lines for PI review
#' res2$lines[res2$lines$line_flag == "possible_consolidation", ]
#' }
#'
#' @seealso \code{\link{tx_intervals}}, \code{\link{tx_duration}},
#'   \code{\link{tx_pooled_analysis}}, \code{\link{dominant_exclusive}}
#'
#' @importFrom dplyr case_when
#' @importFrom stats wilcox.test kruskal.test median quantile
#' @export

tx_lines <- function(
    timeline,
    annotations          = NULL,
    meta                 = NULL,
    group_var            = NULL,
    ann_id_col           = "AvatarKey",
    ann_line_col         = "MedLineRegimen",
    ann_start_col        = "AgeAtMedStart",
    gap_threshold        = 3 / 52,
    continuation_types   = c("IO"),
    continuation_stages  = c("I", "II", "III"),
    stage_col            = NULL,
    exclude_types        = NULL,
    specimen_age_col     = "Age.At.Specimen.Collection",
    specimen_buffer      = 0.25,
    filter_timeline      = TRUE
) {
  
  # ── 0. Input validation ─────────────────────────────────────────────────────
  if (!is.data.frame(timeline))
    stop("[tx_lines] 'timeline' must be a data.frame (tx_intervals() output).")
  
  required_tl <- c("sample", "start_year", "end_year", "type")
  missing_tl  <- setdiff(required_tl, names(timeline))
  if (length(missing_tl) > 0L)
    stop("[tx_lines] 'timeline' missing columns: ", paste(missing_tl, collapse = ", "))
  
  if (!is.null(annotations)) {
    if (!is.data.frame(annotations))
      stop("[tx_lines] 'annotations' must be a data.frame.")
    for (col in c(ann_id_col, ann_line_col, ann_start_col)) {
      if (!col %in% names(annotations))
        stop("[tx_lines] 'annotations' missing column: ", col)
    }
  }
  
  if (!is.numeric(gap_threshold) || gap_threshold <= 0)
    stop("[tx_lines] 'gap_threshold' must be a positive number (years).")
  
  if (!is.numeric(specimen_buffer) || specimen_buffer < 0)
    stop("[tx_lines] 'specimen_buffer' must be >= 0 (years).")
  
  # ── 1. Build crosswalk and specimen ages from meta ──────────────────────────
  crosswalk    <- NULL
  specimen_lut <- NULL   # lookup: sample → specimen age
  
  if (!is.null(meta)) {
    if (!is.data.frame(meta))
      stop("[tx_lines] 'meta' must be a data.frame.")
    
    if (all(c(ann_id_col, "sample") %in% names(meta))) {
      crosswalk <- unique(meta[, c(ann_id_col, "sample"), drop = FALSE])
    }
    
    if (specimen_age_col %in% names(meta)) {
      specimen_lut <- unique(meta[, c("sample", specimen_age_col), drop = FALSE])
    }
  }
  
  # ── 2. Specimen-anchored record filtering — TIMELINE ────────────────────────
  # Skipped when filter_timeline = FALSE (relative time scale from
  # tx_intervals / TimeSinceTreatmentStart cannot be compared to absolute
  # specimen age — would drop all records). Annotations filter always runs.
  
  n_before <- nrow(timeline)
  
  if (isTRUE(filter_timeline) && !is.null(specimen_lut)) {
    timeline <- merge(timeline, specimen_lut, by = "sample", all.x = TRUE)
    keep <- is.na(timeline[[specimen_age_col]]) |
      timeline$start_year >= (timeline[[specimen_age_col]] - specimen_buffer)
    timeline <- timeline[keep, , drop = FALSE]
    timeline[[specimen_age_col]] <- NULL   # drop helper column
    
    n_dropped <- n_before - nrow(timeline)
    message(sprintf(
      "[tx_lines] Specimen filter (timeline): dropped %d of %d records",
      n_dropped, n_before
    ))
  } else if (isFALSE(filter_timeline)) {
    message(
      "[tx_lines] Timeline filter skipped (filter_timeline = FALSE). ",
      "Annotations filter still applied."
    )
  } else {
    message(
      "[tx_lines] Specimen filter skipped for timeline: ",
      "'", specimen_age_col, "' not found in meta."
    )
  }
  
  # ── 3. Exclude treatment types ──────────────────────────────────────────────
  if (!is.null(exclude_types) && length(exclude_types) > 0L) {
    timeline <- timeline[!timeline$type %in% exclude_types, , drop = FALSE]
  }
  
  # ── 4. Specimen-anchored record filtering — ANNOTATIONS ─────────────────────
  # Same logic applied to the Medication file before any annotation join.
  
  ann_mapped <- NULL
  
  if (!is.null(annotations) && !is.null(crosswalk)) {
    
    # Join sample → specimen age via AvatarKey crosswalk
    ann_xwalk <- merge(annotations, crosswalk,     by = ann_id_col, all.x = TRUE)
    ann_xwalk <- merge(ann_xwalk,   specimen_lut,  by = "sample",   all.x = TRUE)
    
    n_ann_before <- nrow(ann_xwalk)
    keep_ann <- is.na(ann_xwalk[[specimen_age_col]]) |
      ann_xwalk[[ann_start_col]] >= (ann_xwalk[[specimen_age_col]] - specimen_buffer)
    ann_xwalk <- ann_xwalk[keep_ann, , drop = FALSE]
    
    message(sprintf(
      "[tx_lines] Specimen filter (annotations): dropped %d of %d records",
      n_ann_before - nrow(ann_xwalk), n_ann_before
    ))
    
    # Map MedLineRegimen → standardised labels; drop Unknown (→ algorithm handles)
    ann_xwalk$line_label_ann <- .map_line_regimen(ann_xwalk[[ann_line_col]])
    ann_xwalk$line_rank_ann  <- .line_label_to_rank(ann_xwalk$line_label_ann)
    
    ann_xwalk <- ann_xwalk[
      !is.na(ann_xwalk$line_label_ann) &
        !ann_xwalk$line_label_ann %in% c("Maintenance", "Palliative"),
    ]
    
    # Per patient: keep the earliest / lowest-ranked annotation to anchor line 1
    if (nrow(ann_xwalk) > 0L) {
      ann_xwalk <- ann_xwalk[order(ann_xwalk$sample, ann_xwalk$line_rank_ann,
                                   ann_xwalk[[ann_start_col]]), ]
      ann_mapped <- ann_xwalk[!duplicated(ann_xwalk$sample),
                              c("sample", "line_label_ann", "line_rank_ann"),
                              drop = FALSE]
    }
  }
  
  # ── 5. Per-patient line detection ───────────────────────────────────────────
  samples   <- unique(timeline$sample)
  all_lines <- vector("list", length(samples))
  
  for (si in seq_along(samples)) {
    s     <- samples[si]
    tl_pt <- timeline[timeline$sample == s, , drop = FALSE]
    if (nrow(tl_pt) == 0L) next
    
    # Stage for consolidation flagging
    stage_val <- NA_character_
    if (!is.null(stage_col) && !is.null(meta) && stage_col %in% names(meta)) {
      sv <- meta[meta$sample == s, stage_col, drop = TRUE]
      if (length(sv) > 0L) stage_val <- as.character(sv[1L])
    }
    
    # Merge overlapping intervals → blocks
    blocks <- .merge_to_blocks(tl_pt)
    if (nrow(blocks) == 0L) next
    
    # Assign line numbers via gap algorithm
    lines_pt <- .assign_lines_from_blocks(blocks, gap_threshold)
    if (nrow(lines_pt) == 0L) next
    
    # Coalesce: override line 1 label with annotation if available
    if (!is.null(ann_mapped)) {
      ann_pt <- ann_mapped[ann_mapped$sample == s, , drop = FALSE]
      if (nrow(ann_pt) > 0L) {
        idx <- lines_pt$line_number == 1L
        lines_pt$line_label[idx]  <- ann_pt$line_label_ann[1L]
        lines_pt$line_source[idx] <- "annotated"
      }
    }
    
    # Flag possible consolidation
    lines_pt <- .flag_consolidation(
      lines_pt, continuation_types, continuation_stages, stage_val
    )
    
    lines_pt$sample <- s
    all_lines[[si]] <- lines_pt
  }
  
  lines_df <- do.call(rbind, all_lines)
  
  if (is.null(lines_df) || nrow(lines_df) == 0L) {
    warning("[tx_lines] No lines detected. Check timeline and filter settings.")
    return(list(
      lines            = data.frame(),
      patient_summary  = NULL,
      group_comparison = NULL,
      params           = .tx_lines_params(match.call())
    ))
  }
  
  # Reorder columns
  col_order <- c("sample", "line_number", "line_label", "line_types",
                 "line_start", "line_end", "line_duration_months",
                 "line_source", "line_flag")
  lines_df <- lines_df[, intersect(col_order, names(lines_df)), drop = FALSE]
  rownames(lines_df) <- NULL
  
  # ── 6. Patient-level summary ────────────────────────────────────────────────
  patient_summary <- do.call(rbind, lapply(unique(lines_df$sample), function(s) {
    pt <- lines_df[lines_df$sample == s, ]
    ln1 <- pt[pt$line_number == min(pt$line_number), ][1L, ]
    data.frame(
      sample                     = s,
      n_lines                    = max(pt$line_number),
      max_line                   = max(pt$line_number),
      first_line_label           = ln1$line_label,
      first_line_types           = ln1$line_types,
      first_line_duration_months = ln1$line_duration_months,
      n_possible_consolidation   = sum(pt$line_flag == "possible_consolidation",
                                       na.rm = TRUE),
      stringsAsFactors           = FALSE
    )
  }))
  
  rownames(patient_summary) <- NULL
  
  # ── 7. Group comparison ─────────────────────────────────────────────────────
  group_comparison <- NULL
  
  if (!is.null(group_var) && !is.null(meta) && group_var %in% names(meta)) {
    
    grp_lut <- unique(meta[, c("sample", group_var), drop = FALSE])
    grp_df  <- merge(patient_summary, grp_lut, by = "sample", all.x = TRUE)
    grp_df  <- grp_df[!is.na(grp_df[[group_var]]), ]
    
    groups  <- unique(grp_df[[group_var]])
    
    .run_test <- function(vals_list) {
      vals_list <- lapply(vals_list, function(v) v[!is.na(v)])
      ns <- sapply(vals_list, length)
      if (any(ns < 3L)) return(list(p = NA_real_, note = "skipped (min_n < 3)"))
      if (length(vals_list) == 2L) {
        p <- tryCatch(
          wilcox.test(vals_list[[1L]], vals_list[[2L]], exact = FALSE)$p.value,
          error = function(e) NA_real_
        )
        list(p = p, note = "Wilcoxon")
      } else {
        p <- tryCatch(
          kruskal.test(vals_list)$p.value,
          error = function(e) NA_real_
        )
        list(p = p, note = "Kruskal-Wallis")
      }
    }
    
    metrics <- c("n_lines", "first_line_duration_months")
    
    rows <- lapply(metrics, function(m) {
      vals <- lapply(groups, function(g)
        grp_df[grp_df[[group_var]] == g, m, drop = TRUE])
      tr   <- .run_test(vals)
      do.call(rbind, lapply(seq_along(groups), function(i) {
        v <- vals[[i]][!is.na(vals[[i]])]
        data.frame(
          metric    = m,
          group     = as.character(groups[i]),
          n         = length(v),
          mean      = round(mean(v),   3L),
          median    = round(median(v), 3L),
          q25       = round(quantile(v, 0.25), 3L),
          q75       = round(quantile(v, 0.75), 3L),
          p_value   = round(tr$p, 4L),
          test_note = tr$note,
          stringsAsFactors = FALSE
        )
      }))
    })
    
    group_comparison <- do.call(rbind, rows)
    rownames(group_comparison) <- NULL
  }
  
  # ── 8. Return ────────────────────────────────────────────────────────────────
  list(
    lines            = lines_df,
    patient_summary  = patient_summary,
    group_comparison = group_comparison,
    params           = list(
      gap_threshold        = gap_threshold,
      specimen_buffer      = specimen_buffer,
      specimen_age_col     = specimen_age_col,
      continuation_types   = continuation_types,
      continuation_stages  = continuation_stages,
      exclude_types        = exclude_types,
      group_var            = group_var,
      ann_id_col           = ann_id_col,
      ann_line_col         = ann_line_col,
      ann_start_col        = ann_start_col,
      filter_timeline      = filter_timeline
    )
  )
}

# ── internal helper ───────────────────────────────────────────────────────────

.tx_lines_params <- function(call) {
  # Captures call arguments for params slot when early-return occurs
  list(call = deparse(call))
}