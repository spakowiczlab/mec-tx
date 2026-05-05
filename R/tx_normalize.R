# =============================================================================
# MEC-TX --- tx_normalize.R
# analysis/tx_normalize.R
#
# Normalises raw medication data into a grid-based treatment timeline.
# Includes:
#   - Treatment group recoding (PI-confirmed mappings)
#   - Grid expansion at configurable resolution (grid_weeks)
#   - Metadata join with automatic status standardisation
#   - Per-patient dominant_regimen column
# =============================================================================

#' Normalise Raw Medication Data into a Grid-Based Treatment Timeline
#'
#' Converts raw per-medication-record data into a long-format treatment
#' timeline on a regular time grid. Handles missing stop dates, clips
#' exposure to specimen collection date, recodes non-canonical treatment
#' group labels (PI-confirmed mappings), and optionally joins patient
#' metadata and computes a per-patient dominant regimen assignment.
#'
#' @param med_data A data frame of raw medication records. Must contain
#'   columns: \code{sample}, \code{AgeAtMedStart}, \code{AgeAtMedStop},
#'   \code{Age.At.Specimen.Collection}, \code{AgeAtLastContact},
#'   \code{Medication}, and \code{treatment_group}.
#' @param metadata A data frame or \code{NULL}. Patient-level metadata
#'   joined to the timeline output. Must contain a \code{sample} column
#'   (case-insensitive). If supplied, \code{stage} and \code{status}
#'   columns are joined where present. \code{status} is automatically
#'   standardised to integer 0/1 via \code{\link{standardise_status}}.
#'   Default \code{NULL}.
#' @param grid_weeks Numeric. Time bin width in weeks. Controls the
#'   resolution of the output timeline grid. Must match the value used
#'   in downstream calls to \code{\link{tx_intervals}} and
#'   \code{\link{tx_cluster_surv}}. Common values: \code{1} (weekly),
#'   \code{2} (biweekly), \code{4} (monthly, default). Default \code{4}.
#' @param dominant_regimen_share Numeric in \code{(0, 1]}. Minimum
#'   treatment duration share threshold passed to
#'   \code{\link{dominant_exclusive}} for per-patient dominant regimen
#'   assignment. Default \code{0.20} (20\%).
#'
#' @return A tibble in long format with one row per patient -- grid time
#'   point -- treatment type. Columns:
#'   \describe{
#'     \item{sample}{Patient identifier.}
#'     \item{AgeGrid}{Numeric. Absolute age at each grid point (years).}
#'     \item{treatment_group}{Character. Recoded canonical treatment type.
#'       One of the eight types in \code{tx_cols} from
#'       \code{constants.R}.}
#'     \item{start_age}{Numeric. Age at first treatment grid point for
#'       this patient (used to compute \code{TimeSinceTreatmentStart}).}
#'     \item{TimeSinceTreatmentStart}{Numeric. Time since first treatment
#'       in years, snapped to the grid. This is the primary time axis for
#'       all downstream analyses.}
#'     \item{stage}{Character or \code{NA}. Joined from \code{metadata}
#'       if supplied.}
#'     \item{status}{Integer 0/1 or \code{NA}. Joined and standardised
#'       from \code{metadata} if supplied.}
#'     \item{end_followup}{Numeric. Age at last contact, joined from
#'       \code{med_data} for downstream truncation.}
#'     \item{dominant_regimen}{Character. Per-patient mutually exclusive
#'       dominant regimen label from \code{\link{dominant_exclusive}}.
#'       Patients with no qualifying non-ancillary treatment are labelled
#'       \code{"Ancillary/Supportive only"}.}
#'   }
#'
#' @details
#' \strong{Grid construction:} Each medication record is expanded into a
#' sequence of grid points from \code{AgeAtTreatmentStart.mod} to
#' \code{AgeAtMedStop} at \code{grid_weeks / 52} year intervals.
#' \code{TimeSinceTreatmentStart} is then computed as the grid-snapped
#' offset from each patient's earliest grid point.
#'
#' \strong{Specimen clipping:} Treatment exposure is clipped to begin no
#' earlier than \code{Age.At.Specimen.Collection}. Records where
#' \code{AgeAtMedStart < Age.At.Specimen.Collection} have their effective
#' start shifted forward. This prevents pre-diagnosis drug exposure from
#' contaminating the timeline.
#'
#' \strong{Missing stop dates:} \code{AgeAtMedStop} missing values are
#' imputed as \code{AgeAtLastContact}, treating the drug as ongoing until
#' last follow-up.
#'
#' \strong{Treatment group recoding (PI-confirmed):}
#' \describe{
#'   \item{\code{NA} / \code{"Other"}}{--- \code{"Others"}}
#'   \item{\code{"Onco_drug"} + Ado-Trastuzumab Emtansine}{---
#'     \code{"Targeted"} (HER2-directed)}
#'   \item{All other \code{"Onco_drug"}}{--- \code{"Others"} (Ramucirumab,
#'     Rituximab, investigational agents)}
#' }
#'
#' \strong{Dominant regimen:} A temporary interval table is built from
#' the cleaned medication data and passed to
#' \code{\link{dominant_exclusive}} to compute the per-patient dominant
#' regimen. This column is used in Cox models stratified by regimen type.
#'
#' \strong{Parameter alignment:} \code{grid_weeks} must match the value
#' used in \code{\link{tx_intervals}} and \code{\link{tx_cluster_surv}}.
#' Mismatches produce incorrect time bin alignment downstream.
#'
#' @examples
#' med_data <- data.frame(
#'   sample                     = rep(c('P01','P02'), each = 3),
#'   Age.At.Specimen.Collection = rep(c(60, 65), each = 3),
#'   AgeAtLastContact           = rep(c(62, 67), each = 3),
#'   diagsurvtime               = rep(c(2, 2), each = 3),
#'   Status                     = rep(c(1L, 0L), each = 3),
#'   Medication                 = rep(c('DrugA','DrugB','DrugC'), 2),
#'   treatment_group            = c('Chemo','IO','Radiation',
#'                                  'Chemo','Targeted','Others'),
#'   AgeAtMedStart              = c(60.1, 60.5, 61.0, 65.1, 65.4, 66.0),
#'   AgeAtMedStop               = c(60.4, 60.9, 61.3, 65.3, 65.8, 66.2),
#'   AgeAtTreatmentStart.mod    = c(60.1, 60.5, 61.0, 65.1, 65.4, 66.0),
#'   stringsAsFactors           = FALSE
#' )
#' norm <- tx_normalize(med_data)
#' head(norm)
#'
#' @seealso \code{\link{tx_intervals}}, \code{\link{tx_cluster_surv}},
#'   \code{\link{dominant_exclusive}}, \code{\link{standardise_status}}
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate filter select left_join group_by summarise rowwise
#'   ungroup transmute distinct case_when if_else rename rename_with
#'   all_of
#' @importFrom tidyr unnest
#' @export
tx_normalize <- function(
    med_data,                              # --- renamed from Modified_medication (Thread 8)
    metadata               = NULL,         # --- renamed from NSCLC_metadata (Thread 8)
    grid_weeks             = 4,
    dominant_regimen_share = 0.20
) {
  
  # ===========================================================================
  # INPUT VALIDATION
  # ===========================================================================
  
  # --- 1. med_data must be a data frame ---
  if (!is.data.frame(med_data)) {
    stop(
      "tx_normalize(): 'med_data' must be a data frame.\n",
      "  -> You passed an object of class: ",
      paste(class(med_data), collapse = ", "), ".\n",
      "  -> Load your data with readRDS() or read.csv() first."
    )
  }
  
  # --- 2. Required columns ---
  required_cols <- c(
    "sample", "AgeAtMedStart", "AgeAtMedStop",
    "Age.At.Specimen.Collection", "AgeAtLastContact",
    "Medication", "treatment_group"
  )
  missing_cols <- setdiff(required_cols, names(med_data))
  if (length(missing_cols) > 0) {
    stop(
      "tx_normalize(): Required column(s) missing from 'med_data':\n",
      paste0("  x ", missing_cols, collapse = "\n"), "\n",
      "  -> These columns must be present in your medication data frame.\n",
      "  -> Check your column names with: names(your_data)"
    )
  }
  
  # --- 3. sample column must not be all NA ---
  if (all(is.na(med_data$sample))) {
    stop(
      "tx_normalize(): 'sample' column is entirely NA.\n",
      "  -> Every row must have a patient identifier in the 'sample' column."
    )
  }
  
  # --- 4. grid_weeks must be a positive number ---
  if (!is.numeric(grid_weeks) || length(grid_weeks) != 1 || grid_weeks <= 0) {
    stop(
      "tx_normalize(): 'grid_weeks' must be a single positive number.\n",
      "  -> You passed: ", deparse(grid_weeks), "\n",
      "  -> Example: grid_weeks = 4 (monthly, default) or grid_weeks = 2 (biweekly)."
    )
  }
  
  # --- 5. dominant_regimen_share must be in (0, 1] ---
  if (!is.numeric(dominant_regimen_share) ||
      length(dominant_regimen_share) != 1 ||
      dominant_regimen_share <= 0 ||
      dominant_regimen_share > 1) {
    stop(
      "tx_normalize(): 'dominant_regimen_share' must be a number in (0, 1].\n",
      "  -> You passed: ", deparse(dominant_regimen_share), "\n",
      "  -> Example: dominant_regimen_share = 0.20 (default, 20% threshold)."
    )
  }
  
  # --- 6. Warn about non-canonical treatment_group values ---
  canonical_types <- c("Ancillary", "Chemo", "Hormone", "IO",
                       "Small_Molecule", "Targeted", "Radiation", "Others")
  auto_recoded    <- c("Onco_drug", "Other", NA)
  raw_types       <- unique(med_data$treatment_group)
  truly_unknown   <- setdiff(
    raw_types[!raw_types %in% c(canonical_types, auto_recoded)],
    canonical_types
  )
  if (length(truly_unknown) > 0) {
    warning(
      "tx_normalize(): Unrecognised treatment_group values --- recoded to 'Others':\n",
      paste0("  -> ", truly_unknown, collapse = "\n"), "\n",
      "  -> If any should map to a canonical type, update the recoding block."
    )
  }
  
  # --- 7. metadata checks (if supplied) ---
  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) {
      stop(
        "tx_normalize(): 'metadata' must be a data frame or NULL.\n",
        "  -> You passed an object of class: ",
        paste(class(metadata), collapse = ", "), "."
      )
    }
    meta_names_lower <- tolower(names(metadata))
    if (!"sample" %in% meta_names_lower) {
      stop(
        "tx_normalize(): 'metadata' must contain a 'sample' column.\n",
        "  -> Columns found: ", paste(names(metadata), collapse = ", ")
      )
    }
    overlap <- intersect(
      unique(as.character(med_data$sample)),
      unique(as.character(metadata[[which(meta_names_lower == "sample")[1]]]))
    )
    pct_overlap <- length(overlap) / length(unique(med_data$sample)) * 100
    if (pct_overlap < 50) {
      warning(sprintf(
        "tx_normalize(): Only %.0f%% of medication samples matched metadata samples.\n",
        pct_overlap
      ),
      "  -> Check that 'sample' IDs use the same format in both data frames."
      )
    }
  }
  
  # ===========================================================================
  # MAIN FUNCTION BODY
  # ===========================================================================
  
  # --- Clean and define treatment start ---
  df_clean <- med_data %>%
    dplyr::mutate(
      AgeAtMedStart              = as.numeric(AgeAtMedStart),
      AgeAtMedStop               = as.numeric(AgeAtMedStop),
      AgeAtLastContact           = as.numeric(AgeAtLastContact),
      Age.At.Specimen.Collection = as.numeric(Age.At.Specimen.Collection),
      AgeAtMedStop               = dplyr::if_else(
        is.na(AgeAtMedStop), AgeAtLastContact, AgeAtMedStop
      ),
      AgeAtTreatmentStart.mod    = pmax(
        Age.At.Specimen.Collection, AgeAtMedStart, na.rm = TRUE
      )
    ) %>%
    dplyr::filter(
      !is.na(AgeAtMedStart),
      !is.na(AgeAtTreatmentStart.mod),
      !is.na(AgeAtMedStop),
      AgeAtMedStop > AgeAtTreatmentStart.mod
    )
  
  if (nrow(df_clean) == 0) {
    stop(
      "tx_normalize(): No valid rows remain after filtering.\n",
      "  -> Check that AgeAtMedStart, AgeAtMedStop, and\n",
      "     Age.At.Specimen.Collection are numeric and that\n",
      "     AgeAtMedStop > AgeAtMedStart for at least some rows."
    )
  }
  
  # --- Recode non-canonical treatment_group labels (PI-confirmed) ---
  df_clean <- df_clean %>%
    dplyr::mutate(treatment_group = dplyr::case_when(
      is.na(treatment_group)                                ~ "Others",
      treatment_group == "Other"                            ~ "Others",
      treatment_group == "Onco_drug" &
        Medication == "Ado-Trastuzumab Emtansine"          ~ "Targeted",
      treatment_group == "Onco_drug"                       ~ "Others",
      TRUE                                                  ~ treatment_group
    ))
  
  # --- Expand to treatment grid ---
  expanded_timeline <- df_clean %>%
    rowwise() %>%
    dplyr::mutate(grid_seq = list(seq(
      AgeAtTreatmentStart.mod, AgeAtMedStop, by = grid_weeks / 52
    ))) %>%
    dplyr::ungroup() %>%
    tidyr::unnest(grid_seq) %>%
    dplyr::mutate(grid_seq = round(grid_seq, 4)) %>%
    dplyr::distinct(sample, treatment_group, grid_seq)
  
  timeline_long_grid <- expanded_timeline %>%
    dplyr::rename(AgeGrid = grid_seq) %>%
    dplyr::select(sample, AgeGrid, treatment_group)
  
  # --- Normalise to time since first treatment ---
  start_age_per_patient <- timeline_long_grid %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(start_age = min(AgeGrid, na.rm = TRUE), .groups = "drop")
  
  timeline_long_norm <- timeline_long_grid %>%
    dplyr::left_join(start_age_per_patient, by = "sample") %>%
    dplyr::mutate(
      TimeSinceTreatmentStart = round(
        (AgeGrid - start_age) * (52 / grid_weeks)
      ) / (52 / grid_weeks)
    )
  
  # --- Join metadata and standardise status ---
  if (!is.null(metadata)) {
    metadata <- metadata %>% dplyr::rename_with(tolower)
    
    if ("status" %in% names(metadata)) {
      metadata <- standardise_status(metadata, status_col = "status")
    }
    
    keep_cols <- intersect(
      c("sample", "stage", "status"), names(metadata)
    )
    timeline_long_norm <- timeline_long_norm %>%
      dplyr::left_join(
        metadata %>% dplyr::select(dplyr::all_of(keep_cols)),
        by = "sample"
      )
  }
  
  # --- Attach end_followup ---
  timeline_long_norm <- timeline_long_norm %>%
    dplyr::left_join(
      med_data %>%
        dplyr::transmute(
          sample       = sample,
          end_followup = as.numeric(AgeAtLastContact)
        ) %>%
        dplyr::distinct(),
      by = "sample"
    )
  
  # --- Compute dominant_regimen via dominant_exclusive() ---
  temp_intervals <- df_clean %>%
    dplyr::transmute(
      sample     = sample,
      type       = treatment_group,
      start_year = AgeAtTreatmentStart.mod,
      end_year   = AgeAtMedStop
    )
  
  dom_regimen <- dominant_exclusive(
    timeline        = temp_intervals,
    ancillary_types = c("Ancillary", "Others"),
    min_share       = dominant_regimen_share
  )
  
  timeline_long_norm <- timeline_long_norm %>%
    dplyr::left_join(dom_regimen, by = "sample") %>%
    dplyr::mutate(
      dominant_regimen = dplyr::if_else(
        is.na(regimen), "Ancillary/Supportive only", regimen
      )
    ) %>%
    dplyr::select(-regimen)
  
  timeline_long_norm
}
