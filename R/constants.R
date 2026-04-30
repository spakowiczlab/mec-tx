# constants.R
# -----------------------------------------------------------------------------
# Package-level constants for MEC-TX.
# Loaded automatically when the package is attached --- do not export.
# -----------------------------------------------------------------------------

#' Treatment type colour palette
#' Named character vector mapping standardised treatment type labels to hex
#' colours. Used by timeline, UMAP, and KM plot functions throughout the
#' package. Labels must match the recoding block in \code{tx_normalize()}.
#' @noRd
tx_cols <- c(
  Ancillary      = "#E1BE6A",  # supportive / ancillary drugs
  Chemo          = "#FDB863",  # cytotoxic chemotherapy
  Hormone        = "#DC267F",  # hormone therapy
  IO             = "#2CA02C",  # immunotherapy / checkpoint inhibitors
  Small_Molecule = "#76B7B2",  # targeted small molecules (EGFR/ALK inhibitors)
  Targeted       = "#4E79A7",  # other targeted agents
  Radiation      = "#6A51A3",  # radiation therapy
  Others         = "#8C8C8C"   # catch-all: unclassified, NA, or recoded Other
)

# All valid treatment types derived from tx_cols keys.
# Used for input validation in .treatment_exposure() and tx_normalize().
valid_types <- names(tx_cols)[names(tx_cols) != "None"]

# CAlevel two-group colour palette.
# Used as the default in km_panel_from_df() and cox_forest_plot_from_df()
# when group_col == "CAlevel" and both High and Low levels are present.
ca_cols <- c(High = "#F28E2B", Low = "#56B4E9")

# Convert internal Cluster_kN label to short display form (e.g. "k3").
# Used in plot titles and axis labels throughout the package.
k_label <- function(kc) sub("^Cluster_k(\\d+)$", "k\\1", kc)


# Suppress R CMD check notes for dplyr non-standard evaluation
utils::globalVariables(c(
  ".", ":=", "%>%",
  "AgeAtLastContact", "AgeAtMedStart", "AgeAtMedStop",
  "Age.At.Specimen.Collection", "AgeAtTreatmentStart.mod", "AgeGrid",
  "bin_idx", "block", "cluster", "cluster_lab", "cluster_num",
  "conf.high", "conf.low", "dev.off", "diagsurvtime", "dom_share",
  "dom_type_tx", "dur", "dur_treated", "duration", "end_year",
  "estimate", "Feature", "first_start", "first_t0", "first_time",
  "focus_share_tx", "gap", "grid_seq", "group", "grp", "has_all",
  "has_all_focus", "has_other", "label", "m", "med_focus",
  "n_focus_types", "n_patients", "n_risk", "new_block", "new_run",
  "only_focus_types", "order_score", "pdf", "regimen", "rowwise",
  "run", "sample_idx", "seq_ok", "share", "start_age", "start_year",
  "status", "t0", "t0_a", "t0_b", "t1", "t1_a", "t1_b",
  "t_year_raw", "TimeBin", "TimeSinceTreatmentStart", "time",
  "total", "total_dur", "treatment_group", "Treatment", "tx_seq",
  "type", "type_idx", "types_first", "Value", "y", "y_mid", "y_pos",
  "combn", "head"
))
