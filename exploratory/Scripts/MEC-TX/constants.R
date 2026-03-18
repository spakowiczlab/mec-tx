# constants.R
# ─────────────────────────────────────────────────────────────────────────────
# Package-level constants for MEC-TX.
# Loaded automatically when the package is attached — do not export.
# ─────────────────────────────────────────────────────────────────────────────

#' @noRd
#' Treatment type colour palette
#'
#' Named character vector mapping standardised treatment type labels to hex
#' colours. Used by timeline, UMAP, and KM plot functions throughout the
#' package. Labels must match the recoding block in \code{tx_normalize()}.
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