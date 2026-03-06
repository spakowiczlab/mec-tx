# ============================================================
# MEC-TX package constants
# helpers/constants.R
# ============================================================

tx_cols <- c(
  Ancillary      = "#E1BE6A",
  Chemo          = "#FDB863",
  Hormone        = "#DC267F",
  IO             = "#2CA02C",
  Small_Molecule = "#76B7B2",
  Targeted       = "#4E79A7",
  Radiation      = "#6A51A3",
  Others         = "#8C8C8C"
)

# All valid treatment types — derived from tx_cols, excluding "None"
valid_types <- names(tx_cols)[names(tx_cols) != "None"]

# Default colour palette for CAlevel (used as fallback in KM/forest)
ca_cols <- c(High = "#F28E2B", Low = "#56B4E9")

# Helper: convert Cluster_kN label to kN for display
k_label <- function(kc) sub("^Cluster_k(\\d+)$", "k\\1", kc)
