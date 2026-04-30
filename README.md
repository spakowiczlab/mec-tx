# mec-tx

MEchanistic Clustering - Treatment eXposure Framework

## Overview

mectx implements the MEC-TX framework for encoding, clustering, and
survival analysis of real-world oncology treatment timelines. It was
developed for registry-based cohorts such as the ORIEN AVATAR dataset.

## Installation

Install the released version from CRAN:

``` r
install.packages("mectx")
```

Or the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("spakowiczlab/mectx")
```

## Core functions

| Function | Description |
|----------|-------------|
| `tx_normalize()` | Normalize raw medication records to a time grid |
| `tx_intervals()` | Compute treatment intervals per patient |
| `tx_cluster_surv()` | K-means clustering in PCA space with survival output |
| `tx_lines()` | Assign line-of-therapy labels |
| `tx_pooled_analysis()` | Compare survival across treatment groups |
| `tx_duration()` | Summarise treatment exposure duration by group |
| `tx_compare_groups()` | Statistical comparison across patient groups |
| `dominant_exclusive()` | Assign mutually exclusive dominant regimen per patient |
| `get_focus_cohort()` | Filter cohort by focus treatment type |
| `tx_focus_dt()` | Build digital-twin timeline for focus treatment |

## Basic usage

``` r
library(mectx)

# Step 1: Normalize raw medication data
norm <- tx_normalize(raw_medication_df)

# Step 2: Compute treatment intervals
intervals <- tx_intervals(norm)

# Step 3: Cluster patients by treatment pattern
clustered <- tx_cluster_surv(norm, meta_df)

# Step 4: Assign line-of-therapy
lines <- tx_lines(intervals)

# Step 5: Pooled survival analysis
results <- tx_pooled_analysis(intervals, meta_df, group_var = "CAlevel")

# Step 6: Compare treatment duration by group
duration <- tx_duration(intervals, meta_df, group_var = "CAlevel")

```

## Pipeline

The canonical MEC-TX pipeline order:

raw data
└─ tx_normalize()
└─ tx_cluster_surv()
└─ tx_intervals()
└─ tx_lines()
└─ tx_pooled_analysis()
└─ tx_compare_groups()
└─ tx_duration()

## Citation

If you use mectx in your research, please cite: need to double  check

Dhrubo and Spakowicz (2026). MEchanistic Clustering - Treatment eXposure
Framework for real-world oncology treatment timeline analysis.

## License

MIT + file LICENSE

