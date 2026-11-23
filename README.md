# kukusafe
Scripts and code for analyzing and visualizing cefotaxime-resistant E. coli strain-sharing in Burma Market, Kenya.

## Repository structure

```
kukusafe/
├── LICENSE            # MIT
├── r-scripts/         # R code (see below)
├── data/
│   ├── raw/           # StrainGE and KMA outputs
│   ├── reference/     # StrainGE reference database metadata; sample type codes
│   ├── derived/       # cleaned lab, survey, and observation data; phylogenetic tree
│   └── processed/     # created by r-scripts/01_data_preprocessing.R
├── docs/
│   └── data_dictionary.md
└── results/
    ├── figures/       # created by the figure scripts
    └── tables/        # created by the table scripts
```

## R scripts

### 1. Setup and data preprocessing

- `00_functions.R`: Shared constants (sample type labels) and helper functions, sourced by every other script.
- `01_data_preprocessing.R`: Builds every processed dataset in `data/processed/`:
  1. Summarizes StrainGST results per sample: *E. coli* relative abundance, strain count, and phylogroups.
  2. Combines StrainGR strain-sharing calls from both batches and splits them into within-vendor, between-vendor, and latrine datasets.
  3. Derives vendor-level covariates for the regression from the vendor survey and stall observations.
  4. Derives slaughter-point-level covariates for the regression from slaughter point observations.
  5. Builds all possible sample pairs with the strain-sharing outcome and pair-level covariates.
  6. Builds the within-vendor regression dataset.
  7. Builds the cross-vendor regression dataset.

### 2. Tables and figures

| Manuscript item | Description | Script |
|---|---|---|
| **Main text** | | |
| Table 1 | Vendor socio-demographics and chicken handling practices | `Table1_vendor_characteristics.R` |
| Table 2 | Risk factors for strain sharing (unadjusted and adjusted ORs) | `Table2_ExtDataTable4_strain_sharing_risk_factors.R` |
| Figure 1 | Sample collection timeline | Schematic; no code |
| Figure 2 | Phylogenetic tree with detection heatmaps by vendor and sample type | `Figure2_phylogenetic_tree.R` |
| Figure 3 | Unique and total strains by sample type | `Figure3_strains_by_sample_type.R` |
| Figure 4 | Strain sharing across sample types (A), reference strains involved (B), cross-vendor network for 2017.04.03CC (C) | `Figure4_strain_sharing_across_sample_types.R` |
| Figure 5 | Cross-vendor strain-sharing rates by vending area | `Figure5_cross_vendor_sharing_rates.R` |
| **Extended Data** | | |
| Extended Data Table 1 | Slaughter point infrastructure | `ExtDataTable1_slaughter_point_infrastructure.R` |
| Extended Data Table 2 | Vendor food safety and hygiene behaviors | `ExtDataTable2_vendor_hygiene_behaviors.R` |
| Extended Data Table 3 | Slaughter and processing practices | `ExtDataTable3_slaughter_processing_practices.R` |
| Extended Data Table 4 | Distribution of risk factors by strain-sharing status | `Table2_ExtDataTable4_strain_sharing_risk_factors.R` |
| Extended Data Figure 1 | Cefotaxime-resistant *E. coli* on chicken samples by sampling day | `ExtDataFigure1_SuppFigureS1_ecoli_load.R` |
| Extended Data Figure 2 | Resistome patterns by sample type | `ExtDataFigure2_ARG_heatmap.R` |
| Extended Data Figure 3 | Within-vendor strain-sharing rates (A) and reference strains (B) | `ExtDataFigure3_within_vendor_sharing.R` |
| Extended Data Figure 4 | Strain sharing between samples from the same chickens | `ExtDataFigure4_same_chicken_sharing.R` |
| Extended Data Figure 5 | Cross-vendor networks for MLI106K1 (A) and WP5-S17-ESBL-11 (B) | `ExtDataFigure5_cross_vendor_networks.R` |
| Extended Data Figure 6 | Timeline of strain sharing involving latrine waste | `ExtDataFigure6_latrine_sharing_timeline.R` |
| **Supplementary Information** | | |
| Supplementary Table S1 | Market spot-check data | No code |
| Supplementary Table S2 | Chicken cage conditions | `SuppTableS2_cage_conditions.R` |
| Supplementary Table S3 | StrainGR statistics for shared strains | `SuppTableS3_straingr_metrics.R` |
| Supplementary Table S4 | Limits of detection and quantification | No code |
| Supplementary Figure S1 | Cefotaxime-resistant *E. coli* load (A) and relative abundance (B) by sample type | `ExtDataFigure1_SuppFigureS1_ecoli_load.R` |
| Supplementary Figure S2 | Proportion of sample pairs with strain sharing, by sample-type dyad | `SuppFigureS2_sharing_by_sample_type_dyad.R` |

Output files in `results/` are named after the manuscript item they belong to. Final figure panels were assembled and annotated in Adobe Illustrator.

### R packages

```r
install.packages(c("tidyverse", "readxl", "broom", "broom.mixed", "lme4", "car", "scales",
                   "ggExtra", "ggh4x", "pheatmap", "circlize", "igraph", "ggraph", "svglite"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("treeio", "ggtree", "ComplexHeatmap"))
```

## Data availability

Raw sequencing reads are available from the NCBI Sequence Read Archive under BioProject [PRJNA1541341](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1541341). The de-identified data needed to run this code are in [`data/`](data/) and described in [`docs/data_dictionary.md`](docs/data_dictionary.md).
