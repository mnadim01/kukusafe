# ExtDataFigure5_cross_vendor_networks.R
# Cross-vendor strain-sharing networks for two reference strains.
#
# Produces: Extended Data Figure 5A (MLST 38, MLI106K1)
#           Extended Data Figure 5B (MLST 58, WP5-S17-ESBL-11)
#
# Inputs:  data/processed/strainge/strain_sharing_within_vendor.csv
#          data/processed/strainge/strain_sharing_between_vendor.csv
#          data/processed/strainge/straingst_ecolidb_summary.csv
#          data/derived/lab/sample_collection_dates.csv
# Outputs: results/figures/ExtDataFigure5A_network_MLI106K1.svg
#          results/figures/ExtDataFigure5B_network_WP5-S17-ESBL-11.svg

source("r-scripts/00_functions.R")
library(igraph)
library(ggraph)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

bwvendor <- load_vendor_sharing()$bwvendor

ggsave("results/figures/ExtDataFigure5A_network_MLI106K1.svg",
       plot = plot_bwvend_network("Esch_coli_MLI106K1.fa.gz", bwvendor), width = 6, height = 6)
ggsave("results/figures/ExtDataFigure5B_network_WP5-S17-ESBL-11.svg",
       plot = plot_bwvend_network("Esch_coli_WP5-S17-ESBL-11.fa.gz", bwvendor), width = 6, height = 6)
