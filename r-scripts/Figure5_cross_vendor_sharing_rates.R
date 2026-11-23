# Figure5_cross_vendor_sharing_rates.R
# Cross-vendor strain-sharing rates for cefotaxime-resistant E. coli, with vendors
# grouped by vending area (block) in Burma Market.
#
# Produces: Figure 5
#
# Inputs:  data/processed/strainge/strain_sharing_within_vendor.csv
#          data/processed/strainge/strain_sharing_between_vendor.csv
#          data/processed/strainge/straingst_ecolidb_summary.csv
#          data/derived/lab/sample_collection_dates.csv
# Output:  results/figures/Figure5_cross_vendor_sharing_rates.svg

source("r-scripts/00_functions.R")
library(ComplexHeatmap)
library(circlize)
library(svglite)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

between_vendor_results <- load_vendor_sharing()$between_vendor_results

# Vendors ordered by vending area
vendor_area  <- c(`6` = "A", `12` = "A", `2` = "B", `4` = "B", `9` = "B", `11` = "B", `5` = "C", `8` = "C")
vendor_order <- c(6, 12, 2, 4, 9, 11, 5, 8)
area_order   <- factor(vendor_area[as.character(vendor_order)], levels = c("A", "B", "C"))

# Symmetric vendor x vendor matrix of sharing rates (diagonal undefined)
mat_bw <- matrix(NA_real_, nrow = length(vendor_order), ncol = length(vendor_order),
                 dimnames = list(vendor_order, vendor_order))
for (i in seq_len(nrow(between_vendor_results))) {
  v1 <- as.character(between_vendor_results$v1[i])
  v2 <- as.character(between_vendor_results$v2[i])
  mat_bw[v1, v2] <- between_vendor_results$rate[i]
  mat_bw[v2, v1] <- between_vendor_results$rate[i]
}

block_colors <- c("A" = "#8DA0CB", "B" = "#FC8D62", "C" = "#66C2A5")

ht <- Heatmap(
  mat_bw,
  name              = "Sharing rate",
  col               = colorRamp2(c(0, max(mat_bw, na.rm = TRUE)), c("white", "darkgoldenrod1")),
  na_col            = "slategray",
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  row_split         = area_order,
  column_split      = area_order,
  row_title         = NULL,
  column_title      = NULL,
  left_annotation   = HeatmapAnnotation(
    block = anno_block(gp = gpar(fill = block_colors[levels(area_order)]),
                       labels = c("Block A", "Block B", "Block C"),
                       labels_gp = gpar(col = "white", fontsize = 12)),
    which = "row", width = unit(5, "mm")),
  row_names_side    = "left",
  column_names_side = "bottom",
  rect_gp           = gpar(col = "grey80"),
  row_gap           = unit(2, "mm"),
  column_gap        = unit(2, "mm")
)

svglite("results/figures/Figure5_cross_vendor_sharing_rates.svg", width = 6, height = 6)
draw(ht)
dev.off()
