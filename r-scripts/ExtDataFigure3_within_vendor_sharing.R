# ExtDataFigure3_within_vendor_sharing.R
# Within-vendor strain sharing of cefotaxime-resistant E. coli.
#
# Produces: Extended Data Figure 3A (within-vendor strain-sharing rates)
#           Extended Data Figure 3B (reference strains among within-vendor sharing events)
#
# Inputs:  data/processed/strainge/strain_sharing_within_vendor.csv
#          data/processed/strainge/strain_sharing_between_vendor.csv
#          data/processed/strainge/straingst_ecolidb_summary.csv
#          data/derived/lab/sample_collection_dates.csv
#          data/reference/strainge_refdb_mlst_phylotype.csv
# Outputs: results/figures/ExtDataFigure3A_within_vendor_rates.svg
#          results/figures/ExtDataFigure3B_within_vendor_reference_strains.svg

source("r-scripts/00_functions.R")
library(ggExtra)
library(scales)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

vendor_sharing  <- load_vendor_sharing()
wivendor        <- vendor_sharing$wivendor
byvendssresults <- vendor_sharing$byvendssresults
clermont        <- read_csv("data/reference/strainge_refdb_mlst_phylotype.csv")

# ==============================================================================
# Extended Data Figure 3A: Within-vendor strain-sharing rates
# ==============================================================================
# Dashed line: median within-vendor sharing rate
median_rate <- median(byvendssresults$rate, na.rm = TRUE)

fig3a <- ggplot(byvendssresults, aes(x = factor(vendor), y = rate)) +
  geom_col(fill = "lightblue", color = "black") +
  geom_text(aes(label = round(rate, digits = 2)), vjust = -0.3, size = 3) +
  geom_hline(yintercept = median_rate, linetype = "dashed", color = "red", linewidth = 1) +
  annotate("text", x = -Inf, y = median_rate, label = "Median", hjust = 1.1, vjust = -0.5, color = "red", size = 4) +
  labs(title = expression(atop("Strain sharing rate within vendors among samples with cefotaxime-resistant",
                               italic("E. coli") ~ "relative abundance above 0.016%")),
       x = "Vendor", y = "Strain sharing rate") +
  scale_y_continuous(limits = c(0, 0.4)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 13), plot.margin = margin(10, 10, 10, 10)) +
  removeGrid()

ggsave("results/figures/ExtDataFigure3A_within_vendor_rates.svg", plot = fig3a, width = 6, height = 6)

# ==============================================================================
# Extended Data Figure 3B: Reference strains among within-vendor sharing events
# ==============================================================================
byref <- wivendor %>%
  count(vendor, ref, name = "count") %>%
  left_join(clermont %>% select(symlink, MLST, phylotype), by = c("ref" = "symlink")) %>%
  mutate(ref_label = ref_label(ref, MLST, phylotype))

# Legend ordered by total number of sharing events across vendors
ref_label_order <- byref %>%
  group_by(ref_label) %>%
  summarise(total = sum(count), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(ref_label)

# Kelly's maximum-contrast palette (white and black dropped)
kelly_colors <- c("#F3C300", "#875692", "#F38400", "#A1CAF1", "#BE0032", "#C2B280", "#848482", "#008856",
                  "#E68FAC", "#0067A5", "#F99379", "#604E97", "#F6A600", "#B3446C", "#DCD300", "#882D17",
                  "#8DB600", "#654522", "#E25822", "#2B3D26")

fig3b <- byref %>%
  mutate(ref_label = factor(ref_label, levels = ref_label_order)) %>%
  ggplot(aes(x = factor(vendor), y = count, fill = ref_label)) +
  geom_col(position = "fill", color = "black", width = 0.7) +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = setNames(rep_len(kelly_colors, length(ref_label_order)), ref_label_order),
                    breaks = ref_label_order) +
  labs(x = "Vendor", y = "Proportion of sharing events", fill = "Reference strain") +
  guides(fill = guide_legend(nrow = 6, byrow = TRUE)) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  removeGrid()

ggsave("results/figures/ExtDataFigure3B_within_vendor_reference_strains.svg", plot = fig3b, width = 6, height = 6)
