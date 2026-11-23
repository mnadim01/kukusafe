# Figure4_strain_sharing_across_sample_types.R
# Cefotaxime-resistant E. coli strain sharing across sample types.
#
# Produces: Figure 4A (chord diagram of strain sharing across sample types)
#           Figure 4B (reference strains involved in strain sharing)
#           Figure 4C (cross-vendor network for reference strain 2017.04.03CC)
#
# Inputs:  data/processed/regression/all_strain_sharing_pairs.csv
#          data/processed/strainge/strain_sharing_by_sample_all.csv
#          data/processed/strainge/strain_sharing_within_vendor.csv
#          data/processed/strainge/strain_sharing_between_vendor.csv
#          data/processed/strainge/straingst_ecolidb_summary.csv
#          data/derived/lab/sample_collection_dates.csv
#          data/reference/strainge_refdb_mlst_phylotype.csv
# Outputs: results/figures/Figure4A_chord_diagram.svg
#          results/figures/Figure4B_reference_strains.svg
#          results/figures/Figure4C_network_2017_04_03CC.svg

source("r-scripts/00_functions.R")
library(circlize)
library(svglite)
library(igraph)
library(ggraph)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

allpairs_raw <- read_csv("data/processed/regression/all_strain_sharing_pairs.csv")
byall_raw    <- read_csv("data/processed/strainge/strain_sharing_by_sample_all.csv")
clermont     <- read_csv("data/reference/strainge_refdb_mlst_phylotype.csv")

# ==============================================================================
# Figure 4A: Chord diagram
# ==============================================================================
# Unique sharing pairs counted once per sample-type combination
chord_df <- allpairs_raw %>%
  filter(strain_sharing == 1) %>%
  mutate(type_smp1 = dplyr::recode(smp_type_smp1, !!!smp_code_map),
         type_smp2 = dplyr::recode(smp_type_smp2, !!!smp_code_map),
         from      = pmin(type_smp1, type_smp2),
         to        = pmax(type_smp1, type_smp2)) %>%
  distinct(sample1, sample2, from, to) %>%
  count(from, to, name = "value")

# Sectors grouped by host category; larger gap after the last sector in each category
category_order <- c("Cage fecal droppings", "Live chicken rectal swab", "Pre-evisceration carcass", "Ready-to-sell carcass",
                    "Food fed to chickens", "Vendor hand rinse", "Slaughterhouse table surface", "Water fed to chickens",
                    "Latrine waste", "Vendor feces")
sector_order <- category_order[category_order %in% unique(c(chord_df$from, chord_df$to))]
gap_vec <- setNames(ifelse(sector_order %in% c("Ready-to-sell carcass", "Water fed to chickens", "Vendor feces"), 18, 3),
                    sector_order)

# Chicken = green, Environment = yellow, Human = orange-red
sector_colors <- c(
  "Cage fecal droppings"         = "#009E73",
  "Live chicken rectal swab"     = "#57C4AD",
  "Pre-evisceration carcass"     = "#1B7F5E",
  "Ready-to-sell carcass"        = "#66C2A5",
  "Food fed to chickens"         = "#E6AB02",
  "Vendor hand rinse"            = "#FDD262",
  "Slaughterhouse table surface" = "#B8A000",
  "Water fed to chickens"        = "#A07800",
  "Latrine waste"                = "#D55E00",
  "Vendor feces"                 = "#A6401A"
)

sector_display <- c(
  "Cage fecal droppings"         = "Cage fecal\ndroppings",
  "Live chicken rectal swab"     = "Live chicken\nrectal swab",
  "Pre-evisceration carcass"     = "Pre-\nevisceration\ncarcass",
  "Ready-to-sell carcass"        = "Ready-to-sell\ncarcass",
  "Food fed to chickens"         = "Food fed\nto chickens",
  "Vendor hand rinse"            = "Vendor\nhand rinse",
  "Slaughterhouse table surface" = "Slaughter\ntable surface",
  "Water fed to chickens"        = "Water fed\nto chickens",
  "Latrine waste"                = "Latrine waste",
  "Vendor feces"                 = "Vendor feces"
)

svglite("results/figures/Figure4A_chord_diagram.svg", width = 4, height = 4)
circos.clear()
circos.par(gap.after = gap_vec)
chordDiagram(chord_df,
             order             = sector_order,
             grid.col          = sector_colors,
             transparency      = 0.4,
             directional       = 0,
             annotationTrack   = "grid",
             preAllocateTracks = list(track.height = mm_h(16)))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(mean(get.cell.meta.data("xlim")), get.cell.meta.data("ylim")[1] + mm_y(1),
              sector_display[get.cell.meta.data("sector.index")],
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.731)
}, bg.border = NA)
circos.clear()
dev.off()

# ==============================================================================
# Figure 4B: Reference strains involved in all strain-sharing events
# ==============================================================================
# Reference strains accounting for < 2% of sharing events are grouped as "Other"
pie_data <- byall_raw %>%
  count(ref) %>%
  left_join(clermont %>% select(symlink, MLST, phylotype), by = c("ref" = "symlink")) %>%
  mutate(pct       = n / sum(n),
         ref_label = if_else(pct >= 0.02, ref_label(ref, MLST, phylotype), "Other (each less <2%)")) %>%
  group_by(ref_label) %>%
  summarise(pct = sum(pct), .groups = "drop") %>%
  mutate(label     = paste0(round(pct * 100, 1), "%"),
         ref_label = fct_relevel(ref_label, "Other (each less <2%)", after = Inf))

# Colorblind-safe categorical hues; "Other" is gray
categorical_hues <- c("#2a78d6", "#1baf7a", "#eda100", "#008300", "#4a3aa7", "#e34948", "#e87ba4", "#eb6834", "#8a5a2b")
n_ref <- n_distinct(pie_data$ref_label) - 1
pie_colors <- setNames(c(categorical_hues[seq_len(n_ref)], "#a2a19f"), levels(pie_data$ref_label))

fig4b <- ggplot(pie_data, aes(x = 2, y = pct, fill = ref_label)) +
  geom_col(width = 1) +
  coord_polar("y") +
  xlim(0, 3) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 5, color = "black", fontface = "bold") +
  scale_fill_manual(values = pie_colors) +
  labs(title = expression(atop("Reference strains involved in cefotaxime-resistant", italic(E.coli) ~ "strain sharing")),
       fill = "Reference strain") +
  theme_void() +
  theme(plot.title.position = "plot",
        plot.title          = element_text(hjust = 0, size = 20),
        legend.position     = "right",
        legend.key.size     = unit(6, "pt"),
        legend.spacing.y    = unit(2, "pt"),
        plot.margin         = margin(t = 10, r = 60, b = 10, l = 10),
        legend.text         = element_text(size = 12),
        legend.title        = element_text(size = 16))

ggsave("results/figures/Figure4B_reference_strains.svg", plot = fig4b, width = 6, height = 6)

# ==============================================================================
# Figure 4C: Cross-vendor network for reference strain 2017.04.03CC
# ==============================================================================
vendor_sharing <- load_vendor_sharing()

ggsave("results/figures/Figure4C_network_2017_04_03CC.svg",
       plot = plot_bwvend_network("Esch_coli_2017_04_03CC.fa.gz", vendor_sharing$bwvendor), width = 6, height = 6)
