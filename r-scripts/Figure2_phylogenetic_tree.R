# Figure2_phylogenetic_tree.R
# Phylogenetic tree of all E. coli reference strains detected at Burma Market,
# with heatmaps of reference strain detection frequency by vendor and by sample type.
#
# Produces: Figure 2 (tree + vendor heatmap + sample type heatmap)
#
# Inputs:  data/derived/phylotree/kukusafe_ecoli.nwk  (JolyTree v2.1 output)
#          data/raw/strainge/straingst/ecoli_straingst_strains.tsv
#          data/reference/strainge_refdb_mlst_phylotype.csv
#          data/derived/lab/lab_data.csv
# Outputs: results/figures/Figure2_phylotree.svg
#          results/figures/Figure2_vendor_heatmap.svg
#          results/figures/Figure2_sample_type_heatmap.svg

source("r-scripts/00_functions.R")
library(treeio)
library(ggtree)
library(ggh4x)
library(pheatmap)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# 1. Load tree ####
tree <- read.newick("data/derived/phylotree/kukusafe_ecoli.nwk")

# 2. Build tip metadata ####
strains <- read_straingst_strains()

clermont <- read_csv("data/reference/strainge_refdb_mlst_phylotype.csv") %>%
  rename(strain = symlink) %>%
  select(strain, phylotype, MLST, original_file) %>%
  mutate(label = gsub(".fna.gz", "", original_file))

lab <- read.csv("data/derived/lab/lab_data.csv")

metadata <- lab %>%
  select(four_dig_id, smp_type, vendor) %>%
  left_join(strains) %>%
  left_join(clermont) %>%
  filter(!is.na(strain)) %>%
  select(strain, smp_type, vendor, phylotype, MLST, label) %>%
  mutate(phylo_cat = case_when(phylotype %in% c("B2", "D") ~ "B2/D",
                               phylotype == "A"            ~ "A",
                               phylotype == "B1"           ~ "B1",
                               TRUE                        ~ "Other"))

# 3. Detection frequency (%) of each reference strain by vendor and by sample type ####
vendor_long <- metadata %>%
  group_by(label, vendor) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(vendor) %>%
  mutate(freq = n / sum(n) * 100) %>%
  ungroup() %>%
  filter(vendor != "LW") %>%
  mutate(vendor = factor(vendor, levels = c("2", "4", "5", "6", "8", "9", "11", "12")))

smp_long <- metadata %>%
  group_by(label, smp_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(smp_type) %>%
  mutate(freq = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(smp_type = factor(smp_type, levels = smp_level_order))

vendor_levels <- levels(vendor_long$vendor)
smp_levels    <- levels(smp_long$smp_type)

# 4. Tree annotated with phylogroup ####
phylo_colors <- c("A" = "#E41A1C", "B1" = "#377EB8", "B2/D" = "#4DAF4A", "Other" = "#984EA3")

p <- ggtree(tree, layout = "roundrect", branch.length = "none") %<+%
  (metadata %>% select(label, phylo_cat, phylotype, MLST)) +
  geom_tippoint(aes(color = phylo_cat), size = 0.8, na.rm = TRUE) +
  scale_color_manual(values = phylo_colors, name = "Phylogroup",
                     guide = guide_legend(override.aes = list(size = 2.4))) +
  geom_treescale(x = 0, y = (Ntip(tree) - 4), fontsize = 4, linesize = 1, offset = 1)

fig_tree <- p +
  geom_facet(panel = "Vendor", data = vendor_long, geom = geom_tile,
             mapping = aes(x = as.numeric(vendor), fill = freq)) +
  geom_facet(panel = "Sample type", data = smp_long, geom = geom_tile,
             mapping = aes(x = as.numeric(smp_type), fill = freq)) +
  scale_fill_gradient(low = "gray90", high = "darkslategrey",
                      name = "Detection\nfrequency (%)", na.value = "white") +
  facetted_pos_scales(x = list(
    NULL,
    scale_x_continuous(breaks = seq_along(vendor_levels), labels = vendor_levels, position = "top"),
    scale_x_continuous(breaks = seq_along(smp_levels),    labels = smp_levels,    position = "top")
  )) +
  theme(strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0, vjust = -1),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2))

ggsave("results/figures/Figure2_phylotree.svg", plot = fig_tree, width = 10, height = 6)

# 5. Heatmaps with hierarchically clustered columns, rows in tree tip order ####
heatmap_colors <- colorRampPalette(c("white", "darkslategrey"))(100)
tip_order <- get_taxa_name(p)

make_heatmap <- function(long_df, col_var) {
  mat <- long_df %>%
    select(label, all_of(col_var), freq) %>%
    pivot_wider(names_from = all_of(col_var), values_from = freq, values_fill = 0) %>%
    column_to_rownames("label") %>%
    as.matrix()
  mat <- mat[tip_order[tip_order %in% rownames(mat)], , drop = FALSE]

  pheatmap::pheatmap(mat,
           color = heatmap_colors,
           cluster_rows = FALSE,
           cluster_cols = TRUE,
           clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           show_rownames = FALSE,
           na_col = "white",
           fontsize_col = 9,
           silent = TRUE)
}

ggsave("results/figures/Figure2_vendor_heatmap.svg",      plot = make_heatmap(vendor_long, "vendor"),   width = 3, height = 6)
ggsave("results/figures/Figure2_sample_type_heatmap.svg", plot = make_heatmap(smp_long,    "smp_type"), width = 3, height = 6)
