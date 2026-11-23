# ExtDataFigure4_same_chicken_sharing.R
# Cefotaxime-resistant E. coli strain sharing between samples taken from the same
# individual chickens (live chicken rectal swab, pre-evisceration carcass, and
# ready-to-sell carcass). Each row is one chicken; each line is a strain shared
# between two of that chicken's samples, colored by reference strain.
#
# Produces: Extended Data Figure 4
#
# Inputs:  data/processed/strainge/strain_sharing_within_vendor.csv
#          data/processed/strainge/strain_sharing_between_vendor.csv
#          data/processed/strainge/straingst_ecolidb_summary.csv
#          data/derived/lab/sample_collection_dates.csv
#          data/derived/lab/lab_data.csv
#          data/reference/strainge_refdb_mlst_phylotype.csv
# Output:  results/figures/ExtDataFigure4_same_chicken_sharing.svg

source("r-scripts/00_functions.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

wivendor <- load_vendor_sharing()$wivendor
lab      <- read_csv("data/derived/lab/lab_data.csv")
clermont <- read_csv("data/reference/strainge_refdb_mlst_phylotype.csv")

# Chicken ID (e.g., V5D2CX1 = vendor 5, day 2, chicken 1) for each sample
chk_lookup <- lab %>% select(four_dig_id, chk_link)

# Reference strain labels, e.g., "Phylo A\ST 206(2017_04_03CC)"
ref_labels <- clermont %>%
  transmute(ref = symlink,
            refabrev = sub("\\.fa\\.gz$", "", sub("^Esch_coli_", "", symlink)),
            mlst_label = if_else(MLST == "-",
                                 paste0("Phylo ", phylotype, " (", refabrev, ")"),
                                 paste0("Phylo ", phylotype, "\\ST ", MLST, "(", refabrev, ")"))) %>%
  select(ref, mlst_label)

# Within-vendor sharing events between two samples from the same chicken
samechick <- wivendor %>%
  left_join(chk_lookup %>% rename(chk1 = chk_link), by = c("sample1" = "four_dig_id")) %>%
  left_join(chk_lookup %>% rename(chk2 = chk_link), by = c("sample2" = "four_dig_id")) %>%
  filter(chk1 == chk2) %>%
  left_join(ref_labels, by = "ref")

# One row per chicken, ordered by vendor
chicken_rows <- samechick %>%
  distinct(vendor_smp1, chk1) %>%
  arrange(vendor_smp1, chk1) %>%
  mutate(y = row_number())

sample_positions <- tibble(sample_type = c("Live chicken rectal swab", "Pre-evisceration carcass", "Ready-to-sell carcass"),
                           x = 1:3)

# Edges between sample types; strains shared by the same chicken are offset vertically
edges <- samechick %>%
  distinct(vendor_smp1, chk1, code_smp1, code_smp2, mlst_label) %>%
  left_join(chicken_rows, by = c("vendor_smp1", "chk1")) %>%
  left_join(sample_positions, by = c("code_smp1" = "sample_type")) %>% rename(x1 = x) %>%
  left_join(sample_positions, by = c("code_smp2" = "sample_type")) %>% rename(x2 = x) %>%
  group_by(vendor_smp1, chk1) %>%
  mutate(y_offset = y + (row_number() - (n() + 1) / 2) * 0.1) %>%
  ungroup()

nodes <- bind_rows(
  edges %>% select(vendor_smp1, chk1, y, sample_type = code_smp1, x = x1),
  edges %>% select(vendor_smp1, chk1, y, sample_type = code_smp2, x = x2)
) %>%
  distinct()

fig <- ggplot() +
  geom_segment(data = edges, aes(x = x1, xend = x2, y = y_offset, yend = y_offset, color = mlst_label),
               linewidth = 1, lineend = "round") +
  geom_point(data = nodes, aes(x = x, y = y, shape = sample_type), size = 4, color = "black") +
  scale_shape_manual(values = c("Live chicken rectal swab" = 16, "Pre-evisceration carcass" = 18,
                                "Ready-to-sell carcass" = 15)) +
  scale_y_reverse(breaks = chicken_rows$y, labels = chicken_rows$chk1) +
  scale_x_continuous(breaks = 1:3,
                     labels = c("Live chicken\nrectal swab", "Pre-evisceration\ncarcass", "Ready-to-sell\ncarcass")) +
  labs(x = NULL, y = NULL, shape = "Sample type", color = "Strain reference") +
  theme_classic() +
  theme(axis.ticks = element_blank(), axis.line = element_blank(), legend.position = "right")

ggsave("results/figures/ExtDataFigure4_same_chicken_sharing.svg", plot = fig, width = 8, height = 6)
