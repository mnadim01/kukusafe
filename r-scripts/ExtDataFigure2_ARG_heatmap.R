# ExtDataFigure2_ARG_heatmap.R
# Resistome patterns among cefotaxime-resistant E. coli by sample type: heatmap of
# the proportion of samples of each type harboring each antibiotic resistance gene (ARG).
#
# Produces: Extended Data Figure 2
#
# Inputs:  data/raw/resfinder/all.kma.res  (KMA alignment to ResFinder, all samples)
#          data/processed/strainge/straingst_ecolidb_summary.csv
# Output:  results/figures/ExtDataFigure2_ARG_heatmap.svg

source("r-scripts/00_functions.R")
library(pheatmap)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# 1. Load data ####
strainGST <- read_csv("data/processed/strainge/straingst_ecolidb_summary.csv")
ecoli <- strainGST %>% filter(!is.na(ecoli_rel_abund))  # samples with E. coli detected

smp_type  <- strainGST %>% select(four_dig_id, smp_type)
smp_count <- smp_type %>% group_by(smp_type) %>% summarize(smp_count = n())

# Keep ARG hits with template coverage > 90% and query identity > 95%
resgene <- read.table("data/raw/resfinder/all.kma.res", header = TRUE) %>%
  filter(Template_Coverage > 90, Query_Identity > 95) %>%
  separate(file_name, into = c(NA, "four_dig_id", NA, NA, NA), sep = "_") %>%
  mutate(four_dig_id = gsub(".kma.res", "", four_dig_id),
         four_dig_id = gsub("ER", "", four_dig_id))

# 2. Reshape ####
# Collapse alleles to genes, replace "-" with "_", and keep
# samples with E. coli detected
first.word <- function(my.string) unlist(strsplit(my.string, "_"))[1]

resgene2 <- resgene %>%
  mutate(Template = as.character(Template)) %>%
  mutate(Gene = sapply(Template, first.word)) %>%
  group_by(four_dig_id, Gene) %>%
  select(-Template) %>%
  summarise_all(sum) %>%
  mutate(Gene = gsub("-", "_", Gene)) %>%
  filter(four_dig_id %in% ecoli$four_dig_id)

# Collapse blaCTX-M genes into groups; collapse other bla genes into gene families
blactam <- resgene2 %>%
  mutate(Gene2 = case_when(Gene %in% c("blaCTX_M_14", "blaCTX_M_121", "blaCTX_M_27") ~ "blaCTX_M_9_group",
                           Gene %in% c("blaCTX_M_15", "blaCTX_M_55")                 ~ "blaCTX_M_1_group",
                           TRUE ~ NA_character_)) %>%
  select(-Gene) %>%
  rename(Gene = Gene2) %>%
  filter(!is.na(Gene))

blactam2 <- resgene2 %>%
  filter(str_detect(Gene, "^bla")) %>%
  mutate(Gene_collapse = stringr::word(Gene, 1, sep = fixed("_"))) %>%
  group_by(four_dig_id, Gene_collapse) %>%
  select(-Gene) %>%
  summarise_all(sum) %>%
  rename(Gene = Gene_collapse) %>%
  filter(!str_detect(Gene, "^blaCTX"))

# Replace bla alleles with bla groups
resgene3 <- resgene2 %>%
  filter(!str_detect(Gene, "^bla")) %>%
  bind_rows(blactam, blactam2) %>%
  distinct()

# 3. Heatmap: proportion of samples of each type harboring each ARG ####
smp_long <- resgene3 %>%
  mutate(four_dig_id = as.numeric(four_dig_id)) %>%
  select(four_dig_id, Gene, Depth) %>%
  left_join(smp_type) %>%
  group_by(Gene, smp_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  left_join(smp_count) %>%
  mutate(freq = n / smp_count * 100) %>%
  mutate(smp_type = factor(smp_type, levels = smp_level_order))

# Remove ARGs present in <= 10% of samples of every sample type
rare_ARG <- smp_long %>% group_by(Gene) %>% summarize(max = max(freq)) %>% filter(max <= 10)

smp_mat <- smp_long %>%
  filter(!Gene %in% rare_ARG$Gene) %>%
  select(Gene, smp_type, freq) %>%
  pivot_wider(names_from = smp_type, values_from = freq, values_fill = 0) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

smp_map <- pheatmap::pheatmap(smp_mat,
                    color = colorRampPalette(c("white", "darkslategrey"))(100),
                    cluster_rows = FALSE,
                    cluster_cols = TRUE,
                    clustering_distance_cols = "euclidean",
                    clustering_method = "complete",
                    show_rownames = TRUE,
                    na_col = "white",
                    fontsize_col = 9,
                    cellheight = 12,
                    silent = TRUE)

ggsave("results/figures/ExtDataFigure2_ARG_heatmap.svg", plot = smp_map, width = 4, height = 30)
