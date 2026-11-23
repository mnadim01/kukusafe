# SuppFigureS2_sharing_by_sample_type_dyad.R
# Proportion of pairwise comparisons in which cefotaxime-resistant E. coli strain
# sharing was detected, by host-category dyad (e.g., Chicken-Human).
#
# Produces: Supplementary Figure S2
#
# Inputs:  data/processed/strainge/straingst_ecolidb_summary.csv
#          data/processed/strainge/strain_sharing_by_sample_all.csv
# Output:  results/figures/SuppFigureS2_sharing_by_sample_type_dyad.svg

source("r-scripts/00_functions.R")
library(scales)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

phylogroups <- read_csv("data/processed/strainge/straingst_ecolidb_summary.csv")
byall_raw   <- read_csv("data/processed/strainge/strain_sharing_by_sample_all.csv")

# Denominator: all possible pairs of samples, by dyad
samples <- phylogroups %>%
  mutate(codecat = factor(dplyr::recode(code, !!!codecat_recode), levels = c("Chicken", "Environment", "Human"))) %>%
  select(four_dig_id, codecat) %>%
  distinct()

pair_df <- combn(samples$four_dig_id, 2) %>% t() %>% as.data.frame()
colnames(pair_df) <- c("s1", "s2")

denom_all <- pair_df %>%
  left_join(samples, by = c("s1" = "four_dig_id")) %>% rename(type1 = codecat) %>%
  left_join(samples, by = c("s2" = "four_dig_id")) %>% rename(type2 = codecat) %>%
  rowwise() %>%
  mutate(dyad = paste(sort(c(type1, type2)), collapse = "-")) %>%
  ungroup() %>%
  count(dyad, name = "n_possible")

# Numerator: unique sample pairs with strain sharing, by dyad
num_all <- byall_raw %>%
  mutate(type1 = dplyr::recode(code_smp1, !!!codecat_recode),
         type2 = dplyr::recode(code_smp2, !!!codecat_recode),
         s1    = pmin(sample1, sample2),
         s2    = pmax(sample1, sample2)) %>%
  rowwise() %>%
  mutate(dyad = paste(sort(c(type1, type2)), collapse = "-")) %>%
  ungroup() %>%
  distinct(s1, s2, dyad) %>%
  count(dyad, name = "n_sharing")

rates_all <- denom_all %>%
  left_join(num_all, by = "dyad") %>%
  mutate(n_sharing = replace_na(n_sharing, 0),
         rate      = n_sharing / n_possible)

fig_s2 <- ggplot(rates_all, aes(x = reorder(dyad, -rate), y = rate)) +
  geom_col(fill = "#9B59F5", alpha = 0.8) +
  geom_text(aes(label = paste0("n = ", n_possible)), vjust = -0.3, size = 3) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(x = "Sample Type Dyads", y = "Proportion (%)", title = "Proportion of dyads with strain sharing") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/figures/SuppFigureS2_sharing_by_sample_type_dyad.svg", plot = fig_s2, width = 5, height = 4)
