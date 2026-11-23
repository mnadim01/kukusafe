# Figure3_strains_by_sample_type.R
# Number of samples with at least one strain, unique strains, and total strains of
# cefotaxime-resistant E. coli detected in each sample type.
#
# Produces: Figure 3
#
# Inputs:  data/processed/strainge/straingst_ecolidb_summary.csv
#          data/raw/strainge/straingst/ecoli_straingst_strains.tsv
#          data/derived/lab/lab_data.csv
# Output:  results/figures/Figure3_strains_by_sample_type.svg

source("r-scripts/00_functions.R")
library(ggExtra)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

lab            <- read_csv("data/derived/lab/lab_data.csv")
phylogroups    <- read_csv("data/processed/strainge/straingst_ecolidb_summary.csv")
refstrains_raw <- read_tsv("data/raw/strainge/straingst/ecoli_straingst_strains.tsv")

newphylogroups <- phylogroups

# Reference strains detected in each sample
refstrains <- refstrains_raw %>%
  mutate(samp = as.numeric(gsub("^A_|ER.*$", "", sample))) %>%
  left_join(lab %>% select(four_dig_id, smp_type), by = c("samp" = "four_dig_id")) %>%
  mutate(type = unname(smp_code_map[smp_type]))

# Total strains (summed across samples) by sample type
samples_summary <- newphylogroups %>%
  group_by(code) %>%
  summarise(n_samples = n(), n_strains = sum(strain_count, na.rm = TRUE), .groups = "drop")

# Unique strains by sample type
unique_strains_by_code <- refstrains %>%
  group_by(type) %>%
  summarise(n_unique_strains = n_distinct(strain), .groups = "drop") %>%
  rename(code = type)

# Samples with at least one strain
num_samples_with_strains <- refstrains %>%
  group_by(sample) %>% slice(1) %>%
  group_by(type) %>%
  summarise(n_with_strain = n()) %>%
  rename(code = type)

samples_long <- samples_summary %>%
  pivot_longer(cols = c(n_samples, n_strains), names_to = "metric", values_to = "count") %>%
  bind_rows(unique_strains_by_code %>%
              pivot_longer(cols = n_unique_strains, names_to = "metric", values_to = "count")) %>%
  mutate(codecat = dplyr::recode(code, !!!codecat_recode)) %>%
  left_join(num_samples_with_strains, by = "code") %>%
  mutate(count  = if_else(metric == "n_samples", as.numeric(n_with_strain), count),
         metric = factor(metric, levels = c("n_strains", "n_unique_strains", "n_samples"))) %>%
  select(-n_with_strain)

# y-axis labels: wrapped sample type names + number of samples analyzed
smp_wrap_map <- c(
  "Cage fecal droppings"         = "Cage fecal\ndroppings",
  "Live chicken rectal swab"     = "Live chicken\nrectal swab",
  "Pre-evisceration carcass"     = "Pre-evisceration\ncarcass",
  "Ready-to-sell carcass"        = "Ready-to-sell\ncarcass",
  "Vendor hand rinse"            = "Vendor\nhand rinse",
  "Food fed to chickens"         = "Food fed\nto chickens",
  "Water fed to chickens"        = "Water fed\nto chickens",
  "Slaughterhouse table surface" = "Slaughterhouse\ntable surface",
  "Latrine waste"                = "Latrine waste",
  "Vendor feces"                 = "Vendor feces"
)

samptype_n <- lab %>%
  filter(!is.na(smp_type)) %>%
  count(smp_type) %>%
  mutate(code = dplyr::recode(smp_type, !!!smp_code_map))

samptype_labels <- setNames(paste0(smp_wrap_map[samptype_n$code], " (n=", samptype_n$n, ")"), samptype_n$code)

fig3 <- ggplot(samples_long, aes(y = factor(code, levels = rev(unname(smp_code_map))), x = count, fill = metric)) +
  facet_grid(codecat ~ ., scales = "free_y", space = "free_y", switch = "y") +
  geom_col(position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_text(aes(label = count), position = position_dodge(width = 0.8), hjust = -0.2, vjust = -0.02, size = 2.7) +
  scale_fill_manual(
    values = c(n_strains = "darkblue", n_unique_strains = "cornflowerblue", n_samples = "lightblue"),
    labels = c(n_strains        = "Total number of strains",
               n_unique_strains = "Number of unique strains",
               n_samples        = "Number of samples with at least 1 strain")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_discrete(labels = samptype_labels) +
  guides(fill = guide_legend(reverse = TRUE)) +
  labs(x = "Count", y = NULL, fill = NULL) +
  theme_minimal() +
  theme(text              = element_text(size = 11),
        strip.placement   = "outside",
        strip.text.y.left = element_text(angle = 0),
        panel.spacing.y   = unit(0.5, "lines"),
        plot.margin       = margin(t = 10, r = 10, b = 10, l = 20)) +
  removeGrid()

ggsave("results/figures/Figure3_strains_by_sample_type.svg", plot = fig3, width = 8, height = 6)
