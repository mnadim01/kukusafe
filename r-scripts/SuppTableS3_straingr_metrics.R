# SuppTableS3_straingr_metrics.R
# StrainGR statistics (commonPct, singleAgreePct, gapJaccardSim) for shared strains,
# comparing sample pairs that share an identical reference genome (batch 1,
# Jaccard = 1) with pairs that share a highly similar reference genome
# (batch 2, Jaccard >= 0.7 and < 1). Reports median, minimum, and maximum.
#
# Produces: Supplementary Table S3
#
# Inputs:  data/raw/strainge/straingr/all_sample_comparisons_filtered.tsv       (batch 1)
#          data/raw/strainge/straingr/6May_all_sample_comparisons_filtered.tsv  (batch 2)
# Output:  results/tables/SuppTableS3_straingr_metrics.csv

source("r-scripts/00_functions.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

straingr <- bind_rows(
  read_straingr("data/raw/strainge/straingr/all_sample_comparisons_filtered.tsv") %>% mutate(batch = 1),
  read_straingr("data/raw/strainge/straingr/6May_all_sample_comparisons_filtered.tsv") %>% mutate(batch = 2)
)

metric_comparison <- straingr %>%
  group_by(batch) %>%
  summarize(across(c(commonPct, singleAgreePct, gapJaccardSim),
                   list(median = median, min = min, max = max),
                   .names = "{.fn}_{.col}")) %>%
  pivot_longer(!batch, names_to = "stat", values_to = "value") %>%
  pivot_wider(names_from = batch, values_from = value, names_prefix = "batch_")

print(metric_comparison)

write_csv(metric_comparison, "results/tables/SuppTableS3_straingr_metrics.csv")
