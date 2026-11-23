# ExtDataFigure1_SuppFigureS1_ecoli_load.R
# Cefotaxime-resistant E. coli load and relative abundance by sample type and sampling day.
#
# Produces: Extended Data Figure 1   (log10 CFU on chicken samples, day 1 vs day 2)
#           Supplementary Figure S1A (log10 CFU by sample type)
#           Supplementary Figure S1B (relative abundance by sample type)
#
# Inputs:  data/derived/lab/lab_data.csv
#          data/reference/sample_type_codes.csv
# Outputs: results/figures/ExtDataFigure1_log10CFU_by_day.svg
#          results/figures/SuppFigureS1A_log10CFU_by_sample_type.svg
#          results/figures/SuppFigureS1B_relative_abundance_by_sample_type.svg

source("r-scripts/00_functions.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# 1. Import ####
lab  <- read.csv("data/derived/lab/lab_data.csv")
code <- read_csv("data/reference/sample_type_codes.csv")

# 2. Prepare ####
# Relative abundance = cefotaxime-resistant E. coli / total E. coli. Samples 6705 and
# 3304 exceeded 100% because one measurement was above the ULOQ; set to 100%.
# Sampling day applies only to the slaughter point samples (CX, CXS, FDCX, SSx).
lab <- lab %>%
  mutate(mc_yes_relab = (10^mc_yes_quant / 10^mc_no_quant) * 100) %>%
  rename(log10mc_yes = mc_yes_quant) %>%
  mutate(mc_yes_relab = case_when(four_dig_id %in% c("6705", "3304") ~ 100, TRUE ~ mc_yes_relab),
         day = as.factor(case_when(smp_type %in% c("CX", "CXS", "FDCX", "SSx") ~ day, TRUE ~ 1)),
         smp_type = as.factor(smp_type)) %>%
  left_join(code)

level_order  <- smp_level_order
level_order2 <- c("Vendor hand rinse", "Food fed to chickens", "Water fed to chickens", "Cage fecal droppings",
                  "Live chicken rectal swab", "Pre-evisceration carcass", "Slaughterpoint table surface",
                  "Ready-to-sell carcass", "Vendor feces", "Latrine waste")
labels <- c("Vendor hand rinse (n=8)", "Food fed to chickens (n=8)", "Water fed to chickens (n=7)",
            "Cage fecal droppings (n=6)", "Live chicken rectal swab (n=28)", "Pre-evisceration carcass (n=28)",
            "Slaughterpoint table surface (n=14)", "Ready-to-sell carcass (n=32)", "Vendor feces (n=5)")
fill_colors <- c("#84969f", "#9BD6ED", "#5e676d", "#86CDE9", "#38383a", "#FFFFFF", "#a2a19f", "#94B7C4", "#9bacb2")

# 3. Extended Data Figure 1: chicken samples by sampling day ####
subset1 <- lab %>% filter(smp_type %in% c("CX", "CXS", "FDCX"))

set.seed(1)  # reproducible point jitter
byday <- ggplot(subset1, aes(x = factor(code, level = level_order2), y = log10mc_yes, fill = day)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.8, size = 1.25, color = "darkslategrey") +
  scale_x_discrete(name = "Sample type",
                   labels = function(x) str_wrap(c("Live chicken rectal swab (n=28)", "Pre-evisceration carcass (n=28)",
                                                   "Ready-to-sell carcass (n=32)"), width = 15)) +
  scale_fill_manual(values = c("white", "skyblue"), name = "Day", labels = c("Day 1", "Day 2")) +
  theme_classic() +
  ylab(str_wrap("Cefotaxime-resistant E. coli (log10 CFU)", width = 25)) +
  theme(text = element_text(size = 16), axis.text = element_text(color = "black"))

ggsave("results/figures/ExtDataFigure1_log10CFU_by_day.svg", plot = byday, width = 9, height = 4)

# 4. Supplementary Figure S1: all market samples (latrine waste not quantified) ####
subset2 <- lab %>% filter(smp_type != "LW", !is.na(mc_yes_relab))

# S1A: log10 CFU
set.seed(1)
smptype <- ggplot(subset2, aes(x = factor(smp_type, level = level_order), y = log10mc_yes,
                               fill = factor(smp_type, level = level_order))) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
             alpha = 0.8, size = 1.25, color = "darkslategrey") +
  scale_x_discrete(name = "Sample type", labels = function(x) str_wrap(labels, width = 5)) +
  scale_fill_manual(values = fill_colors) +
  theme_classic() +
  ylab(str_wrap("Cefotaxime-resistant E. coli (log10 CFU)", width = 25)) +
  theme(text = element_text(size = 16), legend.position = "none", axis.text = element_text(color = "black"))

ggsave("results/figures/SuppFigureS1A_log10CFU_by_sample_type.svg", plot = smptype, width = 13, height = 5)

# S1B: relative abundance (mean ± SE)
rel_abund <- ggplot(subset2, aes(x = factor(smp_type, level = level_order), y = mc_yes_relab,
                                 fill = factor(smp_type, level = level_order))) +
  stat_summary(fun = "mean", geom = "bar", position = "dodge", color = "black") +
  stat_summary(fun.data = "mean_se", geom = "errorbar", position = position_dodge(width = 0.9), width = 0.2) +
  scale_x_discrete(name = "Sample type", labels = function(x) str_wrap(labels, width = 5)) +
  scale_fill_manual(values = fill_colors) +
  theme_classic() +
  ylab(str_wrap("Relative abundance cefotaxime-resistant E. coli (%)", width = 25)) +
  theme(text = element_text(size = 16), legend.position = "none", axis.text = element_text(color = "black"))

ggsave("results/figures/SuppFigureS1B_relative_abundance_by_sample_type.svg", plot = rel_abund, width = 13, height = 5)
