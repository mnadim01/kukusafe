# ExtDataFigure6_latrine_sharing_timeline.R
# Timeline of cefotaxime-resistant E. coli strain-sharing events involving latrine
# waste, for pairs of samples collected on different days.
#
# Produces: Extended Data Figure 6
#
# Inputs:  data/processed/strainge/strain_sharing_by_latrine.csv
#          data/derived/lab/sample_collection_dates.csv
# Output:  results/figures/ExtDataFigure6_latrine_timeline.svg

source("r-scripts/00_functions.R")

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------
bylatrine_raw <- read_csv("data/processed/strainge/strain_sharing_by_latrine.csv")
dates         <- read_csv("data/derived/lab/sample_collection_dates.csv")

bylatrine <- bylatrine_raw %>%
  add_smpday(dates)

# ------------------------------------------------------------------------------
# Extended Data Figure 6: timeline of latrine sharing events on different days
# ------------------------------------------------------------------------------
# Rows are grouped by the other sample's type, with vertical jitter to separate events
set.seed(1)
bylatrine_diff <- bylatrine %>%
  filter(smpday == "Different day") %>%
  mutate(y_jitter = as.numeric(factor(code_smp2)) * 2 + runif(n(), -0.4, 0.4))

fig_timeline <- ggplot(bylatrine_diff) +
  geom_segment(aes(x = date1, xend = date2, y = y_jitter, yend = y_jitter),
               color = "grey", linewidth = 1) +
  geom_point(aes(x = date1, y = y_jitter, color = "Latrine waste sampling date"), size = 1) +
  geom_point(aes(x = date2, y = y_jitter, color = "Other sample type sampling date"), size = 1) +
  scale_color_manual(values = c("Latrine waste sampling date"     = "darkorange",
                                "Other sample type sampling date" = "darkorchid4")) +
  scale_y_continuous(breaks = seq_along(levels(factor(bylatrine_diff$code_smp2))) * 2,
                     labels = levels(factor(bylatrine_diff$code_smp2))) +
  labs(x = "Date", y = "Sample type", color = NULL) +
  theme_minimal() +
  theme(legend.position  = "bottom",
        # 1 pt black border: ggplot linewidth is in mm-based units (1 lwd = 1/96 in = 0.75 pt)
        legend.background = element_rect(color = "black", fill = NA, linewidth = 1 / (.pt * 0.75)),
        axis.text.y      = element_text(size = 12),
        panel.grid.minor = element_blank())

ggsave("results/figures/ExtDataFigure6_latrine_timeline.svg", plot = fig_timeline, width = 8.5, height = 10)
