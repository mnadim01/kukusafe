# SuppTableS2_cage_conditions.R
# Chicken cage conditions among the 8 enrolled poultry vendors at Burma Market
# (vendor stall observations).
#
# Produces: Supplementary Table S2
#
# Feed container questions were only asked when chicken feed was available, and
# water container questions only when water was available; percentages for those
# rows use the number of cages with feed or water as the denominator.
#
# Input:   data/derived/field/vendor_obsv_collapsed.csv
# Output:  results/tables/SuppTableS2_cage_conditions.csv

source("r-scripts/00_functions.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cages <- read_csv("data/derived/field/vendor_obsv_collapsed.csv")
N       <- nrow(cages)
n_feed  <- sum(cages$feed_av == 1, na.rm = TRUE)
n_water <- sum(cages$water_av == 1, na.rm = TRUE)

supp_table2 <- bind_rows(
  tab_binary(cages, "Cage", c("Overcrowded"                         = "cg_ovrcrowd",
                              "Excess fecal waste"                  = "cg_exc_wst",
                              "Placed directly on the ground"       = "cg_on_ground",
                              "Stacked"                             = "cg_stacked",
                              "Objects on top (e.g., utensils)"     = "cg_objs",
                              "Carcasses on top"                    = "cg_carc",
                              "Live chickens on top"                = "cg_live_chik",
                              "Eggs on top"                         = "cg_eggs",
                              "Adjacent to other vendors"           = "cg_adj",
                              "Located near produce or food stalls" = "cg_food")),
  tab_binary(cages, "Feed", c("Chicken feed available" = "feed_av")),
  tab_binary(cages, paste0("Feed (among cages with feed, n=", n_feed, ")"),
             c("Chicken feed container covered" = "feed_cont_cover",
               "Chicken feed container clean"   = "feed_cont_cln",
               "Flies around feed container"    = "feed_cont_fly"),
             total = n_feed),
  tab_binary(cages, "Water", c("Water available" = "water_av")),
  tab_binary(cages, paste0("Water (among cages with water, n=", n_water, ")"),
             c("Water visibly murky"     = "water_murky",
               "Water container clean"   = "water_cont_cln",
               "Water container covered" = "water_cont_cov"),
             total = n_water)
) %>%
  rename(!!paste0("N=", N, " (%)") := `n (%)`)

print(supp_table2, n = Inf)

write_csv(supp_table2, "results/tables/SuppTableS2_cage_conditions.csv")
