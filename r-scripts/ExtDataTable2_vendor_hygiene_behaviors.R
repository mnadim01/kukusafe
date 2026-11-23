# ExtDataTable2_vendor_hygiene_behaviors.R
# Food safety and hygiene behaviors among chicken vendors during 24 customer
# interactions (8 vendors x 3 interactions) in Burma Market.
#
# PPE was recorded once per vendor at the start of the observation period and is
# applied to each of that vendor's three customer interactions.
#
# Produces: Extended Data Table 2
#
# Input:   data/derived/field/vendor_obsv_collapsed.csv
# Output:  results/tables/ExtDataTable2_vendor_hygiene_behaviors.csv

source("r-scripts/00_functions.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

vendor_obsv <- read_csv("data/derived/field/vendor_obsv_collapsed.csv")

# One row per customer interaction
interactions <- vendor_obsv %>%
  select(vd_id, starts_with("vend_ppe"), matches("^hw_.*_rep[123]$")) %>%
  pivot_longer(matches("_rep[123]$"), names_to = c(".value", "rep"), names_pattern = "(.*)_rep([123])")

N <- nrow(interactions)

ext_table2 <- bind_rows(
  tab_binary(interactions, "Use of personal protective equipment",
             c("Gloves worn"    = "vend_ppe-glove",
               "Apron worn"     = "vend_ppe-apron",
               "Face mask worn" = "vend_ppe-fmask",
               "Boots worn"     = "vend_ppe-boot",
               "Hairnet worn"   = "vend_ppe-hnet")),
  tab_binary(interactions, "Hand hygiene behaviors",
             c("Handwashing before handling chicken"    = "hw_chik_touch",
               "Handwashing after customer interaction" = "hw_af_pymt"))
) %>%
  rename(!!paste0("N=", N, " (%)") := `n (%)`)

print(ext_table2, n = Inf)

write_csv(ext_table2, "results/tables/ExtDataTable2_vendor_hygiene_behaviors.csv")
