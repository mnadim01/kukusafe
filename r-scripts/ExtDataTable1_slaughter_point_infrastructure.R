# ExtDataTable1_slaughter_point_infrastructure.R
# Infrastructural conditions of the 5 sampled slaughter points in Burma Market
# (market conditions spot check).
#
# Produces: Extended Data Table 1
#
# Input:   data/derived/field/slaughterhouse_condition.csv
# Output:  results/tables/ExtDataTable1_slaughter_point_infrastructure.csv

source("r-scripts/00_functions.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

sp <- read_csv("data/derived/field/slaughterhouse_condition.csv")
N <- nrow(sp)

# Surface types are multiple choice (e.g., "1 7"); split into one indicator per type
sp <- sp %>%
  mutate(surf_steel = as.integer(str_detect(slaught_surf, "\\b1\\b")),
         surf_wood  = as.integer(str_detect(slaught_surf, "\\b7\\b")))

ext_table1 <- bind_rows(
  tab_categorical(sp$slaught_type, "Type of slaughter point", c("1" = "Dry", "3" = "Mixed")),
  tab_categorical(sp$slaught_water_source, "Type of water source present",
                  c("3" = "Water tank", "5" = "No water source observed")),
  tab_categorical(sp$slaught_drain_type, "Drainage system present", c("3" = "No drainage observed")),
  tab_categorical(sp$slaught_hw_type, "Handwashing facility present", c("3" = "No handwashing place observed")),
  tab_categorical(sp$slaught_waste_systm, "Solid waste disposal system present",
                  c("1" = "Single/central waste container", "2" = "Multiple waste bins/containers",
                    "3" = "No designated waste container")),
  # Slaughter points could have multiple surface types
  tab_binary(sp, "Type of surfaces present", c("Stainless steel countertops" = "surf_steel",
                                               "Wooden tables"               = "surf_wood")),
  tab_categorical(sp$slaught_surf_cond, "Condition of surfaces present",
                  c("1" = "Without visible cracks or damage", "2" = "Minor-damage (some cracks)")),
  tab_categorical(sp$slaught_pest, "Pests", c("3" = "Flies", "5" = "No pests observed"))
) %>%
  rename(!!paste0("N=", N, " (%)") := `n (%)`)

print(ext_table1, n = Inf)

write_csv(ext_table1, "results/tables/ExtDataTable1_slaughter_point_infrastructure.csv")
