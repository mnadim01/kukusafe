# ExtDataTable3_slaughter_processing_practices.R
# Food safety and hygiene practices during the slaughter and processing of the 32
# observed chickens at Burma Market (slaughter point observations).
#
# Produces: Extended Data Table 3
#
# Input:   data/derived/field/slaughterhouse_collapsed.csv
# Output:  results/tables/ExtDataTable3_slaughter_processing_practices.csv

source("r-scripts/00_functions.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

sh <- read_csv("data/derived/field/slaughterhouse_collapsed.csv")
N <- nrow(sh)

# Cleaning materials and disposal locations keep their original codes;
# NA = surface/utensil not cleaned
no_clean <- "No cleaning observed"

ext_table3 <- bind_rows(
  tab_binary(sh, "PPE worn during slaughtering",
             c("Gloves" = "slghtr_glov", "Face mask" = "slghtr_fmask", "Hairnet" = "slghtr_hairnet",
               "Apron" = "slghtr_apron", "Boots" = "slghtr_boots")),
  tab_categorical(sh$slghtr_surf_clnm, "Material used to clean slaughtering surface",
                  c("5" = "Feathers"), missing_label = no_clean),
  tab_binary(sh, "Scalding occurred", c("Scalding occurred" = "scald")),
  tab_categorical(sh$def_surf_clnm, "Material used to clean defeathering surface",
                  c("5" = "Feathers"), missing_label = no_clean),
  tab_binary(sh, "Handwashing practices during slaughter",
             c("Before slaughter" = "slghtr_hw", "Before scalding" = "scald_hw", "Before defeathering" = "def_hw")),
  tab_binary(sh, "PPE worn during processing",
             c("Gloves" = "proc_glov", "Face mask" = "proc_fmask", "Hairnet" = "proc_hairnet",
               "Apron" = "proc_apron", "Boots" = "proc_boots")),
  tab_categorical(sh$evisc_surf_cln_mtrl, "Material used to clean evisceration surface",
                  c("5" = "Feathers"), missing_label = no_clean),
  tab_categorical(sh$evisc_utens_clnm, "Material used to clean evisceration utensil",
                  c("5" = "Feathers", "2" = "Water"), missing_label = no_clean),
  tab_binary(sh, "Processing practices",
             c("Feathers inserted into body cavity for cleaning" = "off_feath_cln",
               "Evisceration finished in water"                  = "evisc_water",
               "Offal rinsed"                                    = "offrins",
               "Carcass rinsed"                                  = "carcrins",
               # Only recorded when the carcass was rinsed (n = 11)
               "Carcass rinsed in clean water"                   = "carcrins_water_chng")),
  tab_categorical(sh$proc_swst_disp, "Solid waste disposal",
                  c("5" = "Slaughter point floor", "2" = "Uncovered waste bin in the slaughter point")),
  tab_binary(sh, "Handwashing practices during processing",
             c("Before evisceration" = "evisc_hw", "Before carcass rinsing" = "carcrins_hw",
               "After solid waste disposal" = "proc_swst_hw")),
  tab_categorical(sh$proc_pest, "Pests observed", c("3" = "Flies", "5" = "None"))
) %>%
  rename(!!paste0("N=", N, " (%)") := `n (%)`)

print(ext_table3, n = Inf)

write_csv(ext_table3, "results/tables/ExtDataTable3_slaughter_processing_practices.csv")
