# Table1_vendor_characteristics.R
# Socio-demographic characteristics and chicken handling practices among the 8
# poultry vendors at Burma Market (vendor survey, day 1 and day 2).
#
# Produces: Table 1
#
# Overnight storage location (chik_strg_d1) and fecal disposal method (wste_disp_d1)
# were multiple-choice questions; answers such as "2 1" mean the vendor selected
# both options. Each option is counted separately, so percentages can sum to more
# than 100%.
#
# Input:   data/derived/field/survey_collapsed.csv
# Output:  results/tables/Table1_vendor_characteristics.csv

source("r-scripts/00_functions.R")

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

survey <- read_csv("data/derived/field/survey_collapsed.csv", col_types = cols(chik_strg_d1 = "c", wste_disp_d1 = "c"))
N <- nrow(survey)

# 1 if a multiple-choice answer (e.g., "2 1") includes the given option code
has_option <- function(x, code) as.integer(str_detect(x, paste0("\\b", code, "\\b")) %in% TRUE)

survey <- survey %>%
  mutate(
    purch_group = cut(qty_purch_tday_d1, breaks = c(-Inf, 99, 299, Inf), labels = c("0-99", "100-299", "300+")),
    # Commercial feed given at any point on day 1
    commercial_feed = as.integer(feed_type_d1 %in% 1 | feed_type_later_d1 %in% 1),
    # Multiple-choice answers split into one indicator per option
    strg_market   = has_option(chik_strg_d1, 2),
    strg_stall    = has_option(chik_strg_d1, 1),
    strg_coop     = has_option(chik_strg_d1, 4),
    disp_dumpster = has_option(wste_disp_d1, 1),
    disp_compost  = has_option(wste_disp_d1, 2),
    disp_drainage = has_option(wste_disp_d1, 3),
    # Poultry health reported on either survey day
    any_sick       = as.integer(sick_chik_tday_d1 %in% 1 | sick_chik_last_d2 %in% 1),
    any_dead       = as.integer(dead_chik_tday_d1 %in% 1 | dead_chik_last_d2 %in% 1),
    gave_away_sick = as.integer(dead_chik_outc_d1 %in% 2),
    any_medication = as.integer(med_tday_d1 %in% 1 | med_last_d2 %in% 1)
  )

table1 <- bind_rows(
  tab_categorical(survey$sex_d1, "Sex", c("1" = "Female", "2" = "Male")),
  tab_categorical(survey$age_group_d1, "Age",
                  c("18-29" = "18-29 years old", "30-39" = "30-39 years old",
                    "40-49" = "40-49 years old", "50-59" = "50-59 years old")),
  tab_categorical(survey$edu_level_d1, "Education level",
                  c("3" = "Completed primary school (8 years)", "4" = "Some secondary school (8-11 years)",
                    "5" = "Completed secondary school (12 years)", "6" = "Vocational or higher education")),
  tab_categorical(survey$prim_income_d1, "Primary source of income",
                  c("1" = "Selling chickens in the market", "2" = "Other")),
  tab_categorical(survey$yrs_in_mkt_d1, "Years working in the market",
                  c("1" = "Less than 1 year", "2" = "1 to 2 years", "3" = "3 to 5 years", "4" = "More than 5 years")),
  tab_categorical(survey$fs_training_d1, "Food safety training", c("1" = "Received", "0" = "Never received")),
  tab_categorical(survey$clean_freq_d1, "Cage cleaning frequency",
                  c("2" = "Every day", "4" = "Twice a week", "5" = "Once a week", "6" = "Less than once a week")),
  tab_categorical(survey$cln_meth_d1, "Cage cleaning methods",
                  c("1" = "Sweep, brush, scrape (dry clean)", "3" = "Rinse with water only",
                    "5" = "None (cages never cleaned)")),
  tab_binary(survey, "Preferred disposal method for chicken fecal droppings (could select more than one)",
             c("Designated waste/dumpster area"       = "disp_dumpster",
               "Composting or fertilizer application" = "disp_compost",
               "Nearby drainage or sewer system"      = "disp_drainage")),
  tab_categorical(survey$commercial_feed, "Food fed to chickens during the study",
                  c("1" = "Commercial feed", "0" = "Food scraps/other")),
  tab_binary(survey, "Location where unsold chickens are stored overnight (could select more than one)",
             c("Communal storage area in the market" = "strg_market",
               "Vending stall"                       = "strg_stall",
               "Vendor's chicken coop in the market" = "strg_coop")),
  tab_binary(survey, "Poultry health during study",
             c("Vendors reported sick chicken(s) during the study"  = "any_sick",
               "Vendors reported chicken(s) dying during the study" = "any_dead",
               "Vendors reported giving away diseased chickens"     = "gave_away_sick")),
  tab_categorical(survey$any_medication, "Frequency of medication administration during the study",
                  c("0" = "Never", "1" = "At least once")),
  tab_categorical(survey$purch_group, "Number of chickens vendors purchased for resale",
                  c("0-99" = "0-99", "100-299" = "100-299", "300+" = "300+")),
  tab_categorical(survey$purch_freq_d1, "How frequently vendors purchase chickens",
                  c("1" = "Every day", "2" = "Every other day", "3" = "Twice a week")),
  tab_categorical(survey$src_d1, "Primary source of chickens",
                  c("1" = "Commercial poultry company", "2" = "A farmer (small scale/backyard)",
                    "3" = "Middleman or informal seller"))
) %>%
  rename(!!paste0("Proportion (N=", N, ")") := `n (%)`)

print(table1, n = Inf)

write_csv(table1, "results/tables/Table1_vendor_characteristics.csv")
