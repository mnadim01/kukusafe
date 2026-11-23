# 01_data_preprocessing.R
# Builds every processed dataset used by the figure and table scripts.
# Run from the repository root. All outputs are written to data/processed/.
#
# Inputs include the cleaned field datasets in data/derived/field/ (vendor survey,
# vendor and cage observations, slaughter point observations and infrastructure).
#
#  1. StrainGST sample summary         -> strainge/straingst_ecolidb_summary.csv
#  2. StrainGR strain-sharing datasets -> strainge/strain_sharing_*.csv
#  3. Vendor-level covariates          -> regression/vendor_variables.csv
#  4. Slaughter-point-level covariates -> regression/slaughterhouse_variables.csv
#  5. All sample pairs                 -> regression/all_strain_sharing_pairs.csv, market.csv, vendor.csv, cross_vendor.csv
#  6. Within-vendor regression dataset -> regression/vendor_log_updated.csv
#  7. Cross-vendor regression dataset  -> regression/cross_vendor_log.csv
#
# Each section runs in its own local() environment so intermediate objects do not carry over.

source("r-scripts/00_functions.R")
for (d in c("strainge", "regression")) {
  dir.create(file.path("data/processed", d), recursive = TRUE, showWarnings = FALSE)
}

# ==============================================================================
# 1. StrainGST sample summary: E. coli relative abundance, strain count, phylogroups
# ==============================================================================
local({
  summary <- read_tsv("data/raw/strainge/straingst/ecoli_straingst_stats.tsv") %>%
    parse_straingst_id() %>%
    select(four_dig_id, ecoli_rel_abund = `pan%`)

  strains  <- read_straingst_strains()
  clermont <- read_csv("data/reference/strainge_refdb_mlst_phylotype.csv") %>%
    rename(strain = symlink) %>%
    select(strain, phylotype, MLST)

  data <- read.csv("data/derived/lab/lab_data.csv")
  code <- read_csv("data/reference/sample_type_codes.csv")

  merge <- data %>%
    select(four_dig_id, smp_type, vendor, day) %>%
    left_join(summary) %>%
    left_join(strains) %>%
    left_join(clermont)

  strain_summary <- merge %>%
    group_by(four_dig_id) %>%
    filter(!is.na(strain)) %>%
    summarize(strain_count = n())

  phylotype_summary <- merge %>%
    group_by(four_dig_id, phylotype) %>%
    summarize(count = n()) %>%
    pivot_wider(names_from = phylotype, values_from = count) %>%
    select(-"NA")

  merge_final <- data %>%
    select(four_dig_id, smp_type, vendor, day) %>%
    left_join(code) %>%
    left_join(summary) %>%
    left_join(strain_summary) %>%
    left_join(phylotype_summary)

  write_csv(merge_final, "data/processed/strainge/straingst_ecolidb_summary.csv")
})

# ==============================================================================
# 2. StrainGR strain-sharing datasets
# ==============================================================================
# Batch 1: pairs of samples assigned the identical reference genome
# Batch 2: pairs of samples assigned highly similar reference genomes (Jaccard >= 0.7)
local({
  straingr <- bind_rows(
    read_straingr("data/raw/strainge/straingr/all_sample_comparisons_filtered.tsv"),
    read_straingr("data/raw/strainge/straingr/6May_all_sample_comparisons_filtered.tsv")
  )

  straingst <- read_csv("data/processed/strainge/straingst_ecolidb_summary.csv")
  lab <- read_csv("data/derived/lab/lab_data.csv") %>% select(four_dig_id, vendor)

  sample1 <- straingst %>% select(four_dig_id, smp_type, day, code) %>% left_join(lab) %>%
    rename_with(~ paste0(., "_smp1")) %>% rename(sample1 = four_dig_id_smp1)
  sample2 <- straingst %>% select(four_dig_id, smp_type, day, code) %>% left_join(lab) %>%
    rename_with(~ paste0(., "_smp2")) %>% rename(sample2 = four_dig_id_smp2)

  straingr_summary <- straingr %>% left_join(sample1) %>% left_join(sample2)

  write_csv(straingr_summary, "data/processed/strainge/strain_sharing_by_sample_all.csv")
  write_csv(straingr_summary %>% filter(vendor_smp1 == vendor_smp2),
            "data/processed/strainge/strain_sharing_within_vendor.csv")
  write_csv(straingr_summary %>% filter(vendor_smp1 != vendor_smp2, smp_type_smp1 != "LW", smp_type_smp2 != "LW"),
            "data/processed/strainge/strain_sharing_between_vendor.csv")
  write_csv(straingr_summary %>% filter(smp_type_smp1 == "LW" | smp_type_smp2 == "LW"),
            "data/processed/strainge/strain_sharing_by_latrine.csv")
})

# ==============================================================================
# 3. Vendor-level covariates (vendor survey and stall observations)
# ==============================================================================
local({
  survey_clean      <- read_csv("data/derived/field/survey_collapsed.csv")
  observation_clean <- read_csv("data/derived/field/vendor_obsv_collapsed.csv")

  survey_variables <- survey_clean %>%
    mutate(
      edu_secondary       = ifelse(edu_level_d1 >= 5, 1, 0),
      work_five_yrs       = ifelse(yrs_in_mkt_d1 >= 4, 1, 0),
      purch_300           = ifelse(qty_purch_tday_d1 >= 300, 1, 0),  # chickens purchased, day 1
      scr_commercial      = ifelse(src_d1 == 1, 1, 0),
      clean_cage_everyday = ifelse(clean_freq_d1 == 2, 1, 0),
      cage_35_chicks      = ifelse(cg_chik_qty_d2 >= 35, 1, 0)       # chickens in cage, day 2
    ) %>%
    rename(fs_training = fs_training_d1) %>%
    select(vd_id, edu_secondary, work_five_yrs, fs_training, purch_300,
           scr_commercial, clean_cage_everyday, cage_35_chicks)

  observation_variables <- observation_clean %>%
    mutate(
      vending_area = case_when(
        vd_id %in% c("V12", "V6")             ~ "Block A",
        vd_id %in% c("V2", "V9", "V11", "V4") ~ "Block B",
        vd_id %in% c("V8", "V5")              ~ "Block C"),
      cagebed_metal     = ifelse(cg_mtrl == 4, 1, 0),
      cagebed_plastic   = ifelse(cg_mtrl == 3, 1, 0),
      cagebed_cardboard = ifelse(cg_mtrl == 5, 1, 0),
      cagebed_cat       = case_when(cg_mtrl == 5 ~ 0, cg_mtrl == 3 ~ 1, TRUE ~ 2)
    ) %>%
    select(vd_id, `vend_ppe-apron`, cagebed_metal, cagebed_plastic, cagebed_cardboard, cagebed_cat, vending_area)

  vendor_variables <- survey_variables %>%
    left_join(observation_variables, by = "vd_id") %>%
    mutate(vd_id = as.integer(gsub("V", "", vd_id)))

  write.csv(vendor_variables, "data/processed/regression/vendor_variables.csv", row.names = FALSE)
})

# ==============================================================================
# 4. Slaughter-point-level covariates (behavioral and infrastructure observations)
# ==============================================================================
local({
  behavior_obsv_clean       <- read_csv("data/derived/field/slaughterhouse_collapsed.csv")
  infrastructure_obsv_clean <- read_csv("data/derived/field/slaughterhouse_condition.csv")

  behavior_obsv <- behavior_obsv_clean %>%
    mutate(
      dry_defeathering            = ifelse(sh_method == 1, 1, 0),
      feather_clean_slaugt_surf   = ifelse(slghtr_surf_clnm == 5, 1, 0),
      flies_present               = ifelse(slghtr_pest == 3, 1, 0),
      feather_clean_evisc_surf    = ifelse(evisc_surf_cln_mtrl == 5, 1, 0),
      feather_clean_evisc_utensil = ifelse(evisc_utens_clnm == 5, 1, 0)
    )

  # Proportion of observations at each slaughter point in which the behavior occurred
  # (missing values kept in the denominator)
  behavior_obsv_var <- behavior_obsv %>%
    select(sh_id, dry_defeathering, slghtr_apron, slghtr_boots, feather_clean_slaugt_surf,
           flies_present, proc_apron, proc_boots, feather_clean_evisc_surf,
           feather_clean_evisc_utensil, carcrins) %>%
    group_by(sh_id) %>%
    summarise(across(everything(), ~ sum(. == 1, na.rm = TRUE) / length(.), .names = "prop_{.col}"))

  behavior_obsv_var_2 <- behavior_obsv %>%
    select(vd_id, sh_id) %>%
    left_join(behavior_obsv_var, by = "sh_id") %>%
    distinct(vd_id, .keep_all = TRUE)

  infrastructure_obv_var <- infrastructure_obsv_clean %>%
    mutate(
      dry_slaughterhouse = ifelse(slaught_type == 1, 1, 0),
      wood_surfaces      = ifelse(slaught_surf == "1 7", 1, 0),
      no_waste_container = ifelse(slaught_waste_systm == 3, 1, 0)
    ) %>%
    select(slaught_id, dry_slaughterhouse, wood_surfaces, no_waste_container) %>%
    group_by(slaught_id) %>%
    summarise(across(everything(), ~ sum(. == 1, na.rm = TRUE) / length(.))) %>%
    rename(sh_id = slaught_id)

  slaughterhouse_variables <- behavior_obsv_var_2 %>%
    left_join(infrastructure_obv_var, by = "sh_id") %>%
    mutate(vd_id = as.numeric(gsub("V", "", vd_id)))

  write.csv(slaughterhouse_variables, "data/processed/regression/slaughterhouse_variables.csv", row.names = FALSE)
})

# ==============================================================================
# 5. All possible sample pairs, strain-sharing outcome, and pair-level covariates
# ==============================================================================
local({
  straingr <- bind_rows(
    read_straingr("data/raw/strainge/straingr/all_sample_comparisons_filtered.tsv"),
    read_straingr("data/raw/strainge/straingr/6May_all_sample_comparisons_filtered.tsv")
  )

  lab <- read_csv("data/derived/lab/lab_data.csv") %>%
    select(four_dig_id, vendor, smp_type, mc_yes_quant)
  clermont <- read_csv("data/reference/strainge_refdb_mlst_phylotype.csv") %>%
    select(symlink, phylotype, MLST) %>%
    rename(ref = symlink)
  dates <- read_csv("data/derived/lab/sample_collection_dates.csv")
  llod <- read_csv("data/derived/lab/lod_calcs.csv") %>%
    filter(media == "mc" & cef == "yes") %>%
    mutate(llod = log10(LLOD)) %>%
    select(smp_type, llod)
  vendor <- read.csv("data/processed/regression/vendor_variables.csv") %>%
    select(vd_id, vending_area) %>%
    rename(vendor = vd_id)

  # Cefotaxime-resistant E. coli load above the within-sample-type median and above the LLOD
  smp_type_medians <- lab %>%
    group_by(smp_type) %>%
    summarise(median_mc = median(mc_yes_quant, na.rm = TRUE), .groups = "drop")

  lab <- lab %>%
    left_join(smp_type_medians, by = "smp_type") %>%
    left_join(llod, by = "smp_type") %>%
    mutate(
      gt_median = case_when(is.na(mc_yes_quant) ~ NA_real_, mc_yes_quant > median_mc ~ 1, TRUE ~ 0),
      gt_llod   = case_when(is.na(mc_yes_quant) ~ NA_real_, mc_yes_quant > llod      ~ 1, TRUE ~ 0)
    ) %>%
    select(-median_mc, -llod, -mc_yes_quant) %>%
    left_join(vendor)

  # All pairwise comparisons, in both orders (sample1/sample2 can appear either way in StrainGR)
  pairs_unique <- as.data.frame(t(combn(lab$four_dig_id, 2))) %>% rename(sample1 = V1, sample2 = V2)
  all_pairs <- bind_rows(pairs_unique, pairs_unique %>% rename(sample1 = sample2, sample2 = sample1))

  # Pairs with strain sharing. When a pair shares more than one strain, ref, ST, and
  # phylogroup are set to "multiple".
  straingr_pairs <- straingr %>%
    select(1:3) %>%
    left_join(clermont) %>%
    mutate(
      ref       = str_remove(str_replace(ref, "Esch_coli_", ""), "\\.fa\\.gz"),
      phylo_new = case_when(phylotype %in% c("C", "E", "G", "F") ~ "Other",
                            phylotype %in% c("B2", "D")          ~ "B2/D",
                            TRUE                                 ~ phylotype),
      MLST      = as.character(MLST)
    ) %>%
    select(-phylotype) %>%
    group_by(sample1, sample2) %>%
    summarise(
      strain_sharing = 1,
      ref       = if_else(n_distinct(ref)       > 1, "multiple", first(ref)),
      phylo_new = if_else(n_distinct(phylo_new) > 1, "multiple", first(phylo_new)),
      MLST      = if_else(n_distinct(MLST)      > 1, "multiple", first(MLST)),
      .groups   = "drop"
    )

  # One row per unique pair, with pair-level covariates
  pairs <- all_pairs %>%
    left_join(straingr_pairs) %>%
    filter(sample1 < sample2) %>%
    left_join(lab %>% rename(sample1 = four_dig_id, vendor_smp1 = vendor, smp_type_smp1 = smp_type,
                             gt_median_smp1 = gt_median, gt_llod_smp1 = gt_llod, block_smp1 = vending_area),
              by = "sample1") %>%
    left_join(lab %>% rename(sample2 = four_dig_id, vendor_smp2 = vendor, smp_type_smp2 = smp_type,
                             gt_median_smp2 = gt_median, gt_llod_smp2 = gt_llod, block_smp2 = vending_area),
              by = "sample2") %>%
    left_join(dates %>% rename(sample1 = four_dig_id, date_smp1 = Today), by = "sample1") %>%
    left_join(dates %>% rename(sample2 = four_dig_id, date_smp2 = Today), by = "sample2") %>%
    mutate(
      strain_sharing = case_when(is.na(strain_sharing) ~ 0, TRUE ~ 1),
      same_vendor    = case_when(is.na(vendor_smp1 | vendor_smp2) ~ NA_real_,
                                 vendor_smp1 == vendor_smp2 ~ 1, TRUE ~ 0),
      same_smp_type  = case_when(smp_type_smp1 == smp_type_smp2 ~ 1, TRUE ~ 0),
      same_block     = case_when((vendor_smp1 == vendor_smp2) | is.na(vendor_smp1 | vendor_smp2) ~ NA_real_,
                                 block_smp1 == block_smp2 ~ 1, TRUE ~ 0),
      same_host      = case_when(
        smp_type_smp1 %in% c("CX", "CXS", "FD", "FDCX") & smp_type_smp2 %in% c("CX", "CXS", "FD", "FDCX") ~ 1,
        smp_type_smp1 %in% c("LW", "VF")                & smp_type_smp2 %in% c("LW", "VF")                ~ 1,
        smp_type_smp1 %in% c("FS", "HR", "WT", "SSx")   & smp_type_smp2 %in% c("FS", "HR", "WT", "SSx")   ~ 1,
        TRUE ~ 0),
      both_gt_median = case_when(is.na(gt_median_smp1) | is.na(gt_median_smp2) ~ NA_real_,
                                 gt_median_smp1 == 1 & gt_median_smp2 == 1     ~ 1, TRUE ~ 0),
      both_gt_llod   = case_when(is.na(gt_llod_smp1) | is.na(gt_llod_smp2)     ~ NA_real_,
                                 gt_llod_smp1 == 1 & gt_llod_smp2 == 1         ~ 1, TRUE ~ 0),
      days_diff = as.numeric(abs(date_smp2 - date_smp1))
    ) %>%
    select(-gt_median_smp1, -gt_median_smp2, -gt_llod_smp1, -gt_llod_smp2,
           -date_smp1, -date_smp2, -block_smp1, -block_smp2)

  write_csv(pairs, "data/processed/regression/all_strain_sharing_pairs.csv")

  # Any sample pair, excluding latrine waste
  write_csv(pairs %>% filter(smp_type_smp1 != "LW", smp_type_smp2 != "LW") %>% select(-vendor_smp1, -vendor_smp2),
            "data/processed/regression/market.csv")

  # Within-vendor pairs
  write_csv(pairs %>% filter(same_vendor == 1) %>% select(-vendor_smp1, -same_vendor, -same_block) %>%
              rename(vd_id = vendor_smp2),
            "data/processed/regression/vendor.csv")

  # Cross-vendor pairs, excluding latrine waste
  write_csv(pairs %>% filter(same_vendor == 0 | is.na(same_vendor)) %>%
              filter(smp_type_smp1 != "LW", smp_type_smp2 != "LW"),
            "data/processed/regression/cross_vendor.csv")
})

# ==============================================================================
# 6. Within-vendor regression dataset
# ==============================================================================
# Slaughter-point covariates are NA when both samples are vendor-level sample types
# (HR, FS, VF, WT), since those samples never passed through a slaughter point.
local({
  vendor_variables         <- read.csv("data/processed/regression/vendor_variables.csv")
  slaughterhouse_variables <- read.csv("data/processed/regression/slaughterhouse_variables.csv")
  vendor                   <- read.csv("data/processed/regression/vendor.csv")

  vendor_types <- c("HR", "FS", "VF", "WT")
  sh_cols <- setdiff(colnames(slaughterhouse_variables), "vd_id")

  vendor_log <- vendor %>%
    left_join(vendor_variables, by = "vd_id") %>%
    mutate(is_vendor_only = smp_type_smp1 %in% vendor_types & smp_type_smp2 %in% vendor_types) %>%
    left_join(slaughterhouse_variables, by = "vd_id") %>%
    mutate(across(all_of(sh_cols), ~ if_else(is_vendor_only, NA, .x))) %>%
    select(-is_vendor_only) %>%
    relocate(vd_id, .after = sample2) %>%
    relocate(sh_id, .after = cagebed_cat)

  write.csv(vendor_log, "data/processed/regression/vendor_log_updated.csv", row.names = FALSE)
})

# ==============================================================================
# 7. Cross-vendor regression dataset
# ==============================================================================
local({
  vendor_variables         <- read.csv("data/processed/regression/vendor_variables.csv")
  slaughterhouse_variables <- read.csv("data/processed/regression/slaughterhouse_variables.csv")
  cross_vendor             <- read.csv("data/processed/regression/cross_vendor.csv")

  vendor_types <- c("HR", "FS", "VF", "WT")
  slaughterhouse_vars <- setdiff(colnames(slaughterhouse_variables), "vd_id")

  cross_vendor_log <- cross_vendor %>%
    # Vendor covariates for each sample (suffix _s1 / _s2)
    left_join(vendor_variables %>% rename_with(~ paste0(., "_s1")), by = c("vendor_smp1" = "vd_id_s1")) %>%
    left_join(vendor_variables %>% rename_with(~ paste0(., "_s2")), by = c("vendor_smp2" = "vd_id_s2")) %>%
    # Slaughter-point covariates for each sample; NA when both samples are vendor-level types
    mutate(is_vendor_only = smp_type_smp1 %in% vendor_types & smp_type_smp2 %in% vendor_types) %>%
    left_join(slaughterhouse_variables %>% rename_with(~ paste0("sh_", ., "_s1"), all_of(slaughterhouse_vars)),
              by = c("vendor_smp1" = "vd_id")) %>%
    left_join(slaughterhouse_variables %>% rename_with(~ paste0("sh_", ., "_s2"), all_of(slaughterhouse_vars)),
              by = c("vendor_smp2" = "vd_id")) %>%
    mutate(across(starts_with("sh_"), ~ if_else(is_vendor_only, NA, .x))) %>%
    select(-is_vendor_only) %>%
    rename(sh_id_s1 = sh_sh_id_s1, sh_id_s2 = sh_sh_id_s2) %>%
    # Keep slaughter-point covariates only when both vendors use the same slaughter point
    mutate(same_sh = ifelse(sh_id_s1 == sh_id_s2, 1, 0)) %>%
    mutate(across(starts_with("sh_"), ~ ifelse(same_sh == 1, ., NA))) %>%
    select(-matches("^sh_.*s2$")) %>%
    rename_with(~ str_remove(., "_s1"), matches("^sh_.*s1$")) %>%
    # Pair-level vendor covariates
    mutate(
      both_secondary_school = ifelse(edu_secondary_s1 == 1 & edu_secondary_s2 == 1, 1, 0),
      both_scr_commercial   = ifelse(scr_commercial_s1 == 1 & scr_commercial_s2 == 1, 1, 0),
      both_clean_everyday   = ifelse(clean_cage_everyday_s1 == 1 & clean_cage_everyday_s2 == 1, 1, 0),
      both_fs_training      = ifelse(fs_training_s1 == 1 & fs_training_s2 == 1, 1, 0),
      both_ppe_apron        = ifelse(vend_ppe.apron_s1 == 1 & vend_ppe.apron_s2 == 1, 1, 0),
      same_block            = ifelse(vending_area_s1 == vending_area_s2, 1, 0)
    ) %>%
    select(-c(ends_with("_s1"), ends_with("_s2"))) %>%
    relocate(c(same_block, same_sh), .before = sh_id)

  write.csv(cross_vendor_log, "data/processed/regression/cross_vendor_log.csv", row.names = FALSE)
})

message("Preprocessing complete. Processed datasets are in data/processed/.")
