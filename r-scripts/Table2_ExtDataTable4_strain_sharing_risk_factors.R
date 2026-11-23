# Table2_ExtDataTable4_strain_sharing_risk_factors.R
# Risk factors for cefotaxime-resistant E. coli strain sharing among (1) any sample
# pair excluding latrine waste, (2) within-vendor pairs, and (3) cross-vendor pairs.
#
# Produces: Table 2                (unadjusted and adjusted odds ratios)
#           Extended Data Table 4  (distribution of risk factors by strain-sharing status)
#
# Approach: univariable logistic regression (p < 0.2 retained) -> pairwise correlation
# screen (|r| >= 0.7) -> multivariable logistic regression -> variables with VIF > 5
# removed and model refit. Variables with <= 10 observations in any outcome cell are
# excluded. Within-vendor models include a random intercept for vendor.
#
# Inputs:  data/processed/regression/market.csv
#          data/processed/regression/vendor_log_updated.csv
#          data/processed/regression/cross_vendor_log.csv
# Outputs: results/tables/Table2_ExtDataTable4_any_pair.csv
#          results/tables/Table2_ExtDataTable4_within_vendor.csv
#          results/tables/Table2_ExtDataTable4_cross_vendor.csv

source("r-scripts/00_functions.R")
library(broom)
library(broom.mixed)
library(lme4)
library(car)

dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Analysis function
# ------------------------------------------------------------------------------
# df            analysis dataset (one row per sample pair; outcome = strain_sharing)
# summary_vars  variables described in Extended Data Table 4
# cont_vars     continuous variables (summarized as median and SD)
# reg_vars      candidate variables for regression
# corr_drop     variables removed after the correlation screen
# random        random-effect term (e.g., "(1 | vd_id)"), or NULL for fixed effects only
# sd_style      "plusminus" -> "median ± SD"; "paren" -> "median (SD: x)"
analyze_sharing <- function(df, summary_vars, cont_vars, reg_vars, corr_drop = character(0),
                            random = NULL, sd_style = "plusminus") {
  cat_vars <- setdiff(summary_vars, cont_vars)
  n_shared_total     <- sum(df$strain_sharing == 1, na.rm = TRUE)
  n_not_shared_total <- sum(df$strain_sharing == 0, na.rm = TRUE)

  # --- Extended Data Table 4: distributions by strain-sharing status ---
  # Categorical: n (%) for level = 1 and level = 0 within each outcome group
  cat_summary <- map_dfr(cat_vars, function(var) {
    map_dfr(c(1, 0), function(lvl) {
      tibble(variable     = var,
             level        = lvl,
             n_shared     = sum(df[[var]] == lvl & df$strain_sharing == 1, na.rm = TRUE),
             n_not_shared = sum(df[[var]] == lvl & df$strain_sharing == 0, na.rm = TRUE))
    })
  }) %>%
    mutate(shared_col     = paste0(n_shared,     " (", round(100 * n_shared / n_shared_total, 1),         "%)"),
           not_shared_col = paste0(n_not_shared, " (", round(100 * n_not_shared / n_not_shared_total, 1), "%)"))

  fmt_cont <- function(x) {
    m <- round(median(x, na.rm = TRUE), 1)
    s <- round(sd(x, na.rm = TRUE), 1)
    if (sd_style == "paren") paste0(m, " (SD: ", s, ")") else paste0(m, " ± ", s)
  }
  cont_summary <- map_dfr(cont_vars, function(var) {
    tibble(variable       = var,
           level          = NA_real_,
           shared_col     = fmt_cont(as.numeric(df[[var]][df$strain_sharing == 1])),
           not_shared_col = fmt_cont(as.numeric(df[[var]][df$strain_sharing == 0])))
  })

  summary_table <- bind_rows(cat_summary, cont_summary) %>%
    select(variable, level, shared_col, not_shared_col) %>%
    rename(!!paste0("Strain sharing (n=",    n_shared_total,     ")") := shared_col,
           !!paste0("No strain sharing (n=", n_not_shared_total, ")") := not_shared_col)

  # --- Exclude categorical variables with perfect or near separation (<= 10 in any cell) ---
  separated <- keep(setdiff(reg_vars, cont_vars), function(var) {
    tab <- table(df$strain_sharing, df[[var]], useNA = "no")
    ncol(tab) == 2 && any(tab <= 10)
  })
  reg_vars <- setdiff(reg_vars, separated)

  # --- Model fitting ---
  fit <- function(vars) {
    f <- as.formula(paste("strain_sharing ~", paste(c(vars, random), collapse = " + ")))
    if (is.null(random)) glm(f, data = df, family = binomial) else glmer(f, data = df, family = binomial)
  }
  tidy_fixed <- function(model) {
    out <- tidy(model, exponentiate = TRUE, conf.int = TRUE)
    if ("effect" %in% names(out)) out <- out %>% filter(effect == "fixed") %>% select(-effect, -group)
    out %>% filter(term != "(Intercept)")
  }

  # Univariable models; retain p < 0.2
  results <- map_dfr(reg_vars, ~ tidy_fixed(fit(.x))) %>% filter(p.value < 0.2)
  final_vars <- unique(sub("[123]$", "", results$term))

  # Pairwise correlations among retained variables
  # (both categorical -> Pearson; both continuous -> Spearman; mixed -> point-biserial)
  strong_corr <- map_dfr(combn(final_vars, 2, simplify = FALSE), function(pair) {
    x <- as.numeric(df[[pair[1]]])
    y <- as.numeric(df[[pair[2]]])
    if (sum(is.finite(x) & is.finite(y)) < 3) return(NULL)
    both_cont <- !(pair[1] %in% cat_vars) && !(pair[2] %in% cat_vars)
    test <- cor.test(x, y, method = if (both_cont) "spearman" else "pearson", use = "complete.obs", exact = FALSE)
    tibble(var1 = pair[1], var2 = pair[2], estimate = unname(test$estimate), p.value = test$p.value)
  }) %>%
    filter(p.value < 0.05, abs(estimate) >= 0.7)
  print(strong_corr)

  # Multivariable model; drop variables with VIF > 5 and refit
  multi_vars  <- setdiff(final_vars, corr_drop)
  vif_values  <- car::vif(fit(multi_vars)) %>% as.data.frame() %>% rownames_to_column("var") %>% filter(. > 5)
  multi_vars  <- setdiff(multi_vars, vif_values$var)
  multi_results <- tidy_fixed(fit(multi_vars))

  # --- Table 2: unadjusted and adjusted ORs (on the first row of each variable) ---
  format_or <- function(res, label) {
    res %>%
      mutate(variable  = sub("1$", "", term),
             !!label  := paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")"),
             first_row = TRUE) %>%
      select(variable, first_row, all_of(label))
  }

  summary_table %>%
    group_by(variable) %>%
    mutate(first_row = row_number() == 1) %>%
    ungroup() %>%
    left_join(format_or(results,       "Unadjusted OR (95% CI)"), by = c("variable", "first_row")) %>%
    left_join(format_or(multi_results, "Adjusted OR (95% CI)"),   by = c("variable", "first_row")) %>%
    select(-first_row)
}

# ------------------------------------------------------------------------------
# 1. Any sample pair, excluding latrine waste
# ------------------------------------------------------------------------------
market <- read_csv("data/processed/regression/market.csv") %>%
  mutate(across(!c(days_diff), as.factor), across(c(days_diff), as.numeric)) %>%
  mutate(sameday = case_when(is.na(days_diff) ~ NA_real_, days_diff == 0 ~ 1, TRUE ~ 0))

market_vars <- setdiff(colnames(market),
                       c("strain_sharing", "sample1", "sample2", "ref", "phylo_new", "MLST",
                         "smp_type_smp1", "smp_type_smp2", "days_diff", "same_block"))

any_pair <- analyze_sharing(market,
                            summary_vars = market_vars,
                            cont_vars    = character(0),
                            reg_vars     = market_vars,
                            corr_drop    = c("both_gt_llod", "same_block"),
                            sd_style     = "paren")

write_csv(any_pair, "results/tables/Table2_ExtDataTable4_any_pair.csv")

# ------------------------------------------------------------------------------
# 2. Within-vendor pairs (random intercept for vendor)
# ------------------------------------------------------------------------------
# Slaughter-point proportions are rescaled to percentages
vendor_log <- read_csv("data/processed/regression/vendor_log_updated.csv") %>%
  mutate(across(!c(starts_with("prop"), days_diff), as.factor),
         across(c(starts_with("prop"), days_diff), as.numeric)) %>%
  mutate(across(starts_with("prop"), ~ . * 100))

within_vars <- setdiff(colnames(vendor_log),
                       c("strain_sharing", "sample1", "sample2", "vd_id", "ref", "phylo_new", "MLST",
                         "smp_type_smp1", "smp_type_smp2", "sh_id", "same_vendor", "cage_35_chicks",
                         "days_diff", "cagebed_metal", "cagebed_plastic", "cagebed_cardboard",
                         "vending_area"))  # vending area is constant within a vendor

# Cage material: the regression uses the 3-level cagebed_cat (0 = cardboard, 1 = plastic,
# 2 = other/metal), but Extended Data Table 4 reports "Vendor stored chickens in metallic cage"
# as a Yes/No row, so the summary uses the binary cagebed_metal. (Summarizing cagebed_cat
# levels 1 and 0 reported plastic and cardboard counts under the metallic-cage label.)
within_summary_vars <- c(setdiff(within_vars, "cagebed_cat"), "cagebed_metal")

# Variables excluded from regression because they do not explain a unique slice of the variance
within_reg_vars <- setdiff(within_vars,
                           c("work_five_yrs", "prop_dry_defeathering", "prop_slghtr_boots", "prop_proc_apron",
                             "prop_proc_boots", "prop_feather_clean_evisc_surf", "prop_carcrins", "prop_slghtr_apron"))

within_vendor <- analyze_sharing(vendor_log,
                                 summary_vars = within_summary_vars,
                                 cont_vars    = within_summary_vars[startsWith(within_summary_vars, "prop")],
                                 reg_vars     = within_reg_vars,
                                 corr_drop    = c("both_gt_llod", "prop_feather_clean_evisc_utensil",
                                                  "prop_feather_clean_slaugt_surf"),
                                 random       = "(1 | vd_id)")

write_csv(within_vendor, "results/tables/Table2_ExtDataTable4_within_vendor.csv")

# ------------------------------------------------------------------------------
# 3. Cross-vendor pairs
# ------------------------------------------------------------------------------
cross_vendor_log <- read_csv("data/processed/regression/cross_vendor_log.csv") %>%
  select(!sh_id) %>%
  mutate(across(!c(starts_with("sh_prop"), days_diff), as.factor),
         across(c(starts_with("sh_prop"), days_diff), as.numeric)) %>%
  mutate(sameday = case_when(is.na(days_diff) ~ NA_real_, days_diff == 0 ~ 1, TRUE ~ 0)) %>%
  mutate(across(starts_with("sh_prop"), ~ . * 100))

cross_vars <- setdiff(colnames(cross_vendor_log),
                      c("strain_sharing", "sample1", "sample2", "vendor_smp1", "vendor_smp2", "ref", "phylo_new",
                        "MLST", "smp_type_smp1", "smp_type_smp2", "same_vendor", "days_diff"))

cross_vendor <- analyze_sharing(cross_vendor_log,
                                summary_vars = cross_vars,
                                cont_vars    = cross_vars[startsWith(cross_vars, "sh_prop")],
                                reg_vars     = cross_vars)

write_csv(cross_vendor, "results/tables/Table2_ExtDataTable4_cross_vendor.csv")
