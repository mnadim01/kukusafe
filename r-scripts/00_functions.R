# 00_functions.R
# Shared constants and helper functions. Sourced at the top of every other script.

library(tidyverse)

# Note: samples excluded from all analyses (associated with a positive field blank,
# vendor 3, duplicates, and an unidentifiable sample) have already been removed from
# the input files in data/. See docs/data_dictionary.md.

# ------------------------------------------------------------------------------
# Sample types
# ------------------------------------------------------------------------------
smp_code_map <- c(
  "FD"   = "Cage fecal droppings",
  "FDCX" = "Live chicken rectal swab",
  "CXS"  = "Pre-evisceration carcass",
  "CX"   = "Ready-to-sell carcass",
  "HR"   = "Vendor hand rinse",
  "FS"   = "Food fed to chickens",
  "WT"   = "Water fed to chickens",
  "SSx"  = "Slaughterhouse table surface",
  "LW"   = "Latrine waste",
  "VF"   = "Vendor feces"
)

# Host category of each sample type
codecat_recode <- c(
  "Cage fecal droppings"         = "Chicken",
  "Live chicken rectal swab"     = "Chicken",
  "Pre-evisceration carcass"     = "Chicken",
  "Ready-to-sell carcass"        = "Chicken",
  "Vendor hand rinse"            = "Environment",
  "Food fed to chickens"         = "Environment",
  "Water fed to chickens"        = "Environment",
  "Slaughterhouse table surface" = "Environment",
  "Latrine waste"                = "Human",
  "Vendor feces"                 = "Human"
)

# Plot order of sample type codes
smp_level_order <- c("HR", "FS", "WT", "FD", "FDCX", "CXS", "SSx", "CX", "VF", "LW")

# ------------------------------------------------------------------------------
# Reading StrainGE outputs
# ------------------------------------------------------------------------------
# StrainGR comparison tables: strip the "ER" suffix from sample IDs
read_straingr <- function(path) {
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(sample1 = as.numeric(gsub("ER", "", sample1)),
           sample2 = as.numeric(gsub("ER", "", sample2)))
}

# StrainGST tables: extract the four-digit sample ID from the sample name
parse_straingst_id <- function(df) {
  df %>%
    separate(sample, into = c(NA, "four_dig_id", NA, NA, NA), sep = "_") %>%
    mutate(four_dig_id = as.numeric(str_extract(four_dig_id, "^.+(?=E)")))
}

read_straingst_strains <- function(path = "data/raw/strainge/straingst/ecoli_straingst_strains.tsv") {
  read_tsv(path, show_col_types = FALSE) %>%
    parse_straingst_id() %>%
    select(1:2) %>%
    mutate(strain = gsub("_new", "", strain))
}

# Label a StrainGE reference genome as "name (ST x/phylogroup)"
ref_label <- function(ref, MLST, phylotype) {
  paste0(str_remove(str_remove(ref, "^Esch_coli_"), "\\.fa\\.gz$"), " (ST ", MLST, "/", phylotype, ")")
}

# ------------------------------------------------------------------------------
# Strain-sharing pairs
# ------------------------------------------------------------------------------
# Flag whether the two samples in each pair were collected on the same day
add_smpday <- function(df, dates) {
  # Coerce collection dates to Date so the difference is numeric even if `Today`
  # was read as character (e.g., dates loaded without readr's type guessing)
  dates <- dates %>% mutate(Today = as.Date(Today))
  df %>%
    left_join(dates %>% select(four_dig_id, date1 = Today), by = c("sample1" = "four_dig_id")) %>%
    left_join(dates %>% select(four_dig_id, date2 = Today), by = c("sample2" = "four_dig_id")) %>%
    mutate(datediff = abs(date2 - date1),
           smpday   = ifelse(datediff == 0, "Same day", "Different day"))
}

# Host category of each sample and transmission pathway of each pair (e.g., "Chicken-Human")
add_transpath <- function(df) {
  df %>%
    mutate(
      type1     = dplyr::recode(code_smp1, !!!codecat_recode),
      type2     = dplyr::recode(code_smp2, !!!codecat_recode),
      transpath = case_when(type1 == type2 ~ paste(type1, type2, sep = "-"),
                            TRUE           ~ paste(pmin(type1, type2), pmax(type1, type2), sep = "-"))
    )
}

# Within- and between-vendor strain-sharing datasets and sharing rates
# (used by Figure 4C, Figure 5, and Extended Data Figures 3 and 5)
load_vendor_sharing <- function() {
  phylogroups  <- read_csv("data/processed/strainge/straingst_ecolidb_summary.csv", show_col_types = FALSE)
  wivendor_raw <- read_csv("data/processed/strainge/strain_sharing_within_vendor.csv", show_col_types = FALSE)
  bwvendor_raw <- read_csv("data/processed/strainge/strain_sharing_between_vendor.csv", show_col_types = FALSE)
  dates        <- read_csv("data/derived/lab/sample_collection_dates.csv", show_col_types = FALSE)

  # Within-vendor sharing; pairs canonicalized so (A,B) and (B,A) are the same pair
  wivendor <- wivendor_raw %>%
    add_transpath() %>%
    add_smpday(dates) %>%
    filter(sample1 != sample2) %>%
    mutate(s1 = pmin(sample1, sample2), s2 = pmax(sample1, sample2), vendor = vendor_smp1)

  # Between-vendor sharing
  bwvendor <- bwvendor_raw %>%
    add_transpath() %>%
    add_smpday(dates)

  # Denominator: number of samples per vendor
  vendsample_counts <- phylogroups %>%
    filter(!is.na(vendor)) %>%
    count(vendor, name = "n_samples")

  # Within-vendor rate = unique sharing pairs / possible pairs. Vendor 6 had no within-vendor sharing.
  byvendssresults <- vendsample_counts %>%
    full_join(wivendor %>% count(vendor, name = "total_sharing_events"), by = "vendor") %>%
    full_join(wivendor %>% distinct(vendor, s1, s2) %>% count(vendor, name = "unique_sharing_pairs"), by = "vendor") %>%
    mutate(possible_pairs       = choose(n_samples, 2),
           rate                 = unique_sharing_pairs / possible_pairs,
           total_sharing_events = case_when(vendor == 6 ~ 0, TRUE ~ total_sharing_events),
           rate                 = case_when(vendor == 6 ~ 0, TRUE ~ rate))

  # Between-vendor rate = unique sharing pairs / (samples from vendor 1 x samples from vendor 2)
  bw_unique_pairs <- bwvendor %>%
    filter(vendor_smp1 != vendor_smp2) %>%
    mutate(v1 = pmin(vendor_smp1, vendor_smp2), v2 = pmax(vendor_smp1, vendor_smp2),
           s1 = pmin(sample1, sample2),         s2 = pmax(sample1, sample2)) %>%
    distinct(v1, v2, s1, s2) %>%
    count(v1, v2, name = "unique_sharing_pairs")

  between_vendor_results <- expand.grid(v1 = vendsample_counts$vendor, v2 = vendsample_counts$vendor) %>%
    as_tibble() %>%
    filter(v1 < v2) %>%
    left_join(vendsample_counts, by = c("v1" = "vendor")) %>% rename(n1 = n_samples) %>%
    left_join(vendsample_counts, by = c("v2" = "vendor")) %>% rename(n2 = n_samples) %>%
    mutate(possible_pairs = n1 * n2) %>%
    left_join(bw_unique_pairs, by = c("v1", "v2")) %>%
    mutate(unique_sharing_pairs = replace_na(unique_sharing_pairs, 0),
           rate                 = unique_sharing_pairs / possible_pairs)

  list(wivendor = wivendor, bwvendor = bwvendor,
       byvendssresults = byvendssresults, between_vendor_results = between_vendor_results)
}

# ------------------------------------------------------------------------------
# Network plots
# ------------------------------------------------------------------------------
shape_map <- c(
  "Live chicken rectal swab"     = 16,
  "Pre-evisceration carcass"     = 18,
  "Ready-to-sell carcass"        = 15,
  "Vendor hand rinse"            = 8,
  "Food fed to chickens"         = 17,
  "Water fed to chickens"        = 25,
  "Slaughterhouse table surface" = 4,
  "Vendor feces"                 = 3
)

# Node sizes (pre-evisceration carcass diamond enlarged relative to other shapes)
size_map_network <- 2 * c(
  "Live chicken rectal swab"     = 2.5,
  "Pre-evisceration carcass"     = 3.5,
  "Ready-to-sell carcass"        = 2.5,
  "Vendor hand rinse"            = 2.5,
  "Food fed to chickens"         = 2.5,
  "Water fed to chickens"        = 2.5,
  "Slaughterhouse table surface" = 2.5,
  "Vendor feces"                 = 2.5
)

# Cross-vendor strain-sharing network for one reference strain.
# Nodes are samples (shape = sample type, color = vendor); edges are sharing events.
plot_bwvend_network <- function(strain_type, data) {
  data_filt <- data %>% filter(ref == strain_type)

  edges <- data_filt %>%
    mutate(s1 = pmin(sample1, sample2), s2 = pmax(sample1, sample2)) %>%
    distinct(s1, s2, smpday)

  g <- igraph::graph_from_data_frame(edges, directed = FALSE)

  sample_types <- data_filt %>%
    select(sample = sample1, smp_type = code_smp1, vendor = vendor_smp1) %>%
    bind_rows(data_filt %>% select(sample = sample2, smp_type = code_smp2, vendor = vendor_smp2)) %>%
    distinct()
  igraph::V(g)$smp_type <- sample_types$smp_type[match(igraph::V(g)$name, sample_types$sample)]
  igraph::V(g)$vendor   <- sample_types$vendor[match(igraph::V(g)$name, sample_types$sample)]

  vendor_colors <- setNames(scales::hue_pal()(n_distinct(sample_types$vendor)), unique(sample_types$vendor))

  set.seed(1)  # reproducible force-directed layout
  ggraph::ggraph(g, layout = "fr") +
    ggraph::geom_edge_link(aes(linetype = smpday), colour = "grey60", alpha = 0.6, width = 0.8) +
    ggraph::geom_node_point(aes(shape = smp_type, color = factor(vendor), size = smp_type)) +
    scale_size_manual(values = size_map_network, guide = "none") +
    ggraph::scale_edge_linetype_manual(values = c("Same day" = "solid", "Different day" = "dashed"),
                                       name = "Sample day") +
    scale_shape_manual(values = shape_map, name = "Sample type") +
    scale_color_manual(values = vendor_colors, name = "Vendor") +
    guides(shape = guide_legend(override.aes = list(size = 4)),
           color = guide_legend(override.aes = list(size = 4))) +
    theme_void() +
    ggtitle(paste("Between-vendor strain sharing:", strain_type))
}

# ------------------------------------------------------------------------------
# Descriptive tables
# ------------------------------------------------------------------------------
# Format a count as "n (%)", rounding percentages half up (e.g., 62.5% -> 63%)
n_pct <- function(n, total) paste0(n, " (", floor(100 * n / total + 0.5), "%)")

# Tabulate one categorical characteristic: one row per category, in the order given.
# `levels` is a named vector mapping coded values to display labels; values not
# listed (including NA) are counted under `missing_label`, if any occur.
tab_categorical <- function(x, characteristic, levels, total = length(x), missing_label = "Missing") {
  x <- as.character(x)
  out <- tibble(Characteristic = characteristic,
                Category       = unname(levels),
                n              = map_int(names(levels), ~ sum(x == .x, na.rm = TRUE)))
  n_other <- sum(is.na(x) | !x %in% names(levels))
  if (n_other > 0) out <- bind_rows(out, tibble(Characteristic = characteristic, Category = missing_label, n = n_other))
  out %>% mutate(`n (%)` = n_pct(n, total)) %>% select(-n)
}

# Tabulate yes/no indicators: one row per indicator, counting values equal to 1
tab_binary <- function(df, characteristic, vars, total = nrow(df)) {
  tibble(Characteristic = characteristic,
         Category       = names(vars),
         `n (%)`        = map_chr(vars, ~ n_pct(sum(df[[.x]] == 1, na.rm = TRUE), total)))
}
