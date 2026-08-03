############################################################
# File: Scripts/04c_apply_proxies_array.R
# Purpose: Apply best LD proxies to one exposure-outcome pair
# Order:
#   1. Read exposure + outcome data
#   2. Load proxy candidates
#   3. Keep proxies present in exposure full GWAS
#   4. Select best proxy by R2
#   5. Format proxy exposure data
#   6. Remove missing SNPs
#   7. Merge proxy SNPs into exposure dataset
#   8. Save proxied exposure file
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript Scripts/04c_apply_proxies_array.R <input_dir> <output_dir> <task_id>")
}

input_dir  <- args[1]
output_dir <- args[2]
task_id    <- as.integer(args[3])

if (is.na(task_id) || task_id < 1) {
  stop("task_id must be a positive integer.")
}

source("Utils/Packages.R")
source("Utils/Config.R")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(output_dir, "trait_manifest_standardised.csv")
if (!file.exists(manifest_path)) {
  manifest_path <- file.path(input_dir, "trait_manifest.csv")
}

manifest <- read.csv(
  manifest_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

plan <- read.csv(
  file.path(input_dir, "analysis_plan.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# -------------------------
# Helper functions
# -------------------------

normalise_run_flag <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1")
}

get_analysis_mode <- function(pair) {
  if ("Analysis" %in% names(pair)) {
    mode <- tolower(trimws(as.character(pair$Analysis[[1]])))
    if (mode %in% c("cis_only", "cis_trans")) {
      return(mode)
    }
  }

  "cis_trans"
}

make_safe_pair <- function(exposure_id, outcome_id, analysis_mode, include_analysis) {
  base <- paste0(make.names(exposure_id), "__", make.names(outcome_id))

  if (include_analysis) {
    return(paste0(base, "__", make.names(analysis_mode)))
  }

  base
}

is_provided <- function(x) {
  !(is.na(x) || trimws(as.character(x)) == "" ||
      trimws(as.character(x)) %in% c("NA", "na", "Na", "n/a", "N/A"))
}

has_protein_coordinates <- function(row) {
  is_provided(row$protein_chr_col[[1]]) &&
    is_provided(row$protein_start_col[[1]]) &&
    is_provided(row$protein_end_col[[1]])
}

filter_exposure_by_analysis <- function(exp_df, exposure_row, analysis_mode) {
  if (analysis_mode == "cis_trans") {
    exp_df$instrument_region <- "cis_trans"
    return(exp_df)
  }

  if (!has_protein_coordinates(exposure_row)) {
    stop("Analysis = cis_only requires protein coordinates for exposure: ", exposure_row$trait_id[[1]])
  }

  needed <- c("chr.exposure", "pos.exposure")
  if (!all(needed %in% names(exp_df))) {
    stop("Exposure data must contain columns for cis filtering: ", paste(needed, collapse = ", "))
  }

  protein_chr <- as.character(exposure_row$protein_chr_col[[1]])
  start <- suppressWarnings(as.numeric(exposure_row$protein_start_col[[1]])) - CIS_WINDOW
  end <- suppressWarnings(as.numeric(exposure_row$protein_end_col[[1]])) + CIS_WINDOW

  exp_df %>%
    dplyr::mutate(
      instrument_region = dplyr::if_else(
        as.character(.data$chr.exposure) == protein_chr &
          suppressWarnings(as.numeric(.data$pos.exposure)) >= start &
          suppressWarnings(as.numeric(.data$pos.exposure)) <= end,
        "cis",
        "trans"
      )
    ) %>%
    dplyr::filter(.data$instrument_region == "cis")
}

get_sep_value <- function(x) {
  if (is.na(x) || x == "") {
    stop("Missing 'sep' value in trait_manifest.csv")
  }
  if (x == "\\t") "\t" else x
}

safe_write_table <- function(dat, file) {
  write.table(
    dat,
    file = file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

write_per_task_summary <- function(summary_row, output_dir, exposure_name, outcome_name) {
  summary_dir <- file.path(output_dir, "proxy_application_summaries")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  file_stem <- paste0(make.names(exposure_name), "__", make.names(outcome_name))
  if (
    exists("include_analysis_in_name", inherits = TRUE) &&
    isTRUE(include_analysis_in_name) &&
    exists("analysis_mode", inherits = TRUE)
  ) {
    file_stem <- paste0(file_stem, "__", make.names(analysis_mode))
  }

  out <- file.path(
    summary_dir,
    paste0(file_stem, ".tsv")
  )

  safe_write_table(summary_row, out)
}

empty_proxy_map <- function() {
  data.frame(
    exposure = character(),
    outcome = character(),
    missing_snp = character(),
    proxy_snp = character(),
    R = numeric(),
    R2 = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_no_proxy <- function() {
  data.frame(
    exposure = character(),
    outcome = character(),
    missing_snp = character(),
    reason = character(),
    stringsAsFactors = FALSE
  )
}

standardise_raw_gwas <- function(dat, row) {

  required_cols <- c(
    row$snp_col[[1]],
    row$effect_col[[1]],
    row$se_col[[1]],
    row$ea_col[[1]],
    row$oa_col[[1]],
    row$pval_col[[1]]
  )

  missing_required <- required_cols[!required_cols %in% names(dat)]

  if (length(missing_required) > 0) {
    stop(
      "Raw exposure GWAS is missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }

  out <- data.frame(
    SNP = dat[[row$snp_col[[1]]]],
    beta = suppressWarnings(as.numeric(dat[[row$effect_col[[1]]]])),
    se = suppressWarnings(as.numeric(dat[[row$se_col[[1]]]])),
    effect_allele = dat[[row$ea_col[[1]]]],
    other_allele = dat[[row$oa_col[[1]]]],
    pval = suppressWarnings(as.numeric(dat[[row$pval_col[[1]]]])),
    stringsAsFactors = FALSE
  )

  if ("eaf_col" %in% names(row) && !is.na(row$eaf_col[[1]]) && row$eaf_col[[1]] != "" && row$eaf_col[[1]] %in% names(dat)) {
    out$eaf <- suppressWarnings(as.numeric(dat[[row$eaf_col[[1]]]]))
  } else {
    out$eaf <- NA_real_
  }

  if ("chromosome_col" %in% names(row) && !is.na(row$chromosome_col[[1]]) && row$chromosome_col[[1]] != "" && row$chromosome_col[[1]] %in% names(dat)) {
    out$chr <- dat[[row$chromosome_col[[1]]]]
  } else {
    out$chr <- NA
  }

  if ("position_col" %in% names(row) && !is.na(row$position_col[[1]]) && row$position_col[[1]] != "" && row$position_col[[1]] %in% names(dat)) {
    out$pos <- dat[[row$position_col[[1]]]]
  } else {
    out$pos <- NA
  }

  if ("effect_type" %in% names(row)) {
    effect_type <- row$effect_type[[1]]
  } else {
    effect_type <- NA_character_
  }

  if (!is.na(effect_type) && toupper(effect_type) == "OR") {
    out$beta <- log(out$beta)
  }

  out <- out[!is.na(out$SNP) & out$SNP != "", , drop = FALSE]
  out <- out[!duplicated(out$SNP), , drop = FALSE]

  out
}

format_proxy_exposure_data <- function(proxy_raw, exposure_id) {

  proxy_raw$exposure <- exposure_id

  TwoSampleMR::format_data(
    proxy_raw,
    type = "exposure",
    phenotype_col = "exposure",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    eaf_col = "eaf",
    pval_col = "pval",
    chr_col = "chr",
    pos_col = "pos"
  )
}

add_f_stat <- function(dat) {
  needed <- c("beta.exposure", "se.exposure")
  if (!all(needed %in% names(dat))) {
    stop("Exposure data must contain columns: beta.exposure, se.exposure")
  }

  dat$F_stat <- suppressWarnings((as.numeric(dat$beta.exposure) / as.numeric(dat$se.exposure))^2)
  dat
}

# -------------------------
# Select exposure-outcome pair for task
# -------------------------

if (!all(c("exposure", "outcome", "run") %in% names(plan))) {
  stop("analysis_plan.csv must contain columns: exposure, outcome, run")
}

if (!all(c("trait_id", "file", "sep") %in% names(manifest))) {
  stop("trait_manifest.csv must contain at least: trait_id, file, sep")
}

plan$run <- normalise_run_flag(plan$run)
plan <- subset(plan, run)

if (nrow(plan) == 0) {
  stop("No rows with run = TRUE found in analysis_plan.csv")
}

if (task_id > nrow(plan)) {
  stop(
    "task_id ", task_id, " is out of range. There are only ",
    nrow(plan), " active exposure-outcome pairs."
  )
}

pair <- plan[task_id, , drop = FALSE]

exposure_id <- pair$exposure[[1]]
outcome_id  <- pair$outcome[[1]]
analysis_mode <- get_analysis_mode(pair)
include_analysis_in_name <- "Analysis" %in% names(plan)

if (!exposure_id %in% manifest$trait_id) {
  stop("Exposure trait_id not found in trait_manifest.csv: ", exposure_id)
}

row <- manifest[match(exposure_id, manifest$trait_id), , drop = FALSE]

safe_pair <- make_safe_pair(exposure_id, outcome_id, analysis_mode, include_analysis_in_name)

message("========================================")
message("Applying proxies")
message("Exposure: ", exposure_id)
message("Outcome : ", outcome_id)
message("Analysis: ", analysis_mode)
message("========================================")

# -------------------------
# Output files
# -------------------------

out_dir <- file.path(output_dir, "Proxied_exposures", "per_pair")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_data  <- file.path(out_dir, paste0(safe_pair, "__exposure_proxied.tsv"))
out_map   <- file.path(out_dir, paste0(safe_pair, "__proxy_map.tsv"))
out_noprx <- file.path(out_dir, paste0(safe_pair, "__no_proxy.tsv"))
done_file <- file.path(out_dir, paste0(safe_pair, ".__DONE"))

result_row <- data.frame(
  exposure = exposure_id,
  outcome = outcome_id,
  analysis = analysis_mode,
  status = "failed",
  err_message = "",
  original_exposure_snps = 0L,
  outcome_snps = 0L,
  missing_in_outcome = 0L,
  original_retained = 0L,
  candidate_proxy_rows = 0L,
  best_proxy_rows = 0L,
  proxy_f_filtered = 0L,
  proxies_added = 0L,
  final_total = 0L,
  no_proxy = 0L,
  stringsAsFactors = FALSE
)

if (file.exists(done_file) && file.exists(out_data)) {
  message("Already completed: ", exposure_id, " -> ", outcome_id)
  quit(save = "no", status = 0)
}

# -------------------------
# Input files
# -------------------------

exposure_file <- file.path(output_dir, "all_exposures.txt")
outcome_file  <- file.path(output_dir, "all_outcomes.txt")
proxy_file    <- file.path(output_dir, "Proxy_candidates", "proxy_candidates.tsv")

if (!file.exists(exposure_file)) {
  result_row$err_message <- paste("Exposure file not found:", exposure_file)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

if (!file.exists(outcome_file)) {
  result_row$err_message <- paste("Outcome file not found:", outcome_file)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

if (!file.exists(proxy_file)) {
  result_row$err_message <- paste("Proxy candidate file not found:", proxy_file)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

exposures <- data.table::fread(
  exposure_file,
  data.table = FALSE,
  showProgress = FALSE
)

outcomes <- data.table::fread(
  outcome_file,
  data.table = FALSE,
  showProgress = FALSE
)

proxies <- data.table::fread(
  proxy_file,
  data.table = FALSE,
  showProgress = FALSE
)

if (!all(c("SNP", "exposure") %in% names(exposures))) {
  result_row$err_message <- "all_exposures.txt must contain columns: SNP, exposure"
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

if (!all(c("SNP", "outcome") %in% names(outcomes))) {
  result_row$err_message <- "all_outcomes.txt must contain columns: SNP, outcome"
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

if (nrow(proxies) > 0 && !all(c("exposure", "outcome", "missing_snp", "proxy_snp", "R2") %in% names(proxies))) {
  result_row$err_message <- "proxy_candidates.tsv must contain columns: exposure, outcome, missing_snp, proxy_snp, R2"
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

exp_df <- exposures %>%
  dplyr::filter(.data$exposure == exposure_id) %>%
  dplyr::distinct(.data$SNP, .keep_all = TRUE)

exp_df <- filter_exposure_by_analysis(exp_df, row, analysis_mode)

out_df <- outcomes %>%
  dplyr::filter(.data$outcome == outcome_id) %>%
  dplyr::distinct(.data$SNP, .keep_all = TRUE)

result_row$original_exposure_snps <- nrow(exp_df)
result_row$outcome_snps <- nrow(out_df)

if (nrow(exp_df) == 0) {
  result_row$err_message <- paste("No exposure rows found for:", exposure_id)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

if (nrow(out_df) == 0) {
  result_row$err_message <- paste("No outcome rows found for:", outcome_id)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

# -------------------------
# Identify missing SNPs
# -------------------------

missing_snps <- setdiff(unique(exp_df$SNP), unique(out_df$SNP))
result_row$missing_in_outcome <- length(missing_snps)

if (length(missing_snps) == 0) {

  result_row$status <- "completed"
  result_row$original_retained <- nrow(exp_df)
  result_row$final_total <- nrow(exp_df)
  result_row$err_message <- ""

  safe_write_table(exp_df, out_data)
  safe_write_table(empty_proxy_map(), out_map)
  safe_write_table(empty_no_proxy(), out_noprx)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)

  file.create(done_file)

  message("No missing SNPs. Original exposure data retained.")
  quit(save = "no", status = 0)
}

# -------------------------
# Read full exposure GWAS
# -------------------------

file_path <- file.path(input_dir, row$file[[1]])
sep_value <- get_sep_value(row$sep[[1]])

if (!file.exists(file_path)) {
  result_row$err_message <- paste("Raw exposure GWAS file not found:", file_path)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

raw_dat <- data.table::fread(
  file_path,
  sep = sep_value,
  data.table = FALSE,
  showProgress = FALSE
)

raw_std <- standardise_raw_gwas(raw_dat, row)

# -------------------------
# Select best proxies
# -------------------------

if (nrow(proxies) == 0) {
  pair_proxies <- data.frame()
} else {
  pair_proxies <- proxies %>%
    dplyr::filter(
      .data$exposure == exposure_id,
      .data$outcome == outcome_id,
      .data$missing_snp %in% missing_snps,
      .data$proxy_snp %in% out_df$SNP,
      .data$proxy_snp %in% raw_std$SNP
    )

  if ("analysis" %in% names(pair_proxies)) {
    pair_proxies <- pair_proxies %>%
      dplyr::filter(.data$analysis == analysis_mode)
  }

  pair_proxies <- pair_proxies %>%
    dplyr::mutate(
      R = if ("R" %in% names(.)) suppressWarnings(as.numeric(.data$R)) else NA_real_,
      R2 = suppressWarnings(as.numeric(.data$R2))
    ) %>%
    dplyr::filter(!is.na(.data$R2)) %>%
    dplyr::group_by(.data$missing_snp) %>%
    dplyr::slice_max(order_by = .data$R2, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
}

result_row$candidate_proxy_rows <- nrow(pair_proxies)

if (nrow(pair_proxies) > 0) {
  best_proxies <- pair_proxies %>%
    dplyr::group_by(.data$proxy_snp) %>%
    dplyr::slice_max(order_by = .data$R2, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()
} else {
  best_proxies <- pair_proxies
}

result_row$best_proxy_rows <- nrow(best_proxies)

lost_due_to_duplicate_proxy <- setdiff(pair_proxies$missing_snp, best_proxies$missing_snp)
no_proxy_snps <- setdiff(missing_snps, best_proxies$missing_snp)

no_proxy_tbl <- data.frame(
  exposure = exposure_id,
  outcome = outcome_id,
  missing_snp = c(no_proxy_snps, lost_due_to_duplicate_proxy),
  reason = c(
    rep("No proxy present in both outcome GWAS and exposure full GWAS", length(no_proxy_snps)),
    rep("Proxy duplicated another missing SNP and lower R2 was dropped", length(lost_due_to_duplicate_proxy))
  ),
  stringsAsFactors = FALSE
)

result_row$no_proxy <- nrow(no_proxy_tbl)

exp_clean <- exp_df %>%
  dplyr::filter(!.data$SNP %in% missing_snps)

result_row$original_retained <- nrow(exp_clean)

if (nrow(best_proxies) == 0) {

  result_row$status <- "completed"
  result_row$final_total <- nrow(exp_clean)
  result_row$err_message <- "No valid proxies found"

  safe_write_table(exp_clean, out_data)
  safe_write_table(empty_proxy_map(), out_map)
  safe_write_table(no_proxy_tbl, out_noprx)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)

  file.create(done_file)

  message("No valid proxies found. Missing SNPs removed.")
  quit(save = "no", status = 0)
}

# -------------------------
# Format proxy rows
# -------------------------

proxy_raw <- raw_std %>%
  dplyr::filter(.data$SNP %in% best_proxies$proxy_snp) %>%
  dplyr::left_join(
    best_proxies %>%
      dplyr::select(
        "missing_snp",
        "proxy_snp",
        "R",
        "R2"
      ),
    by = c("SNP" = "proxy_snp")
  )

proxy_fmt <- format_proxy_exposure_data(proxy_raw, exposure_id)

proxy_fmt <- proxy_fmt %>%
  dplyr::left_join(
    proxy_raw %>%
      dplyr::select("SNP", "missing_snp", "R", "R2"),
    by = "SNP"
  )

proxy_fmt <- add_f_stat(proxy_fmt)

proxy_fmt_strong <- proxy_fmt %>%
  dplyr::filter(!is.na(.data$F_stat), .data$F_stat > EXPOSURE_F_THRESHOLD)

weak_proxy_snps <- setdiff(proxy_fmt$missing_snp, proxy_fmt_strong$missing_snp)

if (length(weak_proxy_snps) > 0) {
  no_proxy_tbl <- dplyr::bind_rows(
    no_proxy_tbl,
    data.frame(
      exposure = exposure_id,
      outcome = outcome_id,
      missing_snp = weak_proxy_snps,
      reason = paste0("Selected proxy failed F-statistic threshold > ", EXPOSURE_F_THRESHOLD),
      stringsAsFactors = FALSE
    )
  ) %>%
    dplyr::distinct(.data$exposure, .data$outcome, .data$missing_snp, .keep_all = TRUE)
}

result_row$proxy_f_filtered <- length(weak_proxy_snps)
result_row$no_proxy <- nrow(no_proxy_tbl)

proxy_fmt_novel <- proxy_fmt_strong %>%
  dplyr::anti_join(exp_clean, by = "SNP")

final_df <- dplyr::bind_rows(exp_clean, proxy_fmt_novel)

mapping_tbl <- proxy_fmt_novel %>%
  dplyr::transmute(
    exposure = exposure_id,
    outcome = outcome_id,
    missing_snp = .data$missing_snp,
    proxy_snp = .data$SNP,
    R = .data$R,
    R2 = .data$R2
  )

result_row$status <- "completed"
result_row$err_message <- ""
result_row$proxies_added <- nrow(proxy_fmt_novel)
result_row$final_total <- nrow(final_df)

safe_write_table(final_df, out_data)
safe_write_table(mapping_tbl, out_map)
safe_write_table(no_proxy_tbl, out_noprx)
write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)

file.create(done_file)

message("Finished applying proxies: ", exposure_id, " -> ", outcome_id)
message("Original SNPs retained: ", nrow(exp_clean))
message("Proxy SNPs added: ", nrow(proxy_fmt_novel))
message("Final SNPs: ", nrow(final_df))
