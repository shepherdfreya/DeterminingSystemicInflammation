############################################################
# File: Scripts/05a_run_mr_array.R
# Purpose: Harmonise one exposure-outcome pair and run MR
# Order:
#   1. Read proxied exposure data
#   2. Read matching outcome data
#   3. Harmonise exposure and outcome
#   4. Apply Steiger filtering where possible
#   5. Run MR methods based on SNP count
#   6. Run heterogeneity and MR-Egger intercept tests where possible
#   7. Save pair-specific MR results
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript Scripts/05a_run_mr_array.R <input_dir> <output_dir> <task_id>")
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

plan <- read.csv(
  file.path(input_dir, "analysis_plan.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)


manifest_path <- file.path(output_dir, "trait_manifest_standardised.csv")
if (!file.exists(manifest_path)) {
  manifest_path <- file.path(input_dir, "trait_manifest.csv")
}

manifest <- read.csv(
    manifest_path,
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

annotate_instrument_region <- function(exp_dat, exposure_row) {
  if (!has_protein_coordinates(exposure_row)) {
    exp_dat$instrument_region <- "unknown"
    return(exp_dat)
  }

  needed <- c("chr.exposure", "pos.exposure")
  if (!all(needed %in% names(exp_dat))) {
    stop("Exposure data must contain columns for cis/trans annotation: ", paste(needed, collapse = ", "))
  }

  protein_chr <- as.character(exposure_row$protein_chr_col[[1]])
  start <- suppressWarnings(as.numeric(exposure_row$protein_start_col[[1]])) - CIS_WINDOW
  end <- suppressWarnings(as.numeric(exposure_row$protein_end_col[[1]])) + CIS_WINDOW

  exp_dat %>%
    dplyr::mutate(
      instrument_region = dplyr::if_else(
        as.character(.data$chr.exposure) == protein_chr &
          suppressWarnings(as.numeric(.data$pos.exposure)) >= start &
          suppressWarnings(as.numeric(.data$pos.exposure)) <= end,
        "cis",
        "trans"
      )
    )
}

filter_exposure_by_analysis <- function(exp_dat, exposure_row, analysis_mode) {
  exp_dat <- annotate_instrument_region(exp_dat, exposure_row)

  if (analysis_mode == "cis_only") {
    if (!has_protein_coordinates(exposure_row)) {
      stop("Analysis = cis_only requires protein coordinates for exposure: ", exposure_row$trait_id[[1]])
    }

    exp_dat <- exp_dat %>%
      dplyr::filter(.data$instrument_region == "cis")
  }

  exp_dat
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
  summary_dir <- file.path(output_dir, "mr_summaries")
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

get_method_list <- function(n_snps) {
  if (n_snps == 1) {
    return(c("mr_wald_ratio"))
  }

  if (n_snps == 2) {
    return(c("mr_ivw"))
  }

  c("mr_ivw", "mr_weighted_median", "mr_egger_regression")
}

run_steiger_safe <- function(hdat) {

  needed <- c("samplesize.exposure", "samplesize.outcome")

  if (!all(needed %in% names(hdat))) {
    hdat$steiger_available <- FALSE
    hdat$steiger_filter_applied <- FALSE
    return(hdat)
  }

  if (all(is.na(hdat$samplesize.exposure)) || all(is.na(hdat$samplesize.outcome))) {
    hdat$steiger_available <- FALSE
    hdat$steiger_filter_applied <- FALSE
    return(hdat)
  }

  out <- tryCatch(
    {
      TwoSampleMR::steiger_filtering(hdat) %>%
        dplyr::filter(.data$mr_keep) %>%
        dplyr::mutate(
          steiger_available = TRUE,
          steiger_filter_applied = TRUE
        )
    },
    error = function(e) {
      warning("Steiger filtering failed: ", conditionMessage(e))
      hdat$steiger_available <- FALSE
      hdat$steiger_filter_applied <- FALSE
      hdat
    }
  )

  out
}

run_mr_for_pair <- function(hdat) {

  hdat <- hdat %>%
    dplyr::filter(
      !is.na(.data$beta.exposure),
      !is.na(.data$se.exposure),
      !is.na(.data$beta.outcome),
      !is.na(.data$se.outcome)
    ) %>%
    dplyr::distinct(.data$id.exposure, .data$id.outcome, .data$SNP, .keep_all = TRUE)

  n_snps <- dplyr::n_distinct(hdat$SNP)

  if (n_snps == 0) {
    return(data.frame())
  }

  method_list <- get_method_list(n_snps)

  mr_res <- tryCatch(
    {
      TwoSampleMR::mr(hdat, method_list = method_list) %>%
        TwoSampleMR::generate_odds_ratios()
    },
    error = function(e) {
      warning("MR failed: ", conditionMessage(e))
      data.frame()
    }
  )

  if (nrow(mr_res) == 0) {
    return(data.frame())
  }

  pleio <- if ("mr_egger_regression" %in% method_list) {
    tryCatch(
      {
        TwoSampleMR::mr_pleiotropy_test(hdat) %>%
          dplyr::transmute(
            id.exposure,
            id.outcome,
            egger_intercept = .data$egger_intercept,
            egger_intercept_se = .data$se,
            egger_intercept_pval = .data$pval
          )
      },
      error = function(e) {
        data.frame(
          id.exposure = unique(hdat$id.exposure),
          id.outcome = unique(hdat$id.outcome),
          egger_intercept = NA_real_,
          egger_intercept_se = NA_real_,
          egger_intercept_pval = NA_real_
        )
      }
    )
  } else {
    data.frame(
      id.exposure = unique(hdat$id.exposure),
      id.outcome = unique(hdat$id.outcome),
      egger_intercept = NA_real_,
      egger_intercept_se = NA_real_,
      egger_intercept_pval = NA_real_
    )
  }

  het <- if (n_snps >= 2 && "mr_ivw" %in% method_list) {
    tryCatch(
      {
        TwoSampleMR::mr_heterogeneity(hdat, method_list = "mr_ivw") %>%
          dplyr::select(
            .data$id.exposure,
            .data$id.outcome,
            .data$method,
            .data$Q,
            .data$Q_df,
            .data$Q_pval
          )
      },
      error = function(e) {
        data.frame()
      }
    )
  } else {
    data.frame()
  }

  out <- mr_res %>%
    dplyr::left_join(pleio, by = c("id.exposure", "id.outcome"))

  if (nrow(het) > 0) {
    out <- out %>%
      dplyr::left_join(het, by = c("id.exposure", "id.outcome", "method"))
  } else {
    out$Q <- NA_real_
    out$Q_df <- NA_real_
    out$Q_pval <- NA_real_
  }

  out$n_snps <- n_snps

  out
}

# -------------------------
# Select exposure-outcome pair for task
# -------------------------

if (!all(c("exposure", "outcome", "run") %in% names(plan))) {
  stop("analysis_plan.csv must contain columns: exposure, outcome, run")
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

safe_pair <- make_safe_pair(exposure_id, outcome_id, analysis_mode, include_analysis_in_name)

message("========================================")
message("Running MR")
message("Exposure: ", exposure_id)
message("Outcome : ", outcome_id)
message("Analysis: ", analysis_mode)
message("========================================")

# -------------------------
# Output files
# -------------------------

mr_dir <- file.path(output_dir, "MR_results", "per_pair")
dir.create(mr_dir, recursive = TRUE, showWarnings = FALSE)

harmonised_dir <- file.path(output_dir, "Harmonised", "per_pair")
dir.create(harmonised_dir, recursive = TRUE, showWarnings = FALSE)

out_mr         <- file.path(mr_dir, paste0(safe_pair, "__mr_results.tsv"))
out_harmonised <- file.path(harmonised_dir, paste0(safe_pair, "__harmonised.tsv"))
done_file      <- file.path(mr_dir, paste0(safe_pair, ".__DONE"))

result_row <- data.frame(
  exposure = exposure_id,
  outcome = outcome_id,
  analysis = analysis_mode,
  status = "failed",
  err_message = "",
  exposure_snps = 0L,
  outcome_snps = 0L,
  harmonised_snps = 0L,
  mr_keep_snps = 0L,
  steiger_filter_applied = FALSE,
  mr_rows = 0L,
  stringsAsFactors = FALSE
)

if (file.exists(done_file) && file.exists(out_mr)) {
  message("Already completed: ", exposure_id, " -> ", outcome_id)
  quit(save = "no", status = 0)
}

# -------------------------
# Input files
# -------------------------

exposure_file <- file.path(
  output_dir,
  "Proxied_exposures",
  "per_pair",
  paste0(safe_pair, "__exposure_proxied.tsv")
)

outcome_file <- file.path(output_dir, "all_outcomes.txt")

if (!file.exists(exposure_file)) {
  result_row$err_message <- paste("Proxied exposure file not found:", exposure_file)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

if (!file.exists(outcome_file)) {
  result_row$err_message <- paste("Outcome file not found:", outcome_file)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

exp_dat <- data.table::fread(
  exposure_file,
  data.table = FALSE,
  showProgress = FALSE
)

exposure_row <- manifest[match(exposure_id, manifest$trait_id), , drop = FALSE]
exp_dat <- filter_exposure_by_analysis(exp_dat, exposure_row, analysis_mode)

# -------------------------
# Standardise exposure IDs
# -------------------------

if ("id.exposure" %in% names(exp_dat)) {

  first_id <- exp_dat$id.exposure[!is.na(exp_dat$id.exposure) & exp_dat$id.exposure != ""][1]

  if (is.na(first_id) || length(first_id) == 0) {
    first_id <- exposure_id
  }

  exp_dat$id.exposure <- first_id

} else {
  exp_dat$id.exposure <- exposure_id
}

exp_dat$exposure <- exposure_id

all_outcomes <- data.table::fread(
  outcome_file,
  data.table = FALSE,
  showProgress = FALSE
)

out_dat <- all_outcomes %>%
  dplyr::filter(.data$outcome == outcome_id)

result_row$exposure_snps <- dplyr::n_distinct(exp_dat$SNP)
result_row$outcome_snps <- dplyr::n_distinct(out_dat$SNP)

if (nrow(exp_dat) == 0) {
  result_row$err_message <- paste("No exposure data found for:", exposure_id)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

if (nrow(out_dat) == 0) {
  result_row$err_message <- paste("No outcome data found for:", outcome_id)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

# -------------------------
# Harmonise
# -------------------------

hdat <- tryCatch(
  {
    TwoSampleMR::harmonise_data(
      exposure_dat = exp_dat,
      outcome_dat = out_dat,
      action = 2
    )
  },
  error = function(e) {
    result_row$err_message <<- paste("Harmonisation failed:", conditionMessage(e))
    data.frame()
  }
)

if (nrow(hdat) == 0) {
  result_row$err_message <- ifelse(
    result_row$err_message == "",
    "No harmonised SNPs",
    result_row$err_message
  )
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

result_row$harmonised_snps <- dplyr::n_distinct(hdat$SNP)

# -------------------------
# Add sample sizes for Steiger filtering
# -------------------------

if (!"samplesize" %in% names(manifest)) {
  warning("No 'samplesize' column found in trait_manifest.csv. Steiger filtering will be skipped.")
} else {

  exposure_n <- manifest %>%
    dplyr::filter(.data$trait_id == exposure_id) %>%
    dplyr::pull(.data$samplesize) %>%
    .[1]

  outcome_n <- manifest %>%
    dplyr::filter(.data$trait_id == outcome_id) %>%
    dplyr::pull(.data$samplesize) %>%
    .[1]

  hdat$samplesize.exposure <- suppressWarnings(as.numeric(exposure_n))
  hdat$samplesize.outcome  <- suppressWarnings(as.numeric(outcome_n))
}




# Keep only MR-valid SNPs
hdat <- hdat %>%
  dplyr::filter(.data$mr_keep)

result_row$mr_keep_snps <- dplyr::n_distinct(hdat$SNP)

if (nrow(hdat) == 0) {
  result_row$status <- "completed"
  result_row$err_message <- "No SNPs retained after harmonisation mr_keep filtering"

  safe_write_table(data.frame(), out_mr)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  file.create(done_file)
  quit(save = "no", status = 0)
}

# -------------------------
# Steiger filtering, if possible
# -------------------------

hdat <- run_steiger_safe(hdat)

result_row$steiger_filter_applied <- any(hdat$steiger_filter_applied, na.rm = TRUE)

if (nrow(hdat) == 0) {
  result_row$status <- "completed"
  result_row$err_message <- "No SNPs retained after Steiger filtering"

  safe_write_table(data.frame(), out_mr)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  file.create(done_file)
  quit(save = "no", status = 0)
}

safe_write_table(hdat, out_harmonised)

# -------------------------
# Run MR
# -------------------------

mr_res <- run_mr_for_pair(hdat)

if (nrow(mr_res) == 0) {
  result_row$status <- "completed"
  result_row$err_message <- "MR returned no results"
  result_row$mr_rows <- 0L

  safe_write_table(data.frame(), out_mr)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  file.create(done_file)
  quit(save = "no", status = 0)
}

mr_res <- mr_res %>%
  dplyr::mutate(
    exposure = exposure_id,
    outcome = outcome_id,
    analysis_set = analysis_mode
  )

result_row$status <- "completed"
result_row$err_message <- ""
result_row$mr_rows <- nrow(mr_res)

safe_write_table(mr_res, out_mr)
write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)

file.create(done_file)

message("MR complete for: ", exposure_id, " -> ", outcome_id)
message("MR rows: ", nrow(mr_res))
