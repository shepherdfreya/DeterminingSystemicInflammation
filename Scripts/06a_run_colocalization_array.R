############################################################
# File: Scripts/06a_run_colocalization_array.R
# Purpose: Run genetic colocalization for one exposure-outcome pair
# Usage:
#   Rscript Scripts/06a_run_colocalization_array.R <input_dir> <output_dir> <task_id>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript Scripts/06a_run_colocalization_array.R <input_dir> <output_dir> <task_id>")
}

input_dir  <- args[1]
output_dir <- args[2]
task_id    <- as.integer(args[3])

if (is.na(task_id) || task_id < 1) {
  stop("task_id must be a positive integer.")
}

source("Utils/Packages.R")
source("Utils/Config.R")

if (!requireNamespace("coloc", quietly = TRUE)) {
  stop(
    "The 'coloc' package is required for colocalization. Install it with: ",
    "install.packages('coloc')"
  )
}

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

empty_coloc_result <- function() {
  data.frame(
    exposure = character(),
    outcome = character(),
    analysis = character(),
    locus_id = character(),
    chr = character(),
    start = numeric(),
    end = numeric(),
    n_overlap = integer(),
    status = character(),
    err_message = character(),
    PP.H0.abf = numeric(),
    PP.H1.abf = numeric(),
    PP.H2.abf = numeric(),
    PP.H3.abf = numeric(),
    PP.H4.abf = numeric(),
    H4_threshold = numeric(),
    H4_flag = character(),
    stringsAsFactors = FALSE
  )
}

standardise_raw_gwas <- function(dat, row) {
  required_cols <- c(
    row$snp_col[[1]],
    row$position_col[[1]],
    row$effect_col[[1]],
    row$se_col[[1]],
    row$pval_col[[1]],
    row$eaf_col[[1]],
    row$chromosome_col[[1]]
  )

  missing_required <- required_cols[!required_cols %in% names(dat)]

  if (length(missing_required) > 0) {
    stop(
      "Raw GWAS is missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }

  out <- data.frame(
    SNP = dat[[row$snp_col[[1]]]],
    chr = as.character(dat[[row$chromosome_col[[1]]]]),
    pos = suppressWarnings(as.numeric(dat[[row$position_col[[1]]]])),
    beta = suppressWarnings(as.numeric(dat[[row$effect_col[[1]]]])),
    se = suppressWarnings(as.numeric(dat[[row$se_col[[1]]]])),
    pval = suppressWarnings(as.numeric(dat[[row$pval_col[[1]]]])),
    eaf = suppressWarnings(as.numeric(dat[[row$eaf_col[[1]]]])),
    stringsAsFactors = FALSE
  )

  out <- out[!is.na(out$SNP) & out$SNP != "", , drop = FALSE]
  out <- out[!is.na(out$chr) & out$chr != "", , drop = FALSE]
  out <- out[!is.na(out$pos), , drop = FALSE]
  out <- out[!is.na(out$beta), , drop = FALSE]
  out <- out[!is.na(out$se) & out$se > 0, , drop = FALSE]
  out <- out[!is.na(out$eaf) & out$eaf > 0 & out$eaf < 1, , drop = FALSE]
  out <- out[!duplicated(out$SNP), , drop = FALSE]

  out
}

get_trait_n <- function(row) {
  if ("samplesize" %in% names(row) && is_provided(row$samplesize[[1]])) {
    return(suppressWarnings(as.numeric(row$samplesize[[1]])))
  }
  NA_real_
}

get_loci_for_pair <- function(exposure_id, exposure_row, output_dir) {
  if (
    is_provided(exposure_row$protein_chr_col[[1]]) &&
    is_provided(exposure_row$protein_start_col[[1]]) &&
    is_provided(exposure_row$protein_end_col[[1]])
  ) {
    start <- suppressWarnings(as.numeric(exposure_row$protein_start_col[[1]])) - cis_window
    end   <- suppressWarnings(as.numeric(exposure_row$protein_end_col[[1]])) + cis_window

    return(data.frame(
      locus_id = paste0(exposure_id, "_cis_region"),
      chr = as.character(exposure_row$protein_chr_col[[1]]),
      start = max(0, start),
      end = end,
      stringsAsFactors = FALSE
    ))
  }

  exposure_file <- file.path(output_dir, "all_exposures.txt")

  if (!file.exists(exposure_file)) {
    stop(
      "No protein coordinates found for exposure and all_exposures.txt is missing: ",
      exposure_file
    )
  }

  exposures <- data.table::fread(
    exposure_file,
    data.table = FALSE,
    showProgress = FALSE
  )

  needed <- c("SNP", "exposure", "chr.exposure", "pos.exposure")
  if (!all(needed %in% names(exposures))) {
    stop("all_exposures.txt must contain columns: ", paste(needed, collapse = ", "))
  }

  exposures %>%
    dplyr::filter(.data$exposure == exposure_id) %>%
    dplyr::distinct(.data$SNP, .keep_all = TRUE) %>%
    dplyr::transmute(
      locus_id = .data$SNP,
      chr = as.character(.data$chr.exposure),
      start = pmax(0, suppressWarnings(as.numeric(.data$pos.exposure)) - cis_window),
      end = suppressWarnings(as.numeric(.data$pos.exposure)) + cis_window
    )
}

run_coloc_for_locus <- function(locus, exp_raw, out_raw, exposure_id, outcome_id, analysis_mode, exposure_n, outcome_n) {
  result <- empty_coloc_result()
  result[1, ] <- NA

  result$exposure <- exposure_id
  result$outcome <- outcome_id
  result$analysis <- analysis_mode
  result$locus_id <- locus$locus_id[[1]]
  result$chr <- locus$chr[[1]]
  result$start <- locus$start[[1]]
  result$end <- locus$end[[1]]
  result$n_overlap <- 0L
  result$status <- "failed"
  result$err_message <- ""
  result$H4_threshold <- COLOC_H4_THRESHOLD
  result$H4_flag <- "Fail"

  if (is.na(exposure_n) || is.na(outcome_n)) {
    result$err_message <- "Missing samplesize for exposure or outcome"
    return(result)
  }

  exp_region <- exp_raw %>%
    dplyr::filter(
      .data$chr == locus$chr[[1]],
      .data$pos >= locus$start[[1]],
      .data$pos <= locus$end[[1]]
    )

  out_region <- out_raw %>%
    dplyr::filter(
      .data$chr == locus$chr[[1]],
      .data$pos >= locus$start[[1]],
      .data$pos <= locus$end[[1]]
    )

  merged <- dplyr::inner_join(
    exp_region,
    out_region,
    by = "SNP",
    suffix = c(".exposure", ".outcome")
  )

  result$n_overlap <- nrow(merged)

  if (nrow(merged) < COLOC_MIN_OVERLAP) {
    result$err_message <- paste0(
      "Fewer than COLOC_MIN_OVERLAP overlapping SNPs: ",
      nrow(merged),
      " < ",
      COLOC_MIN_OVERLAP
    )
    return(result)
  }

  coloc_res <- tryCatch(
    {
      coloc::coloc.abf(
        dataset1 = list(
          beta = merged$beta.exposure,
          varbeta = merged$se.exposure^2,
          snp = merged$SNP,
          MAF = merged$eaf.exposure,
          N = exposure_n,
          type = "quant"
        ),
        dataset2 = list(
          beta = merged$beta.outcome,
          varbeta = merged$se.outcome^2,
          snp = merged$SNP,
          MAF = merged$eaf.outcome,
          N = outcome_n,
          type = "quant"
        )
      )
    },
    error = function(e) {
      result$err_message <<- paste("coloc.abf failed:", conditionMessage(e))
      NULL
    }
  )

  if (is.null(coloc_res)) {
    return(result)
  }

  summary <- coloc_res$summary

  result$status <- "completed"
  result$err_message <- ""
  result$PP.H0.abf <- unname(summary["PP.H0.abf"])
  result$PP.H1.abf <- unname(summary["PP.H1.abf"])
  result$PP.H2.abf <- unname(summary["PP.H2.abf"])
  result$PP.H3.abf <- unname(summary["PP.H3.abf"])
  result$PP.H4.abf <- unname(summary["PP.H4.abf"])
  result$H4_flag <- ifelse(result$PP.H4.abf > COLOC_H4_THRESHOLD, "Pass", "Fail")

  result
}

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

missing_traits <- setdiff(c(exposure_id, outcome_id), manifest$trait_id)
if (length(missing_traits) > 0) {
  stop("Traits missing from trait_manifest.csv: ", paste(missing_traits, collapse = ", "))
}

exposure_row <- manifest[match(exposure_id, manifest$trait_id), , drop = FALSE]
outcome_row  <- manifest[match(outcome_id, manifest$trait_id), , drop = FALSE]

safe_pair <- make_safe_pair(exposure_id, outcome_id, analysis_mode, include_analysis_in_name)

message("========================================")
message("Running colocalization")
message("Exposure: ", exposure_id)
message("Outcome : ", outcome_id)
message("Analysis: ", analysis_mode)
message("H4 pass threshold: ", COLOC_H4_THRESHOLD)
message("========================================")

out_dir <- file.path(output_dir, "Colocalization", "per_pair")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, paste0(safe_pair, "__colocalization.tsv"))
done_file <- file.path(out_dir, paste0(safe_pair, ".__DONE"))

if (file.exists(done_file) && file.exists(out_file)) {
  message("Already completed: ", exposure_id, " -> ", outcome_id)
  quit(save = "no", status = 0)
}

exp_path <- file.path(input_dir, exposure_row$file[[1]])
out_path <- file.path(input_dir, outcome_row$file[[1]])

if (!file.exists(exp_path)) {
  stop("Exposure GWAS file not found: ", exp_path)
}
if (!file.exists(out_path)) {
  stop("Outcome GWAS file not found: ", out_path)
}

exp_dat <- data.table::fread(
  exp_path,
  sep = get_sep_value(exposure_row$sep[[1]]),
  data.table = FALSE,
  showProgress = FALSE
)

out_dat <- data.table::fread(
  out_path,
  sep = get_sep_value(outcome_row$sep[[1]]),
  data.table = FALSE,
  showProgress = FALSE
)

exp_raw <- standardise_raw_gwas(exp_dat, exposure_row)
out_raw <- standardise_raw_gwas(out_dat, outcome_row)

loci <- get_loci_for_pair(exposure_id, exposure_row, output_dir)
loci <- loci[!is.na(loci$chr) & !is.na(loci$start) & !is.na(loci$end), , drop = FALSE]

if (nrow(loci) == 0) {
  stop("No colocalization loci could be defined for exposure: ", exposure_id)
}

exposure_n <- get_trait_n(exposure_row)
outcome_n  <- get_trait_n(outcome_row)

results <- lapply(seq_len(nrow(loci)), function(i) {
  run_coloc_for_locus(
    locus = loci[i, , drop = FALSE],
    exp_raw = exp_raw,
    out_raw = out_raw,
    exposure_id = exposure_id,
    outcome_id = outcome_id,
    analysis_mode = analysis_mode,
    exposure_n = exposure_n,
    outcome_n = outcome_n
  )
})

results <- dplyr::bind_rows(results)

safe_write_table(results, out_file)
file.create(done_file)

message("Colocalization complete for: ", exposure_id, " -> ", outcome_id)
message("Loci tested: ", nrow(results))
message("H4 pass loci: ", sum(results$H4_flag == "Pass", na.rm = TRUE))
