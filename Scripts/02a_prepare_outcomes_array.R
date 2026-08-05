############################################################
# File: Scripts/02_prepare_outcomes_array.R
# Purpose: Prepare ONE outcome GWAS file per array task
# Usage:
#   Rscript Scripts/02_prepare_outcomes_array.R <input_dir> <output_dir> <task_id>
#
# Example:
#   Rscript Scripts/02_prepare_outcomes_array.R data results 3
#
# Notes:
# - task_id should correspond to the row number of the unique active outcomes
# - use a separate combine script afterwards to bind all prepared outcomes
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript Scripts/02_prepare_outcomes_array.R <input_dir> <output_dir> <task_id>")
}

input_dir  <- args[1]
output_dir <- args[2]
task_id    <- as.integer(args[3])

if (is.na(task_id) || task_id < 1) {
  stop("task_id must be a positive integer.")
}

source("Utils/Packages.R")

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

fread_args_for_path <- function(file_path) {
  if (grepl("\\.gz$", file_path)) {
    list(cmd = paste("gzip -dc", shQuote(file_path)))
  } else {
    list(file = file_path)
  }
}

read_gwas_table <- function(file_path, sep_value) {
  do.call(
    data.table::fread,
    c(
      fread_args_for_path(file_path),
      list(
        sep = sep_value,
        data.table = FALSE,
        showProgress = FALSE
      )
    )
  )
}

resolve_input_file <- function(input_dir, file_value) {
  file_value <- trimws(as.character(file_value))
  if (grepl("^(/|~)", file_value)) {
    return(path.expand(file_value))
  }
  file.path(input_dir, file_value)
}

# Keep only active analyses
plan$run <- tolower(as.character(plan$run)) %in% c("true", "t", "1")
plan <- subset(plan, run)

if (nrow(plan) == 0) {
  stop("No rows with run = TRUE found in analysis_plan.csv")
}

# Unique outcomes needed
outcome_traits <- unique(plan$outcome)

# Subset manifest to required outcomes, preserving analysis_plan order
outcome_manifest <- manifest[match(outcome_traits, manifest$trait_id), , drop = FALSE]

# Check for missing outcomes in manifest
missing_outcomes <- outcome_traits[is.na(match(outcome_traits, manifest$trait_id))]
if (length(missing_outcomes) > 0) {
  stop(
    "The following outcome trait_id values were found in analysis_plan.csv but not trait_manifest.csv: ",
    paste(missing_outcomes, collapse = ", ")
  )
}

# Check task_id range
if (task_id > nrow(outcome_manifest)) {
  stop(
    "task_id ", task_id, " is out of range. There are only ",
    nrow(outcome_manifest), " unique active outcomes."
  )
}

# Output folder for per-outcome files
outcome_dir <- file.path(output_dir, "Prepared_outcomes")
dir.create(outcome_dir, recursive = TRUE, showWarnings = FALSE)

prepare_one_outcome <- function(row, input_dir) {
  file_path <- resolve_input_file(input_dir, row$file)

  if (!file.exists(file_path)) {
    stop("GWAS file not found: ", file_path)
  }

  sep_value <- if (row$sep == "\\t") "\t" else row$sep

  message("Preparing outcome: ", row$trait_id)
  message("Using file: ", file_path)

  raw_dat <- read_gwas_table(file_path, sep_value)

  out <- TwoSampleMR::format_data(
    dat = raw_dat,
    type = "outcome",
    snp_col = row$snp_col,
    beta_col = row$effect_col,
    se_col = row$se_col,
    effect_allele_col = row$ea_col,
    other_allele_col = row$oa_col,
    eaf_col = row$eaf_col,
    pval_col = row$pval_col,
    chr_col = row$chromosome_col,
    pos_col = row$position_col
  )

  out$outcome <- row$trait_id
  out
}

row <- outcome_manifest[task_id, , drop = FALSE]
df  <- prepare_one_outcome(row, input_dir)

outfile <- file.path(outcome_dir, paste0(row$trait_id, "_outcome_data.txt"))

write.table(
  df,
  file = outfile,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Outcome preparation complete for: ", row$trait_id)
message("Written to: ", outfile)
