############################################################
# File: Scripts/04b_combine_proxy_candidates.R
# Purpose: Combine proxy candidate files across all exposure-outcome pairs
# Usage:
#   Rscript Scripts/04b_combine_proxy_candidates.R <output_dir>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Scripts/04b_combine_proxy_candidates.R <output_dir>")
}

output_dir <- args[1]

source("Utils/Packages.R")

# -------------------------
# Input/output directories
# -------------------------

proxy_dir   <- file.path(output_dir, "Proxy_candidates", "per_pair")
summary_dir <- file.path(output_dir, "proxy_candidate_summaries")

combined_dir <- file.path(output_dir, "Proxy_candidates")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(proxy_dir)) {
  stop("Proxy candidate directory does not exist: ", proxy_dir)
}

# -------------------------
# Helper functions
# -------------------------

read_many_tables <- function(files) {

  if (length(files) == 0) {
    return(data.frame())
  }

  dat_list <- lapply(files, function(f) {
    utils::read.table(
      f,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE,
      colClasses = "character"
    )
  })

  dplyr::bind_rows(dat_list)
}

write_table <- function(dat, file) {
  write.table(
    dat,
    file = file,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

# -------------------------
# Locate files
# -------------------------

candidate_files <- list.files(
  proxy_dir,
  pattern = "__proxy_candidates\\.tsv$",
  full.names = TRUE
)

missing_files <- list.files(
  proxy_dir,
  pattern = "__missing_snps\\.tsv$",
  full.names = TRUE
)

summary_files <- character()

if (dir.exists(summary_dir)) {
  summary_files <- list.files(
    summary_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
}

if (length(candidate_files) == 0) {
  stop("No proxy candidate files found in: ", proxy_dir)
}

message("Combining ", length(candidate_files), " proxy candidate files...")

# -------------------------
# Combine files
# -------------------------

proxy_candidates <- read_many_tables(candidate_files)
missing_snps     <- read_many_tables(missing_files)
proxy_summary    <- read_many_tables(summary_files)
numeric_cols <- c(
  "R",
  "R2",
  "original_exposure_snps",
  "outcome_snps",
  "missing_in_outcome",
  "raw_proxy_rows",
  "proxy_rows_in_outcome",
  "missing_with_candidate_proxy"
)

for (col in numeric_cols) {
  if (col %in% names(proxy_candidates)) {
    proxy_candidates[[col]] <- suppressWarnings(as.numeric(proxy_candidates[[col]]))
  }
  if (col %in% names(proxy_summary)) {
    proxy_summary[[col]] <- suppressWarnings(as.numeric(proxy_summary[[col]]))
  }
}
if (nrow(proxy_candidates) > 0) {

  required_cols <- c("exposure", "outcome", "missing_snp", "proxy_snp")

  if (!all(required_cols %in% names(proxy_candidates))) {
    stop(
      "Combined proxy candidate table is missing required columns: ",
      paste(setdiff(required_cols, names(proxy_candidates)), collapse = ", ")
    )
  }

  distinct_cols <- c("exposure", "outcome", "missing_snp", "proxy_snp")
  proxy_candidates <- proxy_candidates[
    !duplicated(proxy_candidates[distinct_cols]),
    ,
    drop = FALSE
  ]
}

if (nrow(missing_snps) > 0) {
  required_missing_cols <- c("exposure", "outcome", "missing_snp")

  if (all(required_missing_cols %in% names(missing_snps))) {
    missing_snps <- missing_snps[
      !duplicated(missing_snps[required_missing_cols]),
      ,
      drop = FALSE
    ]
  }
}

if (nrow(proxy_summary) > 0) {
  required_summary_cols <- c("exposure", "outcome")

  if (all(required_summary_cols %in% names(proxy_summary))) {
    proxy_summary <- proxy_summary[
      !duplicated(proxy_summary[required_summary_cols]),
      ,
      drop = FALSE
    ]
  }
}

# -------------------------
# Save combined outputs
# -------------------------

candidate_out <- file.path(combined_dir, "proxy_candidates.tsv")
missing_out   <- file.path(combined_dir, "missing_snps.tsv")
summary_out   <- file.path(combined_dir, "proxy_discovery_summary.tsv")

write_table(proxy_candidates, candidate_out)
write_table(missing_snps, missing_out)
write_table(proxy_summary, summary_out)

message("Combined proxy candidate file written to: ", candidate_out)
message("Combined missing SNP file written to: ", missing_out)
message("Combined proxy discovery summary written to: ", summary_out)

message("Proxy candidate rows: ", nrow(proxy_candidates))
message("Missing SNP rows: ", nrow(missing_snps))
message("Summary rows: ", nrow(proxy_summary))

if (nrow(proxy_summary) > 0 && "status" %in% names(proxy_summary)) {
  completed <- sum(proxy_summary$status == "completed", na.rm = TRUE)
  failed    <- sum(proxy_summary$status == "failed", na.rm = TRUE)

  message("Proxy searches completed: ", completed)
  message("Proxy searches failed: ", failed)
}
