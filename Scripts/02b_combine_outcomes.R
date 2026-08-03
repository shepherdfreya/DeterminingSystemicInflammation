############################################################
# File: Scripts/02b_combine_outcomes.R
# Purpose: Combine per-outcome files produced by the array job
# Usage:
#   Rscript Scripts/02b_combine_outcomes.R <output_dir>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Scripts/02b_combine_outcomes.R <output_dir>")
}

output_dir <- args[1]

source("Utils/Packages.R")

outcome_dir <- file.path(output_dir, "Prepared_outcomes")

if (!dir.exists(outcome_dir)) {
  stop("Prepared_outcomes directory does not exist: ", outcome_dir)
}

files <- list.files(
  outcome_dir,
  pattern = "_outcome_data\\.txt$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No prepared outcome files found in: ", outcome_dir)
}

message("Combining ", length(files), " prepared outcome files...")

outcome_list <- lapply(files, function(f) {
  read.table(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

all_outcomes <- do.call(rbind, outcome_list)

combined_file <- file.path(output_dir, "all_outcomes.txt")

write.table(
  all_outcomes,
  file = combined_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Combined outcome file written to: ", combined_file)