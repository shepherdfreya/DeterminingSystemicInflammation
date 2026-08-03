############################################################
# File: Scripts/06b_combine_colocalization.R
# Purpose: Combine pair-specific colocalization outputs
# Usage:
#   Rscript Scripts/06b_combine_colocalization.R <output_dir>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Scripts/06b_combine_colocalization.R <output_dir>")
}

output_dir <- args[1]

source("Utils/Packages.R")

coloc_dir <- file.path(output_dir, "Colocalization", "per_pair")

if (!dir.exists(coloc_dir)) {
  stop("Colocalization directory does not exist: ", coloc_dir)
}

files <- list.files(
  coloc_dir,
  pattern = "__colocalization\\.tsv$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No colocalization files found in: ", coloc_dir)
}

message("Combining ", length(files), " colocalization files...")

coloc_list <- lapply(files, function(f) {
  utils::read.table(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

coloc_summary <- dplyr::bind_rows(coloc_list)

out_file <- file.path(output_dir, "colocalization_summary.tsv")

write.table(
  coloc_summary,
  file = out_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Combined colocalization summary written to: ", out_file)
message("Colocalization rows: ", nrow(coloc_summary))
message("H4 pass rows: ", sum(coloc_summary$H4_flag == "Pass", na.rm = TRUE))
