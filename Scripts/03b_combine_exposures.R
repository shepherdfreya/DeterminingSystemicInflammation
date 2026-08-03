############################################################
# File: Scripts/03b_combine_exposures.R
# Purpose: Combine per-exposure files produced by the array job
# Usage:
#   Rscript Scripts/03b_combine_exposures.R <output_dir>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Scripts/03b_combine_exposures.R <output_dir>")
}

output_dir <- args[1]

source("Utils/Packages.R")

exposure_dir <- file.path(output_dir, "Prepared_exposures")
summary_dir  <- file.path(output_dir, "per_task_summaries")

if (!dir.exists(exposure_dir)) {
  stop("Prepared_exposures directory does not exist: ", exposure_dir)
}

files <- list.files(
  exposure_dir,
  pattern = "_exposure_data\\.txt$",
  full.names = TRUE
)

if (length(files) == 0) {
  stop("No prepared exposure files found in: ", exposure_dir)
}

message("Combining ", length(files), " prepared exposure files...")

exposure_list <- lapply(files, function(f) {
  utils::read.table(
    f,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})

all_exposures <- dplyr::bind_rows(exposure_list)

combined_file <- file.path(output_dir, "all_exposures.txt")

write.table(
  all_exposures,
  file = combined_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Combined exposure file written to: ", combined_file)

# Optionally combine per-task summaries too
if (dir.exists(summary_dir)) {
  summary_files <- list.files(
    summary_dir,
    pattern = "\\.tsv$",
    full.names = TRUE
  )
  
  if (length(summary_files) > 0) {
    summary_list <- lapply(summary_files, function(f) {
      utils::read.table(
        f,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    })
    
    processing_summary <- dplyr::bind_rows(summary_list)
    
    summary_file <- file.path(output_dir, "processing_summary_exposures.txt")
    
    write.table(
      processing_summary,
      file = summary_file,
      sep = "\t",
      row.names = FALSE,
      quote = FALSE
    )
    
    message("Exposure processing summary written to: ", summary_file)
    
    successful <- sum(processing_summary$status == "completed" & processing_summary$final_snps > 0, na.rm = TRUE)
    failed     <- sum(processing_summary$status == "failed", na.rm = TRUE)
    no_snps    <- sum(processing_summary$status == "completed" & processing_summary$final_snps == 0, na.rm = TRUE)
    
    message("Files successfully processed: ", successful)
    message("Files failed: ", failed)
    message("Files with no significant SNPs: ", no_snps)
    message("Total files attempted: ", nrow(processing_summary))
  }
}