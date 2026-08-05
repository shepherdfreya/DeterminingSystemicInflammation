############################################################
# File: Scripts/05b_combine_mr_outputs.R
# Purpose: Combine pair-specific MR summaries and MR result files
# Usage:
#   Rscript Scripts/05b_combine_mr_outputs.R <output_dir>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Scripts/05b_combine_mr_outputs.R <output_dir>")
}

output_dir <- args[1]

source("Utils/Packages.R")

mr_summary_dir <- file.path(output_dir, "mr_summaries")
mr_results_dir <- file.path(output_dir, "MR_results", "per_pair")

read_many_tables <- function(files) {
  if (length(files) == 0) {
    return(data.frame())
  }

  dat_list <- lapply(files, function(f) {
    if (file.info(f)$size <= 1) {
      return(NULL)
    }

    tryCatch(
      utils::read.table(
        f,
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      error = function(e) {
        warning("Skipping unreadable table: ", f, " (", conditionMessage(e), ")")
        NULL
      }
    )
  })

  dat_list <- Filter(Negate(is.null), dat_list)

  if (length(dat_list) == 0) {
    return(data.frame())
  }

  dplyr::bind_rows(dat_list)
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

if (!dir.exists(mr_summary_dir)) {
  stop("MR summary directory does not exist: ", mr_summary_dir)
}

summary_files <- list.files(
  mr_summary_dir,
  pattern = "\\.tsv$",
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No MR summary files found in: ", mr_summary_dir)
}

message("Combining ", length(summary_files), " MR summary files...")

mr_summary <- read_many_tables(summary_files)
mr_summary_out <- file.path(output_dir, "mr_summary.tsv")
safe_write_table(mr_summary, mr_summary_out)

message("Combined MR summary written to: ", mr_summary_out)
message("MR summary rows: ", nrow(mr_summary))

if (dir.exists(mr_results_dir)) {
  result_files <- list.files(
    mr_results_dir,
    pattern = "__mr_results\\.tsv$",
    full.names = TRUE
  )

  if (length(result_files) > 0) {
    message("Combining ", length(result_files), " MR result files...")

    mr_results <- read_many_tables(result_files)
    mr_results_out <- file.path(output_dir, "mr_results.tsv")
    safe_write_table(mr_results, mr_results_out)

    message("Combined MR results written to: ", mr_results_out)
    message("MR result rows: ", nrow(mr_results))
  }
}
