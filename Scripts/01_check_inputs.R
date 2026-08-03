############################################################
# File: Scripts/01_check_inputs.R
# Purpose: Validate trait manifest, analysis plan, and GWAS file headers.
#          If effect_col == OR, add a beta column to that same GWAS file.
# Usage: Rscript Scripts/01_check_inputs.R <input_dir> <output_dir>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript Scripts/01_check_inputs.R <input_dir> <output_dir>")
}

input_dir  <- args[1]
output_dir <- args[2]

source("Utils/Packages.R")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest_path <- file.path(input_dir, "trait_manifest.csv")
plan_path     <- file.path(input_dir, "analysis_plan.csv")

if (!file.exists(manifest_path)) {
  stop("trait_manifest.csv not found in: ", input_dir)
}

if (!file.exists(plan_path)) {
  stop("analysis_plan.csv not found in: ", input_dir)
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
plan <- read.csv(plan_path, stringsAsFactors = FALSE, check.names = FALSE)

required_manifest_cols <- c(
  "trait_id",
  "file",
  "sep",
  "snp_col",
  "position_col",
  "effect_col",
  "se_col",
  "ea_col",
  "oa_col",
  "pval_col",
  "build",
  "eaf_col",
  "chromosome_col",
  "protein_start_col",
  "protein_end_col",
  "protein_chr_col"
)

required_plan_cols <- c(
  "exposure",
  "outcome",
  "run"
)

missing_manifest_cols <- setdiff(required_manifest_cols, names(manifest))
missing_plan_cols <- setdiff(required_plan_cols, names(plan))

if (length(missing_manifest_cols) > 0) {
  stop(
    "Missing required columns in trait_manifest.csv: ",
    paste(missing_manifest_cols, collapse = ", ")
  )
}

if (length(missing_plan_cols) > 0) {
  stop(
    "Missing required columns in analysis_plan.csv: ",
    paste(missing_plan_cols, collapse = ", ")
  )
}

if (anyDuplicated(manifest$trait_id) > 0) {
  dup_ids <- unique(manifest$trait_id[duplicated(manifest$trait_id)])
  stop("Duplicate trait_id values found in trait_manifest.csv: ",
       paste(dup_ids, collapse = ", "))
}

if (anyDuplicated(manifest$file) > 0) {
  dup_files <- unique(manifest$file[duplicated(manifest$file)])
  stop("Duplicate file values found in trait_manifest.csv: ",
       paste(dup_files, collapse = ", "))
}

plan$run <- as.character(plan$run)
valid_run_values <- c("TRUE", "FALSE", "T", "F", "true", "false", "True", "False")

if (any(!plan$run %in% valid_run_values)) {
  bad_run <- unique(plan$run[!plan$run %in% valid_run_values])
  stop("Invalid values in analysis_plan.csv 'run' column: ",
       paste(bad_run, collapse = ", "),
       ". Use TRUE or FALSE.")
}

plan$run <- tolower(plan$run) %in% c("true", "t")

if ("Analysis" %in% names(plan)) {
  plan$Analysis <- tolower(trimws(as.character(plan$Analysis)))
  valid_analysis_values <- c("cis_only", "cis_trans")

  if (any(!plan$Analysis %in% valid_analysis_values)) {
    bad_analysis <- unique(plan$Analysis[!plan$Analysis %in% valid_analysis_values])
    stop(
      "Invalid values in analysis_plan.csv 'Analysis' column: ",
      paste(bad_analysis, collapse = ", "),
      ". Use cis_only or cis_trans."
    )
  }
}

traits_in_manifest <- unique(manifest$trait_id)
traits_in_plan <- unique(c(plan$exposure, plan$outcome))

missing_traits <- setdiff(traits_in_plan, traits_in_manifest)

if (length(missing_traits) > 0) {
  stop(
    "These traits are listed in analysis_plan.csv but not in trait_manifest.csv: ",
    paste(missing_traits, collapse = ", ")
  )
}

is_provided <- function(x) {
  !(is.na(x) || trimws(x) == "" || trimws(x) %in% c("NA", "na", "Na", "n/a", "N/A"))
}

validation_results <- lapply(seq_len(nrow(manifest)), function(i) {
  row <- manifest[i, ]

  file_path <- file.path(input_dir, row$file)

  if (!file.exists(file_path)) {
    return(data.frame(
      trait_id = row$trait_id,
      file = row$file,
      file_exists = FALSE,
      header_read = FALSE,
      missing_columns = "FILE NOT FOUND",
      status = "FAIL",
      stringsAsFactors = FALSE
    ))
  }

  sep_value <- if (row$sep == "\\t") "\t" else row$sep

  header <- tryCatch(
    names(data.table::fread(
      file_path,
      nrows = 0,
      sep = sep_value,
      data.table = FALSE,
      showProgress = FALSE
    )),
    error = function(e) NULL
  )

  if (is.null(header)) {
    return(data.frame(
      trait_id = row$trait_id,
      file = row$file,
      file_exists = TRUE,
      header_read = FALSE,
      missing_columns = "FAILED TO READ HEADER",
      status = "FAIL",
      stringsAsFactors = FALSE
    ))
  }

  required_mapped_cols <- c(
    row$snp_col,
    row$position_col,
    row$effect_col,
    row$se_col,
    row$ea_col,
    row$oa_col,
    row$pval_col,
    row$eaf_col,
    row$chromosome_col
  )

  required_mapped_cols <- unique(required_mapped_cols[vapply(required_mapped_cols, is_provided, logical(1))])
  missing_file_cols <- setdiff(required_mapped_cols, header)
  status <- if (length(missing_file_cols) == 0) "OK" else "FAIL"

  data.frame(
    trait_id = row$trait_id,
    file = row$file,
    file_exists = TRUE,
    header_read = TRUE,
    missing_columns = if (length(missing_file_cols) == 0) "" else paste(missing_file_cols, collapse = "; "),
    status = status,
    stringsAsFactors = FALSE
  )
})

validation_results <- do.call(rbind, validation_results)

summary_path <- file.path(output_dir, "input_validation_summary.csv")
write.csv(validation_results, summary_path, row.names = FALSE)

active_plan <- subset(plan, run)
active_plan_path <- file.path(output_dir, "active_analysis_plan.csv")
write.csv(active_plan, active_plan_path, row.names = FALSE)

message("Input validation complete.")
message("Validation summary written to: ", summary_path)
message("Active analysis plan written to: ", active_plan_path)

if (any(validation_results$status == "FAIL")) {
  failed_traits <- validation_results$trait_id[validation_results$status == "FAIL"]
  message("Validation failed for: ", paste(failed_traits, collapse = ", "))
  stop("One or more input files failed validation. See input_validation_summary.csv for details.")
} else {
  message("All listed GWAS files passed header validation.")
}

# ---------------------------------------------------------
# Add beta column to raw GWAS files where effect_col == OR
# ---------------------------------------------------------

add_beta_if_or <- function(manifest_row, input_dir) {
  effect_name <- trimws(manifest_row$effect_col)

  if (tolower(effect_name) != "or") {
    return(
      data.frame(
        trait_id = manifest_row$trait_id,
        file = manifest_row$file,
        action = "No change",
        details = paste0("effect_col = ", effect_name),
        stringsAsFactors = FALSE
      )
    )
  }

  file_path <- file.path(input_dir, manifest_row$file)
  sep_value <- if (manifest_row$sep == "\\t") "\t" else manifest_row$sep

  message("Checking OR-to-beta conversion for: ", manifest_row$trait_id)

  dt <- data.table::fread(
    file_path,
    sep = sep_value,
    data.table = FALSE,
    showProgress = FALSE
  )

  if ("beta" %in% names(dt)) {
    return(
      data.frame(
        trait_id = manifest_row$trait_id,
        file = manifest_row$file,
        action = "Skipped",
        details = "beta column already present; standardised manifest will use effect_col = beta",
        stringsAsFactors = FALSE
      )
    )
  }

  if (!(effect_name %in% names(dt))) {
    stop("Effect column '", effect_name, "' not found in file: ", manifest_row$file)
  }

  # Create backup
  backup_path <- paste0(file_path, ".bak")
  if (!file.exists(backup_path)) {
    file.copy(file_path, backup_path, overwrite = FALSE)
  }

  dt$beta <- log(as.numeric(dt[[effect_name]]))

  # Overwrite original file
  data.table::fwrite(
    dt,
    file = file_path,
    sep = sep_value,
    quote = FALSE,
    na = "NA"
  )

  data.frame(
    trait_id = manifest_row$trait_id,
    file = manifest_row$file,
    action = "Updated",
    details = "Added beta column as log(OR)",
    stringsAsFactors = FALSE
  )
}

beta_update_results <- do.call(
  rbind,
  lapply(seq_len(nrow(manifest)), function(i) add_beta_if_or(manifest[i, ], input_dir))
)

beta_update_path <- file.path(output_dir, "beta_column_updates.csv")
write.csv(beta_update_results, beta_update_path, row.names = FALSE)

standardised_manifest <- manifest
or_rows <- tolower(trimws(standardised_manifest$effect_col)) == "or"
standardised_manifest$effect_col[or_rows] <- "beta"

standardised_manifest_path <- file.path(output_dir, "trait_manifest_standardised.csv")
write.csv(standardised_manifest, standardised_manifest_path, row.names = FALSE)

message("Beta-column update summary written to: ", beta_update_path)
message("Standardised trait manifest written to: ", standardised_manifest_path)
