############################################################
# File: Scripts/04a_find_proxy_candidates_array.R
# Purpose: Find LD proxy candidates for one exposure-outcome pair
# Order:
#   1. Read exposure + outcome data
#   2. Identify missing SNPs
#   3. Select correct LD reference panel
#   4. Find LD proxies
#   5. Keep proxies present in outcome GWAS
#   6. Save proxy candidate table
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript Scripts/04a_find_proxy_candidates_array.R <input_dir> <output_dir> <task_id>")
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

manifest <- read.csv(
  file.path(input_dir, "trait_manifest.csv"),
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

get_ref_panel_prefix <- function(build) {
  build <- trimws(as.character(build))
  
  if (!build %in% names(REF_PANEL_PREFIXES)) {
    stop(
      "Unsupported genome build '", build, "'. Supported builds are: ",
      paste(names(REF_PANEL_PREFIXES), collapse = ", ")
    )
  }
  
  ref_prefix <- REF_PANEL_PREFIXES[[build]]
  
  needed <- paste0(ref_prefix, c(".bed", ".bim", ".fam"))
  missing <- needed[!file.exists(needed)]
  
  if (length(missing) > 0) {
    stop(
      "Reference panel files not found for build '", build, "'. Missing: ",
      paste(basename(missing), collapse = ", "),
      ". Expected prefix: ", ref_prefix
    )
  }
  
  ref_prefix
}

write_per_task_summary <- function(summary_row, output_dir, exposure_name, outcome_name) {
  summary_dir <- file.path(output_dir, "proxy_candidate_summaries")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  
  out <- file.path(
    summary_dir,
    paste0(make.names(exposure_name), "__", make.names(outcome_name), ".tsv")
  )
  
  write.table(
    summary_row,
    out,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

empty_proxy_candidates <- function() {
  data.frame(
    exposure = character(),
    outcome = character(),
    build = character(),
    missing_snp = character(),
    proxy_snp = character(),
    R = numeric(),
    R2 = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_missing_snps <- function() {
  data.frame(
    exposure = character(),
    outcome = character(),
    build = character(),
    missing_snp = character(),
    stringsAsFactors = FALSE
  )
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

# -------------------------
# Select exposure-outcome pair for task
# -------------------------

if (!all(c("exposure", "outcome", "run") %in% names(plan))) {
  stop("analysis_plan.csv must contain columns: exposure, outcome, run")
}

if (!all(c("trait_id", "build") %in% names(manifest))) {
  stop("trait_manifest.csv must contain columns: trait_id, build")
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

missing_traits <- setdiff(
  unique(c(exposure_id, outcome_id)),
  manifest$trait_id
)

if (length(missing_traits) > 0) {
  stop(
    "The following trait_id values were found in analysis_plan.csv but not trait_manifest.csv: ",
    paste(missing_traits, collapse = ", ")
  )
}

exposure_manifest <- manifest[match(exposure_id, manifest$trait_id), , drop = FALSE]
build_val <- exposure_manifest$build[[1]]

safe_pair <- paste0(make.names(exposure_id), "__", make.names(outcome_id))

proxy_dir <- file.path(output_dir, "Proxy_candidates", "per_pair")
dir.create(proxy_dir, recursive = TRUE, showWarnings = FALSE)

out_candidates <- file.path(proxy_dir, paste0(safe_pair, "__proxy_candidates.tsv"))
out_missing    <- file.path(proxy_dir, paste0(safe_pair, "__missing_snps.tsv"))
done_file      <- file.path(proxy_dir, paste0(safe_pair, ".__DONE"))

result_row <- data.frame(
  exposure = exposure_id,
  outcome = outcome_id,
  build = build_val,
  status = "failed",
  err_message = "",
  original_exposure_snps = 0L,
  outcome_snps = 0L,
  missing_in_outcome = 0L,
  raw_proxy_rows = 0L,
  proxy_rows_in_outcome = 0L,
  missing_with_candidate_proxy = 0L,
  stringsAsFactors = FALSE
)

if (file.exists(done_file) && file.exists(out_candidates)) {
  message("Already completed: ", exposure_id, " -> ", outcome_id)
  quit(save = "no", status = 0)
}

message("========================================")
message("Finding proxy candidates")
message("Exposure    : ", exposure_id)
message("Outcome     : ", outcome_id)
message("Genome build: ", build_val)
message("========================================")

# -------------------------
# Input files
# -------------------------

exposure_file <- file.path(output_dir, "all_exposures.txt")
outcome_file  <- file.path(output_dir, "all_outcomes.txt")

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

exp_df <- exposures %>%
  dplyr::filter(.data$exposure == exposure_id) %>%
  dplyr::distinct(.data$SNP, .keep_all = TRUE)

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
  missing_tbl <- empty_missing_snps()
  candidates <- empty_proxy_candidates()
  
  result_row$status <- "completed"
  result_row$err_message <- ""
  
  safe_write_table(missing_tbl, out_missing)
  safe_write_table(candidates, out_candidates)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  
  file.create(done_file)
  
  message("No missing SNPs for: ", exposure_id, " -> ", outcome_id)
  quit(save = "no", status = 0)
}

missing_tbl <- data.frame(
  exposure = exposure_id,
  outcome = outcome_id,
  build = build_val,
  missing_snp = missing_snps,
  stringsAsFactors = FALSE
)

safe_write_table(missing_tbl, out_missing)

# -------------------------
# Find LD proxies
# -------------------------

ref_prefix <- get_ref_panel_prefix(build_val)

gwasvcf::set_plink(PLINK_BIN)

ld_raw <- tryCatch(
  {
    gwasvcf::get_ld_proxies(
      rsid = missing_snps,
      bfile = ref_prefix,
      tag_kb = PROXY_TAG_KB,
      tag_nsnp = PROXY_TAG_NSNP,
      tag_r2 = PROXY_TAG_R2,
      threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")),
      out = tempfile()
    )
  },
  error = function(e) {
    result_row$err_message <<- paste("LD proxy search failed:", conditionMessage(e))
    warning(result_row$err_message)
    data.frame()
  }
)

result_row$raw_proxy_rows <- nrow(ld_raw)

if (nrow(ld_raw) == 0) {
  candidates <- empty_proxy_candidates()
  
  result_row$status <- "completed"
  if (result_row$err_message == "") {
    result_row$err_message <- "No LD proxies found"
  }
  
  safe_write_table(candidates, out_candidates)
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  
  file.create(done_file)
  
  message("No LD proxies found for: ", exposure_id, " -> ", outcome_id)
  quit(save = "no", status = 0)
}

if (!all(c("SNP_A", "SNP_B") %in% names(ld_raw))) {
  result_row$err_message <- paste(
    "LD proxy output is missing expected columns:",
    paste(setdiff(c("SNP_A", "SNP_B"), names(ld_raw)), collapse = ", ")
  )
  write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)
  stop(result_row$err_message)
}

ld_raw <- as.data.frame(ld_raw)

if (!"R2" %in% names(ld_raw)) {
  if ("R" %in% names(ld_raw)) {
    ld_raw$R2 <- suppressWarnings(as.numeric(ld_raw$R)^2)
  } else {
    ld_raw$R2 <- NA_real_
  }
}

if (!"R" %in% names(ld_raw)) {
  ld_raw$R <- NA_real_
}

candidates <- ld_raw %>%
  dplyr::rename(
    missing_snp = "SNP_A",
    proxy_snp = "SNP_B"
  ) %>%
  dplyr::mutate(
    exposure = exposure_id,
    outcome = outcome_id,
    build = build_val,
    R = suppressWarnings(as.numeric(.data$R)),
    R2 = suppressWarnings(as.numeric(.data$R2))
  ) %>%
  dplyr::filter(.data$proxy_snp %in% out_df$SNP) %>%
  dplyr::select(
    exposure,
    outcome,
    build,
    missing_snp,
    proxy_snp,
    R,
    R2,
    dplyr::everything()
  )

result_row$proxy_rows_in_outcome <- nrow(candidates)
result_row$missing_with_candidate_proxy <- length(unique(candidates$missing_snp))
result_row$status <- "completed"
result_row$err_message <- ""

safe_write_table(candidates, out_candidates)
write_per_task_summary(result_row, output_dir, exposure_id, outcome_id)

file.create(done_file)

message("Proxy candidate search complete for: ", exposure_id, " -> ", outcome_id)
message("Missing SNPs: ", length(missing_snps))
message("Proxy rows in outcome: ", nrow(candidates))