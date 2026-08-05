############################################################
# File: Scripts/03a_prepare_exposures_array.R
# Purpose: Prepare ONE exposure GWAS file per array task
# Order:
#   1. Read raw GWAS
#   2. MAF filter
#   3. P-value filter
#   4. Format as exposure using TwoSampleMR
#   5. LD clump
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop("Usage: Rscript Scripts/03a_prepare_exposures_array.R <input_dir> <output_dir> <task_id>")
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

# -------------------------
# Helper functions
# -------------------------

normalise_run_flag <- function(x) {
  tolower(as.character(x)) %in% c("true", "t", "1")
}

get_sep_value <- function(x) {
  if (is.na(x) || x == "") {
    stop("Missing 'sep' value in trait_manifest.csv")
  }
  if (x == "\\t") "\t" else x
}

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

calculate_maf_and_filter <- function(dat, eaf_col, maf_threshold = 0.01) {

  if (!eaf_col %in% names(dat)) {
    stop("Column '", eaf_col, "' not found in raw exposure data")
  }

  dat[[eaf_col]] <- suppressWarnings(as.numeric(dat[[eaf_col]]))
  dat$maf <- pmin(dat[[eaf_col]], 1 - dat[[eaf_col]])

  n_before <- nrow(dat)
  dat <- dat[!is.na(dat$maf) & dat$maf >= maf_threshold, , drop = FALSE]
  n_after <- nrow(dat)

  list(
    data = dat,
    removed = n_before - n_after
  )
}

filter_significant_snps <- function(dat, pval_col, p_threshold = EXPOSURE_PVAL_THRESHOLD) {

  if (!pval_col %in% names(dat)) {
    stop("Column '", pval_col, "' not found in raw exposure data")
  }

  dat[[pval_col]] <- suppressWarnings(as.numeric(dat[[pval_col]]))
  dat <- dat[!is.na(dat[[pval_col]]), , drop = FALSE]
  dat[dat[[pval_col]] < p_threshold, , drop = FALSE]
}

calculate_f_stat_and_filter <- function(dat, f_threshold = EXPOSURE_F_THRESHOLD) {

  needed <- c("beta.exposure", "se.exposure")
  if (!all(needed %in% names(dat))) {
    stop("Formatted exposure data must contain columns: beta.exposure, se.exposure")
  }

  dat$F_stat <- suppressWarnings((as.numeric(dat$beta.exposure) / as.numeric(dat$se.exposure))^2)

  n_before <- nrow(dat)
  dat <- dat[!is.na(dat$F_stat) & dat$F_stat > f_threshold, , drop = FALSE]

  list(
    data = dat,
    removed = n_before - nrow(dat)
  )
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

perform_ld_clumping_safe <- function(df, trait_name, build, snp_col, pval_col) {

  if (nrow(df) == 0) {
    return(list(data = df, note = "No significant SNPs to clump"))
  }

  ref_prefix <- get_ref_panel_prefix(build)

  clump_input <- data.frame(
    rsid = df[[snp_col]],
    pval = df[[pval_col]],
    id   = trait_name,
    stringsAsFactors = FALSE
  )

  clumped_snps <- tryCatch(
    ieugwasr::ld_clump(
      dat = clump_input,
      plink_bin = PLINK_BIN,
      bfile = ref_prefix,
      clump_kb = EXPOSURE_CLUMP_KB,
      clump_r2 = EXPOSURE_CLUMP_R2
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      no_clumps <- grepl("No significant --clump results", msg, ignore.case = TRUE) ||
        grepl("\\.clumped", msg, ignore.case = TRUE) ||
        grepl("cannot open the connection", msg, ignore.case = TRUE)

      if (!no_clumps) {
        stop(e)
      }

      return(clump_input[0, , drop = FALSE])
    }
  )

  out <- dplyr::semi_join(df, clumped_snps, by = setNames("rsid", snp_col))

  list(data = out, note = NULL)
}
write_per_task_summary <- function(summary_row, output_dir, exposure_name) {
  summary_dir <- file.path(output_dir, "per_task_summaries")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)

  out <- file.path(summary_dir, paste0(exposure_name, ".tsv"))

  write.table(
    summary_row,
    out,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
}

# -------------------------
# Select exposure for task
# -------------------------

plan$run <- normalise_run_flag(plan$run)
plan <- subset(plan, run)

if (nrow(plan) == 0) {
  stop("No rows with run = TRUE found in analysis_plan.csv")
}

exposure_traits <- unique(plan$exposure)

missing_exposures <- exposure_traits[!exposure_traits %in% manifest$trait_id]

if (length(missing_exposures) > 0) {
  stop(
    "The following exposure trait_id values were found in analysis_plan.csv but not trait_manifest.csv: ",
    paste(missing_exposures, collapse = ", ")
  )
}

exposure_manifest <- manifest[match(exposure_traits, manifest$trait_id), , drop = FALSE]

if (task_id > nrow(exposure_manifest)) {
  stop(
    "task_id ", task_id, " is out of range. There are only ",
    nrow(exposure_manifest), " unique active exposures."
  )
}

row <- exposure_manifest[task_id, , drop = FALSE]

trait_id  <- row$trait_id[[1]]
file_path <- resolve_input_file(input_dir, row$file[[1]])
sep_value <- get_sep_value(row$sep[[1]])
build_val <- row$build[[1]]

exposure_dir <- file.path(output_dir, "Prepared_exposures")
dir.create(exposure_dir, recursive = TRUE, showWarnings = FALSE)

result_row <- data.frame(
  trait_id = trait_id,
  file_name = row$file[[1]],
  status = "failed",
  err_message = "",
  rows_processed = 0L,
  maf_filtered_snps = 0L,
  significant_snps = 0L,
  clumped_snps = 0L,
  f_filtered_snps = 0L,
  final_snps = 0L,
  stringsAsFactors = FALSE
)

message("========================================")
message("Preparing exposure: ", trait_id)
message("Input file        : ", file_path)
message("Genome build      : ", build_val)
message("========================================")

if (!file.exists(file_path)) {
  result_row$err_message <- paste("GWAS file not found:", file_path)
  write_per_task_summary(result_row, output_dir, trait_id)
  stop(result_row$err_message)
}

# 1. Read raw GWAS
raw_dat <- read_gwas_table(file_path, sep_value)
result_row$rows_processed <- nrow(raw_dat)

# 2. MAF filter
maf_res <- calculate_maf_and_filter(
  dat = raw_dat,
  eaf_col = row$eaf_col[[1]],
  maf_threshold = EXPOSURE_MAF_THRESHOLD
)

raw_dat_maf <- maf_res$data
result_row$maf_filtered_snps <- maf_res$removed

# 3. P-value filter
sig_raw_dat <- filter_significant_snps(
  dat = raw_dat_maf,
  pval_col = row$pval_col[[1]],
  p_threshold = EXPOSURE_PVAL_THRESHOLD
)
result_row$significant_snps <- nrow(sig_raw_dat)

# 4. LD clumping on raw filtered data
clump_res <- perform_ld_clumping_safe(
  df = sig_raw_dat,
  trait_name = trait_id,
  build = build_val,
  snp_col = row$snp_col[[1]],
  pval_col = row$pval_col[[1]]
)

clumped_raw_dat <- clump_res$data
result_row$clumped_snps <- nrow(clumped_raw_dat)

outfile <- file.path(exposure_dir, paste0(trait_id, "_exposure_data.txt"))

if (nrow(clumped_raw_dat) == 0) {
  empty_exp <- data.frame(
    SNP = character(),
    `chr.exposure` = character(),
    `pos.exposure` = numeric(),
    `effect_allele.exposure` = character(),
    `other_allele.exposure` = character(),
    `eaf.exposure` = numeric(),
    `beta.exposure` = numeric(),
    `se.exposure` = numeric(),
    `pval.exposure` = numeric(),
    exposure = character(),
    `mr_keep.exposure` = logical(),
    `pval_origin.exposure` = character(),
    `id.exposure` = character(),
    `data_source.exposure` = character(),
    F_stat = numeric(),
    check.names = FALSE
  )

  write.table(
    empty_exp,
    file = outfile,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )

  result_row$status <- "completed"
  result_row$err_message <- "No SNPs remained after LD clumping"
  result_row$final_snps <- 0L
  write_per_task_summary(result_row, output_dir, trait_id)

  message("No SNPs remained after LD clumping for: ", trait_id)
  message("Written empty exposure file to: ", outfile)
  quit(save = "no", status = 0)
}

# 5. Format final clumped data with TwoSampleMR
tmp_clumped_file <- file.path(
  tempdir(),
  paste0(trait_id, "_clumped_raw.txt")
)

data.table::fwrite(
  clumped_raw_dat,
  file = tmp_clumped_file,
  sep = sep_value,
  quote = FALSE,
  na = "NA"
)

exp_dat <- TwoSampleMR::read_exposure_data(
  filename = tmp_clumped_file,
  sep = sep_value,
  snp_col = row$snp_col[[1]],
  beta_col = row$effect_col[[1]],
  se_col = row$se_col[[1]],
  effect_allele_col = row$ea_col[[1]],
  other_allele_col = row$oa_col[[1]],
  eaf_col = row$eaf_col[[1]],
  pval_col = row$pval_col[[1]],
  chr_col = row$chromosome_col[[1]],
  pos_col = row$position_col[[1]],
  phenotype_col = "exposure",
  clump = FALSE
)

exp_dat$exposure <- trait_id
exp_dat$id.exposure <- trait_id

f_res <- calculate_f_stat_and_filter(exp_dat)
exp_dat <- f_res$data
result_row$f_filtered_snps <- f_res$removed

write.table(
  exp_dat,
  file = outfile,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
result_row$final_snps <- nrow(exp_dat)
result_row$status <- "completed"

result_row$err_message <- ""

write_per_task_summary(result_row, output_dir, trait_id)

message("Exposure preparation complete for: ", trait_id)
message("Written to: ", outfile)
