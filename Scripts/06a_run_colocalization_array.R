############################################################
# File: Scripts/06a_run_colocalization_array.R
# Purpose: Run PWCoCo for one exposure-outcome-analysis task
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

get_trait_n <- function(row) {
  if ("samplesize" %in% names(row) && is_provided(row$samplesize[[1]])) {
    return(suppressWarnings(as.numeric(row$samplesize[[1]])))
  }
  NA_real_
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
      list(sep = sep_value, data.table = FALSE, showProgress = FALSE)
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

resolve_binary <- function(bin) {
  bin <- path.expand(bin)
  if (file.exists(bin)) {
    return(bin)
  }

  resolved <- Sys.which(bin)
  if (nzchar(resolved)) {
    return(unname(resolved))
  }

  stop(
    "PWCoCo binary not found. Set PWCOCO_BIN in Utils/Config.R or export PWCOCO_BIN. Current value: ",
    bin
  )
}

safe_write_table <- function(dat, file) {
  write.table(dat, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
}

empty_coloc_result <- function() {
  data.frame(
    exposure = character(),
    outcome = character(),
    analysis = character(),
    locus_id = character(),
    top_snp = character(),
    instrument_region = character(),
    MHC = character(),
    chr = character(),
    start = numeric(),
    end = numeric(),
    selected_instrument_count = integer(),
    n_overlap = integer(),
    status = character(),
    err_message = character(),
    Dataset1 = character(),
    Dataset2 = character(),
    SNP1 = character(),
    SNP2 = character(),
    nsnps = integer(),
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

has_protein_coordinates <- function(row) {
  is_provided(row$protein_chr_col[[1]]) &&
    is_provided(row$protein_start_col[[1]]) &&
    is_provided(row$protein_end_col[[1]])
}

annotate_instrument_region <- function(exp_dat, exposure_row) {
  if (!has_protein_coordinates(exposure_row)) {
    exp_dat$instrument_region <- "unknown"
    return(exp_dat)
  }

  needed <- c("chr.exposure", "pos.exposure")
  if (!all(needed %in% names(exp_dat))) {
    stop("Exposure data must contain columns: ", paste(needed, collapse = ", "))
  }

  protein_chr <- as.character(exposure_row$protein_chr_col[[1]])
  start <- suppressWarnings(as.numeric(exposure_row$protein_start_col[[1]])) - PWCOCO_WINDOW
  end <- suppressWarnings(as.numeric(exposure_row$protein_end_col[[1]])) + PWCOCO_WINDOW

  exp_dat %>%
    dplyr::mutate(
      instrument_region = dplyr::if_else(
        as.character(.data$chr.exposure) == protein_chr &
          suppressWarnings(as.numeric(.data$pos.exposure)) >= start &
          suppressWarnings(as.numeric(.data$pos.exposure)) <= end,
        "C",
        "T"
      )
    )
}

select_pwcoco_instruments <- function(exp_dat, exposure_row, analysis_mode) {
  exp_dat <- exp_dat %>% dplyr::distinct(.data$SNP, .keep_all = TRUE)
  exp_dat <- annotate_instrument_region(exp_dat, exposure_row)

  if (analysis_mode == "cis_only") {
    if (!has_protein_coordinates(exposure_row)) {
      stop("Analysis = cis_only requires protein coordinates for exposure: ", exposure_row$trait_id[[1]])
    }
    exp_dat <- exp_dat %>% dplyr::filter(.data$instrument_region == "C")
  }

  if (nrow(exp_dat) == 0) {
    return(exp_dat)
  }

  cis_instruments <- exp_dat %>% dplyr::filter(.data$instrument_region == "C")
  if (nrow(cis_instruments) > 0) {
    return(cis_instruments)
  }

  if (nrow(exp_dat) == 1) {
    return(exp_dat)
  }

  exp_dat %>%
    dplyr::mutate(p_for_selection = suppressWarnings(as.numeric(.data$pval.exposure))) %>%
    dplyr::arrange(.data$p_for_selection) %>%
    dplyr::slice(1) %>%
    dplyr::select(-"p_for_selection")
}

standardise_raw_gwas <- function(dat, row) {
  required_cols <- c(
    row$snp_col[[1]],
    row$position_col[[1]],
    row$effect_col[[1]],
    row$se_col[[1]],
    row$pval_col[[1]],
    row$eaf_col[[1]],
    row$chromosome_col[[1]],
    row$ea_col[[1]],
    row$oa_col[[1]]
  )

  missing_required <- required_cols[!required_cols %in% names(dat)]
  if (length(missing_required) > 0) {
    stop("Raw GWAS is missing required columns: ", paste(missing_required, collapse = ", "))
  }

  out <- data.frame(
    SNP = dat[[row$snp_col[[1]]]],
    chr = as.character(dat[[row$chromosome_col[[1]]]]),
    pos = suppressWarnings(as.numeric(dat[[row$position_col[[1]]]])),
    A1 = dat[[row$ea_col[[1]]]],
    A2 = dat[[row$oa_col[[1]]]],
    A1_freq = suppressWarnings(as.numeric(dat[[row$eaf_col[[1]]]])),
    beta = suppressWarnings(as.numeric(dat[[row$effect_col[[1]]]])),
    se = suppressWarnings(as.numeric(dat[[row$se_col[[1]]]])),
    p = suppressWarnings(as.numeric(dat[[row$pval_col[[1]]]])),
    stringsAsFactors = FALSE
  )

  sample_n <- get_trait_n(row)
  if (!is.na(sample_n)) {
    out$n <- sample_n
  }

  case_col <- NA_character_
  if ("ncase_col" %in% names(row) && is_provided(row$ncase_col[[1]])) {
    case_col <- row$ncase_col[[1]]
  } else if ("case_col" %in% names(row) && is_provided(row$case_col[[1]])) {
    case_col <- row$case_col[[1]]
  }

  if (!is.na(case_col) && case_col %in% names(dat)) {
    out$case <- suppressWarnings(as.numeric(dat[[case_col]]))
  }

  out <- out[!is.na(out$SNP) & out$SNP != "", , drop = FALSE]
  out <- out[!is.na(out$chr) & out$chr != "", , drop = FALSE]
  out <- out[!is.na(out$pos), , drop = FALSE]
  out <- out[!is.na(out$A1) & out$A1 != "", , drop = FALSE]
  out <- out[!is.na(out$A2) & out$A2 != "", , drop = FALSE]
  out <- out[!is.na(out$beta), , drop = FALSE]
  out <- out[!is.na(out$se) & out$se > 0, , drop = FALSE]
  out <- out[!is.na(out$p), , drop = FALSE]
  out <- out[!is.na(out$A1_freq) & out$A1_freq > 0 & out$A1_freq < 1, , drop = FALSE]
  out <- out[!duplicated(out$SNP), , drop = FALSE]

  out
}

get_pwcoco_ref_panel_prefix <- function(build, chr) {
  build <- trimws(as.character(build))

  if (!build %in% names(PWCOCO_REF_PANEL_PREFIXES)) {
    stop(
      "Unsupported genome build '", build, "'. Supported PWCoCo builds are: ",
      paste(names(PWCOCO_REF_PANEL_PREFIXES), collapse = ", ")
    )
  }

  ref_prefix <- PWCOCO_REF_PANEL_PREFIXES[[build]]
  ref_prefix <- gsub("\\{chr\\}", as.character(chr), ref_prefix)
  ref_prefix <- gsub("\\{CHR\\}", as.character(chr), ref_prefix)

  needed <- paste0(ref_prefix, c(".bed", ".bim", ".fam"))
  missing <- needed[!file.exists(needed)]
  if (length(missing) > 0) {
    stop(
      "PWCoCo reference panel files not found for build '", build, "' chromosome ",
      chr, ". Missing: ", paste(basename(missing), collapse = ", "),
      ". Expected prefix: ", ref_prefix
    )
  }

  ref_prefix
}

make_loci <- function(selected_instruments) {
  selected_instruments %>%
    dplyr::mutate(
      chr = as.character(.data$chr.exposure),
      top_pos = suppressWarnings(as.numeric(.data$pos.exposure)),
      start = pmax(0, .data$top_pos - PWCOCO_WINDOW),
      end = .data$top_pos + PWCOCO_WINDOW,
      MHC = dplyr::if_else(
        .data$chr == MHC_CHR & .data$top_pos > MHC_START & .data$top_pos < MHC_END,
        "caution_MHC",
        "no_MHC"
      ),
      locus_id = paste0(.data$exposure, ":", .data$SNP, ":", .data$chr, ":", .data$top_pos)
    ) %>%
    dplyr::select(
      "locus_id",
      top_snp = "SNP",
      "instrument_region",
      "MHC",
      "chr",
      "start",
      "end"
    )
}

write_pwcoco_sumstats <- function(dat, out_file) {
  cols <- c("SNP", "A1", "A2", "A1_freq", "beta", "se", "p")
  if ("n" %in% names(dat)) {
    cols <- c(cols, "n")
  }
  if ("case" %in% names(dat)) {
    cols <- c(cols, "case")
  }
  safe_write_table(dat[, cols, drop = FALSE], out_file)
}

read_pwcoco_output <- function(coloc_file) {
  if (!file.exists(coloc_file) || file.info(coloc_file)$size <= 1) {
    return(data.frame())
  }

  utils::read.table(
    coloc_file,
    header = TRUE,
    sep = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

normalise_pp <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  ifelse(!is.na(x) & x > 1, x / PWCOCO_PPH4_SCALE, x)
}

run_pwcoco_for_locus <- function(locus, exp_raw, out_raw, exposure_id, outcome_id, analysis_mode, build, pwcoco_bin, work_dir) {
  result <- empty_coloc_result()
  result[1, ] <- NA

  result$exposure <- exposure_id
  result$outcome <- outcome_id
  result$analysis <- analysis_mode
  result$locus_id <- locus$locus_id[[1]]
  result$top_snp <- locus$top_snp[[1]]
  result$instrument_region <- locus$instrument_region[[1]]
  result$MHC <- locus$MHC[[1]]
  result$chr <- locus$chr[[1]]
  result$start <- locus$start[[1]]
  result$end <- locus$end[[1]]
  result$selected_instrument_count <- 1L
  result$n_overlap <- 0L
  result$status <- "failed"
  result$err_message <- ""
  result$H4_threshold <- PWCOCO_H4_THRESHOLD
  result$H4_flag <- "Fail"

  exp_region <- exp_raw %>%
    dplyr::filter(.data$chr == locus$chr[[1]], .data$pos >= locus$start[[1]], .data$pos <= locus$end[[1]])

  out_region <- out_raw %>%
    dplyr::filter(.data$chr == locus$chr[[1]], .data$pos >= locus$start[[1]], .data$pos <= locus$end[[1]])

  overlapping_snps <- intersect(exp_region$SNP, out_region$SNP)
  result$n_overlap <- length(overlapping_snps)

  if (length(overlapping_snps) < COLOC_MIN_OVERLAP) {
    result$err_message <- paste0(
      "Fewer than COLOC_MIN_OVERLAP overlapping SNPs: ",
      length(overlapping_snps),
      " < ",
      COLOC_MIN_OVERLAP
    )
    return(result)
  }

  exp_region <- exp_region[exp_region$SNP %in% overlapping_snps, , drop = FALSE]
  out_region <- out_region[out_region$SNP %in% overlapping_snps, , drop = FALSE]

  locus_safe <- make.names(paste0(exposure_id, "__", outcome_id, "__", analysis_mode, "__", locus$top_snp[[1]]))
  exp_file <- file.path(work_dir, paste0(locus_safe, "__exposure.txt"))
  out_file <- file.path(work_dir, paste0(locus_safe, "__outcome.txt"))
  out_prefix <- file.path(work_dir, paste0(locus_safe, "__pwcoco"))
  coloc_file <- paste0(out_prefix, ".coloc")

  write_pwcoco_sumstats(exp_region, exp_file)
  write_pwcoco_sumstats(out_region, out_file)

  if (file.exists(coloc_file)) {
    file.remove(coloc_file)
  }

  ref_prefix <- get_pwcoco_ref_panel_prefix(build, locus$chr[[1]])
  threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
  if (is.na(threads) || threads < 1) {
    threads <- 1L
  }

  cmd_args <- c(
    "--bfile", ref_prefix,
    "--sum_stats1", exp_file,
    "--sum_stats2", out_file,
    "--out", out_prefix,
    "--chr", as.character(locus$chr[[1]]),
    "--maf", as.character(PWCOCO_MAF),
    "--init_h4", as.character(PWCOCO_H4_THRESHOLD * PWCOCO_PPH4_SCALE),
    "--threads", as.character(threads)
  )

  pwcoco_status <- system2(
    command = pwcoco_bin,
    args = cmd_args,
    stdout = file.path(work_dir, paste0(locus_safe, "__pwcoco.stdout")),
    stderr = file.path(work_dir, paste0(locus_safe, "__pwcoco.stderr"))
  )

  if (!identical(pwcoco_status, 0L)) {
    result$err_message <- paste("PWCoCo failed with exit status", pwcoco_status)
    return(result)
  }

  pwcoco_res <- read_pwcoco_output(coloc_file)
  if (nrow(pwcoco_res) == 0) {
    result$err_message <- paste("PWCoCo output was empty or missing:", coloc_file)
    return(result)
  }

  required_cols <- c("Dataset1", "Dataset2", "SNP1", "SNP2", "nsnps", "H0", "H1", "H2", "H3", "H4")
  missing_cols <- setdiff(required_cols, names(pwcoco_res))
  if (length(missing_cols) > 0) {
    result$err_message <- paste("PWCoCo output missing columns:", paste(missing_cols, collapse = ", "))
    return(result)
  }

  pwcoco_res <- pwcoco_res %>%
    dplyr::mutate(H4_probability = normalise_pp(.data$H4)) %>%
    dplyr::arrange(dplyr::desc(.data$H4_probability))

  best <- pwcoco_res[1, , drop = FALSE]

  result$status <- "completed"
  result$err_message <- ""
  result$Dataset1 <- best$Dataset1[[1]]
  result$Dataset2 <- best$Dataset2[[1]]
  result$SNP1 <- best$SNP1[[1]]
  result$SNP2 <- best$SNP2[[1]]
  result$nsnps <- suppressWarnings(as.integer(best$nsnps[[1]]))
  result$PP.H0.abf <- normalise_pp(best$H0[[1]])
  result$PP.H1.abf <- normalise_pp(best$H1[[1]])
  result$PP.H2.abf <- normalise_pp(best$H2[[1]])
  result$PP.H3.abf <- normalise_pp(best$H3[[1]])
  result$PP.H4.abf <- best$H4_probability[[1]]
  result$H4_flag <- ifelse(result$PP.H4.abf >= PWCOCO_H4_THRESHOLD, "Pass", "Fail")

  result
}

manifest_path <- file.path(output_dir, "trait_manifest_standardised.csv")
if (!file.exists(manifest_path)) {
  manifest_path <- file.path(input_dir, "trait_manifest.csv")
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
plan <- read.csv(file.path(input_dir, "analysis_plan.csv"), stringsAsFactors = FALSE, check.names = FALSE)

plan$run <- normalise_run_flag(plan$run)
plan <- subset(plan, run)

if (nrow(plan) == 0) {
  stop("No rows with run = TRUE found in analysis_plan.csv")
}

if (task_id > nrow(plan)) {
  stop("task_id ", task_id, " is out of range. There are only ", nrow(plan), " active exposure-outcome pairs.")
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
message("Running PWCoCo colocalization")
message("Exposure: ", exposure_id)
message("Outcome : ", outcome_id)
message("Analysis: ", analysis_mode)
message("Window  : +/- ", PWCOCO_WINDOW, " bp around selected instruments")
message("H4 pass threshold: ", PWCOCO_H4_THRESHOLD)
message("========================================")

out_dir <- file.path(output_dir, "Colocalization", "per_pair")
work_dir <- file.path(output_dir, "Colocalization", "pwcoco_work", safe_pair)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, paste0(safe_pair, "__colocalization.tsv"))
done_file <- file.path(out_dir, paste0(safe_pair, ".__DONE"))

if (file.exists(done_file) && file.exists(out_file)) {
  message("Already completed: ", exposure_id, " -> ", outcome_id)
  quit(save = "no", status = 0)
}

pwcoco_bin <- resolve_binary(PWCOCO_BIN)

exposure_file <- file.path(output_dir, "all_exposures.txt")
if (!file.exists(exposure_file)) {
  stop("all_exposures.txt not found: ", exposure_file)
}

all_exposures <- data.table::fread(exposure_file, data.table = FALSE, showProgress = FALSE)
exp_instruments <- all_exposures %>% dplyr::filter(.data$exposure == exposure_id)
selected_instruments <- select_pwcoco_instruments(exp_instruments, exposure_row, analysis_mode)

if (nrow(selected_instruments) == 0) {
  result <- empty_coloc_result()
  result[1, ] <- NA
  result$exposure <- exposure_id
  result$outcome <- outcome_id
  result$analysis <- analysis_mode
  result$status <- "completed"
  result$err_message <- paste("No instruments selected for PWCoCo after", analysis_mode, "filtering")
  result$H4_threshold <- PWCOCO_H4_THRESHOLD
  result$H4_flag <- "Fail"
  safe_write_table(result, out_file)
  file.create(done_file)
  quit(save = "no", status = 0)
}

loci <- make_loci(selected_instruments)
loci <- loci[!is.na(loci$chr) & !is.na(loci$start) & !is.na(loci$end), , drop = FALSE]

if (nrow(loci) == 0) {
  stop("No PWCoCo loci could be defined for exposure: ", exposure_id)
}

exp_path <- resolve_input_file(input_dir, exposure_row$file[[1]])
out_path <- resolve_input_file(input_dir, outcome_row$file[[1]])

if (!file.exists(exp_path)) {
  stop("Exposure GWAS file not found: ", exp_path)
}
if (!file.exists(out_path)) {
  stop("Outcome GWAS file not found: ", out_path)
}

exp_dat <- read_gwas_table(exp_path, get_sep_value(exposure_row$sep[[1]]))
out_dat <- read_gwas_table(out_path, get_sep_value(outcome_row$sep[[1]]))

exp_raw <- standardise_raw_gwas(exp_dat, exposure_row)
out_raw <- standardise_raw_gwas(out_dat, outcome_row)

build_val <- exposure_row$build[[1]]

results <- lapply(seq_len(nrow(loci)), function(i) {
  run_pwcoco_for_locus(
    locus = loci[i, , drop = FALSE],
    exp_raw = exp_raw,
    out_raw = out_raw,
    exposure_id = exposure_id,
    outcome_id = outcome_id,
    analysis_mode = analysis_mode,
    build = build_val,
    pwcoco_bin = pwcoco_bin,
    work_dir = work_dir
  )
})

results <- dplyr::bind_rows(results)
results$selected_instrument_count <- nrow(loci)

safe_write_table(results, out_file)
file.create(done_file)

message("PWCoCo colocalization complete for: ", exposure_id, " -> ", outcome_id)
message("Loci tested: ", nrow(results))
message("H4 pass loci: ", sum(results$H4_flag == "Pass", na.rm = TRUE))
