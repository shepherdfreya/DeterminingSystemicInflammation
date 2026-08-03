############################################################
# File: Scripts/07_evidence_appraisal.R
# Purpose: Appraise MR evidence tiers using MR, diagnostic, and coloc outputs
# Usage:
#   Rscript Scripts/07_evidence_appraisal.R <output_dir>
############################################################

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript Scripts/07_evidence_appraisal.R <output_dir>")
}

output_dir <- args[1]

source("Utils/Packages.R")
source("Utils/Config.R")

mr_results_file <- file.path(output_dir, "mr_results.tsv")
mr_summary_file <- file.path(output_dir, "mr_summary.tsv")
coloc_file <- file.path(output_dir, "colocalization_summary.tsv")

if (!file.exists(mr_results_file)) {
  stop("Combined MR results file not found: ", mr_results_file)
}

normalise_analysis <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x[x == "" | is.na(x)] <- "cis_trans"
  x
}

yes_no <- function(x) {
  ifelse(is.na(x), "Not available", ifelse(x, "Yes", "No"))
}

diagnostic_label <- function(pval, applicable) {
  dplyr::case_when(
    !applicable ~ "Not applicable",
    is.na(pval) ~ "Not available",
    pval > EVIDENCE_ALPHA ~ "Yes",
    TRUE ~ "No"
  )
}

pass_or_not_applicable <- function(label) {
  label %in% c("Yes", "Not applicable")
}

safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  min(x, na.rm = TRUE)
}

safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (length(x) == 0 || all(is.na(x))) {
    return(NA_real_)
  }
  max(x, na.rm = TRUE)
}

select_primary_mr <- function(dat) {
  dat <- dat %>%
    dplyr::mutate(
      method_rank = dplyr::case_when(
        .data$method == "Wald ratio" ~ 1L,
        .data$method == "Inverse variance weighted" ~ 2L,
        TRUE ~ 99L
      )
    ) %>%
    dplyr::arrange(.data$method_rank)

  dat[1, , drop = FALSE]
}

mr_results <- utils::read.table(
  mr_results_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_mr_cols <- c("exposure", "outcome", "method", "pval")
missing_mr_cols <- setdiff(required_mr_cols, names(mr_results))

if (length(missing_mr_cols) > 0) {
  stop("mr_results.tsv is missing required columns: ", paste(missing_mr_cols, collapse = ", "))
}

if (!"analysis_set" %in% names(mr_results)) {
  mr_results$analysis_set <- "cis_trans"
}

mr_results$analysis_set <- normalise_analysis(mr_results$analysis_set)

primary <- mr_results %>%
  dplyr::group_by(.data$exposure, .data$outcome, .data$analysis_set) %>%
  dplyr::group_modify(~ select_primary_mr(.x)) %>%
  dplyr::ungroup()

primary$primary_pval <- suppressWarnings(as.numeric(primary$pval))
primary$fdr_pval <- stats::p.adjust(primary$primary_pval, method = "BH")

bonferroni_n_tests <- suppressWarnings(as.numeric(EVIDENCE_BONFERRONI_N_TESTS))
if (is.na(bonferroni_n_tests) || bonferroni_n_tests <= 0) {
  bonferroni_n_tests <- nrow(primary)
}

bonferroni_threshold <- EVIDENCE_ALPHA / bonferroni_n_tests

diagnostics <- mr_results %>%
  dplyr::group_by(.data$exposure, .data$outcome, .data$analysis_set) %>%
  dplyr::summarise(
    min_Q_pval = safe_min(.data$Q_pval),
    min_egger_intercept_pval = safe_min(.data$egger_intercept_pval),
    .groups = "drop"
  )

appraisal <- primary %>%
  dplyr::select(
    "exposure",
    "outcome",
    analysis = "analysis_set",
    primary_method = "method",
    "nsnp",
    "b",
    "se",
    "primary_pval",
    "fdr_pval"
  ) %>%
  dplyr::left_join(
    diagnostics,
    by = c("exposure", "outcome", "analysis" = "analysis_set")
  )

if (file.exists(mr_summary_file)) {
  mr_summary <- utils::read.table(
    mr_summary_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"analysis" %in% names(mr_summary)) {
    mr_summary$analysis <- "cis_trans"
  }

  mr_summary$analysis <- normalise_analysis(mr_summary$analysis)

  if ("steiger_filter_applied" %in% names(mr_summary)) {
    direction_summary <- mr_summary %>%
      dplyr::group_by(.data$exposure, .data$outcome, .data$analysis) %>%
      dplyr::summarise(
        correct_causal_direction = any(.data$steiger_filter_applied %in% TRUE, na.rm = TRUE),
        .groups = "drop"
      )

    appraisal <- appraisal %>%
      dplyr::left_join(
        direction_summary,
        by = c("exposure", "outcome", "analysis")
      )
  }
}

if (!"correct_causal_direction" %in% names(appraisal)) {
  appraisal$correct_causal_direction <- NA
}

if (file.exists(coloc_file)) {
  coloc <- utils::read.table(
    coloc_file,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"analysis" %in% names(coloc)) {
    coloc$analysis <- "cis_trans"
  }

  coloc$analysis <- normalise_analysis(coloc$analysis)

  coloc_summary <- coloc %>%
    dplyr::group_by(.data$exposure, .data$outcome, .data$analysis) %>%
    dplyr::summarise(
      max_PP.H4.abf = safe_max(.data$PP.H4.abf),
      coloc_H4_pass = any(.data$H4_flag == "Pass", na.rm = TRUE),
      .groups = "drop"
    )

  appraisal <- appraisal %>%
    dplyr::left_join(
      coloc_summary,
      by = c("exposure", "outcome", "analysis")
    )
} else {
  appraisal$max_PP.H4.abf <- NA_real_
  appraisal$coloc_H4_pass <- NA
}

appraisal <- appraisal %>%
  dplyr::mutate(
    pQTLs = dplyr::if_else(.data$analysis == "cis_only", "Cis", "Cis/Trans"),
    bonferroni_threshold = bonferroni_threshold,
    bonferroni_n_tests = bonferroni_n_tests,
    significant_bonferroni = !is.na(.data$primary_pval) & .data$primary_pval <= .data$bonferroni_threshold,
    significant_fdr = !is.na(.data$fdr_pval) & .data$fdr_pval <= EVIDENCE_ALPHA,
    significant_bonferroni_label = yes_no(.data$significant_bonferroni),
    significant_fdr_label = yes_no(.data$significant_fdr),
    coloc_H4_label = yes_no(.data$coloc_H4_pass),
    heterogeneity_label = diagnostic_label(.data$min_Q_pval, .data$nsnp >= 2),
    pleiotropy_label = diagnostic_label(.data$min_egger_intercept_pval, .data$nsnp >= 3),
    correct_causal_direction_label = yes_no(.data$correct_causal_direction),
    diagnostics_ok = pass_or_not_applicable(.data$heterogeneity_label) &
      pass_or_not_applicable(.data$pleiotropy_label),
    coloc_ok = .data$coloc_H4_label == "Yes",
    causal_direction_ok = .data$correct_causal_direction_label == "Yes",
    evidence_tier = dplyr::case_when(
      .data$significant_fdr &
      .data$pQTLs == "Cis" &
        .data$significant_bonferroni &
        .data$diagnostics_ok &
        .data$coloc_ok &
        .data$causal_direction_ok ~ "Tier A",
      .data$significant_fdr &
        .data$pQTLs == "Cis" &
        !.data$significant_bonferroni &
        .data$diagnostics_ok &
        .data$coloc_ok &
        .data$causal_direction_ok ~ "Tier B",
      .data$significant_fdr &
        .data$pQTLs == "Cis/Trans" &
        .data$diagnostics_ok &
        .data$coloc_ok &
        .data$causal_direction_ok ~ "Tier C",
      TRUE ~ "No Tier"
    )
  ) %>%
  dplyr::select(
    "exposure",
    "outcome",
    "analysis",
    "pQTLs",
    "evidence_tier",
    "primary_method",
    "nsnp",
    "b",
    "se",
    "primary_pval",
    "fdr_pval",
    "bonferroni_n_tests",
    "bonferroni_threshold",
    "significant_bonferroni_label",
    "significant_fdr_label",
    "max_PP.H4.abf",
    "coloc_H4_label",
    "min_Q_pval",
    "heterogeneity_label",
    "min_egger_intercept_pval",
    "pleiotropy_label",
    "correct_causal_direction_label"
  )

out_file <- file.path(output_dir, "evidence_appraisal.tsv")

write.table(
  appraisal,
  file = out_file,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

message("Evidence appraisal written to: ", out_file)
message("Rows appraised: ", nrow(appraisal))
message("Tier counts:")
print(table(appraisal$evidence_tier, useNA = "ifany"))
