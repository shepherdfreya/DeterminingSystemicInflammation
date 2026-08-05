# Determining Systemic Inflammation Pipeline

This repository contains a reproducible pipeline for analysing GWAS summary statistics and performing Mendelian randomisation (MR) analyses.

This pipeline was generated for use in:

Shepherd, Freya, Chloe Slaney, Hannah J. Jones, et al. 2026. *Prioritising Determinants of Systemic Inflammation across Molecular, Physiological and Disease Phenotypes.* Preprint, Epidemiology, April 14. <https://doi.org/10.64898/2026.04.10.26350510>.

Please cite appropriately.

------------------------------------------------------------------------

# Overview

The pipeline allows users to:

- Input a folder of GWAS summary statistics
- Run a standardised MR workflow
- Generate harmonised files, results tables, and plots

A key feature of this pipeline is that it allows multiple pairwise analyses of traits, where traits can be analysed flexibly as both exposures and outcomes. This enables systematic MR analyses across many trait pairs.

At the current stage, the pipeline includes:

- input checking and GWAS standardisation
- outcome preparation
- exposure preparation, including significance filtering and LD clumping
- LD proxy identification for exposure SNPs missing from outcome GWAS files
- application of best available proxies to exposure datasets before harmonisation
- MR with Steiger filtering, heterogeneity tests, and MR-Egger intercept tests where applicable
- optional genetic colocalization
- evidence appraisal using tiered criteria

------------------------------------------------------------------------

# Computing Requirements

This pipeline requires:

- R (\>= 4.0)
- PLINK 1.9 or later
- Access to build-specific LD reference panels

Required R packages are checked and loaded by `Utils/Packages.R`.

Packages are **not automatically downloaded**. If packages are missing, the setup script reports which packages need to be installed.

Core packages used include:

- data.table
- dplyr
- TwoSampleMR
- ieugwasr
- gwasvcf

Some `gwasvcf` dependencies are installed through Bioconductor. If `gwasvcf` is installed but fails to load, check that its Bioconductor dependencies are also installed.

The colocalization step additionally requires:

- coloc

Install it with:

``` r
install.packages("coloc")
```

If you install packages into a project-local R library, set `R_LIBS_USER` before running the pipeline so that each `Rscript` call can find them. For example:

``` bash
mkdir -p .r-lib
R_LIBS_USER="$(pwd)/.r-lib" Rscript -e "install.packages(c('data.table', 'dplyr', 'TwoSampleMR', 'gwasvcf'), repos = c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))"
R_LIBS_USER="$(pwd)/.r-lib" Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager', repos = 'https://cloud.r-project.org'); BiocManager::install(c('BiocGenerics', 'Biostrings', 'GenomeInfoDb', 'GenomicRanges', 'IRanges', 'Rsamtools', 'S4Vectors', 'SummarizedExperiment', 'VariantAnnotation'), ask = FALSE, update = FALSE)"
```

------------------------------------------------------------------------

# Input Requirements

The input directory must contain:

- GWAS summary statistic files
- `trait_manifest.csv`
- `analysis_plan.csv`

------------------------------------------------------------------------

# trait_manifest.csv

This file tells the pipeline which GWAS files are available and how to interpret their columns.

It must contain the following columns:

| Column | Description |
|----|----|
| `trait_id` | Short unique identifier for the trait. Must match traits in `analysis_plan.csv`. |
| `file` | GWAS filename located in the input directory. |
| `sep` | File delimiter, usually `\t` for tab-delimited files. |
| `snp_col` | Column containing SNP identifiers. |
| `position_col` | Column containing base-pair position. |
| `effect_col` | Column containing the effect estimate, e.g. beta or OR. |
| `se_col` | Column containing the standard error. |
| `ea_col` | Column containing the effect allele. |
| `oa_col` | Column containing the other allele. |
| `pval_col` | Column containing the p-value. |
| `build` | Genome build for the file, e.g. `GRCh37` or `GRCh38`. |
| `eaf_col` | Column containing effect allele frequency. |
| `chromosome_col` | Column containing chromosome. |
| `protein_start_col` | Protein start position, if relevant. Otherwise use `NA`. |
| `protein_end_col` | Protein end position, if relevant. Otherwise use `NA`. |
| `protein_chr_col` | Protein chromosome, if relevant. Otherwise use `NA`. |

Important: empty cells in `trait_manifest.csv` should be filled with `NA`.

------------------------------------------------------------------------

# Important note on effect_col

The `effect_col` can contain either:

- beta values, or
- odds ratio values

If `effect_col = OR`, the input-checking step appends a `beta` column to that GWAS file using:

``` r
beta = log(OR)
```

The original effect column is retained in the GWAS file.

The pipeline then writes an output-side standardised manifest:

``` text
trait_manifest_standardised.csv
```

Downstream scripts automatically use this standardised manifest when it is present, so the original `trait_manifest.csv` does not need to be edited manually between steps.

------------------------------------------------------------------------

# analysis_plan.csv

This file specifies which analyses to run.

It must contain:

| exposure   | outcome | run   | Analysis  |
|------------|---------|-------|-----------|
| CRP        | IL6     | FALSE | cis_trans |
| IL6        | CRP     | TRUE  | cis_only  |
| IL6        | CRP     | TRUE  | cis_trans |
| Depression | CRP     | TRUE  | cis_trans |
| Depression | IL6     | TRUE  | cis_trans |

Only rows with `run = TRUE` are processed.

The `Analysis` column is optional for backwards compatibility. If provided, values must be:

``` text
cis_only
cis_trans
```

Use `cis_only` for protein exposures where instruments should be restricted to the protein coordinate window. Use `cis_trans` to include all retained instruments. The cis window is defined in `Utils/Config.R` as `CIS_WINDOW <- 1e6`, meaning plus or minus 1 Mb around the protein start/end coordinates supplied in `trait_manifest.csv`.

The `Analysis` value is carried through MR for naming and summaries, but proxy candidate search and proxy application are de-duplicated across unique exposure-outcome pairs. The cis-only instrument restriction is applied immediately before MR in `Scripts/05a_run_mr_array.R`. This allows proxy recovery to use the full prepared exposure set before the final cis/trans MR split is made.

------------------------------------------------------------------------

# Configuration and Settings

User-defined settings are stored in:

``` text
Utils/Config.R
```

These include:

- `PLINK_BIN`: path to the PLINK executable
- `AUX_DIR`: folder containing LD reference panels
- `REF_PANEL_PREFIXES`: build-specific LD reference panel prefixes
- exposure MAF threshold
- exposure p-value threshold
- exposure F-statistic threshold
- exposure LD clumping r2 threshold
- exposure LD clumping window size
- cis analysis window
- proxy search r2 threshold
- proxy search window size
- maximum number of proxy SNPs to return
- colocalization H4 pass threshold
- minimum overlapping SNPs required for colocalization
- evidence appraisal alpha threshold

Example proxy settings:

``` r
PROXY_TAG_R2   <- 0.8
PROXY_TAG_KB   <- 5000
PROXY_TAG_NSNP <- 5000
```

Example colocalization settings:

``` r
COLOC_H4_THRESHOLD <- 0.8
COLOC_MIN_OVERLAP  <- 10
```

Example cis/trans analysis setting:

``` r
CIS_WINDOW <- 1e6
```

Example evidence appraisal setting:

``` r
EVIDENCE_ALPHA <- 0.05
EVIDENCE_BONFERRONI_N_TESTS <- NA
```

Example reference panel settings:

``` r
PLINK_BIN <- "/software/local/apps/plink1.9-2025.01/plink"

REF_PANEL_PREFIXES <- list(
  GRCh37 = "Data/Aux/EUR_HG19GRCh37/EUR_HG19GRCh37",
  GRCh38 = "Data/Aux/EUR_HG38GRCh38/EUR_HG38GRCh38"
)
```

These paths should be checked and edited by the user if local paths differ.

------------------------------------------------------------------------

# Running the Pipeline

The pipeline is designed as a mixture of serial and array-based steps. Array-based steps are used for computationally heavier tasks and are suitable for HPC/SLURM execution.

Current structure:

``` text
01_check_inputs.R                 serial
02a_prepare_outcomes_array.R      array
02b_combine_outcomes.R            serial
03a_prepare_exposures_array.R     array
03b_combine_exposures.R           serial
04a_find_proxy_candidates_array.R array
04b_combine_proxy_candidates.R    serial
04c_apply_proxies_array.R         array
05a_run_mr_array.R                array
05b_combine_mr_outputs.R          serial
06a_run_colocalization_array.R    array
06b_combine_colocalization.R      serial
07_evidence_appraisal.R           serial
```

For a local end-to-end demonstration using the bundled trimmed example data, run:

``` bash
bash Shell/run_pipeline_local.sh example_data Results/ExampleRun
```

If using a project-local R library, prefix the command with `R_LIBS_USER`:

``` bash
R_LIBS_USER="$(pwd)/.r-lib" bash Shell/run_pipeline_local.sh example_data Results/ExampleRun
```

The local wrapper runs outcome preparation, exposure preparation, proxy finding, proxy application, MR, MR output combination, and evidence appraisal. Colocalization is optional because it can be expensive:

``` bash
RUN_COLOC=TRUE bash Shell/run_pipeline_local.sh example_data Results/ExampleRun
```

The second argument names the output folder. For example, `Results/ExampleRun`, `Results/TestRun`, and `Results/SmokeTest` are separate runs created by passing different output directories.

For SLURM execution from the project root, run:

``` bash
bash Shell/run_pipeline.sh <input_dir> <output_dir>
```

The SLURM wrapper submits the full dependency chain through outcome/exposure preparation, proxy finding, proxy application, MR, MR summary combination, and evidence appraisal. Colocalization is optional because it can be expensive at scale:

``` bash
RUN_COLOC=TRUE bash Shell/run_pipeline.sh <input_dir> <output_dir>
```

If using a project-local R library on the cluster, export or prefix `R_LIBS_USER` so each submitted job can find the same packages:

``` bash
R_LIBS_USER="/path/to/DeterminingSystemicInflammation/.r-lib" bash Shell/run_pipeline.sh <input_dir> <output_dir>
```

The wrapper inherits `R_LIBS_USER` inside each SLURM job. It also allows cluster-specific resource settings without editing the script:

``` bash
R_MODULE=R/4.3.1 \
EXPOSURE_MEM=16G EXPOSURE_TIME=08:00:00 \
PROXY_MEM=32G PROXY_TIME=12:00:00 \
MR_MEM=16G MR_TIME=08:00:00 \
bash Shell/run_pipeline.sh <input_dir> <output_dir>
```

Useful flags:

| Setting | Default | Purpose |
|----|----|----|
| `RUN_COLOC` | `FALSE` | Run colocalization and combine colocalization outputs |
| `RUN_APPRAISAL` | `TRUE` | Run `Scripts/07_evidence_appraisal.R` after MR, or after colocalization if enabled |
| `R_MODULE` | `R` | Module name loaded inside each SLURM job |
| `OUTCOME_MEM`, `OUTCOME_TIME` | `8G`, `02:00:00` | Resources for outcome preparation |
| `EXPOSURE_MEM`, `EXPOSURE_TIME` | `8G`, `04:00:00` | Resources for exposure preparation and clumping |
| `PROXY_MEM`, `PROXY_TIME` | `16G`, `08:00:00` | Resources for proxy search and proxy application |
| `MR_MEM`, `MR_TIME` | `8G`, `04:00:00` | Resources for MR array jobs |
| `COLOC_MEM`, `COLOC_TIME` | `16G`, `08:00:00` | Resources for colocalization array jobs |
| `COMBINE_MEM`, `COMBINE_TIME` | `4G`, `01:00:00` | Resources for serial combine/appraisal jobs |

For real runs with thousands of proteins, start with a pilot manifest containing a small active subset, check `Results/<run>/pipeline_log.txt`, the files under `Results/<run>/logs/`, `mr_summary.tsv`, `mr_results.tsv`, and `evidence_appraisal.tsv`, then scale up by increasing the number of active rows in `analysis_plan.csv`.

The evidence appraisal criteria are also summarised in:

``` text
docs/evidence_appraisal_criteria.md
```

The example files are small extracts of redistributable GWAS summary statistics. They are intended for checking installation and demonstrating the workflow, not for scientific interpretation.

The example data directory may also contain `example_subset_summary.tsv`, a static summary of the number of rows requested and written when the bundled example files were created. It is not generated by the pipeline.

------------------------------------------------------------------------

# GWAS Standardisation

The input-checking step checks that required columns are present in each GWAS file.

If `effect_col` contains OR values, the pipeline can append a beta column using:

``` r
beta = log(OR)
```

The original OR column is retained.

The input-checking step writes `trait_manifest_standardised.csv` to the output directory and changes converted traits to `effect_col = beta` in that file. Downstream scripts automatically use this standardised manifest when it is present.

------------------------------------------------------------------------

# Outcome Preparation

The outcome preparation step reads all traits listed as outcomes in `analysis_plan.csv` where `run = TRUE`.

Column mappings are taken from `trait_manifest.csv`, allowing different GWAS files to use different header names.

This step writes:

``` text
Prepared_outcomes/
all_outcomes.txt
```

------------------------------------------------------------------------

# Exposure Preparation

The exposure preparation step reads all traits listed as exposures in `analysis_plan.csv` where `run = TRUE`.

This step includes:

- reading exposure GWAS data using manifest-defined columns
- minor allele frequency filtering
- genome-wide significance filtering
- LD clumping using PLINK
- F-statistic filtering after exposure formatting, using `(beta / se)^2`
- build-specific reference panel selection based on the `build` column in `trait_manifest.csv`

Current exposure filtering uses:

``` text
MAF threshold: > 0.01
p-value threshold: p < 5e-8
F-statistic threshold: F > 10
```

After this, variants are LD clumped using:

- PLINK
- a build-matched LD reference panel
- thresholds defined in `Utils/Config.R`

This step writes:

``` text
Prepared_exposures/
all_exposures.txt
per_task_summaries/
processing_summary_exposures.txt
```

------------------------------------------------------------------------

# Proxy Candidate Search

Some exposure SNPs may be absent from a given outcome GWAS. The proxy search step identifies LD proxies for these missing exposure SNPs.

This is performed by:

``` text
Scripts/04a_find_proxy_candidates_array.R
```

For each exposure-outcome pair, this script:

1.  Reads `all_exposures.txt` and `all_outcomes.txt`
2.  Identifies exposure SNPs missing from the outcome GWAS
3.  Selects the correct LD reference panel using the exposure genome build
4.  Finds LD proxies using `gwasvcf::get_ld_proxies`
5.  Keeps only proxies present in the outcome GWAS
6.  Writes pair-specific proxy candidate files

Proxy candidate search is performed once per unique exposure-outcome pair, even when `analysis_plan.csv` contains both `cis_only` and `cis_trans` rows for the same pair. It does not apply the cis-only filter; cis-only filtering is applied later, immediately before MR.

Outputs are written to:

``` text
Proxy_candidates/per_pair/
proxy_candidate_summaries/
```

Pair-specific files include:

``` text
<exposure>__<outcome>__proxy_candidates.tsv
<exposure>__<outcome>__missing_snps.tsv
```

The per-task summaries record:

- original exposure SNP count
- outcome SNP count
- number of missing SNPs
- number of proxy rows found
- number of proxy rows present in the outcome GWAS
- number of missing SNPs with at least one candidate proxy

------------------------------------------------------------------------

# Combining Proxy Candidates

Proxy candidate files are combined using:

``` text
Scripts/04b_combine_proxy_candidates.R
```

This serial step combines all pair-specific proxy outputs into:

``` text
Proxy_candidates/proxy_candidates.tsv
Proxy_candidates/missing_snps.tsv
Proxy_candidates/proxy_discovery_summary.tsv
```

These combined files are used by the next proxy application step.

------------------------------------------------------------------------

# Applying Proxies

Best available proxies are applied using:

``` text
Scripts/04c_apply_proxies_array.R
```

For each exposure-outcome pair, this script:

1.  Reads `all_exposures.txt` and `all_outcomes.txt`
2.  Loads combined proxy candidates from `Proxy_candidates/proxy_candidates.tsv`
3.  Reads the full raw exposure GWAS from the input directory
4.  Keeps proxy SNPs that are present in:
    - the outcome GWAS
    - the full exposure GWAS
5.  Selects the best proxy for each missing SNP based on highest `R2`
6.  Removes exposure SNPs missing from the outcome GWAS
7.  Formats selected proxy SNPs as exposure data
8.  Applies the same F-statistic filter to proxy rows
9.  Merges proxy rows into the exposure dataset
10. Writes pair-specific proxied exposure files

Proxy application is performed once per unique exposure-outcome pair, even when `analysis_plan.csv` contains both `cis_only` and `cis_trans` rows for the same pair. It reuses the shared exposure-outcome proxy candidates for all analysis modes and does not apply the cis-only filter; the final cis-only restriction is applied in the MR step.

Outputs are written to:

``` text
Proxied_exposures/per_pair/
proxy_application_summaries/
```

Pair-specific files include:

``` text
<exposure>__<outcome>__exposure_proxied.tsv
<exposure>__<outcome>__proxy_map.tsv
<exposure>__<outcome>__no_proxy.tsv
```

The proxy map records:

``` text
exposure
outcome
missing_snp
proxy_snp
R
R2
```

The no-proxy file records missing SNPs that could not be replaced. Proxy rows that fail the F-statistic threshold are also recorded there and are not added to the proxied exposure dataset.

------------------------------------------------------------------------

# Pipeline Structure

Scripts completed so far:

| Script | Type | Purpose |
|----|----|----|
| `01_check_inputs.R` | Serial | Input checking and GWAS standardisation |
| `02a_prepare_outcomes_array.R` | Array | Prepare one outcome per task |
| `02b_combine_outcomes.R` | Serial | Combine prepared outcomes into `all_outcomes.txt` |
| `03a_prepare_exposures_array.R` | Array | Prepare one exposure per task, including MAF filtering, p-value filtering, and LD clumping |
| `03b_combine_exposures.R` | Serial | Combine prepared exposures into `all_exposures.txt` |
| `04a_find_proxy_candidates_array.R` | Array | Identify LD proxy candidates for SNPs missing from each outcome GWAS |
| `04b_combine_proxy_candidates.R` | Serial | Combine proxy candidate files across all exposure-outcome pairs |
| `04c_apply_proxies_array.R` | Array | Apply best available proxies to exposure datasets for each exposure-outcome pair |
| `05a_run_mr_array.R` | Array | Harmonise one exposure-outcome pair and run MR |
| `05b_combine_mr_outputs.R` | Serial | Combine pair-specific MR summaries and MR results for the run |
| `06a_run_colocalization_array.R` | Array | Run genetic colocalization for one exposure-outcome pair |
| `06b_combine_colocalization.R` | Serial | Combine pair-specific colocalization outputs for the run |
| `07_evidence_appraisal.R` | Serial | Assign evidence tiers using MR significance, pQTL class, colocalization, heterogeneity, and pleiotropy diagnostics |

------------------------------------------------------------------------

# MR Outputs

The MR step writes one summary file per exposure-outcome pair to:

``` text
mr_summaries/
```

These files are created by `Scripts/05a_run_mr_array.R`. The local wrapper then runs `Scripts/05b_combine_mr_outputs.R`, which combines them into:

``` text
mr_summary.tsv
```

The same combine step also combines pair-specific MR result files from:

``` text
MR_results/per_pair/
```

into:

``` text
mr_results.tsv
```

For rows where `Analysis = cis_only`, `Scripts/05a_run_mr_array.R` annotates instruments as cis or trans using the exposure protein coordinates and applies the cis-only filter immediately before harmonisation/MR. Rows where `Analysis = cis_trans` retain all pair-specific instruments.

------------------------------------------------------------------------

# Genetic Colocalization

Genetic colocalization can be run after MR as a sensitivity analysis using:

``` text
Scripts/06a_run_colocalization_array.R
```

For each exposure-outcome pair, this script:

1.  Reads the full raw exposure and outcome GWAS files using the manifest.
2.  Defines colocalization loci.
3.  Runs `coloc::coloc.abf` for each locus with overlapping SNPs.
4.  Reports posterior probabilities for H0-H4.
5.  Flags loci as `Pass` when `PP.H4.abf > COLOC_H4_THRESHOLD`.

By default:

``` text
H4 > 0.8 = Pass
```

For exposures with protein coordinates in `trait_manifest.csv`, the script uses the protein region plus `cis_window` from `Utils/Config.R`. For exposures without protein coordinates, it falls back to instrument-centred loci using `all_exposures.txt`.

Pair-specific outputs are written to:

``` text
Colocalization/per_pair/
```

with files named:

``` text
<exposure>__<outcome>__colocalization.tsv
```

The combined colocalization summary is created with:

``` text
Scripts/06b_combine_colocalization.R
```

and written to:

``` text
colocalization_summary.tsv
```

To run locally for all active exposure-outcome pairs after MR:

``` bash
N_PAIRS=$(Rscript -e "plan <- read.csv(file.path('example_data','analysis_plan.csv'), stringsAsFactors=FALSE); plan\$run <- tolower(as.character(plan\$run)) %in% c('true','t','1'); cat(sum(plan\$run))")
for TASK_ID in $(seq 1 "$N_PAIRS"); do
  Rscript Scripts/06a_run_colocalization_array.R example_data Results/ExampleRun "$TASK_ID"
done
Rscript Scripts/06b_combine_colocalization.R Results/ExampleRun
```

If using a project-local R library, prefix these commands with `R_LIBS_USER="$(pwd)/.r-lib"`.

------------------------------------------------------------------------

# Evidence Appraisal

Evidence appraisal is performed after MR outputs have been combined using:

``` text
Scripts/07_evidence_appraisal.R
```

This script reads:

``` text
mr_results.tsv
colocalization_summary.tsv
```

The colocalization file is optional for running the script. If it is not present, the colocalization criterion is reported as `Not available`, and rows cannot be assigned to Tier A, B, or C.

The script selects one primary MR result per exposure-outcome-analysis set:

- `Wald ratio` for single-SNP analyses
- `Inverse variance weighted` for analyses with two or more SNPs

It then calculates:

- Bonferroni significance using `EVIDENCE_ALPHA / number_of_primary_tests`
- FDR significance using Benjamini-Hochberg correction
- colocalization pass status using `PP.H4.abf > COLOC_H4_THRESHOLD`
- heterogeneity pass status using Cochran's Q `p > EVIDENCE_ALPHA`
- pleiotropy pass status using MR-Egger intercept `p > EVIDENCE_ALPHA`
- causal direction status using Steiger filtering from `mr_summary.tsv`

Heterogeneity and pleiotropy diagnostics can be `Not applicable` when there are too few SNPs for the test. `Not applicable` is allowed for tier assignment and does not count as a failed diagnostic.

By default, the Bonferroni denominator is the number of primary MR tests in the run. To use a fixed denominator, set `EVIDENCE_BONFERRONI_N_TESTS` in `Utils/Config.R`.

Evidence tiers are assigned as:

| Tier | pQTLs | FDR | Bonferroni | Colocalization | Heterogeneity | Pleiotropy | Causal direction |
|----|----|----|----|----|----|----|----|
| Tier A | Cis | Pass | Pass | H4 pass | Pass or not applicable | Pass or not applicable | Pass |
| Tier B | Cis | Pass | No pass | H4 pass | Pass or not applicable | Pass or not applicable | Pass |
| Tier C | Cis/Trans | Pass | Not required | H4 pass | Pass or not applicable | Pass or not applicable | Pass |
| No Tier | Any | Criteria not met | Criteria not met | Criteria not met | Criteria not met | Criteria not met | Criteria not met |

The output is:

``` text
evidence_appraisal.tsv
```

Diagnostic columns are retained in the output so each tier assignment can be audited.

------------------------------------------------------------------------

# Notes

- Column names can vary between GWAS files. This is handled via `trait_manifest.csv`.
- GWAS files do not need to be manually renamed to a standard set of headers.
- Empty cells in `trait_manifest.csv` should be filled with `NA`.
- Always run the input-checking step before downstream analyses.
- If `effect_col = OR`, run the input-checking step first. It will write `trait_manifest_standardised.csv`, and downstream scripts will use that file automatically.
- Exposure preparation requires a valid PLINK path and appropriate build-specific LD reference panel files.
- Exposure instruments and proxy instruments are filtered using `EXPOSURE_F_THRESHOLD` from `Utils/Config.R`.
- Proxy search also requires build-specific LD reference panels.
- The `trait_id` values in `trait_manifest.csv` must exactly match the exposure and outcome names used in `analysis_plan.csv`.
- Proxy SNPs are selected using LD information and then reformatted using the full exposure GWAS.
- Colocalization is a sensitivity analysis and requires dense regional GWAS data, not only clumped instruments.
- Allele harmonisation is still performed downstream after proxy application.

------------------------------------------------------------------------

# Next Steps

Potential extensions include:

- additional sensitivity analyses
- plotting and report generation
