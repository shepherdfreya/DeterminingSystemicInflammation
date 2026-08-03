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

If OR values are supplied, they are converted within the pipeline to beta estimates using:

``` r
beta = log(OR)
```

The original effect column is retained in the GWAS file.

To keep downstream analyses simple, the manifest should then be manually updated so that:

``` text
effect_col = beta
```

------------------------------------------------------------------------

# analysis_plan.csv

This file specifies which analyses to run.

It must contain:

| exposure   | outcome | run   |
|------------|---------|-------|
| CRP        | IL6     | FALSE |
| IL6        | CRP     | TRUE  |
| Obesity    | CRP     | TRUE  |
| Obesity    | IL6     | TRUE  |
| Depression | CRP     | TRUE  |
| Depression | IL6     | TRUE  |

Only rows with `run = TRUE` are processed.

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
- exposure LD clumping r2 threshold
- exposure LD clumping window size
- proxy search r2 threshold
- proxy search window size
- maximum number of proxy SNPs to return

Example proxy settings:

``` r
PROXY_TAG_R2   <- 0.8
PROXY_TAG_KB   <- 5000
PROXY_TAG_NSNP <- 5000
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
```

From the project root, run:

``` bash
bash Shell/run_pipeline.sh <input_dir> <output_dir>
```

The wrapper script runs/submits each step in order.

For a local demonstration using the bundled trimmed example data, run:

``` bash
bash Shell/run_pipeline_local.sh example_data Results/ExampleRun
```

The example files are small extracts of redistributable GWAS summary statistics. They are intended for checking installation and demonstrating the workflow, not for scientific interpretation.

------------------------------------------------------------------------

# GWAS Standardisation

The input-checking step checks that required columns are present in each GWAS file.

If `effect_col` contains OR values, the pipeline can append a beta column using:

``` r
beta = log(OR)
```

The original OR column is retained.

Before continuing downstream, update `trait_manifest.csv` so that `effect_col` points to the beta column.

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
- build-specific reference panel selection based on the `build` column in `trait_manifest.csv`

Current exposure filtering uses:

``` text
MAF threshold: > 0.01
p-value threshold: p < 5e-8
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
8.  Merges proxy rows into the exposure dataset
9.  Writes pair-specific proxied exposure files

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

The no-proxy file records missing SNPs that could not be replaced.

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

------------------------------------------------------------------------

# Notes

- Column names can vary between GWAS files. This is handled via `trait_manifest.csv`.
- GWAS files do not need to be manually renamed to a standard set of headers.
- Empty cells in `trait_manifest.csv` should be filled with `NA`.
- Always run the input-checking step before downstream analyses.
- Ensure `effect_col` is set to beta before proceeding beyond standardisation.
- Exposure preparation requires a valid PLINK path and appropriate build-specific LD reference panel files.
- Proxy search also requires build-specific LD reference panels.
- The `trait_id` values in `trait_manifest.csv` must exactly match the exposure and outcome names used in `analysis_plan.csv`.
- Proxy SNPs are selected using LD information and then reformatted using the full exposure GWAS.
- Allele harmonisation is still performed downstream after proxy application.

------------------------------------------------------------------------

# Next Steps

Downstream steps to be added include:

- harmonisation of proxied exposures and outcomes
- Mendelian randomisation analyses
- sensitivity analyses
- results summarisation
- plotting
