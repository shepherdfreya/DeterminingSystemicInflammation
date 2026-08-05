#!/bin/bash

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: bash Shell/run_pipeline.sh <input_dir> <output_dir>"
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
LOG_FILE="$OUTPUT_DIR/pipeline_log.txt"

R_MODULE=${R_MODULE:-R}
RUN_COLOC=${RUN_COLOC:-FALSE}
RUN_APPRAISAL=${RUN_APPRAISAL:-TRUE}

OUTCOME_MEM=${OUTCOME_MEM:-8G}
OUTCOME_TIME=${OUTCOME_TIME:-02:00:00}
EXPOSURE_MEM=${EXPOSURE_MEM:-8G}
EXPOSURE_TIME=${EXPOSURE_TIME:-04:00:00}
PROXY_MEM=${PROXY_MEM:-16G}
PROXY_TIME=${PROXY_TIME:-08:00:00}
MR_MEM=${MR_MEM:-8G}
MR_TIME=${MR_TIME:-04:00:00}
COLOC_MEM=${COLOC_MEM:-16G}
COLOC_TIME=${COLOC_TIME:-08:00:00}
COMBINE_MEM=${COMBINE_MEM:-4G}
COMBINE_TIME=${COMBINE_TIME:-01:00:00}

normalise_bool() {
  case "$(printf "%s" "$1" | tr "[:lower:]" "[:upper:]")" in
    TRUE|T|1|YES|Y) printf "TRUE" ;;
    *) printf "FALSE" ;;
  esac
}

RUN_COLOC=$(normalise_bool "$RUN_COLOC")
RUN_APPRAISAL=$(normalise_bool "$RUN_APPRAISAL")

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================"
echo "Submitting SLURM pipeline"
echo "Date             : $(date)"
echo "Input directory  : $INPUT_DIR"
echo "Output directory : $OUTPUT_DIR"
echo "Log file         : $LOG_FILE"
echo "R module         : $R_MODULE"
echo "R_LIBS_USER      : ${R_LIBS_USER:-<not set>}"
echo "Run coloc        : $RUN_COLOC"
echo "Run appraisal    : $RUN_APPRAISAL"
echo "========================================"

count_outcomes() {
  Rscript -e 'args <- commandArgs(TRUE); plan <- read.csv(file.path(args[1], "analysis_plan.csv"), stringsAsFactors=FALSE); plan$run <- tolower(as.character(plan$run)) %in% c("true","t","1"); plan <- subset(plan, run); cat(length(unique(plan$outcome)))' "$INPUT_DIR"
}

count_exposures() {
  Rscript -e 'args <- commandArgs(TRUE); plan <- read.csv(file.path(args[1], "analysis_plan.csv"), stringsAsFactors=FALSE); plan$run <- tolower(as.character(plan$run)) %in% c("true","t","1"); plan <- subset(plan, run); cat(length(unique(plan$exposure)))' "$INPUT_DIR"
}

count_pairs() {
  Rscript -e 'args <- commandArgs(TRUE); plan <- read.csv(file.path(args[1], "analysis_plan.csv"), stringsAsFactors=FALSE); plan$run <- tolower(as.character(plan$run)) %in% c("true","t","1"); plan <- subset(plan, run); cat(nrow(plan))' "$INPUT_DIR"
}

count_proxy_pairs() {
  Rscript -e 'args <- commandArgs(TRUE); plan <- read.csv(file.path(args[1], "analysis_plan.csv"), stringsAsFactors=FALSE); plan$run <- tolower(as.character(plan$run)) %in% c("true","t","1"); plan <- subset(plan, run); cat(nrow(unique(plan[, c("exposure", "outcome"), drop=FALSE])))' "$INPUT_DIR"
}

############################################################
# STEP 1: CHECK INPUTS (SERIAL)
############################################################

echo ""
echo "Step 1: Checking inputs"
Rscript Scripts/01_check_inputs.R "$INPUT_DIR" "$OUTPUT_DIR"

N_OUTCOMES=$(count_outcomes)
N_EXPOSURES=$(count_exposures)
N_PAIRS=$(count_pairs)
N_PROXY_PAIRS=$(count_proxy_pairs)

echo "Number of outcome tasks       : $N_OUTCOMES"
echo "Number of exposure tasks      : $N_EXPOSURES"
echo "Number of unique proxy-search pairs: $N_PROXY_PAIRS"
echo "Number of exposure-outcome-analysis tasks: $N_PAIRS"

if [ "$N_OUTCOMES" -eq 0 ]; then
  echo "No outcomes to process. Exiting."
  exit 1
fi

if [ "$N_EXPOSURES" -eq 0 ]; then
  echo "No exposures to process. Exiting."
  exit 1
fi

if [ "$N_PAIRS" -eq 0 ]; then
  echo "No exposure-outcome pairs to process. Exiting."
  exit 1
fi

if [ "$N_PROXY_PAIRS" -eq 0 ]; then
  echo "No unique exposure-outcome pairs to process for proxy search. Exiting."
  exit 1
fi

############################################################
# STEP 2: PREPARE OUTCOMES (ARRAY)
############################################################

echo ""
echo "Step 2: Preparing outcomes (array job)"

OUTCOME_JOB_ID=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=prep_outcomes
#SBATCH --array=1-${N_OUTCOMES}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${OUTCOME_MEM}
#SBATCH --time=${OUTCOME_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/prep_outcomes_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/logs/prep_outcomes_%A_%a.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/02a_prepare_outcomes_array.R "$INPUT_DIR" "$OUTPUT_DIR" \${SLURM_ARRAY_TASK_ID}
EOF
)

echo "Submitted outcome array job with ID: $OUTCOME_JOB_ID"

############################################################
# STEP 3: PREPARE EXPOSURES (ARRAY)
############################################################

echo ""
echo "Step 3: Preparing exposures (array job)"

EXPOSURE_JOB_ID=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=prep_exposures
#SBATCH --array=1-${N_EXPOSURES}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${EXPOSURE_MEM}
#SBATCH --time=${EXPOSURE_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/prep_exposures_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/logs/prep_exposures_%A_%a.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/03a_prepare_exposures_array.R "$INPUT_DIR" "$OUTPUT_DIR" \${SLURM_ARRAY_TASK_ID}
EOF
)

echo "Submitted exposure array job with ID: $EXPOSURE_JOB_ID"

############################################################
# STEP 4: COMBINE OUTCOMES AND EXPOSURES
############################################################

echo ""
echo "Step 4: Combining prepared outcomes and exposures"

COMBINE_OUTCOME_JOB_ID=$(sbatch --parsable --dependency=afterok:${OUTCOME_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=combine_outcomes
#SBATCH --cpus-per-task=1
#SBATCH --mem=${COMBINE_MEM}
#SBATCH --time=${COMBINE_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/combine_outcomes_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/combine_outcomes_%j.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/02b_combine_outcomes.R "$OUTPUT_DIR"
EOF
)

COMBINE_EXPOSURE_JOB_ID=$(sbatch --parsable --dependency=afterok:${EXPOSURE_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=combine_exposures
#SBATCH --cpus-per-task=1
#SBATCH --mem=${COMBINE_MEM}
#SBATCH --time=${COMBINE_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/combine_exposures_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/combine_exposures_%j.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/03b_combine_exposures.R "$OUTPUT_DIR"
EOF
)

echo "Submitted combine outcomes job with ID : $COMBINE_OUTCOME_JOB_ID"
echo "Submitted combine exposures job with ID: $COMBINE_EXPOSURE_JOB_ID"

############################################################
# STEP 5: FIND AND COMBINE PROXY CANDIDATES
############################################################

echo ""
echo "Step 5: Finding proxy candidates after outcomes and exposures are ready"

PROXY_CANDIDATES_JOB_ID=$(sbatch --parsable --dependency=afterok:${COMBINE_OUTCOME_JOB_ID}:${COMBINE_EXPOSURE_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=find_proxies
#SBATCH --array=1-${N_PROXY_PAIRS}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${PROXY_MEM}
#SBATCH --time=${PROXY_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/find_proxies_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/logs/find_proxies_%A_%a.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/04a_find_proxy_candidates_array.R "$INPUT_DIR" "$OUTPUT_DIR" \${SLURM_ARRAY_TASK_ID}
EOF
)

COMBINE_PROXY_JOB_ID=$(sbatch --parsable --dependency=afterok:${PROXY_CANDIDATES_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=combine_proxies
#SBATCH --cpus-per-task=1
#SBATCH --mem=${COMBINE_MEM}
#SBATCH --time=${COMBINE_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/combine_proxies_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/combine_proxies_%j.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/04b_combine_proxy_candidates.R "$OUTPUT_DIR"
EOF
)

echo "Submitted proxy candidate array job with ID: $PROXY_CANDIDATES_JOB_ID"
echo "Submitted combine proxy job with ID         : $COMBINE_PROXY_JOB_ID"

############################################################
# STEP 6: APPLY PROXIES
############################################################

echo ""
echo "Step 6: Applying proxies"

APPLY_PROXY_JOB_ID=$(sbatch --parsable --dependency=afterok:${COMBINE_PROXY_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=apply_proxies
#SBATCH --array=1-${N_PROXY_PAIRS}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${PROXY_MEM}
#SBATCH --time=${PROXY_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/apply_proxies_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/logs/apply_proxies_%A_%a.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/04c_apply_proxies_array.R "$INPUT_DIR" "$OUTPUT_DIR" \${SLURM_ARRAY_TASK_ID}
EOF
)

echo "Submitted apply proxy array job with ID: $APPLY_PROXY_JOB_ID"

############################################################
# STEP 7: RUN AND COMBINE MR
############################################################

echo ""
echo "Step 7: Running MR"

MR_JOB_ID=$(sbatch --parsable --dependency=afterok:${APPLY_PROXY_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=run_mr
#SBATCH --array=1-${N_PAIRS}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${MR_MEM}
#SBATCH --time=${MR_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/run_mr_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/logs/run_mr_%A_%a.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/05a_run_mr_array.R "$INPUT_DIR" "$OUTPUT_DIR" \${SLURM_ARRAY_TASK_ID}
EOF
)

COMBINE_MR_JOB_ID=$(sbatch --parsable --dependency=afterok:${MR_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=combine_mr
#SBATCH --cpus-per-task=1
#SBATCH --mem=${COMBINE_MEM}
#SBATCH --time=${COMBINE_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/combine_mr_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/combine_mr_%j.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/05b_combine_mr_outputs.R "$OUTPUT_DIR"
EOF
)

echo "Submitted MR array job with ID: $MR_JOB_ID"
echo "Submitted combine MR job with ID: $COMBINE_MR_JOB_ID"

############################################################
# STEP 8: OPTIONAL COLOCALIZATION
############################################################

COLOC_JOB_ID=""
COMBINE_COLOC_JOB_ID=""
APPRAISAL_DEPENDENCY=$COMBINE_MR_JOB_ID

if [ "$RUN_COLOC" = "TRUE" ]; then
  echo ""
  echo "Step 8: Running colocalization"

  COLOC_JOB_ID=$(sbatch --parsable --dependency=afterok:${COMBINE_MR_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=run_coloc
#SBATCH --array=1-${N_PAIRS}
#SBATCH --cpus-per-task=1
#SBATCH --mem=${COLOC_MEM}
#SBATCH --time=${COLOC_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/run_coloc_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/logs/run_coloc_%A_%a.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/06a_run_colocalization_array.R "$INPUT_DIR" "$OUTPUT_DIR" \${SLURM_ARRAY_TASK_ID}
EOF
)

  COMBINE_COLOC_JOB_ID=$(sbatch --parsable --dependency=afterok:${COLOC_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=combine_coloc
#SBATCH --cpus-per-task=1
#SBATCH --mem=${COMBINE_MEM}
#SBATCH --time=${COMBINE_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/combine_coloc_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/combine_coloc_%j.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/06b_combine_colocalization.R "$OUTPUT_DIR"
EOF
)

  APPRAISAL_DEPENDENCY=$COMBINE_COLOC_JOB_ID
  echo "Submitted colocalization array job with ID: $COLOC_JOB_ID"
  echo "Submitted combine colocalization job with ID: $COMBINE_COLOC_JOB_ID"
else
  echo ""
  echo "Step 8: Skipping colocalization. Set RUN_COLOC=TRUE to enable it."
fi

############################################################
# STEP 9: EVIDENCE APPRAISAL
############################################################

APPRAISAL_JOB_ID=""

if [ "$RUN_APPRAISAL" = "TRUE" ]; then
  echo ""
  echo "Step 9: Appraising evidence"

  APPRAISAL_JOB_ID=$(sbatch --parsable --dependency=afterok:${APPRAISAL_DEPENDENCY} <<EOF
#!/bin/bash
#SBATCH --job-name=appraise_evidence
#SBATCH --cpus-per-task=1
#SBATCH --mem=${COMBINE_MEM}
#SBATCH --time=${COMBINE_TIME}
#SBATCH --output=${OUTPUT_DIR}/logs/appraise_evidence_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/appraise_evidence_%j.err

module load "$R_MODULE"
if [ -n "${R_LIBS_USER:-}" ]; then
  export R_LIBS_USER="${R_LIBS_USER:-}"
fi

Rscript Scripts/07_evidence_appraisal.R "$OUTPUT_DIR"
EOF
)

  echo "Submitted evidence appraisal job with ID: $APPRAISAL_JOB_ID"
else
  echo ""
  echo "Step 9: Skipping evidence appraisal. Set RUN_APPRAISAL=TRUE to enable it."
fi

echo ""
echo "========================================"
echo "Pipeline submitted successfully"
echo "Date: $(date)"
echo "Outcome array job       : $OUTCOME_JOB_ID"
echo "Exposure array job      : $EXPOSURE_JOB_ID"
echo "Combine outcomes job    : $COMBINE_OUTCOME_JOB_ID"
echo "Combine exposures job   : $COMBINE_EXPOSURE_JOB_ID"
echo "Proxy candidate job     : $PROXY_CANDIDATES_JOB_ID"
echo "Combine proxy job       : $COMBINE_PROXY_JOB_ID"
echo "Apply proxy job         : $APPLY_PROXY_JOB_ID"
echo "MR array job            : $MR_JOB_ID"
echo "Combine MR job          : $COMBINE_MR_JOB_ID"
echo "Colocalization job      : ${COLOC_JOB_ID:-<skipped>}"
echo "Combine colocalization  : ${COMBINE_COLOC_JOB_ID:-<skipped>}"
echo "Evidence appraisal job  : ${APPRAISAL_JOB_ID:-<skipped>}"
echo "========================================"
