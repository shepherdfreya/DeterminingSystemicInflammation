#!/bin/bash

set -euo pipefail

if [ "$#" -lt 2 ]; then
  echo "Usage: bash Shell/run_pipeline_local.sh <input_dir> <output_dir>"
  exit 1
fi

INPUT_DIR=$1
OUTPUT_DIR=$2
LOG_FILE="$OUTPUT_DIR/pipeline_log.txt"

RUN_COLOC=${RUN_COLOC:-FALSE}
RUN_APPRAISAL=${RUN_APPRAISAL:-TRUE}

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
echo "Running local pipeline"
echo "Date            : $(date)"
echo "Input directory : $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log file        : $LOG_FILE"
echo "R_LIBS_USER     : ${R_LIBS_USER:-<not set>}"
echo "Run coloc       : $RUN_COLOC"
echo "Run appraisal   : $RUN_APPRAISAL"
echo "========================================"

############################################################
# STEP 1: CHECK INPUTS
############################################################

echo ""
echo "Step 1: Checking inputs"
Rscript Scripts/01_check_inputs.R "$INPUT_DIR" "$OUTPUT_DIR"

############################################################
# STEP 2: PREPARE OUTCOMES LOCALLY
############################################################

echo ""
echo "Step 2: Preparing outcomes locally"

N_OUTCOMES=$(Rscript -e "
plan <- read.csv(file.path('$INPUT_DIR','analysis_plan.csv'), stringsAsFactors=FALSE)
plan\$run <- tolower(as.character(plan\$run)) %in% c('true','t','1')
plan <- subset(plan, run)
cat(length(unique(plan\$outcome)))
")

echo "Number of outcome tasks: $N_OUTCOMES"

if [ "$N_OUTCOMES" -eq 0 ]; then
  echo "No outcomes to process. Exiting."
  exit 1
fi

for TASK_ID in $(seq 1 "$N_OUTCOMES"); do
  echo ""
  echo "Preparing outcome task $TASK_ID of $N_OUTCOMES"
  Rscript Scripts/02a_prepare_outcomes_array.R "$INPUT_DIR" "$OUTPUT_DIR" "$TASK_ID"
done

############################################################
# STEP 3: COMBINE OUTCOMES
############################################################

echo ""
echo "Step 3: Combining outcomes"
Rscript Scripts/02b_combine_outcomes.R "$OUTPUT_DIR"

############################################################
# STEP 4: PREPARE EXPOSURES LOCALLY
############################################################

echo ""
echo "Step 4: Preparing exposures locally"

N_EXPOSURES=$(Rscript -e "
plan <- read.csv(file.path('$INPUT_DIR','analysis_plan.csv'), stringsAsFactors=FALSE)
plan\$run <- tolower(as.character(plan\$run)) %in% c('true','t','1')
plan <- subset(plan, run)
cat(length(unique(plan\$exposure)))
")

echo "Number of exposure tasks: $N_EXPOSURES"

if [ "$N_EXPOSURES" -eq 0 ]; then
  echo "No exposures to process. Exiting."
  exit 1
fi

for TASK_ID in $(seq 1 "$N_EXPOSURES"); do
  echo ""
  echo "Preparing exposure task $TASK_ID of $N_EXPOSURES"
  Rscript Scripts/03a_prepare_exposures_array.R "$INPUT_DIR" "$OUTPUT_DIR" "$TASK_ID"
done

############################################################
# STEP 5: COMBINE EXPOSURES
############################################################

echo ""
echo "Step 5: Combining exposures"
Rscript Scripts/03b_combine_exposures.R "$OUTPUT_DIR"

############################################################
# STEP 6: FIND PROXY CANDIDATES LOCALLY
############################################################

echo ""
echo "Step 6: Finding proxy candidates locally"

N_PROXY_PAIRS=$(Rscript -e "
plan <- read.csv(file.path('$INPUT_DIR','analysis_plan.csv'), stringsAsFactors=FALSE)
plan\$run <- tolower(as.character(plan\$run)) %in% c('true','t','1')
plan <- subset(plan, run)
cat(nrow(unique(plan[, c('exposure', 'outcome'), drop=FALSE])))
")

N_PAIRS=$(Rscript -e "
plan <- read.csv(file.path('$INPUT_DIR','analysis_plan.csv'), stringsAsFactors=FALSE)
plan\$run <- tolower(as.character(plan\$run)) %in% c('true','t','1')
plan <- subset(plan, run)
cat(nrow(plan))
")

echo "Number of unique proxy-search pairs: $N_PROXY_PAIRS"
echo "Number of exposure-outcome-analysis tasks: $N_PAIRS"

if [ "$N_PROXY_PAIRS" -eq 0 ]; then
  echo "No exposure-outcome pairs to process for proxy search. Exiting."
  exit 1
fi

for TASK_ID in $(seq 1 "$N_PROXY_PAIRS"); do
  echo ""
  echo "Finding proxy candidates for unique pair task $TASK_ID of $N_PROXY_PAIRS"
  Rscript Scripts/04a_find_proxy_candidates_array.R "$INPUT_DIR" "$OUTPUT_DIR" "$TASK_ID"
done

############################################################
# STEP 7: COMBINE PROXY CANDIDATES
############################################################

echo ""
echo "Step 7: Combining proxy candidates"
Rscript Scripts/04b_combine_proxy_candidates.R "$OUTPUT_DIR"

############################################################
# STEP 8: APPLY PROXIES LOCALLY
############################################################

echo ""
echo "Step 8: Applying proxies locally"

for TASK_ID in $(seq 1 "$N_PROXY_PAIRS"); do
  echo ""
  echo "Applying proxies for unique pair task $TASK_ID of $N_PROXY_PAIRS"
  Rscript Scripts/04c_apply_proxies_array.R "$INPUT_DIR" "$OUTPUT_DIR" "$TASK_ID"
done

############################################################
# STEP 9: RUN MR LOCALLY
############################################################

echo ""
echo "Step 9: Running MR locally"

for TASK_ID in $(seq 1 "$N_PAIRS"); do
  echo ""
  echo "Running MR for pair task $TASK_ID of $N_PAIRS"
  Rscript Scripts/05a_run_mr_array.R "$INPUT_DIR" "$OUTPUT_DIR" "$TASK_ID"
done

############################################################
# STEP 10: COMBINE MR OUTPUTS
############################################################

echo ""
echo "Step 10: Combining MR outputs"
Rscript Scripts/05b_combine_mr_outputs.R "$OUTPUT_DIR"

############################################################
# STEP 11: OPTIONAL COLOCALIZATION
############################################################

if [ "$RUN_COLOC" = "TRUE" ]; then
  echo ""
  echo "Step 11: Running colocalization locally"

  for TASK_ID in $(seq 1 "$N_PAIRS"); do
    echo ""
    echo "Running colocalization for pair task $TASK_ID of $N_PAIRS"
    Rscript Scripts/06a_run_colocalization_array.R "$INPUT_DIR" "$OUTPUT_DIR" "$TASK_ID"
  done

  echo ""
  echo "Step 12: Combining colocalization outputs"
  Rscript Scripts/06b_combine_colocalization.R "$OUTPUT_DIR"
else
  echo ""
  echo "Step 11: Skipping colocalization. Set RUN_COLOC=TRUE to enable it."
fi

############################################################
# STEP 12/13: EVIDENCE APPRAISAL
############################################################

if [ "$RUN_APPRAISAL" = "TRUE" ]; then
  if [ "$RUN_COLOC" = "TRUE" ]; then
    echo ""
    echo "Step 13: Appraising evidence"
  else
    echo ""
    echo "Step 12: Appraising evidence"
  fi

  Rscript Scripts/07_evidence_appraisal.R "$OUTPUT_DIR"
else
  if [ "$RUN_COLOC" = "TRUE" ]; then
    echo ""
    echo "Step 13: Skipping evidence appraisal. Set RUN_APPRAISAL=TRUE to enable it."
  else
    echo ""
    echo "Step 12: Skipping evidence appraisal. Set RUN_APPRAISAL=TRUE to enable it."
  fi
fi

echo ""
echo "========================================"
echo "Local pipeline completed successfully"
echo "Date: $(date)"
echo "========================================"
