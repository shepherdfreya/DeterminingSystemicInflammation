#!/bin/bash

set -e

INPUT_DIR=$1
OUTPUT_DIR=$2
LOG_FILE="$OUTPUT_DIR/pipeline_log.txt"

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
  echo "Usage: bash Shell/run_pipeline_local.sh <input_dir> <output_dir>"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================"
echo "Running local pipeline"
echo "Date            : $(date)"
echo "Input directory : $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log file        : $LOG_FILE"
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

N_PAIRS=$(Rscript -e "
plan <- read.csv(file.path('$INPUT_DIR','analysis_plan.csv'), stringsAsFactors=FALSE)
plan\$run <- tolower(as.character(plan\$run)) %in% c('true','t','1')
plan <- subset(plan, run)
cat(nrow(plan))
")

echo "Number of exposure-outcome pairs: $N_PAIRS"

if [ "$N_PAIRS" -eq 0 ]; then
  echo "No exposure-outcome pairs to process. Exiting."
  exit 1
fi

for TASK_ID in $(seq 1 "$N_PAIRS"); do
  echo ""
  echo "Finding proxy candidates for pair task $TASK_ID of $N_PAIRS"
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

for TASK_ID in $(seq 1 "$N_PAIRS"); do
  echo ""
  echo "Applying proxies for pair task $TASK_ID of $N_PAIRS"
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

echo ""
echo "========================================"
echo "Local pipeline completed successfully"
echo "Date: $(date)"
echo "========================================"
