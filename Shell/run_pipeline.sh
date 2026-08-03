#!/bin/bash

set -e

INPUT_DIR=$1
OUTPUT_DIR=$2
LOG_FILE="$OUTPUT_DIR/pipeline_log.txt"

if [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
  echo "Usage: bash Shell/run_pipeline.sh <input_dir> <output_dir>"
  exit 1
fi

mkdir -p "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR/logs"

exec > >(tee -a "$LOG_FILE") 2>&1

echo "========================================"
echo "Running pipeline"
echo "Date            : $(date)"
echo "Input directory : $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log file        : $LOG_FILE"
echo "========================================"

############################################################
# STEP 1: CHECK INPUTS (SERIAL)
############################################################

echo ""
echo "Step 1: Checking inputs"
Rscript Scripts/01_check_inputs.R "$INPUT_DIR" "$OUTPUT_DIR"

############################################################
# STEP 2: PREPARE OUTCOMES (ARRAY)
############################################################

echo ""
echo "Step 2: Preparing outcomes (array job)"

# Count number of unique active outcomes
N_OUTCOMES=$(Rscript -e "
plan <- read.csv(file.path('$INPUT_DIR','analysis_plan.csv'), stringsAsFactors=FALSE)
plan\$run <- tolower(as.character(plan\$run)) %in% c('true','t','1')
plan <- subset(plan, run)
length(unique(plan\$outcome))
")

echo "Number of outcome tasks: $N_OUTCOMES"

if [ "$N_OUTCOMES" -eq 0 ]; then
  echo "No outcomes to process. Exiting."
  exit 1
fi

# Submit array job
JOB_ID=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=prep_outcomes
#SBATCH --array=1-${N_OUTCOMES}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=${OUTPUT_DIR}/logs/prep_outcomes_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/logs/prep_outcomes_%A_%a.err

module load R

Rscript Scripts/02a_prepare_outcomes_array.R "$INPUT_DIR" "$OUTPUT_DIR" \${SLURM_ARRAY_TASK_ID}
EOF
)

echo "Submitted array job with ID: $JOB_ID"

############################################################
# STEP 3: WAIT FOR ARRAY + COMBINE (SERIAL)
############################################################

echo ""
echo "Step 3: Combining outcomes (after array completes)"

# Submit dependent job for combining
COMBINE_JOB_ID=$(sbatch --parsable --dependency=afterok:${JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=combine_outcomes
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=${OUTPUT_DIR}/logs/combine_outcomes_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/combine_outcomes_%j.err

module load R

Rscript Scripts/02b_combine_outcomes.R "$OUTPUT_DIR"
EOF
)

echo "Submitted combine job with ID: $COMBINE_JOB_ID"

echo ""
echo "========================================"
echo "Pipeline submitted successfully"
echo "Date: $(date)"
echo "Array job ID   : $JOB_ID"
echo "Combine job ID : $COMBINE_JOB_ID"
echo "========================================"

############################################################
# STEP 4: PREPARE EXPOSURES (ARRAY)
############################################################

echo ""
echo "Step 4: Preparing exposures (array job)"

mkdir -p "$OUTPUT_DIR/logs"

# Count number of unique active exposures
N_EXPOSURES=$(Rscript -e "
plan <- read.csv(file.path('$INPUT_DIR','analysis_plan.csv'), stringsAsFactors=FALSE)
plan\$run <- tolower(as.character(plan\$run)) %in% c('true','t','1')
plan <- subset(plan, run)
length(unique(plan\$exposure))
")

echo "Number of exposure tasks: $N_EXPOSURES"

if [ "$N_EXPOSURES" -eq 0 ]; then
  echo "No exposures to process. Exiting."
  exit 1
fi

# Submit exposure array job
EXPOSURE_JOB_ID=$(sbatch --parsable <<EOF
#!/bin/bash
#SBATCH --job-name=prep_exposures
#SBATCH --array=1-${N_EXPOSURES}
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=${OUTPUT_DIR}/logs/prep_exposures_%A_%a.out
#SBATCH --error=${OUTPUT_DIR}/logs/prep_exposures_%A_%a.err

module load R

Rscript Scripts/03a_prepare_exposures_array.R "$INPUT_DIR" "$OUTPUT_DIR" \${SLURM_ARRAY_TASK_ID}
EOF
)

echo "Submitted exposure array job with ID: $EXPOSURE_JOB_ID"

############################################################
# STEP 5: COMBINE EXPOSURES (SERIAL, AFTER ARRAY)
############################################################

echo ""
echo "Step 5: Combining exposures (after array completes)"

COMBINE_EXPOSURE_JOB_ID=$(sbatch --parsable --dependency=afterok:${EXPOSURE_JOB_ID} <<EOF
#!/bin/bash
#SBATCH --job-name=combine_exposures
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=${OUTPUT_DIR}/logs/combine_exposures_%j.out
#SBATCH --error=${OUTPUT_DIR}/logs/combine_exposures_%j.err

module load R

Rscript Scripts/03b_combine_exposures.R "$OUTPUT_DIR"
EOF
)

echo "Submitted combine exposures job with ID: $COMBINE_EXPOSURE_JOB_ID"
