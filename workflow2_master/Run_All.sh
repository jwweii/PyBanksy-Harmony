#!/bin/bash

# Master script to run the complete BANKSY pipeline
set -e  # Exit on any error

# =============================================================================
# GLOBAL CONFIGURATION PARAMETERS
# =============================================================================

# Main directories and files
PROJECT_DIR="/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-master/TEST/" # Set this to your project directory
CSV_FILE="$PROJECT_DIR/Xenium_tracking.csv"
BASE_OUTPUT_DIR="${PROJECT_DIR}/out"
mkdir -p "$BASE_OUTPUT_DIR"

# Set up logging
LOG_FILE="$BASE_OUTPUT_DIR/run.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "Starting BANKSY Pipeline Master Script"
echo "======================================="
echo "Log file: $LOG_FILE"
echo "Timestamp: $(date)"

# Project identifier (used in output naming)
PROJECT_NAME="Banksy_Run"

# =============================================================================
# CONSTANT PARAMETERS (NO NEED TO CHANGE)
# =============================================================================

# Workflow directory (where this script is located) and environment 
WORKFLOW_DIR="/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-master/workflow2_subworkflow"
CONDA_ENV="/diskmnt/Users2/simonmo/Software/miniforge3/envs/banksy"


# =============================================================================
# STEP 1: ANNDATA GENERATION PARAMETERS
# =============================================================================

# Input/Output
ANNDATA_OUTPUT_DIR="$BASE_OUTPUT_DIR/1-make_anndata"
XENIUM_PATH_COLUMN="katmai_path"
FILTER_COLUMN="Use"

# Processing parameters
PARALLEL_JOBS=8

# File paths for filtering
FILTER_CELLS_FILE="/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-master/resource/5K_gene_to_remove.tsv"
PYTHON_SCRIPT_ANNDATA="/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-master/script/Xenium_anndata_generator_one_sample.py"

# =============================================================================
# STEP 2: MERGE ANNDATA PARAMETERS
# =============================================================================

# Input/Output
MERGED_OUTPUT_DIR="$BASE_OUTPUT_DIR/2-Xenium_merged_anndata"
MERGED_FILE_PREFIX="$MERGED_OUTPUT_DIR/2-Xenium_merged_anndata"
MERGED_ANNDATA_FILE="${MERGED_FILE_PREFIX}_merged.h5ad"

# Merge parameters. Change these is sample overlaps
SAMPLES_PER_ROW=6
GRID_WIDTH=20000
GRID_HEIGHT=20000

# Gene filtering
REMOVE_GENES_TSV="/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-master/resource/5K_gene_to_remove.tsv"
PYTHON_SCRIPT_MERGE="/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-master/script/Xenium_anndata_merge_v2.py"

# =============================================================================
# STEP 3: BANKSY ANALYSIS PARAMETERS
# =============================================================================

# Input/Output
BANKSY_OUTPUT_DIR="$BASE_OUTPUT_DIR/3-banksy"
OUTPUT_PREFIX="$PROJECT_NAME"

# Core Banksy parameters
N_TOP_GENES=2000
K_GEOM=15
MAX_M=0
NBR_WEIGHT_DECAY="scaled_gaussian"
PCA_DIMS=20
LAMBDA_VALUE=0.8

# Harmony integration parameters
HARMONY_BATCH_KEY="dataset"
SAMPLE_ID_COLUMN="dataset"

# Clustering parameters
RUN_CLUSTERING="both"                           # Options: "leiden", "secuer", or "both"
LEIDEN_RESOLUTION=0.5
SECUER_RESOLUTION=1.2

# Visualization
PLOT="true"                                     # Set to "true" to generate spatial scatter plot

# Python script path
PYTHON_SCRIPT_BANKSY="/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-master/script/Xenium_Banksy_Harmony_v1.py"

# =============================================================================
# VALIDATION
# =============================================================================

echo "Validating configuration..."

# Check if required files exist
if [ ! -f "$CSV_FILE" ]; then
    echo "ERROR: Xenium tracking CSV file not found: $CSV_FILE"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT_ANNDATA" ]; then
    echo "ERROR: AnnData generation Python script not found: $PYTHON_SCRIPT_ANNDATA"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT_MERGE" ]; then
    echo "ERROR: Merge Python script not found: $PYTHON_SCRIPT_MERGE"
    exit 1
fi

if [ ! -f "$PYTHON_SCRIPT_BANKSY" ]; then
    echo "ERROR: BANKSY Python script not found: $PYTHON_SCRIPT_BANKSY"
    exit 1
fi

# Check conda environment
if ! conda info --envs | grep -q "$CONDA_ENV"; then
    echo "WARNING: Conda environment '$CONDA_ENV' may not exist"
fi

echo "Configuration validated successfully."

# =============================================================================
# CREATE OUTPUT DIRECTORIES
# =============================================================================

echo "Creating output directories..."
mkdir -p "$ANNDATA_OUTPUT_DIR"
mkdir -p "$MERGED_OUTPUT_DIR"
mkdir -p "$BANKSY_OUTPUT_DIR"

# =============================================================================
# STEP 1: GENERATE ANNDATA FILES
# =============================================================================

echo ""
echo "Step 1: Generating AnnData files from Xenium data..."
echo "====================================================="

# Export parameters for Step 1
export CSV_FILE="$CSV_FILE"
export OUTPUT_DIR="$ANNDATA_OUTPUT_DIR"
export XENIUM_PATH_COLUMN="$XENIUM_PATH_COLUMN"
export FILTER_COLUMN="$FILTER_COLUMN"
export PARALLEL_JOBS="$PARALLEL_JOBS"
export FILTER_CELLS="$FILTER_CELLS_FILE"
export PYTHON_SCRIPT="$PYTHON_SCRIPT_ANNDATA"
export CONDA_ENV="$CONDA_ENV"

bash "$WORKFLOW_DIR/1-run_making_anndata_parallel_csv.sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 1 failed. Exiting."
    exit 1
fi

echo "Step 1 completed successfully."

# =============================================================================
# STEP 2: MERGE ANNDATA FILES
# =============================================================================

echo ""
echo "Step 2: Merging AnnData files..."
echo "================================="

# Export parameters for Step 2
export INPUT_DIR="$ANNDATA_OUTPUT_DIR"
export OUTPUT_DIR="$MERGED_OUTPUT_DIR"
export SAMPLES_PER_ROW="$SAMPLES_PER_ROW"
export GRID_WIDTH="$GRID_WIDTH"
export GRID_HEIGHT="$GRID_HEIGHT"
export REMOVE_GENES_TSV="$REMOVE_GENES_TSV"
export PYTHON_SCRIPT="$PYTHON_SCRIPT_MERGE"
export CONDA_ENV="$CONDA_ENV"

bash "$WORKFLOW_DIR/2-run_merge_anndata.sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 2 failed. Exiting."
    exit 1
fi

# Check if merged file was created
if [ ! -f "$MERGED_ANNDATA_FILE" ]; then
    echo "ERROR: Merged AnnData file not found at $MERGED_ANNDATA_FILE"
    exit 1
fi

echo "Step 2 completed successfully."

# =============================================================================
# STEP 3: RUN BANKSY ANALYSIS
# =============================================================================

echo ""
echo "Step 3: Running BANKSY analysis..."
echo "=================================="

# Export parameters for Step 3
export INPUT_FILE="$MERGED_ANNDATA_FILE"
export OUTPUT_DIR="$BANKSY_OUTPUT_DIR"
export OUTPUT_PREFIX="$OUTPUT_PREFIX"
export N_TOP_GENES="$N_TOP_GENES"
export K_GEOM="$K_GEOM"
export MAX_M="$MAX_M"
export NBR_WEIGHT_DECAY="$NBR_WEIGHT_DECAY"
export PCA_DIMS="$PCA_DIMS"
export LAMBDA_VALUE="$LAMBDA_VALUE"
export HARMONY_BATCH_KEY="$HARMONY_BATCH_KEY"
export SAMPLE_ID_COLUMN="$SAMPLE_ID_COLUMN"
export RUN_CLUSTERING="$RUN_CLUSTERING"
export LEIDEN_RESOLUTION="$LEIDEN_RESOLUTION"
export SECUER_RESOLUTION="$SECUER_RESOLUTION"
export PLOT="$PLOT"
export PYTHON_SCRIPT="$PYTHON_SCRIPT_BANKSY"
export CONDA_ENV="$CONDA_ENV"

bash "$WORKFLOW_DIR/3-run_Banksy.sh"

if [ $? -ne 0 ]; then
    echo "ERROR: Step 3 failed. Exiting."
    exit 1
fi

echo "Step 3 completed successfully."

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "BANKSY Pipeline completed successfully!"
echo "======================================="
echo "Project: $PROJECT_NAME"
echo "CSV file used: $CSV_FILE"
echo "Final output directory: $BANKSY_OUTPUT_DIR"
echo "Merged AnnData file: $MERGED_ANNDATA_FILE"
echo "Log file: $LOG_FILE"
echo ""
echo "Pipeline summary:"
echo "- Processed $(tail -n +2 "$CSV_FILE" | wc -l) samples from CSV"
echo "- AnnData files: $ANNDATA_OUTPUT_DIR"
echo "- Merged file: $MERGED_ANNDATA_FILE"
echo "- BANKSY results: $BANKSY_OUTPUT_DIR"
echo ""
echo "Check individual log files in each output directory for details."
echo "Full pipeline log available at: $LOG_FILE"
