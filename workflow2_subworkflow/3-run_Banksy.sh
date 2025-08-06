#!/bin/bash

# Define paths and parameters - use environment variables if provided
INPUT_FILE="${INPUT_FILE}"
OUTPUT_DIR="${OUTPUT_DIR}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-banksy}"
LOG_FILE="${OUTPUT_DIR}/banksy_pipeline.log"

# BANKSY algorithm parameters
N_TOP_GENES="${N_TOP_GENES:-2000}"
K_GEOM="${K_GEOM:-15}"
MAX_M="${MAX_M:-0}"
NBR_WEIGHT_DECAY="${NBR_WEIGHT_DECAY:-scaled_gaussian}"
PCA_DIMS="${PCA_DIMS:-20}"
LAMBDA_VALUE="${LAMBDA_VALUE:-0.8}"

# Integration and clustering parameters
HARMONY_BATCH_KEY="${HARMONY_BATCH_KEY:-dataset}"
SAMPLE_ID_COLUMN="${SAMPLE_ID_COLUMN:-dataset}"
RUN_CLUSTERING="${RUN_CLUSTERING:-both}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.5}"
SECUER_RESOLUTION="${SECUER_RESOLUTION:-1.2}"

# Visualization
PLOT="${PLOT:-true}"

# File paths
PYTHON_SCRIPT="${PYTHON_SCRIPT:-/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony/src/Xenium_Banksy_Harmony_v1.py}"
CONDA_ENV="${CONDA_ENV:-/diskmnt/Users2/simonmo/Software/miniforge3/envs/banksy}"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Activate Conda environment
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

# Run the Python script
python "$PYTHON_SCRIPT" \
    --input_merged_anndata "$INPUT_FILE" \
    --output_dir "$OUTPUT_DIR" \
    --output_prefix "$OUTPUT_PREFIX" \
    --n_top_genes "$N_TOP_GENES" \
    --k_geom "$K_GEOM" \
    --max_m "$MAX_M" \
    --nbr_weight_decay "$NBR_WEIGHT_DECAY" \
    --pca_dims "$PCA_DIMS" \
    --lambda_list "$LAMBDA_VALUE" \
    --harmony_batch_key "$HARMONY_BATCH_KEY" \
    --run_clustering "$RUN_CLUSTERING" \
    --leiden_resolution "$LEIDEN_RESOLUTION" \
    --secuer_resolution "$SECUER_RESOLUTION" \
    --sample_id_column "$SAMPLE_ID_COLUMN" \
    $( [ "$PLOT" = "true" ] && echo "--plot" ) \
    &> "$LOG_FILE"

# Deactivate Conda environment
conda deactivate

# Print completion message
echo "Banksy pipeline completed. Check the log file at: $LOG_FILE"
