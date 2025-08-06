#!/bin/bash

# Define parameters - use environment variables if provided, otherwise defaults
INPUT_DIR="${INPUT_DIR}"
OUTPUT_DIR="${OUTPUT_DIR}" # "/diskmnt/Projects/HTAN_prostate_analysis/Analysis/Xenium/2-integrate/Primary5K/4b-BanksyPy/out/2-Xenium_merged_anndata/"
OUTPUT_FILE="${OUTPUT_DIR}/2-Xenium_merged_anndata"  # output file name is "${OUTPUT_FILE}_merged.h5ad"
LOG_FILE="${OUTPUT_DIR}/merge_and_stagger.log"

mkdir -p "${OUTPUT_DIR}"

# Merge parameters
SAMPLES_PER_ROW="${SAMPLES_PER_ROW:-6}"
GRID_WIDTH="${GRID_WIDTH:-15000}"
GRID_HEIGHT="${GRID_HEIGHT:-15000}"

# File paths
REMOVE_GENES_TSV="${REMOVE_GENES_TSV:-/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony/5K_gene_to_remove.tsv}"
PYTHON_SCRIPT="${PYTHON_SCRIPT:-/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony/src/Xenium_anndata_merge_v2.py}"

# Conda environment
CONDA_ENV="${CONDA_ENV:-/diskmnt/Users2/simonmo/Software/miniforge3/envs/banksy}"

# Activate the Python environment
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

# Run the Python script
python "$PYTHON_SCRIPT" \
    --input_dir "$INPUT_DIR" \
    --output_prefix "$OUTPUT_FILE" \
    --remove_genes_tsv "$REMOVE_GENES_TSV" \
    --samples_per_row "$SAMPLES_PER_ROW" \
    --grid_width "$GRID_WIDTH" \
    --grid_height "$GRID_HEIGHT" \
    &> "$LOG_FILE"

# Deactivate the Python environment
conda deactivate
