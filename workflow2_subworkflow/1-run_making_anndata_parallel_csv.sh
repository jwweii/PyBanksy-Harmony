#!/bin/bash

# Set the input arguments - use environment variables if provided, otherwise defaults
CSV_FILE="${CSV_FILE}"
XENIUM_PATH_COLUMN="${XENIUM_PATH_COLUMN:-katmai_path}"
OUTPUT_DIR="${OUTPUT_DIR}"
FILTER_COLUMN="${FILTER_COLUMN:-Use}"
FILTER_CELLS="${FILTER_CELLS:-/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-parallel/resource/5K_gene_to_remove.tsv}"
PYTHON_SCRIPT="${PYTHON_SCRIPT:-/diskmnt/Users2/simonmo/Tools/BANKSY/PyBanksy-Harmony-parallel/script/Xenium_anndata_generator_one_sample.py}"

# Number of parallel jobs (adjust based on your system)
PARALLEL_JOBS="${PARALLEL_JOBS:-8}"

# Conda environment
CONDA_ENV="${CONDA_ENV:-/diskmnt/Users2/simonmo/Software/miniforge3/envs/banksy}"

mkdir -p "$OUTPUT_DIR"

# Activate the Python environment
eval "$(conda shell.bash hook)"
conda activate "$CONDA_ENV"

# Function to process a single sample
process_sample() {
    local section_id="$1"
    local casette_id="$2"
    local xenium_path="$3"
    local use_flag="$4"
    local removal_csv="$5"
    local output_dir="$6"
    local python_script="$7"
    
    # Skip if not marked for use
    if [[ "$use_flag" != "Yes" ]]; then
        echo "Skipping sample $section_id (Use = $use_flag)"
        return
    fi
    
    # Check for cell removal file
    local cell_removal_csv=""
    if [[ -f "$removal_csv" ]]; then
        cell_removal_csv=$(grep -E "^${section_id}," "$removal_csv" | cut -d, -f2)
    fi
    
    # Build command arguments
    local cmd_args="--xenium_path \"$xenium_path\" --sample_name \"$section_id\" --casette_id \"$casette_id\" --output_dir \"$output_dir\""
    if [[ -n "$cell_removal_csv" && -f "$cell_removal_csv" ]]; then
        cmd_args="$cmd_args --cell_removal_csv \"$cell_removal_csv\""
    fi
    
    # Run the Python script
    echo "Processing sample: $section_id"
    eval "python \"$python_script\" $cmd_args"
}

# Export the function for parallel
export -f process_sample

# Read CSV file and process with GNU parallel
echo "Starting parallel processing with $PARALLEL_JOBS jobs..."

# Debug: First let's see what we're actually getting
echo "Debugging CSV content:"
tail -n +2 "$CSV_FILE" | head -5 | while IFS=',' read -r col1 col2 col3 col4; do
    echo "Col1: '$col1', Col2: '$col2', Col3: '$col3', Col4: '$col4'"
done

# Skip header, filter by Use=Yes, and process with parallel using column splitting
tail -n +2 "$CSV_FILE" | awk -F',' '$4~/^Yes/ {print $1","$2","$3","$4}' | \
    parallel -j "$PARALLEL_JOBS" --line-buffer --colsep ',' \
    process_sample {1} {2} {3} {4} "$FILTER_CELLS" "$OUTPUT_DIR" "$PYTHON_SCRIPT" \
    2>&1 | tee "$OUTPUT_DIR/parallel_run.log"

echo "All samples processed. Check individual log files in $OUTPUT_DIR for details."

# Deactivate the Python environment
conda deactivate
