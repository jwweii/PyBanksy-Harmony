#!/bin/bash

# Set the input arguments
CSV_FILE="/diskmnt/Projects/MetNet_analysis_2/Colorectal/Xenium/Tracking_sheet/mCRC_tracking_Xenium_20241212.csv"
XENIUM_PATH_COLUMN="Output file path"
OUTPUT_DIR="/diskmnt/Projects/MetNet_analysis/Colorectal/epeng/Xenium/Xenium_anndata"
FILTER_COLUMN="Panel_5k"

mkdir -p "$OUTPUT_DIR"

# Activate the Python environment
eval "$(conda shell.bash hook)"
conda activate scanpy_env

# Run the Python script
python /diskmnt/Projects/Users/Evan.p/scripts/Python/Banksy/src/Xenium_anndata_generator_v1.py \
    --csv_file "$CSV_FILE" \
    --Xenium_path "$XENIUM_PATH_COLUMN" \
    --output_dir "$OUTPUT_DIR" \
    --filter_column "$FILTER_COLUMN" 

# Deactivate the Python environment
conda deactivate
