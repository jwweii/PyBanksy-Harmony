import argparse
import pandas as pd
import spatialdata as sd
from spatialdata_io import xenium
import anndata as ad
import scanpy as sc
import logging
import os

def setup_logging(output_dir, sample_name):
    log_file = f"{output_dir}/{sample_name}_xenium_anndata.log"
    logging.basicConfig(
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ],
        format="%(asctime)s - %(levelname)s - %(message)s",
        level=logging.INFO
    )
    logging.info(f"Logging initialized for sample {sample_name}.")

def remove_cells(adata, cell_removal_csv):
    """Remove specific cells from an AnnData object based on a CSV file."""
    try:
        removal_data = pd.read_csv(cell_removal_csv, skiprows=3, header=None)
        cells_to_remove = removal_data.iloc[:, 0].str.strip().tolist()
        valid_cells_to_remove = [cell for cell in cells_to_remove if cell in adata.obs_names]

        if len(valid_cells_to_remove) > 0:
            adata = adata[~adata.obs_names.isin(valid_cells_to_remove)].copy()
            logging.info(f"Removed {len(valid_cells_to_remove)} cells from the AnnData object.")
        else:
            logging.info("No matching cells found to remove.")

    except Exception as e:
        logging.error(f"Error removing cells: {e}")

    return adata

def process_single_sample(xenium_path, sample_name, casette_id, output_dir, cell_removal_csv=None):
    setup_logging(output_dir, sample_name)
    
    try:
        zarr_path = f"{output_dir}/{sample_name}.zarr"
        h5ad_path = f"{output_dir}/{sample_name}.h5ad"

        logging.info(f"Processing sample: {sample_name}")

        # Generate SpatialData object from Xenium path
        sdata = xenium(xenium_path)

        if not os.path.exists(zarr_path):
            logging.info(f"Writing Zarr output to {zarr_path}.")
            sdata.write(zarr_path)
        else:
            logging.info(f"Zarr output already exists at {zarr_path}. Skipping.")

        # Extract the AnnData table
        adata = sdata.tables["table"]
        adata.obs["dataset"] = sample_name
        adata.obs["casette"] = casette_id
        adata.obs["barcode"] = adata.obs["cell_id"]
        adata.obs.set_index('barcode', inplace=True)
        logging.info(f"Sample {sample_name} AnnData shape: {adata.shape}")

        # Perform QC and filter cells and genes
        sc.pp.filter_cells(adata, min_counts=20)
        sc.pp.filter_cells(adata, min_genes=10)
        logging.info(f"Sample {sample_name} post-QC AnnData shape: {adata.shape}")

        # Check for cell removal
        if cell_removal_csv and os.path.exists(cell_removal_csv):
            logging.info(f"Removing cells for sample {sample_name} using {cell_removal_csv}.")
            adata = remove_cells(adata, cell_removal_csv)
        elif cell_removal_csv:
            logging.warning(f"Cell removal CSV file not found: {cell_removal_csv}")

        # Save the AnnData table to H5AD format
        adata.write(h5ad_path)
        logging.info(f"Sample {sample_name} saved as {h5ad_path}")

    except Exception as e:
        logging.error(f"Error processing sample {sample_name}: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Generate Xenium AnnData for a single sample.")
    parser.add_argument("--xenium_path", required=True, help="Path to the Xenium output directory.")
    parser.add_argument("--sample_name", required=True, help="Sample name (Section ID).")
    parser.add_argument("--casette_id", required=True, help="Casette ID.")
    parser.add_argument("--output_dir", required=True, help="Directory to save the outputs.")
    parser.add_argument("--cell_removal_csv", help="Path to the CSV file specifying cells to remove.")

    args = parser.parse_args()

    print("Input Arguments:")
    print(f"Xenium Path: {args.xenium_path}")
    print(f"Sample Name: {args.sample_name}")
    print(f"Casette ID: {args.casette_id}")
    print(f"Output Directory: {args.output_dir}")
    print(f"Cell Removal CSV: {args.cell_removal_csv}")

    # Call the function with the parsed arguments
    process_single_sample(args.xenium_path, args.sample_name, args.casette_id, args.output_dir, args.cell_removal_csv)

if __name__ == "__main__":
    main()
