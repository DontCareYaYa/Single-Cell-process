import argparse
import os
import scanpy as sc
import pandas as pd
import anndata as ad
import numpy as np
import bbknn

sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
# scanpy==1.9.9.dev1+g596ed01e anndata==0.10.2 umap==0.5.7 numpy==1.26.4 scipy==1.14.1 pandas==2.2.3 scikit-learn==1.5.2 statsmodels==0.14.4 pynndescent==0.5.13

def check_file_exists(file_path, file_type):
    """Check if a required file exists."""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"{file_type} file '{file_path}' not found.")

def read_meta_data(meta_file):
    """Read the meta data file."""
    check_file_exists(meta_file, "Meta data")
    meta_data = pd.read_csv(meta_file, sep="\t", header=None)
    meta_data.columns = ["file_name"] + [f"meta_{i}" for i in range(1, meta_data.shape[1])]
    return meta_data

def load_sample(input_dir, sample_name):
    """Load 10X data and add sample-specific prefixes."""
    sample_dir = os.path.join(input_dir, sample_name)
    if not os.path.isdir(sample_dir):
        raise FileNotFoundError(f"Sample directory '{sample_dir}' not found.")
    print(f"Loading 10X data from '{sample_dir}'...")
    adata = sc.read_10x_mtx(sample_dir, var_names="gene_symbols", cache=True) 
    adata.var_names_make_unique()
    
    # Add sample prefix to cell names
    print(f"Adding sample prefix '{sample_name}_' to cell names...")
    adata.obs_names = [f"{sample_name}_{barcode}" for barcode in adata.obs_names]

    # Ensure unique cell names
    if not adata.obs_names.is_unique:
        print("Non-unique cell names detected. Making them unique...")
        adata.obs_names_make_unique() 
        
    # Replace underscores with dashes in variable names
#   if any('_' in name for name in adata.var_names):
#        print("Replacing underscores ('_') with dashes ('-') in variable names...")
#        adata.var.index = adata.var.index.str.replace('_', '-', regex=False)
#        adata.var_names_make_unique()
    
    # Add sample prefix to cell names    
    # adata.obs_names = [f"{sample_name}_{barcode}" for barcode in adata.obs_names]
    return adata
    
def add_meta_data(adata, meta_row):
    """Add meta data information to AnnData.obs."""
    print("Adding meta data to the sample...")
    for col_name, value in meta_row.items():
        if col_name != "file_name":
            adata.obs[col_name] = value
    return adata

def add_copykat_predictions(adata, sample_name, copykat_dir):
    """Add CopyKAT predictions if available."""
    if not copykat_dir:
        print("CopyKAT directory not provided. Skipping...")
        return adata

    copykat_file = os.path.join(copykat_dir, f"{sample_name}_copykat_prediction.txt")
    if os.path.exists(copykat_file):
        print(f"Adding CopyKAT predictions for sample '{sample_name}'...")
        copykat_data = pd.read_csv(copykat_file, sep="\t", header=0)
        
        # Add sample name prefix to the first column (barcodes)
        copykat_data.iloc[:, 0] = copykat_data.iloc[:, 0].apply(lambda x: f"{sample_name}_{x}")
        
        # Set the first column as index and remove it from the data
        copykat_data.index = copykat_data.iloc[:, 0] # First column as index
        copykat_data = copykat_data.iloc[:, 1:] # Remove first column from datas
        copykat_data = copykat_data.applymap(str)  # Convert all values to string
        adata.obs = adata.obs.merge(copykat_data, left_index=True, right_index=True, how="left")
    else:
        print(f"Warning: CopyKAT file '{copykat_file}' not found. Skipping.")
    return adata

def perform_qc(adata, qc_params):
    """Perform quality control filtering."""
    print("Performing QC filtering...")
    sc.pp.filter_cells(adata, min_genes=qc_params["min_genes"])
    sc.pp.filter_genes(adata, min_cells=qc_params["min_cells"])
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    adata = adata[
        (adata.obs.n_genes_by_counts > qc_params["min_n_genes"]) &
        (adata.obs.n_genes_by_counts < qc_params["max_n_genes"]) &
        (adata.obs.total_counts > qc_params["min_total"]) &
        (adata.obs.total_counts < qc_params["max_total"]) &
        (adata.obs.pct_counts_mt < qc_params["max_mt"]), :
    ]
    return adata

def ensure_categorical_and_clean_feature_names(adata):
    """Ensure categorical columns in obs and var are properly set, and clean feature names."""
    # Replace underscores with dashes in feature names
    if any('_' in name for name in adata.var_names):
        print("Replacing underscores ('_') with dashes ('-') in feature names...")
        adata.var.index = adata.var.index.str.replace('_', '-', regex=False)
        adata.var_names_make_unique()

    # Ensure categorical data in obs
    for col in adata.obs.columns:
        if adata.obs[col].dtype == "object":
            print(f"Converting '{col}' in adata.obs to categorical...")
            adata.obs[col] = adata.obs[col].astype("category")

    # Ensure categorical data in var
    for col in adata.var.columns:
        if adata.var[col].dtype == "object":
            print(f"Converting '{col}' in adata.var to categorical...")
            adata.var[col] = adata.var[col].astype("category")

    return adata

# Load and QC filter
# Add meta data to each sample
# Add CopyKAT predictions if available
def process_samples(input_dir, meta_data, qc_params, copykat_dir):
    """Process each sample and return a list of AnnData objects."""
    adata_list = []
    for _, meta_row in meta_data.iterrows():
        sample_name = meta_row["file_name"]
        try:
            adata = load_sample(input_dir, sample_name)
            adata = add_meta_data(adata, meta_row)
            adata = add_copykat_predictions(adata, sample_name, copykat_dir)
            adata = perform_qc(adata, qc_params)
            adata_list.append(adata)
        except Exception as e:
            print(f"Error processing sample '{sample_name}': {e}")
            continue
    if not adata_list:
        raise ValueError("No valid samples were processed.")
    return adata_list

def main(args):
    # Validate paths
    check_file_exists(args.meta_data, "Meta data file")
    if not os.path.isdir(args.input_dir):
        raise FileNotFoundError(f"Input directory '{args.input_dir}' not found.")
    
    # Load meta data
    meta_data = read_meta_data(args.meta_data)

    # Define QC parameters 
    qc_params = {
        "min_genes": args.min_genes,
        "min_cells": args.min_cells,
        "min_n_genes": args.min_n_genes_by_counts,
        "max_n_genes": args.max_n_genes_by_counts,
        "min_total": args.min_total_counts,
        "max_total": args.max_total_counts,
        "max_mt": args.max_mt,
    }
    
    # Process each sample
    adata_list = process_samples(args.input_dir, meta_data, qc_params, args.copykat) 
    # combined_adata = adata_list[0] if len(adata_list) == 1 else ad.concat(adata_list, join="outer", label="batch")
    
    # Combine all samples
    if len(adata_list) == 1:
        combined_adata = adata_list[0]
        print("Only one sample found. No merging required.")
    else:
        sample_names = [meta_row["file_name"] for _, meta_row in meta_data.iterrows()]
        combined_adata = ad.concat(
        adata_list,
        join="outer",
        label="orig.ident",  # Set orig.ident as label
        keys=sample_names,   # Use sample names as keys
        index_unique=None    # Keep original cell names unique
    )
    
    # Ensure categorical and string data in obs and var
    # Ensure meta data is exported correctly and clean feature names
    
    combined_adata = ensure_categorical_and_clean_feature_names(combined_adata)
    
    # Ensure all columns in obs and var are strings
    combined_adata.obs = combined_adata.obs.astype(str)
    combined_adata.var = combined_adata.var.astype(str)
    print(f"Combined {len(adata_list)} samples.")  
    
    # Preprocessing
    sc.pp.normalize_total(combined_adata, target_sum=1e4)
    sc.pp.log1p(combined_adata)
    sc.pp.highly_variable_genes(
        combined_adata,
        min_mean=0.0125,
        max_mean=3,
        min_disp=0.5,
        n_top_genes=args.n_top_genes
    )
    
    combined_adata = combined_adata[:, combined_adata.var.highly_variable]
    #sc.pp.regress_out(combined_adata, ["total_counts", "pct_counts_mt"])
    sc.pp.scale(combined_adata, max_value=10)

    # Dimensionality reduction and clustering
    sc.tl.pca(combined_adata, svd_solver="arpack")
    sc.external.pp.bbknn(combined_adata, batch_key=args.batch_key, computation="annoy" ) # or "umap" 
    sc.pp.neighbors(combined_adata, n_neighbors=args.n_neighbors, n_pcs=args.n_pcs)
    sc.tl.umap(combined_adata) # init_pos="paga"
    sc.tl.leiden(
        combined_adata,
        resolution=args.resolution,
        #random_state=0,
        #flavor="igraph",
        n_iterations=-1,
        #directed=False
    )
    
    # Save the combined AnnData object 
    output_file = os.path.join(args.output_dir, f"{args.sample}.h5ad")
    os.makedirs(args.output_dir, exist_ok=True)
    combined_adata.write(output_file)
    #combined_adata = sc.read_h5ad(output_file)
    print(f"Combined dataset saved to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process, filter and merge scRNA-Seq data with Scanpy.")
    parser.add_argument("-i", "--input_dir", required=True, help="Input directory containing sample subdirectories. Required!")
    parser.add_argument("-m", "--meta_data", required=True, help="Meta data file. Required!")
    parser.add_argument("-o", "--output_dir", default=".", help="Output directory. Required!")
    parser.add_argument("-s", "--sample", required=True, help="Sample name for the combined .h5ad file. Required!")
    parser.add_argument("-c", "--copykat", help="Directory containing CopyKAT prediction files. Default: NULL.")
    parser.add_argument("--min_genes", type=int, default=200, help="Minimum number of genes per cell. Default: 200.")
    parser.add_argument("--min_cells", type=int, default=3, help="Minimum number of cells per gene. Default: 3.")
    parser.add_argument("--min_total_counts", type=int, default=200, help="Minimum total counts per cell. Default: 200.")
    parser.add_argument("--max_total_counts", type=int, default=30000, help="Maximum total counts per cell. Default: 30000.")
    parser.add_argument("--min_n_genes_by_counts", type=int, default=200, help="Minimum genes by counts. Default: 200.")
    parser.add_argument("--max_n_genes_by_counts", type=int, default=6000, help="Maximum genes by counts. Default: 6000.")
    parser.add_argument("--max_mt", type=float, default=10, help="Maximum mitochondrial percentage. Default: 10.")
    parser.add_argument("-g", "--n_top_genes", type=int, default=2000, help="Number of highly variable genes to select. Default: 2000.")
    parser.add_argument("-b", "--batch_key", default="orig.ident", help="Batch key for BBKNN. Default: orig.ident.")
    parser.add_argument("-n", "--n_neighbors", type=int, default=10, help="Number of neighbors for sc.pp.neighbors. Default: 10.")
    parser.add_argument("-p", "--n_pcs", type=int, default=30, help="Number of principal components for sc.pp.neighbors. Default: 30.")
    parser.add_argument("-r", "--resolution", type=float, default=0.5, help="Resolution for Leiden clustering. Default: 0.5.")
    args = parser.parse_args()
    main(args)
