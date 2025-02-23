import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData

def preprocess_data(file_path, file_type='h5ad', min_genes=200, min_cells=3, 
                    target_sum=1e4, n_top_genes=2000):
    """
    Comprehensive preprocessing pipeline for spatial transcriptomics data.

    Parameters:
    - file_path (str): Path to the input data file.
    - file_type (str): Type of the input file ('h5ad' or 'csv').
    - min_genes (int): Minimum number of genes expressed in a cell to retain it.
    - min_cells (int): Minimum number of cells expressing a gene to retain it.
    - target_sum (float): The total count for each cell after normalization.
    - n_top_genes (int): Number of top highly variable genes to select.

    Returns:
    - adata (AnnData): Preprocessed AnnData object.
    """
    # Load data
    if file_type == 'h5ad':
        adata = sc.read_h5ad(file_path)
    elif file_type == 'csv':
        df = pd.read_csv(file_path, index_col=0)
        adata = AnnData(df)
    else:
        raise ValueError("Unsupported file type. Use 'h5ad' or 'csv'.")

    # Quality control: Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)

    # Filter cells and genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Normalize data
    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # Identify highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes)

    # Filter to keep only highly variable genes
    adata = adata[:, adata.var['highly_variable']]

    return adata
