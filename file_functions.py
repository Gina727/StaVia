import pandas as pd
import scanpy as sc
import os

def read_10x_mtx(file_path):
    import scipy.io
    import pandas as pd
    import scanpy as sc
    import os

    matrix_path = os.path.join(file_path, 'matrix.mtx.gz')
    print(f"Reading matrix from: {matrix_path}")
    matrix = scipy.io.mmread(matrix_path)
    print(f"Matrix shape: {matrix.shape}")
    
    # Read barcodes
    barcodes_path = os.path.join(file_path, 'barcodes.tsv.gz')
    print(f"Reading barcodes from: {barcodes_path}")
    barcodes = pd.read_csv(barcodes_path, header=None, sep='\t')[0].values
    print(f"Barcodes count: {len(barcodes)}")
    
    # Read features
    features_path = os.path.join(file_path, 'features.tsv.gz')
    print(f"Reading features from: {features_path}")
    features = pd.read_csv(features_path, header=None, sep='\t')
    print(f"Features shape: {features.shape}")
    print(f"Features columns: {features.columns.tolist()}")
    
    # Create AnnData with proper structure
    adata = sc.AnnData(X=matrix.T.tocsr())  # Transpose to cells x genes
    
    # Set observation names (barcodes)
    adata.obs_names = barcodes
    adata.obs['barcode'] = barcodes  # Add as column as well
    
    # Set variable names and metadata
    if features.shape[1] >= 2:
        # Use gene symbols as variable names (standard in Scanpy)
        adata.var_names = features[1].values
        adata.var['gene_ids'] = features[0].values
    else:
        adata.var_names = features[0].values
    
    # Set feature types if available
    if features.shape[1] >= 3:
        adata.var['feature_types'] = features[2].values
    else:
        adata.var['feature_types'] = 'Gene Expression'  # Default assumption
    
    # Add gene symbols as a separate column for clarity
    if features.shape[1] >= 2:
        adata.var['gene_symbol'] = features[1].values
    
    print(f"Final AnnData shape: {adata.shape}")
    print(f"AnnData .obs columns: {adata.obs.columns.tolist()}")
    print(f"AnnData .var columns: {adata.var.columns.tolist()}")
    print(adata.var_names[:10])
    
    return adata