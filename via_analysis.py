import pyVIA.core as via #core module of pyVIA (a trajectory inference tool for single-cell data)
import pandas as pd #pandas library (for data manipulation and analysis)
import numpy as np #numpy library (for numerical computing)
import scanpy as sc #scanpy library (a toolkit for analyzing single-cell gene expression data)
import scanpy.external as sce #external extensions of scanpy (additional community-contributed tools)
import anndata as ad #anndata library (used for handling annotated data matrices, commonly used with scanpy)
import umap #umap library (for dimensionality reduction using UMAP)
import phate #phate library (a tool for visualizing high-dimensional data)
import matplotlib.pyplot as plt #pyplot module from matplotlib (for plotting)
import matplotlib as mpl #base matplotlib library (for additional customization)
mpl.use('Agg')
from matplotlib.pyplot import rc_context #rc_context from matplotlib.pyplot, which allows temporary runtime configuration changes.
import seaborn as sns #seaborn library (for statistical data visualization) 
import warnings #handling warnings in Python
import sys, os, glob #system-specific functions, operating system interactions, file path pattern matching
from sklearn.manifold import TSNE #TSNE (t-Distributed Stochastic Neighbor Embedding) function
warnings.filterwarnings('ignore') #Suppresses all warning messages
from importlib import reload #dynamically reload modules during development
from datetime import datetime
import random
from collections import defaultdict
import scvelo as scv
import csv

from flask import jsonify
from io import BytesIO
import base64
import json

def run_via_analysis(adata, params, file_data = None):
    try:
        knn = int(params.get('knn', 30))
        cluster_graph_pruning = float(params.get('cluster_graph_pruning', 0.15))
        edgepruning_clustering_resolution = float(params.get('edgepruning_clustering_resolution', 0.15))
        edgebundle_pruning = float(params.get('edgebundle_pruning', 0.15))
        true_label = params.get('true_label', None)
        root_user = params.get('root_user', None)
        time_series_labels = params.get('time_series_labels', None)
        adata_obs = params.get('adata_obs', None)
      
        data_categories = params.get('par_option', [])
        time_series = 'time-series' in data_categories
        use_velocity = 'rna-velocity' in data_categories
        do_spatial = 'spatial-temporal' in data_categories
        do_cytometry = 'cytometry' in data_categories

        if file_data is not None: 
            time_series_file = file_data.get('time-upload')
            velocity_matrix_file = file_data.get('velocity-matrix-upload')
            gene_matrix_file = file_data.get('gene-matrix-upload')
            root_upload_file = file_data.get('root-upload')
            true_label_file = file_data.get('csv-upload')
            spatial_coords_file = file_data.get('coords-upload')
            cytometry_file = file_data.get('cytometry-upload')

        results = {}

        true_label = None
        true_label = None
        if true_label_file:
            try:
                true_label = []
                reader = csv.reader(true_label_file)
                for row in reader:
                    if row:  
                        true_label.append(row[0])  
                if all(item.lstrip('-').isdigit() for item in true_label):
                    true_label = [int(item) for item in true_label]
            except Exception as e:
                print(f"Error processing true_label CSV: {e}")
                true_label = None
        else: 
            if true_label and isinstance(true_label, str):
                if true_label.lower() == 'none' and adata_obs.lower() == 'none':
                    true_label = None
                elif true_label.lower() != 'none' and adata_obs.lower() == 'none':
                    try:
                        true_label = [item.strip() for item in true_label.split(',')]
                        if all(item.lstrip('-').isdigit() for item in true_label):
                            true_label = [int(item) for item in true_label]
                    except Exception as e:
                        print(f"Error processing true_label: {e}")
                        true_label = None
                else:
                    if "annotation" in adata.obs:
                        true_label = adata.obs["annotation"]
                    if "PARC" in adata.obs:
                        true_label = adata.obs["PARC"]
                    else: 
                        true_label = adata.obs[adata_obs]
        if time_series_file:
            try:
                time_series_labels = []
                reader = csv.reader(time_series_file)
                for row in reader:
                    if row:  
                        time_series_labels.append(row[0])  
                if all(item.lstrip('-').isdigit() for item in time_series_labels):
                    time_series_labels = [int(item) for item in time_series_labels]
            except Exception as e:
                print(f"Error processing time series CSV: {e}")
                time_series_labels = None
        else: 
            if time_series_labels and isinstance(time_series_labels, str):
                if time_series_labels.lower() == 'none':
                    time_series_labels = None
                else:
                    try:
                        time_series_labels = [int(i.strip()) for i in time_series_labels.split(',') 
                                        if i.strip().isdigit()]
                    except Exception as e:
                        print(f"Error processing time_series_labels: {e}")
                        time_series_labels = None
        
        if root_upload_file:
            try:
                root_user = []
                reader = csv.reader(root_upload_file)
                for row in reader:
                    if row:  
                        root_user.append(row[0])  
                if all(item.lstrip('-').isdigit() for item in root_user):
                    root_user = [int(item) for item in root_user]
            except Exception as e:
                print(f"Error processing root file CSV: {e}")
                root_user = None
        else: 
            # Set root user if not provided
            if str(root_user).lower() != 'none': 
                root_user = [i.strip() for i in root_user.split(',')]
            elif not root_user or str(root_user).lower() == 'none':
                # Random the gene if None
                gene = random.choice(adata.var_names.tolist())
                print(f"DEBUG: Selected random gene: '{gene}'")
                print(f"DEBUG: Gene type: {type(gene)}")
                print(f"DEBUG: Gene in var_names? {gene in adata.var_names}")
                root_user = [adata[:, gene].X.argmax()]

        # INITIALIZE PARAMETERS
        n_pcs = 50
        ncomp = 50 
        random_seed = 0 
        small_pop = 2 
        too_big_factor = 0.3 
        memory = 50
        random_seed = 0

        if time_series: 
            time_series_labels = time_series_labels
        else:
            time_series_labels = None

        velocity_matrix = None
        gene_matrix = adata.X.todense() if hasattr(adata.X, 'todense') else adata.X
        velo_weight = 0

        if use_velocity:
            velo_weight = 0.5
            if velocity_matrix_file is not None:
                try: 
                    velocity_df = velocity_matrix_file 
                    
                    # Clean velocity cell names
                    # velocity_df.index = [
                    #     name.replace('Het_CR_outs:', '').replace('x', '-1') 
                    #     for name in velocity_df.index
                    # ]
                    
                    # Find common cells
                    common_cells = velocity_df.index.intersection(adata.obs_names)
                    print(f"Common cells: {len(common_cells)}")
                    
                    if len(common_cells) == 0:
                        print("✗ No common cells - disabling velocity")
                        raise ValueError("No common cells between velocity and adata")
                    
                    # Filter to common cells
                    velocity_df = velocity_df.loc[common_cells]
                    adata = adata[common_cells]
                    
                    # Match genes
                    print(f"Velocity shape: {velocity_df.shape}")
                    print(f"Adata shape: {adata.shape}")
                    
                    # Get the CURRENT adata gene names (after preprocessing filtering)
                    adata_genes = set(adata.var_names)
                    
                    # Check if velocity has column names
                    if hasattr(velocity_df, 'columns') and len(velocity_df.columns) > 0:
                        velocity_genes = set(velocity_df.columns)
                        
                        print(f"Adata genes: {len(adata_genes)}")
                        print(f"Velocity genes: {len(velocity_genes)}")
                        
                        # Find intersection - only genes that exist in BOTH
                        common_genes = list(adata_genes.intersection(velocity_genes))
                        print(f"Common genes: {len(common_genes)}")
                        
                        # Find what's different
                        missing_from_adata = velocity_genes - adata_genes
                        missing_from_velocity = adata_genes - velocity_genes
                        
                        if missing_from_adata:
                            print(f"Genes in velocity but NOT in filtered adata: {len(missing_from_adata)}")
                            print(f"  Sample: {list(missing_from_adata)[:5]}")
                        
                        if missing_from_velocity:
                            print(f"Genes in filtered adata but NOT in velocity: {len(missing_from_velocity)}")
                            print(f"  Sample: {list(missing_from_velocity)[:5]}")
                        
                        if len(common_genes) == 0:
                            print("No common genes - disabling velocity")
                            raise ValueError("No common genes between velocity and adata")
                        
                        # Filter BOTH to common genes
                        velocity_df = velocity_df[common_genes]
                        adata = adata[:, common_genes]  # This is important!
                        
                        print(f"After gene filtering: velocity {velocity_df.shape}, adata {adata.shape}")
                        
                    else:
                        print("ERROR: Velocity DataFrame has no column names!")
                        print("Cannot match genes without column names")
                        raise ValueError("Velocity matrix missing gene names")
                    
                    # Convert to numpy arrays for VIA
                    velocity_matrix = velocity_df.values
                    gene_matrix = adata.X.todense() if hasattr(adata.X, 'todense') else adata.X
                    
                    # Final validation
                    print(f"\n=== FINAL CHECK ===")
                    print(f"Velocity matrix shape: {velocity_matrix.shape}")
                    print(f"Gene matrix shape: {gene_matrix.shape}")
                    print(f"adata.X shape: {adata.X.shape}")
                    
                    if velocity_matrix.shape != gene_matrix.shape:
                        print("Final shape mismatch!")
                        print("  Forcing alignment by taking first n genes...")
                        
                        min_genes = min(velocity_matrix.shape[1], gene_matrix.shape[1])
                        velocity_matrix = velocity_matrix[:, :min_genes]
                        gene_matrix = gene_matrix[:, :min_genes]
                        adata = adata[:, :min_genes]
                        
                        print(f"  After emergency fix: velocity {velocity_matrix.shape}, gene {gene_matrix.shape}")
                        velo_weight = 0.5
                    else:
                        print("All dimensions match!")
                        velo_weight = 0.5
                        
                except Exception as e: 
                    print(f"Error loading velocity data: {e}")
                    import traceback
                    traceback.print_exc()
                    velocity_matrix = None
                    gene_matrix = None
                    velo_weight = 0
            if velocity_matrix_file:
                try: 
                    velocity_matrix = pd.read_csv(velocity_matrix_file)
                    velocity_matrix = velocity_matrix_file.values

                    if not gene_matrix_file:
                        gene_matrix=adata.X.todense()
                    else: 
                        gene_matrix = pd.read_csv(gene_matrix_file)
                        adata.X = gene_matrix

                    velo_weight=0.5
                except Exception as e: 
                    print(f"Error loading velocity data: {e}")
            else: 
                print('Please upload velocity matrix file')
        else: 
            gene_matrix = None
            velocity_matrix = None
            velo_weight=0

        if do_spatial:
            
            # Add text input?
            if spatial_coords_file:
                coords = pd.read_csv(spatial_coords_file) 
            else:
                coords=adata.obsm['X_pca'] 

            spatial_knn_input = 6 
            spatial_weight = 0.3
            spatial_knn_trajectory =6
            X_spatial_exp = via.spatial_input(X_genes = adata.X, spatial_coords=coords, knn_spatial=spatial_knn_input, spatial_weight=spatial_weight)

            adata.obsm['X_spatial_adjusted'] = X_spatial_exp
            adata.obsm["spatial_pca"] = sc.tl.pca(adata.obsm['X_spatial_adjusted'],n_comps=n_pcs)

            print(f'end X_spatial input')
        else:
            spatial_knn_input = 0
            spatial_knn_trajectory = 0
            coords = None  
            spatial_weight = 0
        
        if do_cytometry:
            print(f"Type of cytometry_file: {type(cytometry_file)}")
            try:
                if 'Unnamed: 0' in cytometry_file.columns:
                    cytometry = cytometry_file.drop('Unnamed: 0', axis=1)
                else:
                    cytometry = cytometry_file
                
                cytometry = cytometry.dropna()
                print(f'Loaded cytometry file with shape: {cytometry.shape}')
                
                # Debug: Check what we're working with
                print(f"Columns type: {type(cytometry.columns)}")
                print(f"First few columns: {cytometry.columns[:5].tolist()}")
                
                # Create AnnData object with explicit parameters
                ad = sc.AnnData(X=cytometry.values, dtype=cytometry.values.dtype)
                print(f'Created AnnData with shape: {ad.shape}')
                
                # Set variable names - ensure it's a list
                if hasattr(cytometry.columns, 'tolist'):
                    var_names = cytometry.columns.tolist()
                else:
                    var_names = list(cytometry.columns)
                
                ad.var_names = var_names
                print(f'Set var_names, type: {type(ad.var_names)}')
                
                # Set observation names
                obs_names = [f"cell_{i}" for i in range(cytometry.shape[0])]
                ad.obs_names = obs_names
                
                # Verify the object
                print(f'AnnData X shape: {ad.X.shape}')
                print(f'AnnData var_names length: {len(ad.var_names)}')
                print(f'AnnData obs_names length: {len(ad.obs_names)}')
                
                # Continue with processing
                sc.pp.scale(ad)
                sc.tl.pca(ad, svd_solver='arpack')
                print(f'PCA is finished')
                
                X_in = ad.X
                cytometry_X = pd.DataFrame(X_in, columns=var_names)
                
                X_in = cytometry_X.values
                
                # Create new AnnData for processed data
                ad_processed = sc.AnnData(X=cytometry_X.values, dtype=cytometry_X.values.dtype)
                ad_processed.var_names = var_names
                ad_processed.obs_names = obs_names[:cytometry_X.shape[0]]  # Adjust if needed
                
                sc.tl.pca(ad_processed, svd_solver='arpack')
                print(f'End cytometry input')

            except Exception as e:
                print(f"Error processing cytometry CSV: {e}")
                import traceback
                traceback.print_exc()  # This will show in server logs
                knn = 20
                jac_std_global = 0.5
                random_seed = 1
                root_user = None
        


        print('RUN VIA')
        v0 = via.VIA(adata.obsm['X_pca'][:,:ncomp], true_label = true_label, memory = memory,
                    edgepruning_clustering_resolution=edgepruning_clustering_resolution, 
                    edgepruning_clustering_resolution_local=1, knn=knn,
                    too_big_factor=too_big_factor, root_user=root_user,
                    cluster_graph_pruning=cluster_graph_pruning, 
                    edgebundle_pruning_twice = False, time_series=time_series, time_series_labels=time_series_labels,
                    edgebundle_pruning = edgebundle_pruning,
                    small_pop = small_pop, velo_weight=velo_weight, velocity_matrix=velocity_matrix, gene_matrix=gene_matrix,
                    piegraph_arrow_head_width=0.07,
                    random_seed=random_seed,
                    resolution_parameter=1,
                    is_coarse=True, 
                    x_lazy=0.99, alpha_teleport=0.99, 
                    viagraph_decay = 1.0, 
                    preserve_disconnected=False,
                    do_spatial_knn=do_spatial, do_spatial_layout= do_spatial, spatial_coords = coords, spatial_knn=spatial_knn_trajectory)
        
        v0.run_VIA()
        if 'X_umap' in adata.obsm:
            v0.embedding = adata.obsm['X_umap'][:,:2]
        elif 'X_pca' in adata.obsm:
            v0.embedding = adata.obsm['X_pca'][:,:2]
        results['via_obj'] = v0

        results['adata'] = adata

        return results
    
    except Exception as e:
        return {'error': str(e)}