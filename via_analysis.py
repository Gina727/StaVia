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

        if file_data is not None: 
            time_series_file = file_data.get('time-upload')
            velocity_matrix_file = file_data.get('velocity-matrix-upload')
            gene_matrix_file = file_data.get('gene-matrix-upload')
            root_upload_file = file_data.get('root-upload')
            true_label_file = file_data.get('csv-upload')
            spatial_coords_file = file_data.get('coords-upload')

        results = {}

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
                if true_label.lower() == 'none':
                    true_label = None
                else:
                    try:
                        true_label = [item.strip() for item in true_label.split(',')]
                        if all(item.lstrip('-').isdigit() for item in true_label):
                            true_label = [int(item) for item in true_label]
                    except Exception as e:
                        print(f"Error processing true_label: {e}")
                        true_label = None
        
        if adata_obs and isinstance(adata_obs, str) and true_label is None:
            if adata_obs.lower() == 'none':
                true_label = None
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
            time_series_labels=time_series_labels
        else:
            time_series_labels=None

        if use_velocity:
            velocity_matrix = pd.read_csv(velocity_matrix_file)
            velocity_matrix = velocity_matrix_file.values
            gene_matrix = pd.read_csv(gene_matrix_file)
            adata.X = gene_matrix.values
            velo_weight=0.5
        else: 
            gene_matrix =None
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

        return results
    
    except Exception as e:

        return {'error': str(e)}
