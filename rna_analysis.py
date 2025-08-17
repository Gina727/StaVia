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
import scvelo as scv
from importlib import reload #dynamically reload modules during development

from flask import jsonify
from io import BytesIO
import base64
import json


def run_rna_velocity(filepath, request):
    try:
        # Input validation and extraction
        if 'file' not in request.files:
            return {'error': 'No file uploaded'}, 400
            
        file = request.files['file']
        if not file.filename.endswith('.h5ad'):
            return {'error': 'Invalid file type. Please upload .h5ad file'}, 400

        choice = json.loads(request.form.get('em', '["umap"]'))
        knn = int(request.form.get('knn', 30))
        cluster_graph_pruning = float(request.form.get('cluster_graph_pruning', 0.15))
        edgepruning_clustering_resolution = float(request.form.get('edgepruning_clustering_resolution', 0.15))
        edgebundle_pruning = float(request.form.get('edgebundle_pruning', 0.15))
        true_label = request.form.get('true_label', 'None')
        dpi = int(request.form.get('dpi', 120))
        var_names = request.form.get('var_names', 'None')
        root_user = request.form.get('root_user', 'None')

        csv_file = request.files.get('true_label_csv')
        if csv_file and csv_file.filename.endswith('.csv'):
            try:
                # Read CSV into list
                df = pd.read_csv(csv_file, header=None)
                true_label = df[0].tolist()
            except Exception as e:
                print(f"Error processing CSV file: {e}")
                # Fall back to text input if CSV processing fails
                if true_label == 'None':
                    true_label = None
        elif true_label == 'None':
            true_label = None
        else:
            try:
                true_label = [item.strip() for item in true_label.split(',')]
                if all(item.lstrip('-').isdigit() for item in true_label):
                    true_label = [int(item) for item in true_label]
            except Exception as e:
                print(f"Error processing true_label: {e}")
                true_label = None

        results = {
                'success': True,
                'filename': file.filename,
                'adata_info': None
            }
        
        adata = sc.read(filepath)
        if var_names == 'None':
            var_names = np.random.choice(adata.var_names)
            print(f"Randomly selected gene: {var_names}")
        else:
            if var_names not in adata.var_names:
                return {'error': f"Gene '{var_names}' not found in dataset"}, 400
        print(adata)
        results['adata_info'] = str(adata)
        print(adata.obs.columns.tolist())

        if 'pca' in choice: 
            sc.pl.pca_variance_ratio(adata, log=False, n_pcs=50, show=False)
            pca_img = BytesIO()
            plt.savefig(pca_img, format='png', bbox_inches='tight', dpi=120)
            plt.close()
            results['pca_plot'] = "data:image/png;base64," + base64.b64encode(pca_img.getvalue()).decode('utf-8')

        if 'umap' in choice:
            adata.obsm['X_umap'] = umap.UMAP(n_neighbors=20, min_dist=0.2, spread=5, init='pca').fit_transform(adata.obsm['X_pca'])
            sc.pl.embedding(adata, basis='X_umap', color = ['cluster'], size = 200)
            umap_img = BytesIO()
            plt.savefig(umap_img, format='png', bbox_inches='tight', dpi=120)
            plt.close()
            results['umap_plot'] = "data:image/png;base64," + base64.b64encode(umap_img.getvalue()).decode('utf-8')

        # Plot percentage of spliced and non-spliced mRNA
        scv.pl.proportions(adata, groupby='clusters', show=False)
        per_img = BytesIO()
        plt.savefig(per_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        results['per_plot'] = "data:image/png;base64," + base64.b64encode(per_img.getvalue()).decode('utf-8')

        # Normalize the data + first and second order moments
        scv.pp.filter_and_normalize(adata)
        scv.pp.moments(adata)
        
        # CORE ANALYSIS (computing velocities)
        # scv.tl.recover_dynamics(adata) 
        # scv.tl.velocity(adata, mode='dynamical')
        scv.tl.velocity_graph(adata, n_pcs=30, n_neighbors=20)

        n_pcs = 30
        knn = 20
        scv.pp.moments(adata, n_pcs=None, n_neighbors=None)

        v0_toobig = .3 #(clusters that form more than 30% of entire population will be re-clustered)
        dataset = '' 
        random_seed = 42
        velo_weight=0.5 #weight given to velocity matrix from scvelo. 1-velo_weight is assigned to gene-gene kernel
        embedding = adata.obsm['X_umap'][:, 0:2]
        true_label = adata.obs['clusters'].tolist()
        pca_loadings = adata.varm['PCs'] # this is optional and offers slight adjustment of the locations of cells based on velocity

        if 'velocity_matrix' in request.files and request.files['velocity_matrix'].filename:
            velocity_matrix = pd.read_csv(request.files['velocity_matrix'])
            velocity_matrix = velocity_matrix.values
        else:
            velocity_matrix=adata.layers['velocity']

        # Handle gene matrix upload if provided
        if 'gene_matrix' in request.files and request.files['gene_matrix'].filename:
            gene_matrix = pd.read_csv(request.files['gene_matrix'])
            adata.X = gene_matrix.values
        else:
           gene_matrix=adata.X.todense()

        if root_user != 'None': 
            if all(part.strip().isdigit() for part in root_user.split(',')):  # Comma-separated numbers
                root_user = [int(i.strip()) for i in root_user.split(',')]
        else:
            root_user = None

        # Impute genes we want to use in gene-trends later 
        df_ = pd.DataFrame(adata.X.todense())
        df_.columns = [i for i in adata.var_names]

        v0 = via.VIA(adata.obsm['X_pca'][:, 0:n_pcs], true_label=true_label, knn=knn,
                    too_big_factor=v0_toobig, root_user=root_user, dataset=dataset, random_seed=random_seed,
                is_coarse=True, preserve_disconnected=True, pseudotime_threshold_TS=50, edgepruning_clustering_resolution=edgepruning_clustering_resolution,
                    piegraph_arrow_head_width=0.15, cluster_graph_pruning=cluster_graph_pruning,
                    piegraph_edgeweight_scalingfactor=2.5, velocity_matrix=velocity_matrix,
                            gene_matrix=gene_matrix, velo_weight=velo_weight,  edgebundle_pruning_twice=False, edgebundle_pruning=edgebundle_pruning, pca_loadings = adata.varm['PCs']) # pca_loadings is optional #edge_bundling_twice = True would remove more of the edges
        v0.run_VIA()

        # PROBLEM HERE 
        if 'atlas' in choice:    
            # 1. First ensure projections exist
            if not hasattr(v0, 'projections') or v0.projections is None:
                v0.projections = adata.obsm['X_pca'][:, :2]  # Fallback to PCA
            
            # 2. Compute embedding manually
            print("Computing UMAP embedding manually")
            reducer = umap.UMAP(
                n_components=2,
                n_neighbors=15,
                min_dist=0.1,
                random_state=42
            )
            v0.embedding = reducer.fit_transform(v0.projections)

            try:
                via.via_atlas_emb(via_object=v0, n_components=2)  # Minimal params
            except Exception as e:
                print(f"via_atlas_emb failed (but we have manual embedding): {str(e)}")
            
            # 4. Plot the embedding we're certain exists
            plt.figure(figsize=(10, 8))
            plt.scatter(v0.embedding[:, 0], v0.embedding[:, 1], 
                    c=v0.labels, cmap='viridis', s=5, alpha=0.6)
            plt.title("VIA Atlas Embedding (Manual UMAP)")
            atlas_img = BytesIO()
            plt.savefig(atlas_img, format='png', bbox_inches='tight', dpi=120)
            plt.close()
            results['atlas_plot'] = "data:image/png;base64," + base64.b64encode(atlas_img.getvalue()).decode('utf-8')
                    
        via.plot_piechart_viagraph(via_object = v0)
        via_img = BytesIO()
        plt.savefig(via_img, format='png', bbox_inches='tight', dpi=dpi)
        plt.close()
        results['via_plot'] = "data:image/png;base64," + base64.b64encode(via_img.getvalue()).decode('utf-8')

        # Embedding on single-cell level 
        scv.pl.velocity_embedding(adata, basis='umap', color='clusters', dpi=dpi, show=False)
        em_img = BytesIO()
        plt.savefig(em_img, format='png', bbox_inches='tight', dpi=dpi)
        plt.close()
        results['em_plot'] = "data:image/png;base64," + base64.b64encode(em_img.getvalue()).decode('utf-8')

        scv.pl.velocity_embedding_grid(adata, basis='umap', color='clusters', scale=0.25, show=False)
        grid_img = BytesIO()
        plt.savefig(grid_img, format='png', bbox_inches='tight', dpi=dpi)
        plt.close()
        results['grid_plot'] = "data:image/png;base64," + base64.b64encode(grid_img.getvalue()).decode('utf-8')

        via.via_streamplot(via_object=v0, embedding=embedding, scatter_size=150, scatter_alpha=0.2,
                   marker_edgewidth=0.05, density_stream=1.0, density_grid=0.5, smooth_transition=2,
                   smooth_grid=0.5, color_scheme='annotation', add_outline_clusters=False,
                   cluster_outline_edgewidth=0.001)
        stream_img = BytesIO()
        plt.savefig(stream_img, format='png', bbox_inches='tight', dpi=dpi)
        plt.close()
        results['stream_plot'] = "data:image/png;base64," + base64.b64encode(stream_img.getvalue()).decode('utf-8')

        # Connection between the clusters
        scv.pl.velocity_graph(adata, color='clusters', threshold=.3, show=False)
        cl_img = BytesIO()
        plt.savefig(cl_img, format='png', bbox_inches='tight', dpi=dpi)
        plt.close()
        results['cl_plot'] = "data:image/png;base64," + base64.b64encode(cl_img.getvalue()).decode('utf-8')

        # RNA velocity and genetic expression level
        scv.pl.velocity(adata, var_names=[var_names], color='clusters', show=False)
        vel_img = BytesIO()
        plt.savefig(vel_img, format='png', bbox_inches='tight', dpi=dpi)
        plt.close()
        results['vel_plot'] = "data:image/png;base64," + base64.b64encode(vel_img.getvalue()).decode('utf-8')
        

        return results
    
    except Exception as e:
        return {'error': str(e)}, 500