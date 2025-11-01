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
import matplotlib.cm as cm

from flask import jsonify
from io import BytesIO
import base64
import json

# Plotting Lieneage Probability, Gene Trend Heatmaps, and Gene Expression (needs v0 object!!!)
def more_plot(lineages, genes, v0, adata):
    plots = {}
    try:
        if not lineages:
            lineages = []

        if not genes: 
            genes = random.sample(adata.var_names.tolist(), 3)

        fig, ax = plt.subplots()
        via.plot_sc_lineage_probability(via_object=v0, embedding=adata.obsm['X_pca'], idx=None, cmap_name='plasma', dpi=150, 
                                        marker_lineages=lineages, fontsize=8, alpha_factor=0.9, 
                                        majority_cluster_population_dict=None, cmap_sankey='rainbow', do_sankey=False)
        plt.gcf().set_size_inches(16,8)
        lineage_img = BytesIO()
        plt.savefig(lineage_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        plots['lineage'] = "data:image/png;base64," + base64.b64encode(lineage_img.getvalue()).decode('utf-8')


        df_genes_exp = pd.DataFrame(adata[:, genes].X.todense(), columns=genes)

        fig, ax = plt.subplots()
        via.plot_gene_trend_heatmaps(via_object=v0, df_gene_exp=df_genes_exp, marker_lineages=[], 
                                    fontsize=8, cmap='viridis', normalize=True, ytick_labelrotation=0, fig_width=7)
        heat_img = BytesIO()
        plt.savefig(heat_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        plots['heat'] = "data:image/png;base64," + base64.b64encode(heat_img.getvalue()).decode('utf-8')
        
        fig, ax = plt.subplots()
        via.get_gene_expression(via_object=v0, gene_exp=df_genes_exp, cmap='jet', dpi=150, marker_genes=[], 
                                linewidth=2.0, n_splines=10, spline_order=4, fontsize_=8, marker_lineages=[], 
                                optional_title_text='', cmap_dict=None, conf_int=0.95, driver_genes=False, driver_lineage=None)
        plt.gcf().set_size_inches(16,8)
        exp_img = BytesIO()
        plt.savefig(exp_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        plots['exp'] = "data:image/png;base64," + base64.b64encode(exp_img.getvalue()).decode('utf-8')

        return plots

    except Exception as e:
        return {'error': str(e)}

# General VIA Plots
def via_plot(params, v0, file_data):
    try: 
        # Get all the parameters and files from Flask 
        var_names = params.get('var_names', None)
        dpi = int(params.get('dpi', 180))
        time_series_labels = params.get('time_series_labels', None)
        knn = int(params.get('knn', 30))
        true_label = params.get('true_label', None)
        
        data_categories = params.get('par_option', [])
        use_velocity = 'rna-velocity' in data_categories
        do_spatial = 'spatial-temporal' in data_categories
        do_cytomtetry = 'cytometry' in data_categories
        # Create an empty plot object 
        plots = {}

        if file_data is not None: 
            time_series_file = file_data.get('time-upload')
            velocity_matrix_file = file_data.get('velocity-matrix-upload')
            true_label_file = file_data.get('csv-upload')

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
                print(f"Error processing true_label CSV: {e}")
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

        # ATLAS PLOT
        via.via_atlas_emb(via_object=v0, n_components=2, alpha=1.0, negative_sample_rate=5, 
                  gamma=1.0, min_dist=0.1, random_state=0, n_epochs=100)
                # Now plot the embedding
        plt.figure(figsize=(10, 8))
        plt.scatter(v0.embedding[:, 0], v0.embedding[:, 1], 
                    c=v0.labels, cmap='viridis', s=5, alpha=0.6)
        plt.title("VIA Atlas Embedding")
        plt.xlabel("Component 1")
        plt.ylabel("Component 2")
        plt.colorbar(label="Cluster")
        atlas_img = BytesIO()
        plt.savefig(atlas_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        plots['atlas'] = "data:image/png;base64," + base64.b64encode(atlas_img.getvalue()).decode('utf-8')

        # VIA GRAPH
        print(f'{datetime.now()}\t plotting piechart_graph')
        fig,ax1, ax2 = via.plot_piechart_viagraph(via_object=v0, reference_labels=time_series_labels, show_legend=False, ax_text=False,
                                headwidth_arrow=0.8, highlight_terminal_clusters=False, cmap_piechart='plasma', cmap='viridis',
                                pie_size_scale=0.6,size_node_notpiechart=1)
        fig.set_size_inches(40, 15)
        via_img = BytesIO()
        plt.savefig(via_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        plots['via'] = "data:image/png;base64," + base64.b64encode(via_img.getvalue()).decode('utf-8')

        # SC LINEAGES
        # ecto_lineages=[1,71,89]#,136,10]#,83,115] #CNS

        # fig, ax = via.plot_sc_lineage_probability(via_object=v0, alpha_factor=0.6, scatter_size=4,
        #                                 marker_lineages=ecto_lineages)
        # fig.set_size_inches(40 , 15)     
        # line_img = BytesIO()
        # plt.savefig(line_img, format='png', bbox_inches='tight', dpi=120)
        # plt.close()
        # plots['line_plot'] = "data:image/png;base64," + base64.b64encode(line_img.getvalue()).decode('utf-8')
        
        # Plot this if the user ticks 'RNA Velocity'
        if use_velocity: 
            adata = sc.read(velocity_matrix_file)
            if var_names == 'None':
                var_names = np.random.choice(adata.var_names)
                print(f"Randomly selected gene: {var_names}")
            else:
                if var_names not in adata.var_names:
                    return {'error': f"Gene '{var_names}' not found in dataset"}, 400
            print(adata)
            print(adata.obs.columns.tolist())

            # Embedding on single-cell level 
            scv.pl.velocity_embedding(adata, basis='umap', color='clusters', dpi=dpi, show=False)
            em_img = BytesIO()
            plt.savefig(em_img, format='png', bbox_inches='tight', dpi=dpi)
            plt.close()
            plots['em'] = "data:image/png;base64," + base64.b64encode(em_img.getvalue()).decode('utf-8')

            scv.pl.velocity_embedding_grid(adata, basis='umap', color='clusters', scale=0.25, show=False)
            grid_img = BytesIO()
            plt.savefig(grid_img, format='png', bbox_inches='tight', dpi=dpi)
            plt.close()
            plots['grid'] = "data:image/png;base64," + base64.b64encode(grid_img.getvalue()).decode('utf-8')

            via.via_streamplot(via_object=v0, scatter_size=150, scatter_alpha=0.2,
                    marker_edgewidth=0.05, density_stream=1.0, density_grid=0.5, smooth_transition=2,
                    smooth_grid=0.5, color_scheme='annotation', add_outline_clusters=False,
                    cluster_outline_edgewidth=0.001)
            stream_img = BytesIO()
            plt.savefig(stream_img, format='png', bbox_inches='tight', dpi=dpi)
            plt.close()
            plots['stream'] = "data:image/png;base64," + base64.b64encode(stream_img.getvalue()).decode('utf-8')

            # Connection between the clusters
            scv.pl.velocity_graph(adata, color='clusters', threshold=.3, show=False)
            cl_img = BytesIO()
            plt.savefig(cl_img, format='png', bbox_inches='tight', dpi=dpi)
            plt.close()
            plots['cl'] = "data:image/png;base64," + base64.b64encode(cl_img.getvalue()).decode('utf-8')

            # RNA velocity and genetic expression level
            scv.pl.velocity(adata, var_names=[var_names], color='clusters', show=False)
            vel_img = BytesIO()
            plt.savefig(vel_img, format='png', bbox_inches='tight', dpi=dpi)
            plt.close()
            plots['vel'] = "data:image/png;base64," + base64.b64encode(vel_img.getvalue()).decode('utf-8')

        # Plot this if the user ticks "Spatio-Temporal"
        if do_spatial:
            # Cluster on tissue slice
            # fig,ax1, ax2 = via.plot_clusters_spatial(spatial_coords=coords, clusters = [7,11,27,30], color='black', 
            #                       via_labels=v0.labels, title_sup='OD types', alpha = 0.8, s=5)
            # fig.set_size_inches(40, 15)
            # spa_img = BytesIO()
            # plt.savefig(spa_img, format='png', bbox_inches='tight', dpi=120)
            # plt.close()
            # plots['spa'] = "data:image/png;base64," + base64.b64encode(spa_img.getvalue()).decode('utf-8')

            if true_label:
                color_dict = {}
                set_labels = list(set(true_label))
                set_labels.sort(reverse=False)#True)
                palette = cm.get_cmap('rainbow', len(set_labels))
                cmap_ = palette(range(len(set_labels)))
                for index, value in enumerate(set_labels):
                    color_dict[value] = cmap_[index]

                mds_title = 'mds n_pcs ='+str(50) +'_knn='+str(knn)
                emb_mds = via.via_mds(via_object=v0, n_milestones=5000)
                f1, ax1 = via.plot_scatter(embedding=emb_mds, labels=true_label, title=mds_title, alpha=0.5,s=12,show_legend=True,color_dict=color_dict, text_labels=False)
                ax1.get_legend().set_bbox_to_anchor((1.05, 1.05))
                mds_img = BytesIO()
                plt.savefig(mds_img, format='png', bbox_inches='tight', dpi=120)
                plt.close()
                plots['mds'] = "data:image/png;base64," + base64.b64encode(mds_img.getvalue()).decode('utf-8')

        # Return all the plots to the Flask app to be displayed on the frontend 
        return plots
    
    except Exception as e:
        return {'error': str(e)}, 500