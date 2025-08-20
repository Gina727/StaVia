import pandas as pd
import scanpy as sc
from via_analysis import run_via_analysis
from rna_analysis import run_rna_velocity
from plotting import via_plot, more_plot
import pyVIA.core as via
from flask import Flask, render_template, request, jsonify, send_file, session, url_for, redirect
from flask_sqlalchemy import SQLAlchemy
import json
import os
from io import BytesIO
import base64
import zipfile
import atexit
import shutil
from datetime import datetime
import tempfile
import matplotlib.pyplot as plt
import phate  
import umap
import pprint
import numpy as np 
from functools import wraps
import gc
import scipy
from datetime import datetime
import uuid
import gzip
import shutil

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///test.db'
db = SQLAlchemy(app)
app.secret_key = os.urandom(24)

class Todo(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    content = db.Column(db.String(200), nullable=False)
    parameters = db.Column(db.JSON)
    date_created = db.Column(db.DateTime, default=datetime.utcnow)

    def __repr__(self):
        return f'<Task {self.id}>'
    
UPLOAD_FOLDER = tempfile.mkdtemp()
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
PREPROCESS_CACHE = {}
VIA_CACHE = {}

def cleanup():
    try:
        shutil.rmtree(UPLOAD_FOLDER)
    except Exception as e:
        print(f"Error cleaning up upload folder: {e}")

atexit.register(cleanup)

def safe_json(obj):
    """Recursively replace NaN and inf with None for JSON serialization"""
    if isinstance(obj, dict):
        return {k: safe_json(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [safe_json(v) for v in obj]
    elif isinstance(obj, float):
        if np.isnan(obj) or np.isinf(obj):
            return None
        return obj
    return obj

@app.template_filter('format_params')
def format_params(params):
    if not params:
        return "No parameters"
    result = []
    for key, value in params.items():
        if isinstance(value, list):
            value = ", ".join(value)
        result.append(f"{key}: {value}")
    return "\n".join(result)  # or use <br> for HTML line breaks
app.jinja_env.filters['format_params'] = format_params

@app.route('/')
def home():
    return render_template('home.html')

@app.route('/upload_page', methods=['GET'])
def upload_page():
    tasks = Todo.query.order_by(Todo.date_created).all()
    return render_template('index.html', tasks=tasks)   

@app.route('/add_info', methods=['POST'])
def add_info():
    try:
        data = request.get_json()
        print("Received data:", data)  
        
        if not data or 'content' not in data:
            return jsonify({'error': 'Invalid data'}), 400
            
        new_task = Todo(
            content=data['content'],
            parameters=data.get('parameters', {}),  
            date_created=datetime.utcnow()
        )
        
        db.session.add(new_task)
        db.session.commit()
        return jsonify({'success': True, 'id': new_task.id})
    
    except Exception as e:
        print("Error:", str(e))
        return jsonify({'error': str(e)}), 500
    
@app.route('/delete/<int:id>')
def delete(id):
    task_to_delete = Todo.query.get_or_404(id)

    try:
        db.session.delete(task_to_delete)
        db.session.commit()
        return redirect('/upload_page')
    except:
        return 'There was a problem deleting that task'


@app.route('/about_page')
def about_page():
    return render_template('about.html')

@app.route('/upload', methods=['POST'])
def upload():
    try:
        PREPROCESS_CACHE.clear()
        VIA_CACHE.clear()
        session.clear()
        unique_id = str(uuid.uuid4())[:8]

        form_data = request.form
        files = request.files
        
        # Determine upload type
        if 'file' in files:  # H5AD file
            file = files['file']
            if not file.filename.lower().endswith('.h5ad'):
                return jsonify({'error': 'Only .h5ad files supported'}), 400
            
            # Save h5ad file
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], f'uploaded_data_{unique_id}.h5ad')
            file.save(file_path)
            session['file_type'] = 'h5ad'
            
        else:  # 10X files
            
            has_matrix = any(f for f in files if 'matrix' in f.lower())
            has_barcodes = any(f for f in files if 'barcodes' in f.lower())
            has_features = any(f for f in files if 'features' in f.lower() or 'genes' in f.lower())
            
            if has_matrix and has_barcodes and has_features:
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], f'10x_data_{unique_id}')
                os.makedirs(file_path, exist_ok=True)

                # Find the actual files
                mtx_file = next(f for f in files.values() if 'matrix' in f.filename.lower())
                barcodes_file = next(f for f in files.values() if 'barcodes' in f.filename.lower())
                features_file = next(f for f in files.values() if 'features' in f.filename.lower() or 'gene' in f.filename.lower())

            # Save original files first
            mtx_path = os.path.join(file_path, 'matrix.mtx')
            barcodes_path = os.path.join(file_path, 'barcodes.tsv')
            features_path = os.path.join(file_path, 'features.tsv')

            mtx_file.save(mtx_path)
            barcodes_file.save(barcodes_path)
            features_file.save(features_path)

            for path in [mtx_path, barcodes_path, features_path]:
                    if path.endswith('.gz'):
                        continue
                    if os.path.exists(path) and not os.path.exists(path + '.gz'):
                        with open(path, 'rb') as f_in:
                            with gzip.open(path + '.gz', 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        os.remove(path)

            session['file_type'] = '10x'
            print("Session file_type set to:", session['file_type'])
        
        session['file_path'] = file_path
        return jsonify({'success': True, 'message': 'Files uploaded successfully'})
    
    except Exception as e:
        return jsonify({'error': str(e)}), 500

def cache_preprocessed_data(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        print("\n=== CACHE CHECK ===")
        print("Session file_path:", session.get('file_path'))
        print("Session file_type:", session.get('file_type'))

        file_path = session.get('file_path')
        file_type = session.get('file_type')
        
        if not file_path:
            print("ERROR: No file_path in session")
            return jsonify({'error': 'No file uploaded'}), 400

        if file_path in PREPROCESS_CACHE:
            print("Using cached data")
            return func(PREPROCESS_CACHE[file_path], *args, **kwargs)

        try:
            if file_type == '10x':
                print('Loading 10X data from:', file_path)
                mtx_gz_path = os.path.join(file_path, 'matrix.mtx.gz')
                print("Checking for gzipped matrix at:", mtx_gz_path)
                print("Exists:", os.path.exists(mtx_gz_path))

                print("Directory listing:", os.listdir(file_path))

                for f in ["barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"]:
                    fpath = os.path.join(file_path, f)
                    print(f"Checking {fpath} ... exists:", os.path.exists(fpath))
                    try:
                        with gzip.open(fpath, "rt") as g:
                            for i, line in enumerate(g):
                                print(f"First line of {f}: {line.strip()}")
                                break
                    except Exception as e:
                        print(f"Failed to open {f}: {e}")
                
                if os.path.exists(mtx_gz_path):
                    adata = sc.read_10x_mtx(file_path, var_names='gene_symbols')
                else:
                    raise FileNotFoundError("Could not find matrix.mtx or matrix.mtx.gz in the 10X directory")
            else:  
                print('Loading h5ad file:', file_path)
                adata = sc.read_h5ad(file_path)
            
            if adata.X.shape[0] == 0:
                raise ValueError("Empty matrix after load")

            if scipy.sparse.issparse(adata.X) and (adata.n_obs * adata.n_vars < 1e7):
                adata.X = adata.X.toarray()

            print('Preprocessing Data')
            sc.pp.filter_cells(adata, min_genes=100)
            sc.pp.filter_genes(adata, min_cells=10)
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, n_comps=100)

            PREPROCESS_CACHE[file_path] = adata
            print('Cached preprocessed data')
            return func(adata, *args, **kwargs)

        except Exception as e:
            return jsonify({'error': f"Processing failed: {str(e)}"}), 500
    return wrapper

@app.route('/preview', methods=['POST'])
@cache_preprocessed_data
def preview(adata):
    choice = request.form.getlist('em') 
    color = request.form.get('color_umap', 'parc_cluster')
    color_scheme = request.form.get('color_scheme', 'viridis')

    valid = ['viridis', 'rainbow', 'paired', 'plasma', 'inferno']
    if color_scheme not in valid:
        color_scheme = 'viridis'

    pca_plot = umap_plot = phate_plot = None

    if 'pca' in choice:
        sc.pl.pca_variance_ratio(adata, log=False, n_pcs=50, show=False)
        pca_img = BytesIO()
        plt.savefig(pca_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        pca_plot = "data:image/png;base64," + base64.b64encode(pca_img.getvalue()).decode('utf-8')

    if 'umap' in choice:
        adata.obsm['X_umap'] = umap.UMAP(n_neighbors=20, min_dist=0.2, spread=5, init='pca').fit_transform(adata.obsm['X_pca'])
        sc.pl.embedding(adata, basis='X_umap', color=[color], palette=color_scheme, size=200, show=False, return_fig=True)
        umap_img = BytesIO()
        plt.savefig(umap_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        umap_plot = "data:image/png;base64," + base64.b64encode(umap_img.getvalue()).decode('utf-8')

    if 'phate' in choice: 
        if 'X_phate' in adata.obsm:
            fig = sc.pl.embedding(adata, basis='phate', color=[color], palette=color_scheme, show=False, return_fig=True)
            phate_img = BytesIO()
            fig.savefig(phate_img, format='png', bbox_inches='tight', dpi=120)
            plt.close(fig)
            phate_plot = "data:image/png;base64," + base64.b64encode(phate_img.getvalue()).decode('utf-8')
    
    preview_data = {
        'adata_info': {
            'dimensions': f"{adata.n_obs} cells × {adata.n_vars} genes",
            'obs_keys': list(adata.obs.keys()),
            'var_keys': list(adata.var.keys()),
            'layers': list(adata.layers.keys()),
            'uns_keys': list(adata.uns.keys()),
            'obsm_keys': list(adata.obsm.keys()),
            'varm_keys': list(adata.varm.keys())
        },
        'plots': {
            'pca': pca_plot if pca_plot else None,
            'umap': umap_plot if umap_plot else None,
            'phate': phate_plot if phate_plot else None
        },
    }
    
    return jsonify(safe_json(preview_data))

@app.route('/analyze', methods=['POST'])
@cache_preprocessed_data
def analyze(adata):
    try:
        try:
            import psutil
            if psutil.virtual_memory().available < 4 * 1024**3:  # <4GB
                gc.collect()
                if psutil.virtual_memory().available < 4 * 1024**3:
                    raise MemoryError("Insufficient memory for analysis")
        except ImportError:
            pass

        params = {
            'var_names': request.form.get('var_names'),
            'time_series_labels': request.form.get('time_series_labels'),
            'true_label': request.form.get('true_label'),
            'knn': int(request.form.get('knn', 30)),
            'cluster_graph_pruning': float(request.form.get('cluster_graph_pruning', 0.9)),
            'edgebundle_pruning': float(request.form.get('edgebundle_pruning', 0.9)),
            'edgepruning_clustering_resolution': float(request.form.get('edgepruning_clustering_resolution', 1.0)),
            'root_user': request.form.get('root_user'),
            'dpi': int(request.form.get('dpi', 120)),
            'adata_obs': request.form.get('obs'),
            'par_option': request.form.getlist('par-option')
        }

        file_data = {}
        for file_type in ['time-upload', 'velocity-matrix-upload', 'gene-matrix-upload', 'root-upload', 'csv-upload']:
            if file_type in request.files:
                file = request.files[file_type]
                if file.filename != '':
                    if file_type == 'csv-upload':
                        file_data['true_labels'] = pd.read_csv(BytesIO(file.read()))
                    else:
                        file_data[file_type] = pd.read_csv(BytesIO(file.read()))

        results = run_via_analysis(adata=adata, params=params, file_data=file_data)
        if 'error' in results:
            return jsonify(results), 500
            
        v0 = results['via_obj']

        global VIA_CACHE        
        VIA_CACHE['via_obj'] = v0

        plots = via_plot(params=params, v0=v0, file_data=file_data)
        
        return jsonify({'success': True, 'plots': plots})
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    
@app.route('/add_plots', methods=['POST'])
@cache_preprocessed_data
def add_plots(adata):
    try:
        data = request.get_json()
        lineages = data.get('marker_lineages', [])  
        genes = data.get('genes_selected') 
        
        global VIA_CACHE
        v0 = VIA_CACHE.get('via_obj')
        if v0 is None:
            return jsonify({'error': 'VIA object not found. Run analysis first.'}), 400

        plots = more_plot(lineages=lineages, genes=genes, adata=adata, v0=v0)

        return jsonify({'success': True, 'plots': plots})
    
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500


@app.route('/download_all', methods=['POST'])
def download_all():
    try:
        plot_data = request.get_json()
        mem_zip = BytesIO()
        
        with zipfile.ZipFile(mem_zip, mode='w') as zf:
            for plot_type, plot_b64 in plot_data.items():
                if plot_b64 and plot_b64.startswith('data:image/png;base64,'):
                    img_data = base64.b64decode(plot_b64.split(',')[1])
                    zf.writestr(f"{plot_type}_plot.png", img_data)
        
        mem_zip.seek(0)
        return send_file(
            mem_zip,
            as_attachment=True,
            download_name="all_plots.zip",
            mimetype='application/zip'
        )
    except Exception as e:
        return jsonify({'error': str(e)}), 500

if __name__ == '__main__':
    app.run((host='0.0.0.0', debug=True, use_reloader=False)
