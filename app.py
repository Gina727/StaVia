import pandas as pd
import scanpy as sc

from via_analysis import run_via_analysis
from plotting import via_plot, more_plot
from file_functions import read_10x_mtx

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
import umap
import numpy as np 
from functools import wraps
import gc
import scipy
from datetime import datetime
import uuid
import gzip
import shutil
import parc

# Initializing Flask App and the SQL database 
app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///test.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['UPLOAD_FOLDER'] = '/app/instance/uploads'
app.config['MAX_CONTENT_LENGTH'] = 50 * 1024 * 1024
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
INITIAL_ADATA_CACHE = {} 

@app.before_request
def clear_trailing():
    session.modified = True

# Automatically cleans up temporary files when the program exits (reload page)
def cleanup():
    try:
        shutil.rmtree(UPLOAD_FOLDER)
    except Exception as e:
        print(f"Error cleaning up upload folder: {e}")

atexit.register(cleanup)

# Just for handling edge cases 
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

# Formatting parameters stored in the database (since they're stored in json)
@app.template_filter('format_params')
def format_params(params):
    if not params:
        return "No parameters"
    result = []
    for key, value in params.items():
        if isinstance(value, list):
            value = ", ".join(value)
        result.append(f"{key}: {value}")
    return "\n".join(result)  
app.jinja_env.filters['format_params'] = format_params

# Set home route 
@app.route('/')
def home():
    return render_template('home.html')

# Reload page when tasks are added 
@app.route('/upload_page', methods=['GET'])
def upload_page():
    tasks = Todo.query.order_by(Todo.date_created).all()
    return render_template('index.html', tasks=tasks)   

# Adding tasks to the database 
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

# Removing tasks from the database
@app.route('/delete/<int:id>')
def delete(id):
    task_to_delete = Todo.query.get_or_404(id)

    try:
        db.session.delete(task_to_delete)
        db.session.commit()
        return redirect('/upload_page')
    except:
        return 'There was a problem deleting that task'

# Just a route to the About page 
@app.route('/about_page')
def about_page():
    return render_template('about.html')

# Handles main file upload (the two main sections at the top)
@app.route('/upload', methods=['POST'])
def upload():
    try:
        # Reset the cache to accept new files 
        PREPROCESS_CACHE.clear()
        VIA_CACHE.clear()
        session.clear()
        unique_id = str(uuid.uuid4())[:8]

        files = request.files
        
        # Saving these files in temporary directories 

        # For .h5ad files 
        if 'file' in files:  
            file = files['file']
            if not file.filename.lower().endswith('.h5ad'):
                return jsonify({'error': 'Only .h5ad files supported'}), 400
            
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], f'uploaded_data_{unique_id}.h5ad')
            file.save(file_path)
            session['file_type'] = 'h5ad'
        # For matrix.mtx, barcodes, and features files (names MUST be as specified below (convention))
        else:  
            has_matrix = any(f for f in files if 'matrix' in f.lower())
            has_barcodes = any(f for f in files if 'barcodes' in f.lower())
            has_features = any(f for f in files if 'features' in f.lower() or 'genes' in f.lower())

            if has_matrix and has_barcodes and has_features:
                file_path = os.path.join(app.config['UPLOAD_FOLDER'], f'10x_data_{unique_id}')
                os.makedirs(file_path, exist_ok=True)

                mtx_file = next(f for f in files.values() if 'matrix' in f.filename.lower())
                barcodes_file = next(f for f in files.values() if 'barcodes' in f.filename.lower())
                features_file = next(f for f in files.values() if 'features' in f.filename.lower() or 'gene' in f.filename.lower())

            mtx_path = os.path.join(file_path, 'matrix.mtx')
            barcodes_path = os.path.join(file_path, 'barcodes.tsv')
            features_path = os.path.join(file_path, 'features.tsv')

            mtx_file.save(mtx_path)
            barcodes_file.save(barcodes_path)
            features_file.save(features_path)

            # Change all file types into .gz to be processed 
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
        
        # Check if the annotation file is uploaded or not 
        annotation_data = None
        if 'anno' in files: 
            annotation_file = files['anno']
            if annotation_file.filename != '':  
                if not any(annotation_file.filename.lower().endswith(ext) for ext in ['.txt', '.csv', '.tsv']):
                    return jsonify({'error': 'Annotation file must be .txt, .csv, or .tsv'}), 400
                
                file_extension = os.path.splitext(annotation_file.filename)[1]
                annotation_path = os.path.join(app.config['UPLOAD_FOLDER'], f'annotation_{unique_id}{file_extension}')
                annotation_file.save(annotation_path)
                session['annotation_path'] = annotation_path

                try:
                    if annotation_path.endswith('.csv'):
                        annotation_data = pd.read_csv(annotation_path)
                    else: 
                        annotation_data = pd.read_csv(annotation_path, sep='\t')
                    session['annotation_data'] = annotation_data.to_json()
                    
                except Exception as e:
                    print(f"Error reading annotation file: {e}")
                    session['annotation_error'] = str(e)

        # Storing the file path in the user's session 
        session['file_path'] = file_path
        
        # Retrieving Session data 
        file_path = session.get('file_path')
        file_type = session.get('file_type')

        # Loading data from the main files 

        # I added functions for listing directories just for sanity check 
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
                try:
                    print("Using scanpy 10x data reader")
                    print(f"Reading from: {file_path}")
                    print(f"Files in directory: {os.listdir(file_path)}")
        
                    adata = sc.read_10x_mtx(file_path, var_names='gene_symbols', cache=True)
                    print("Scanpy reader successful")
                except Exception as e:
                    print(f"Error with scanpy 10x reader: {str(e)}")
                    print(f"Error type: {type(e).__name__}")
                    import traceback
                    print(f"Traceback: {traceback.format_exc()}")
                    print("Falling back to manual 10x data reader")
    
                try:
                    adata = read_10x_mtx(file_path)
                    print("Manual reader successful")
                except Exception as e2:
                    print(f"Error with manual reader: {e2}")
                    raise Exception(f"Both readers failed. Final error: {e2}")
            
            # if os.path.exists(mtx_gz_path): 
            #     try:
            #         print("Using scanpy 10x data reader")
            #         adata = sc.read_10x_mtx(file_path, var_names='gene_symbols', cache=True)
            #     except Exception as e:
            #         print(f"Error making adata: {e}")
                # except:
                    # print("Using manual 10x data reader")
                    # adata = read_10x_mtx(file_path)
            # raise FileNotFoundError("Could not find matrix.mtx or matrix.mtx.gz in the 10X directory")
        else:  
            print('Loading h5ad file:', file_path)
            adata = sc.read_h5ad(file_path)

        if annotation_data is not None:
            try:
                if len(annotation_data) == adata.n_obs:
                    for col_name in annotation_data.columns:
                        if col_name not in adata.obs.columns:  
                            adata.obs[col_name] = annotation_data[col_name].values
                    print(f"Added {len(annotation_data.columns)} annotation columns")
                else:
                    print(f"Warning: Annotation rows ({len(annotation_data)}) don't match cells ({adata.n_obs})")
            except Exception as e:
                print(f"Error applying annotation: {e}")

        INITIAL_ADATA_CACHE[file_path] = adata.copy()

        # Return Anndata object information read from the main files 
        return jsonify({'success': True, 
                        'message': 'Files uploaded successfully',
                        'adata_info': {
                            'dimensions': f"{adata.n_obs} cells × {adata.n_vars} genes",
                            'obs_keys': list(adata.obs.keys()),
                            'var_keys': list(adata.var.keys()),
                            'layers': list(adata.layers.keys()),
                            'uns_keys': list(adata.uns.keys()),
                            'obsm_keys': list(adata.obsm.keys()),
                            'varm_keys': list(adata.varm.keys())
                            }})
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Preprocessing the data and caching it 
def cache_preprocessed_data(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        print("\n=== CACHE CHECK ===")
        print("Session file_path:", session.get('file_path'))
        print("Session file_type:", session.get('file_type'))

        file_path = session.get('file_path')
        if not file_path:
            print("ERROR: No file_path in session")
            return jsonify({'error': 'No file uploaded'}), 400

        # Using cached data from either this function or initial cache (data before preprocessing)
        if file_path in PREPROCESS_CACHE:
            print("Using cached preprocessed data")
            adata = PREPROCESS_CACHE[file_path]
        else:
            if file_path not in INITIAL_ADATA_CACHE:
                return jsonify({'error': 'Data not found. Please upload again.'}), 400

            print("Processing data from initial cache")
            adata = INITIAL_ADATA_CACHE[file_path].copy()
            if scipy.sparse.issparse(adata.X) and (adata.n_obs * adata.n_vars < 1e7):
                adata.X = adata.X.toarray()

            # Standard procedures of preprocessing data ***WHAT NEEDS TO BE CHANGED***
            print('Preprocessing Data')
            sc.pp.filter_cells(adata, min_genes=100)
            sc.pp.filter_genes(adata, min_cells=10)
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            sc.pp.pca(adata, n_comps=100)

            # Caching the preprocessed data 
            PREPROCESS_CACHE[file_path] = adata
            print('Cached preprocessed data')

        return func(adata, *args, **kwargs)
        
    return wrapper

# Previewing data with PCA, UMAP, or PHATE (Data preprocessing automatically runs before it in the first run)
@app.route('/preview', methods=['POST'])
@cache_preprocessed_data
def preview(adata):
    # Get choices (from checkboxes) and colors 
    choice = request.form.getlist('em') 
    color = request.form.get('color_umap', 'parc_cluster')
    color_scheme = request.form.get('color_scheme', 'viridis')

    valid = ['viridis', 'rainbow', 'paired', 'plasma', 'inferno']
    if color_scheme not in valid:
        color_scheme = 'viridis'

    pca_plot = umap_plot = parc_plot = None

    if 'pca' in choice:
        sc.pl.pca_variance_ratio(adata, log=False, n_pcs=50, show=False)
        pca_img = BytesIO()
        plt.savefig(pca_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        pca_plot = "data:image/png;base64," + base64.b64encode(pca_img.getvalue()).decode('utf-8')

    if 'umap' in choice:
        adata.obs['X_umap'] = umap.UMAP(n_neighbors=20, min_dist=0.2, spread=5, init='pca').fit_transform(adata.obsm['X_pca'])
        sc.pl.embedding(adata, basis='X_umap', color=[color], palette=color_scheme, size=200, show=False, return_fig=True)
        umap_img = BytesIO()
        plt.savefig(umap_img, format='png', bbox_inches='tight', dpi=120)
        plt.close()
        umap_plot = "data:image/png;base64," + base64.b64encode(umap_img.getvalue()).decode('utf-8')

    if 'parc' in choice: 
        parc1 = parc.PARC(adata.obsm['X_pca'], jac_std_global=0.15, random_seed =1, small_pop = 50)
        parc1.run_PARC()
        parc_labels = parc1.labels
        adata.obs["PARC"] = pd.Categorical(parc_labels)
        PREPROCESS_CACHE[session['file_path']] = adata

        sc.settings.n_jobs=4
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata)

        plt.figure(figsize=(8,6))
        sc.pl.umap(adata, color='PARC')
        parc_img = BytesIO()
        plt.savefig(parc_img, format='png', bbox_inches='tight', dpi=120)
        plt.close('all')
        parc_plot = "data:image/png;base64," + base64.b64encode(parc_img.getvalue()).decode('utf-8')
        print(f"PARC plot size: {len(parc_plot)//1024} KB")

    # Also returns Anndata information (after preprocessing)
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
            'parc': parc_plot if parc_plot else None
        },
    }
    
    return jsonify(safe_json(preview_data))

# Main analysis pipeline 
@app.route('/analyze', methods=['POST'])
@cache_preprocessed_data
def analyze(adata):
    try:
        # Checking memory threshold (again, for sanity check because most files are very big)
        try:
            import psutil
            if psutil.virtual_memory().available < 4 * 1024**3:  # <4GB
                gc.collect()
                if psutil.virtual_memory().available < 4 * 1024**3:
                    raise MemoryError("Insufficient memory for analysis")
        except ImportError:
            pass
        
        # Getting all of the parameters from frontend 
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

        # Saving files other than the main oens 
        file_data = {}
        for file_type in ['time-upload', 'velocity-matrix-upload', 'gene-matrix-upload', 'root-upload', 'csv-upload']:
            if file_type in request.files:
                file = request.files[file_type]
                if file.filename != '':
                    # Added this because of naming conflict, you can change the name later to make it neater
                    if file_type == 'csv-upload':
                        file_data['true_labels'] = pd.read_csv(BytesIO(file.read()))
                    else:
                        file_data[file_type] = pd.read_csv(BytesIO(file.read()))

        # RUN VIA ANALYSIS (See via_analysis.py)
        results = run_via_analysis(adata=adata, params=params, file_data=file_data)
        if 'error' in results:
            return jsonify(results), 500
            
        v0 = results['via_obj']

        # Store the via object for further plotting 
        global VIA_CACHE        
        VIA_CACHE['via_obj'] = v0

        # Starts plotting (See plotting.py)
        plots = via_plot(params=params, v0=v0, file_data=file_data)
        
        return jsonify({'success': True, 'plots': plots})
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

# Generates additional plots     
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

        # See plotting.py 
        plots = more_plot(lineages=lineages, genes=genes, adata=adata, v0=v0)

        return jsonify({'success': True, 'plots': plots})
    
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)}), 500

# Downloads all plots into a zip file 
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

# When the script is run directly, run the app 
if __name__ == '__main__':
    app.run(host='0.0.0.0', use_reloader=False, port=10000, debug=True)
# When another script imports this file as a module, create a database 
# Actually problematic. You should move this to the top (below database initialization), 
# but it didn't work for me when exporting the app. You can try.
else: 
    with app.app_context():
        print("Creating database tables...")
        db.create_all()
        print("Database tables created!")