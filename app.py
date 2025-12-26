import pandas as pd
import scanpy as sc
from via_analysis import run_via_analysis
from plotting import via_plot, more_plot
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
import loompy
import phate
import anndata as ad
from scipy import sparse
from google.cloud import storage  # <-- ADD THIS
from google.auth import default   # <-- ADD THIS

app = Flask(__name__)
app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///test.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['UPLOAD_FOLDER'] = '/app/instance/uploads'
app.config['MAX_CONTENT_LENGTH'] = 1024 * 1024 * 1024 * 20
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

        files = request.files
        
        if 'file' in files:  
            file = files['file']
            if not file.filename.lower().endswith('.h5ad'):
                return jsonify({'error': 'Only .h5ad files supported'}), 400
            
            file_path = os.path.join(app.config['UPLOAD_FOLDER'], f'uploaded_data_{unique_id}.h5ad')
            file.save(file_path)
            session['file_type'] = 'h5ad'
            
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

        session['file_path'] = file_path
        
        file_path = session.get('file_path')
        file_type = session.get('file_type')

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
    
def velocity_df_to_adata(velocity_df, adata_reference, common_genes=None, layer_name='velocity'):
    
    if velocity_df is None: 
        return adata_reference
    
    print(f"=== velocity_df_to_adata START ===")
    print(f"Input velocity_df shape: {velocity_df.shape}")
    print(f"Input adata shape: {adata_reference.shape}")
    
    # Clean velocity cell names
    velocity_df_cleaned = velocity_df.copy()
    # velocity_df_cleaned.index = [
    #     str(name).replace('Het_CR_outs:', '').replace('x', '-1') 
    #     for name in velocity_df_cleaned.index
    # ]
    
    # Clean adata cell names
    adata_cleaned = adata_reference.copy()
    # adata_cleaned.obs_names = [
    #     str(name).replace('Het_CR_outs:', '').replace('x', '-1')
    #     for name in adata_cleaned.obs_names
    # ]
    
    # Find common cells
    common_cells = list(velocity_df_cleaned.index.intersection(adata_cleaned.obs_names))
    
    if len(common_cells) == 0:
        raise ValueError(f"No common cells found")
    
    # Use provided common_genes or calculate intersection
    if common_genes is not None:
        print(f"Using provided common_genes: {len(common_genes)} genes")
        genes_to_use = list(set(common_genes).intersection(set(velocity_df_cleaned.columns)))
    else:
        genes_to_use = list(velocity_df_cleaned.columns.intersection(adata_cleaned.var_names))
        print(f"Calculated common genes: {len(genes_to_use)} genes")
    
    # CHECK FOR DUPLICATE GENE NAMES IN VELOCITY
    duplicate_genes = velocity_df_cleaned.columns[velocity_df_cleaned.columns.duplicated()].tolist()
    if duplicate_genes:
        print(f"WARNING: Velocity DataFrame has {len(duplicate_genes)} duplicate gene names")
        print(f"Sample duplicates: {duplicate_genes[:10]}")
        
        # Remove duplicates - keep first occurrence
        velocity_df_cleaned = velocity_df_cleaned.loc[:, ~velocity_df_cleaned.columns.duplicated()]
        print(f"After removing duplicates: {velocity_df_cleaned.shape}")
        
        # Recalculate genes_to_use after removing duplicates
        genes_to_use = [gene for gene in genes_to_use if gene in velocity_df_cleaned.columns]
        print(f"Genes after removing duplicates: {len(genes_to_use)}")
    
    if len(genes_to_use) == 0:
        raise ValueError(f"No common genes found")
    
    print(f"Filtering to {len(common_cells)} cells and {len(genes_to_use)} genes")
    print(f"genes_to_use sample: {genes_to_use[:10]}")
    
    # DEBUG: Check if genes_to_use are actually in velocity_df_cleaned
    missing_genes = [gene for gene in genes_to_use if gene not in velocity_df_cleaned.columns]
    if missing_genes:
        print(f"WARNING: {len(missing_genes)} genes in genes_to_use not in velocity columns")
        print(f"Sample missing: {missing_genes[:10]}")
    
    # Filter velocity DataFrame - be explicit about columns
    velocity_filtered = velocity_df_cleaned[genes_to_use]  # Filter columns first
    velocity_filtered = velocity_filtered.loc[common_cells]  # Then filter rows
    
    print(f"After column filtering: {velocity_filtered.shape}")
    
    # Filter adata
    # Get original cell names for subsetting
    cell_mapping = dict(zip(adata_cleaned.obs_names, adata_reference.obs_names))
    common_cells_original = [cell_mapping[cell] for cell in common_cells]
    
    adata_filtered = adata_reference[common_cells_original, genes_to_use]
    adata_filtered.obs_names = common_cells
    
    print(f"Final velocity shape: {velocity_filtered.shape}")
    print(f"Final adata shape: {adata_filtered.shape}")
    
    if velocity_filtered.shape != adata_filtered.shape:
        print(f"ERROR: Shape mismatch after filtering!")
        print(f"  Velocity columns: {velocity_filtered.columns.tolist()[:10]}")
        print(f"  Adata var_names: {adata_filtered.var_names.tolist()[:10]}")
        raise ValueError(f"Shape mismatch: velocity {velocity_filtered.shape} vs adata {adata_filtered.shape}")
    
    # Convert to sparse matrix
    velocity_sparse = sparse.csr_matrix(velocity_filtered.values)
    adata_filtered.layers[layer_name] = velocity_sparse
    
    print(f"✓ Successfully added velocity to adata layers")
    
    return adata_filtered

def velocity_adata_to_df(adata, layer_name = 'velocity'):
    if layer_name in adata.layers:
        velocity_data = adata.layers[layer_name]
        
        # Convert sparse to dense if needed
        if sparse.issparse(velocity_data):
            velocity_array = velocity_data.toarray()
        else:
            velocity_array = velocity_data
            
        # Create DataFrame with proper index and columns
        df = pd.DataFrame(
            velocity_array,
            index=adata.obs_names,
            columns=adata.var_names
        )
        return df
    
    elif 'velocity' in adata.obsm:
        # Handle velocity in obsm
        velocity_array = adata.obsm['velocity']
        df = pd.DataFrame(
            velocity_array,
            index=adata.obs_names,
            columns=[f"vel_{i}" for i in range(velocity_array.shape[1])]
        )
        return df
    
    else:
        raise KeyError(f"Velocity data not found. Available layers: {list(adata.layers.keys())}")

def get_storage_client():
    """Get Google Cloud Storage client with fallback authentication"""
    try:
        client = storage.Client()
        return client
    except Exception as auth_error:
        print(f"Auth error, using default credentials: {auth_error}")
        credentials, project = default()
        client = storage.Client(credentials=credentials, project=project)
        return client
    
@app.route('/generate-signed-url', methods=['POST'])
def generate_signed_url():
    """Generate a signed URL for direct upload to Cloud Storage"""
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No JSON data received'}), 400
        
        filename = data.get('filename')
        file_type = data.get('file_type', 'h5ad')  
        
        if not filename:
            return jsonify({'error': 'Filename required'}), 400
        
        try:
            client = storage.Client()
        except Exception as auth_error:
            print(f"Auth error: {auth_error}")
            # Fallback: use default credentials
            credentials, project = default()
            client = storage.Client(credentials=credentials, project=project)
        
        bucket_name = 'stavia-uploads'  # Make sure this bucket exists!
        bucket = client.bucket(bucket_name)
        
        # Create unique filename
        unique_filename = f"uploads/{uuid.uuid4()}_{filename}"
        blob = bucket.blob(unique_filename)

        # Generate signed URL (valid for 1 hour)
        url = blob.generate_signed_url(
            version='v4',
            expiration=3600,  # 1 hour
            method='PUT',
            content_type='application/octet-stream'
        )
        
        return jsonify({
            'signed_url': url,
            'file_path': blob.name,
            'filename': filename
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500
    
@app.route('/process-uploaded-file', methods=['POST'])
def process_uploaded_file():
    """Process file after it's uploaded to Cloud Storage"""
    try:
        data = request.get_json()
        file_path = data.get('file_path')
        file_type = data.get('file_type')  # h5ad, mtx, barcodes, genes, anno
        
        if not file_path:
            return jsonify({'error': 'File path required'}), 400
        
        client = get_storage_client()
        bucket = client.bucket('stavia-uploads')
        blob = bucket.blob(file_path)

        # Download file to temporary location
        temp_dir = tempfile.mkdtemp()
        local_path = os.path.join(temp_dir, os.path.basename(file_path))
        blob.download_to_filename(local_path)
        
        # Process based on file type
        if file_type == 'h5ad':
            adata = sc.read_h5ad(local_path)
        elif file_type == 'mtx':
            # Handle 10X files
            pass
        
        # Store in your cache and process...
        INITIAL_ADATA_CACHE[local_path] = adata.copy()
        
        # Clean up
        os.remove(local_path)
        os.rmdir(temp_dir)

        return jsonify({
            'success': True,
            'message': 'File processed successfully',
            'dimensions': f"{adata.n_obs} cells × {adata.n_vars} genes"
        })
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500

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

            print('Preprocessing Data')
            sc.pp.filter_cells(adata, min_genes=100)
            print("Pass cell filtering")
            sc.pp.filter_genes(adata, min_cells=10)
            print("Pass gene filtering")
            sc.pp.normalize_total(adata)
            print("Pass data normalization")
            sc.pp.log1p(adata)
            print("Pass log normalization")
            sc.pp.pca(adata, n_comps=100)
            print("Pass PCA")

            PREPROCESS_CACHE[file_path] = adata
            print('Cached preprocessed data')

        return func(adata, *args, **kwargs)
        
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
        phate_op = phate.PHATE(n_components=2, random_state=42)
        X_phate = phate_op.fit_transform(adata.X.T)
        adata.obsm['X_phate'] = X_phate

        plt.figure(figsize=(10, 8))
        sc.pl.umap(adata, color='PHATE')
        phate_img = BytesIO()
        plt.savefig(phate_img, format='png', bbox_inches='tight', dpi=120)
        plt.close('all')
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

def handle_uploaded_file(file):
    """Handle both CSV and loom files, always return DataFrame"""
    if file.filename.endswith('.loom'):
        return loom_to_dataframe(file)
    
    elif file.filename.endswith('.h5ad'):
        file_content = BytesIO(file.read())
        
        with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp_file:
            file_content.seek(0)
            tmp_file.write(file_content.read())
            tmp_path = tmp_file.name
        
        try:
            # Just read it as AnnData
            adata = sc.read_h5ad(tmp_path)
            print(f"Loaded h5ad: shape={adata.shape}")
            print(f"Layers: {list(adata.layers.keys())}")
            return adata  # Return AnnData object
            
        finally:
            os.unlink(tmp_path)
    else:  
        return pd.read_csv(BytesIO(file.read()))
    
def loom_to_dataframe(file):
    """Convert loom file to pandas DataFrame"""
    file_content = BytesIO(file.read())
    
    with tempfile.NamedTemporaryFile(delete=False, suffix='.loom') as tmp_file:
        file_content.seek(0)
        tmp_file.write(file_content.read())
        tmp_path = tmp_file.name
    
    try:
        with loompy.connect(tmp_path) as ds:
            # Use main matrix
            matrix = ds[:,:].T
            
            # Try to get cell and gene names
            index = ds.ca.CellID if 'CellID' in ds.ca else range(matrix.shape[0])
            columns = ds.ra.Gene if 'Gene' in ds.ra else range(matrix.shape[1])
            
            df = pd.DataFrame(matrix, index=index, columns=columns)
        return df
    finally:
        os.unlink(tmp_path)

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
        velocity_adata = None

        for file_type in ['time-upload', 'velocity-matrix-upload', 'gene-matrix-upload', 'root-upload', 'csv-upload']:
            if file_type in request.files:
                file = request.files[file_type]
                if file.filename != '':
                    try:
                        if file_type == 'csv-upload':
                            file_data['true_labels'] = handle_uploaded_file(file)
                        
                        elif file_type == 'velocity-matrix-upload':
                            vel_file = request.files['velocity-matrix-upload']
                            velocity_data = handle_uploaded_file(vel_file)
                        
                            if isinstance(velocity_data, ad.AnnData):
                                print("Got AnnData from h5ad file")
                                velocity_adata = velocity_data
                                file_data['velocity_adata'] = velocity_adata
                                file_data['velocity-matrix-upload'] = velocity_adata_to_df(adata=velocity_adata, layer_name = 'velocity')
                            else:
                                print("Got DataFrame from file")
                                file_data['velocity-matrix-upload'] = velocity_data
                                file_data['velocity_adata'] = velocity_df_to_adata(velocity_df=velocity_data, adata_reference=adata, layer_name='velocity')
                        else:
                            file_data[file_type] = handle_uploaded_file(file)
                    except Exception as e:
                        print(f"Error processing {file.filename}: {e}")
                        file_data[file_type] = None

        results = run_via_analysis(adata=adata, params=params, file_data=file_data)
        if 'error' in results:
            return jsonify(results), 500
            
        v0 = results['via_obj']
        adata = results['adata']

        global VIA_CACHE        
        VIA_CACHE['via_obj'] = v0

        plots = via_plot(params=params, v0=v0, file_data=file_data, adata=adata)
        
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
    app.run(host='0.0.0.0', use_reloader=False, port=10000, debug=True, reloader_type='watchdog', threaded=True)
else: 
    with app.app_context():
        print("Creating database tables...")
        db.create_all()
        print("Database tables created!")