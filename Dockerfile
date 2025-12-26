FROM python:3.9-slim

# Install ALL system dependencies needed
RUN apt-get update && apt-get install -y \
    build-essential \
    libopenblas-dev \
    libcairo2-dev \
    libffi-dev \
    pkg-config \
    python3-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy requirements first
COPY requirements.txt .

# Upgrade pip and set up environment
RUN pip install --no-cache-dir -r requirements.txt

# Install base packages first (ones that don't need compilation)
RUN pip install --no-cache-dir \
    flask==2.3.3 \
    gunicorn==21.2.0 \
    pandas==2.0.3 \
    numpy==1.24.3 \
    scipy==1.11.1 \
    scikit-learn==1.3.0 \
    flask-sqlalchemy==3.0.5

# Now install the problematic packages with specific flags
RUN pip install --no-cache-dir \
    matplotlib==3.7.2 \
    umap-learn==0.5.3

# Try to install scanpy without some optional dependencies
RUN pip install --no-cache-dir "scanpy[leiden]==1.9.6" --no-deps

# Add SQLite support
RUN apt-get update && apt-get install -y sqlite3 libsqlite3-dev

# Copy application code
COPY . .

# Set port and run
ENV PORT 8080
CMD exec gunicorn --bind :$PORT --workers 1 --threads 8 --timeout 0 app:app