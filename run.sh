#!/bin/bash
echo "Building scRNA-seq analysis container..."
docker build -t scrna-app .

echo "Starting application..."
echo "Open http://localhost:8080 in your browser when ready"
echo "Press Ctrl+C to stop the application"

docker run -p 8080:8080 -v $(pwd)/data:/app/instance scrna-app