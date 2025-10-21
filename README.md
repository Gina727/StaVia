============ WELCOME TO STAVIA ============

## Prerequisites

1. **Docker Desktop**: [Download here](https://www.docker.com/products/docker-desktop)
2. **Git** (optional): To clone this repository

## Quick Start

### For Windows Users:
1. Double-click on `run.bat`
2. Wait for the container to build (first time only)
3. Open your browser to `http://localhost:8080`

### For Mac/Linux Users:
1. Open Terminal
2. Navigate to the application folder: `cd path/to/scrna-web-app`
3. Run: `./run.sh`
4. Open your browser to `http://localhost:8080`

## Using the Application

1. **Upload your data**: Use the file upload button in the web interface
2. **Run analysis**: Follow the steps in the web interface
3. **Download results**: Use the download button to download the plots generated

## Important Notes

- The first run will take some time to download and build the environment
- All uploaded files and results are saved in the `data` folder
- To stop the application, press `Ctrl+C` in the terminal window
- Your data persists in the `data` folder between sessions

## Troubleshooting

**Problem**: "Docker command not found"
**Solution**: Install Docker Desktop and make sure it's running

**Problem**: Port 8080 is already in use
**Solution**: Stop other applications using port 8080 or modify the port in run.sh/run.bat and Dockerfile

**Problem**: Out of memory errors
**Solution**: Close other applications or use a computer with more RAM