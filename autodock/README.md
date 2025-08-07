# AutoDock-GPU API

A FastAPI-based web API for AutoDock-GPU molecular docking.

## Features

- Submit protein and ligand files for molecular docking
- Asynchronous docking using AutoDock-GPU for fast results
- Download docking results as a ZIP archive
- Configurable docking parameters

## Docker Base Image

**Dockerfile**: Uses nvcr.io/hpc/autodock:2020.06

### GPU Compatibility

**Supported NVIDIA GPUs**: The base image supports compute architectures 5.2, 6.0, 6.1, and 7.0 (Tesla K80, P100, V100, GTX 1060+, RTX 20-series).

**A100/H100 GPUs**: Not supported by the 2020.06 base image (requires compute architecture 8.0). If using A100/H100, you'll need to recompile AutoDock-GPU with `TARGETS="80"` or use the OpenCL version.
- NVIDIA HPC container with AutoDock-GPU pre-installed
- Includes CUDA libraries and dependencies
- Supports Meeko, RDKit, and Python 3.9

### Docker Configuration

```bash
docker-compose up --build
```

## Setup

### Prerequisites

- Docker
- NVIDIA GPU with CUDA support (for AutoDock-GPU)
- nvidia-docker runtime
- NGC API key for NVIDIA Container Registry access

### NVIDIA Container Registry Authentication

The Docker build uses NVIDIA's official AutoDock container which requires authentication:

1. Go to https://ngc.nvidia.com/signin and create a free account
2. Once logged in, go to Setup → Generate API Key
3. Copy your API key and login to the container registry:
   ```bash
   docker login nvcr.io
   # Username: $oauthtoken
   # Password: <your-ngc-api-key>
   ```

### Building and Running

```bash
docker-compose up --build
```

## API Usage

### Test the ping endpoint
```bash
curl http://localhost:8080/ping
```

### Submit a docking job
```bash
curl -X POST http://localhost:8080/inference \
  -F "protein_file=@protein.pdb" \
  -F "ligand_file=@ligand.sdf"
```

### Check status
```bash
curl http://localhost:8080/task_status/{task_id}
```

### Download results
```bash
curl -O http://localhost:8080/download_result/{task_id}
```
