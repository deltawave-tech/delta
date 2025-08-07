from fastapi import FastAPI, UploadFile, File, Form, HTTPException
from fastapi import Path as FastAPIPath
from fastapi.responses import JSONResponse, StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
import logging
import uuid
from pathlib import Path
import io
import json
from zipfile import ZipFile
from typing import Optional, List

from autodock_task import tasks, process_and_run_autodock
from models import AutoDockConfig

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Verify Meeko is available at startup
try:
    from meeko.preparation import PDBQTWriterLegacy
    from meeko.preparation import MoleculePreparation
    logger.info("✅ Meeko successfully imported on startup")
except ImportError as e:
    logger.error(f"❌ Failed to import Meeko: {str(e)}")
    logger.error("Please ensure Meeko is properly installed and in the Python path")
    import sys
    logger.error(f"Current Python path: {sys.path}")

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=['*'],
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*']
)

@app.get('/ping')
def ping():
    return JSONResponse(status_code=200, content={'message': 'pong'})

from fastapi import Depends, Request
from fastapi.datastructures import FormData

async def collect_ligand_files(request: Request):
    """
    Dependency to collect all ligand files from the form
    Clients should submit ligand files as individual fields:
    ligand_files_0, ligand_files_1, etc.
    """
    form_data = await request.form()
    ligand_files = []
    
    # Log all form keys for debugging
    logger.debug(f"Form keys: {list(form_data.keys())}")
    
    # Extract and sort all ligand_files_X entries
    indices = []
    for key in form_data.keys():
        if key.startswith("ligand_files_") and key[13:].isdigit():
            index = int(key[13:])
            indices.append((index, key))
            logger.debug(f"Found ligand file key: {key}, index: {index}")
    
    # Sort by index
    indices.sort()
    for _, key in indices:
        ligand_files.append(form_data[key])
    
    logger.debug(f"Collected {len(ligand_files)} ligand files")
    return ligand_files

@app.post('/inference')
async def run_inference(
    request: Request,
    protein_file: UploadFile = File(...),
    config: str = Form(None)
):
    # Log the request details
    logger.debug(f"Received inference request with protein file: {protein_file.filename}")
    logger.debug(f"Content-Type: {request.headers.get('content-type', 'Not provided')}")
    
    # Collect ligand files
    ligand_files = await collect_ligand_files(request)
    logger.info(f"Processing request with {len(ligand_files)} ligand files")
    
    try:
        # Verify meeko is available before processing
        try:
            from meeko.preparation import PDBQTWriterLegacy
            from meeko.preparation import MoleculePreparation
            logger.info("✅ Meeko available for this request")
        except ImportError as e:
            logger.error(f"❌ Meeko not available: {str(e)}")
            import sys
            logger.error(f"Python path: {sys.path}")
            return JSONResponse(
                status_code=500,
                content={'error': "Meeko package not available. Please check server configuration.", 
                         'detail': str(e)}
            )

        task_id = str(uuid.uuid4())
        tasks[task_id] = {'task_id': task_id, 'status': 'queued'}

        # Parse config JSON string if provided
        config_obj = None
        if config:
            try:
                config_dict = json.loads(config)
                config_obj = AutoDockConfig(**config_dict)
            except json.JSONDecodeError as e:
                raise HTTPException(status_code=400, detail=f"Invalid JSON in config: {str(e)}")
        else:
            config_obj = AutoDockConfig()

        # Read protein file content
        protein_content = await protein_file.read()
        if isinstance(protein_content, str):
            protein_content = protein_content.encode('utf-8')
        
        if not ligand_files:
            raise HTTPException(status_code=400, detail="No ligand files provided. Please provide at least one ligand file.")
        
        # Read all ligand file contents
        # If only one ligand is provided, use single ligand mode, otherwise use batch mode
        if len(ligand_files) == 1:
            # Single ligand mode
            ligand_content = await ligand_files[0].read()
            if isinstance(ligand_content, str):
                ligand_content = ligand_content.encode('utf-8')
                
            logger.debug(f"InferenceRequest for task {task_id} (single ligand)")
            logger.debug(f"Config: {config_obj}")
        else:
            # Batch mode - read all ligand files
            ligand_contents = []
            for i, ligand_file in enumerate(ligand_files):
                content = await ligand_file.read()
                if isinstance(content, str):
                    content = content.encode('utf-8')
                ligand_contents.append(content)
                
            ligand_content = ligand_contents  # Assign the list to ligand_content for the batch mode
            
            logger.debug(f"Batch InferenceRequest for task {task_id} with {len(ligand_files)} ligands")
            logger.debug(f"Config: {config_obj}")

        # Start the docking task - process_and_run_autodock handles both single and batch modes
        await process_and_run_autodock(
            task_id,
            protein_content,
            ligand_content,
            config_obj,
            Path("storage/results") / task_id
        )

        return JSONResponse(status_code=202, content=tasks[task_id])

    except Exception as e:
        logger.exception(f'An exception occurred: {str(e)}')
        # Log detailed environment information to help diagnose the issue
        try:
            import sys
            logger.error(f"Python version: {sys.version}")
            logger.error(f"Python path: {sys.path}")
            logger.error(f"Environment variables:")
            import os
            for key, value in os.environ.items():
                if key.startswith('PYTHON') or key.startswith('PATH'):
                    logger.error(f"  {key}: {value}")
        except Exception as env_error:
            logger.error(f"Error fetching environment info: {str(env_error)}")
            
        return JSONResponse(
            status_code=500,
            content={'error': "An error occurred while processing your request", "detail": str(e)}
        )

@app.get('/task_status/{task_id}')
def get_task_status(task_id: str = FastAPIPath(...)):
    task = tasks.get(task_id, None)
    if not task:
        raise HTTPException(status_code=404, detail='Task not found')
    return JSONResponse(status_code=200, content=task)

@app.get("/download_result/{task_id}")
async def download_result(task_id: str):
    if task_id not in tasks:
        raise HTTPException(status_code=404, detail="Task not found")

    task = tasks[task_id]
    if task['status'] != 'completed':
        raise HTTPException(status_code=400, detail="Task not completed")

    results_dir = Path(task['output_dir'])
    if not results_dir.exists():
        raise HTTPException(status_code=404, detail="Results not found")

    zip_buffer = io.BytesIO()
    with ZipFile(zip_buffer, 'w') as zip_file:
        # Recursively add all files and directories
        # But exclude large working directories that aren't needed
        excluded_dirs = ["work"]
        
        # Function to recursively add files
        def add_directory_to_zip(directory, parent_path=""):
            for path in directory.iterdir():
                # Skip excluded directories
                if path.is_dir() and path.name in excluded_dirs:
                    logger.debug(f"Skipping excluded directory: {path.name}")
                    continue
                    
                # Calculate the path within the zip file
                rel_path = f"{parent_path}/{path.name}" if parent_path else path.name
                
                if path.is_file():
                    # Add the file to the zip
                    zip_file.write(path, rel_path)
                    logger.debug(f"Adding file to zip: {rel_path}")
                elif path.is_dir():
                    # Recursively add the directory
                    add_directory_to_zip(path, rel_path)
        
        # Start the recursive process from the results directory
        add_directory_to_zip(results_dir)

    # Now that we've packaged everything, clean up the work directory
    # to free resources and disk space
    work_dir = results_dir / "work"
    if work_dir.exists():
        try:
            import shutil
            logger.info(f"Cleaning up work directory after zip creation: {work_dir}")
            shutil.rmtree(work_dir)
            logger.info("Work directory successfully removed")
        except Exception as e:
            logger.warning(f"Failed to remove work directory: {str(e)}")
    
    zip_buffer.seek(0)
    return StreamingResponse(
        iter([zip_buffer.getvalue()]),
        media_type="application/zip",
        headers={
            "Content-Disposition": f"attachment; filename=autodock_results_{task_id}.zip"
        }
    )
