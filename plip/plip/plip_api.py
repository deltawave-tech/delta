from fastapi import FastAPI, HTTPException, Path as FastAPIPath, File, UploadFile, Form
from fastapi.responses import JSONResponse, StreamingResponse
from fastapi.middleware.cors import CORSMiddleware
from typing import Optional, List
import logging
import sys
from pathlib import Path
import io
import json
from zipfile import ZipFile

from plip_task import process_task, get_task_status, list_tasks

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    force=True
)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

# Disable other verbose loggers
logging.getLogger("plip").setLevel(logging.WARNING)
logging.getLogger("openbabel").setLevel(logging.WARNING)
logging.getLogger("urllib3").setLevel(logging.WARNING)

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
    """Health check endpoint"""
    return JSONResponse(status_code=200, content={'message': 'pong'})

@app.post('/inference')
async def inference(
    file: UploadFile | None = File(None),
    body: str = Form(default="{}")
):
    """Start PLIP analysis and return task ID"""
    try:
        logger.info("Received inference request")

        # Create request data with default output formats
        request_data = {
            "output_format": ["xml", "txt"]
        }

        # Handle file upload or file content from body
        if file:
            content = await file.read()
            request_data["file_content"] = content.decode()
        else:
            # Parse body for PDB ID or file content
            try:
                body_data = json.loads(body)
                if 'file_content' in body_data:
                    request_data["file_content"] = body_data['file_content']
                elif 'pdb_id' in body_data:
                    request_data["pdb_id"] = body_data['pdb_id']
            except json.JSONDecodeError:
                pass

        # Validate input
        if "file_content" not in request_data and "pdb_id" not in request_data:
            raise HTTPException(
                status_code=400,
                detail="Either file or pdb_id in body must be provided"
            )

        task_id = await process_task(request_data)
        logger.info(f"Created task: {task_id}")

        return JSONResponse(
            status_code=202,
            content={'task_id': task_id}
        )

    except Exception as e:
        logger.error(f"Error submitting inference task: {str(e)}")
        raise HTTPException(status_code=500, detail=str(e))

@app.get('/task_status/{task_id}')
async def check_task_status(task_id: str):
    """Get status of a specific task"""
    status = await get_task_status(task_id)
    if status == "not_found":
        raise HTTPException(status_code=404, detail="Task not found")
    return status

@app.get('/tasks')
def get_tasks():
    """List all tasks"""
    return JSONResponse(status_code=200, content=list_tasks())

@app.get('/download/{task_id}')
async def download_results(task_id: str = FastAPIPath(...)):
    """Download task results as txt file"""
    status = await get_task_status(task_id)
    if status == "not_found":
        return JSONResponse(
            content={"error": "Task not found"},
            status_code=404
        )

    if status != "completed":
        return JSONResponse(
            content={"error": "Task not completed"},
            status_code=400
        )

    task_dir = Path("storage") / task_id
    txt_file = task_dir / f"report.txt"  # Assuming the txt file is named report.txt
    
    if not txt_file.exists():
        logger.error(f"Text report not found for task {task_id}")
        return JSONResponse(
            content={"error": "Text report not found"},
            status_code=404
        )

    with open(txt_file, 'rb') as f:
        content = f.read()

    return StreamingResponse(
        io.BytesIO(content),
        media_type="text/plain",
        headers={"Content-Disposition": f"attachment; filename=plip_report_{task_id}.txt"}
    )

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="127.0.0.1", port=8000)
