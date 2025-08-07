import asyncio
import os
import uuid
import logging
from pathlib import Path
from typing import Dict

from plip_inference import PLIPInference

# Configure logging
logger = logging.getLogger(__name__)

# Task storage
tasks: Dict[str, str] = {}  # task_id -> status string

async def process_task(request_data: dict) -> str:
    """Create and manage task for PLIP analysis"""
    task_id = str(uuid.uuid4())
    tasks[task_id] = "queued"

    # Start inference in background
    asyncio.create_task(
        run_inference_task(task_id, request_data)
    )

    return task_id

async def run_inference_task(task_id: str, request_data: dict):
    """Run inference in background"""
    try:
        tasks[task_id] = "running"

        # Prepare input
        input_file = await prepare_input(task_id, request_data)
        logger.info(f"Prepared input file: {input_file}")

        # Run inference
        inference = PLIPInference(task_id)
        await inference.run(
            pdb_file=input_file,
            output_format=request_data.get('output_format', ['xml', 'txt'])
        )

        tasks[task_id] = "completed"
        logger.info(f"Task {task_id} completed successfully")

    except Exception as e:
        logger.error(f"Task {task_id} failed: {str(e)}")
        tasks[task_id] = "failed"

async def prepare_input(task_id: str, request_data: dict) -> str:
    """Prepare input file for inference"""
    if request_data.get('file_content'):
        input_dir = Path(f"storage/{task_id}")
        input_dir.mkdir(parents=True, exist_ok=True)
        input_file = input_dir / "input.pdb"
        content = request_data['file_content']
        with open(input_file, 'w') as f:  # Changed to text mode since we decoded in API
            f.write(content)
        logger.info(f"Saved input file to: {input_file}")
        return str(input_file)
    else:
        return f"pdb:{request_data['pdb_id']}"

async def get_task_status(task_id: str) -> str:
    """Get current status of a task"""
    return tasks.get(task_id, "not_found")

def list_tasks() -> list:
    """Get list of all task IDs"""
    return list(tasks.keys())
