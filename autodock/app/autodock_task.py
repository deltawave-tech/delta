import asyncio
import logging
from pathlib import Path
from typing import List, Union

# Import project modules
from models import AutoDockConfig
from inference import run_autodock_docking, run_batch_autodock_docking, run_sequential_autodock_docking

logger = logging.getLogger(__name__)
tasks = {}

async def process_and_run_autodock(
    task_id: str,
    protein_content: bytes,
    ligand_content: Union[bytes, List[bytes]],
    config: AutoDockConfig,
    output_dir: Path
):
    try:
        tasks[task_id]['status'] = 'running'
        
        # Determine if this is a batch job based on type of ligand_content
        is_batch = isinstance(ligand_content, list)
        if is_batch:
            logger.info(f"Starting processing for task {task_id} with {len(ligand_content)} ligands")
            tasks[task_id]['num_ligands'] = len(ligand_content)
        else:
            logger.info(f"Starting single ligand AutoDock docking for task {task_id}")
            tasks[task_id]['num_ligands'] = 1

        loop = asyncio.get_event_loop()
        
        # Call the appropriate function based on batch status
        if is_batch:
            # Using batch mode for multiple ligands (native AutoDock-GPU batch functionality)
            logger.info("Using batch mode for multiple ligands (native AutoDock-GPU batch functionality)")
            result_dir = await loop.run_in_executor(
                None,
                run_batch_autodock_docking,
                protein_content,
                ligand_content,
                config,
                output_dir
            )
            
            # Sequential processing alternative - commented out but preserved for fallback
            # logger.info("Using sequential processing for multiple ligands")
            # result_dir = await loop.run_in_executor(
            #     None,
            #     run_sequential_autodock_docking,
            #     protein_content,
            #     ligand_content,
            #     config,
            #     output_dir
            # )
        else:
            result_dir = await loop.run_in_executor(
                None,
                run_autodock_docking,
                protein_content,
                ligand_content,
                config,
                output_dir
            )

        if result_dir:
            tasks[task_id]['status'] = 'completed'
            tasks[task_id]['output_dir'] = str(result_dir)
            logger.info(f"AutoDock docking completed for task {task_id}")
        else:
            tasks[task_id]['status'] = 'failed'
            tasks[task_id]['error'] = "AutoDock docking failed"
            logger.error("AutoDock docking failed")

    except Exception as e:
        logger.exception(f'An exception occurred during AutoDock processing: {str(e)}')
        tasks[task_id]['status'] = 'failed'
        tasks[task_id]['error'] = str(e)