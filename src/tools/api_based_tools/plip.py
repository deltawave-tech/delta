import requests
import time
import os
from pathlib import Path
import logging
from typing import Union, Optional, List, Dict
import asyncio
import json
from src.tools.api_based_tools.plip_rep import summarize_plip_output
from datetime import datetime
logger = logging.getLogger(__name__)
class PLIP:
    def __init__(self, api_url: str = "http://127.0.0.1:8000"):
        self.api_url = api_url
        self.output_dir = Path("plip_reports")
        self.output_dir.mkdir(exist_ok=True)

    def analyze_structure(
        self, 
        input_structure: Union[str, Path, List[Union[str, Path]]], 
        pdb_id: Optional[Union[str, List[str]]] = None,
        wait_for_result: bool = True,
        max_wait_time: int = 300,  # 5 minutes timeout
        output_dir: Optional[Path] = None
    ) -> Union[dict, Dict[str, dict]]:
        """
        Analyze one or multiple protein structures using PLIP.
        
        Args:
            input_structure: Path to a structure file or list of paths
            pdb_id: Optional PDB ID or list of PDB IDs (must match length of input_structure if both are lists)
            wait_for_result: Whether to wait for analysis to complete
            max_wait_time: Maximum time to wait for results in seconds
            output_dir: Optional directory to save results to
            
        Returns:
            A single result dict if input was a single structure, or
            a dictionary mapping structure paths/IDs to their result dicts
        """
        # Handle single structure case for backward compatibility
        if not isinstance(input_structure, list):
            # Original function behavior for a single structure
            return self._analyze_single_structure(
                input_structure=input_structure,
                pdb_id=pdb_id,
                wait_for_result=wait_for_result,
                max_wait_time=max_wait_time,
                output_dir=output_dir
            )
        
        # Handle multiple structures
        results = {}
        
        # Convert pdb_id to list if it's not None and not already a list
        pdb_ids = pdb_id if isinstance(pdb_id, list) or pdb_id is None else [pdb_id]
        
        # Validate inputs if both input_structure and pdb_ids are lists
        if pdb_ids is not None and len(pdb_ids) != len(input_structure):
            raise ValueError("If both input_structure and pdb_id are lists, they must have the same length")
        
        # Process each structure

        for i, structure in enumerate(input_structure):
            structure_id = pdb_ids[i] if pdb_ids is not None else None
            try:
                result = self._analyze_single_structure(
                    input_structure=structure,
                    pdb_id=structure_id,
                    wait_for_result=wait_for_result,
                    max_wait_time=max_wait_time,
                    output_dir=output_dir
                )
                
                # Use the structure path or pdb_id as the key
                key = str(structure) if structure_id is None else structure_id
                results[key] = result
                
            except Exception as e:
                logger.error(f"Error analyzing structure {structure}: {str(e)}")
                # Store the error in the results
                key = str(structure) if structure_id is None else structure_id
                results[key] = {"status": "failed", "error": str(e)}
        
        return results

    def _analyze_single_structure(
        self, 
        input_structure: Union[str, Path], 
        pdb_id: Optional[str] = None,
        wait_for_result: bool = True,
        max_wait_time: int = 300,
        output_dir: Optional[Path] = None
    ) -> dict:
        """Implementation of the original analyze_structure method for a single structure"""
        try:
            #print(f"[AGENT] Starting analysis for structure: {input_structure}")

            # Process input as a file path
            with open(input_structure, 'r') as f:
                pdb_content = f.read()
            payload = {
                "file_content": pdb_content,
                "output_format": ["txt"],
                "outpath": str(self.output_dir)
            }
            files = None
            body = json.dumps(payload)
            
            # Use provided pdb_id or extract from filename
            if not pdb_id:
                pdb_id = Path(input_structure).stem
            
            # Submit analysis request
            #print(f"[AGENT] Submitting analysis request to {self.api_url}/inference")
            response = requests.post(
                f"{self.api_url}/inference",
                data={"body": body},
                files=files
            )
            #print(f"[AGENT] Response status code: {response.status_code}")
            
            if response.status_code != 202:
                raise Exception(f"Failed to submit analysis: {response.text}")

            task_id = response.json()['task_id']
            #print(f"[AGENT] Analysis submitted for {pdb_id} (Task ID: {task_id})")

            if not wait_for_result:
                return {"task_id": task_id}

            # Wait for results
            #print(f"[AGENT] Waiting for results (timeout: {max_wait_time}s)")
            start_time = time.time()
            check_count = 0
            
            while time.time() - start_time < max_wait_time:
                check_count += 1
                #print(f"[AGENT] Checking status (attempt {check_count})")

                status_response = requests.get(f"{self.api_url}/task_status/{task_id}")
                status = status_response.json()
                #print(f"[AGENT] Current status: {status['status']}")
   
                if status == 'completed':
                    # Download text results when analysis is complete
                    download_response = requests.get(f"{self.api_url}/download/{task_id}")
                    
                    if download_response.status_code == 200:
                        # Create plip_results directory inside the output_dir
                        plip_results_dir = output_dir / "plip_results"
                        plip_results_dir.mkdir(exist_ok=True, parents=True)
                        
                        # Save the text results in the plip_results directory
                        text_path = plip_results_dir / f"{pdb_id}_{self.generate_unique_id()}_plip.txt"
                        csv_path = plip_results_dir / f"{pdb_id}_{self.generate_unique_id()}_plip.csv"
                        with open(text_path, 'w') as f:
                            f.write(download_response.text)
              
                        summarize_plip_output(str(text_path), output_file=csv_path)

                        with open(csv_path, 'r') as f:
                            plip_output = f.read()
                        return {
                            "status": "completed", 
                            "task_id": task_id, 
                            "text_path": str(csv_path),
                            "text": plip_output
                        }
                    else:
                        raise Exception(f"Failed to download results: {plip_output}")
                elif status == 'failed':
                    #print(f"[AGENT] Analysis failed with error: {status.get('error')}")
                    raise Exception(f"Analysis failed: {status.get('error')}")

                #print(f"[AGENT] Waiting 2 seconds before next check")
                time.sleep(2)


            raise TimeoutError(f"Analysis timed out after {max_wait_time} seconds")

        except Exception as e:

            logger.error(f"Error analyzing structure {pdb_id}: {str(e)}")
            raise
    
    def generate_unique_id(self) -> str:
        #return hour minute second
        return datetime.now().strftime("%H%M%S")

from ..tool_definitions import tool_registry

@tool_registry.register("get_plip_report")
def get_plip_report(arguments: dict, output_dir_id: str) -> dict:
    print("PLIP is being used")
    # Initialize PLIP
    plip = PLIP()

    # Get output directory path
    output_dir = tool_registry.output_dir/output_dir_id
    # Extract parameters from arguments
    input_structure = arguments.get("input_structure")
    pdb_id = arguments.get("pdb_id", None)
    wait_for_result = arguments.get("wait_for_result", True)
    max_wait_time = arguments.get("max_wait_time", 300)

    # Run PLIP analysis
    results = plip.analyze_structure(
        input_structure=input_structure,
        pdb_id=pdb_id,
        wait_for_result=wait_for_result,
        max_wait_time=max_wait_time,
        output_dir=output_dir
    )

    return {"results": results}