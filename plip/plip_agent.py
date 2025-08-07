import requests
import time
import os
from pathlib import Path
import logging
from typing import Union, Optional
import asyncio

import json
logger = logging.getLogger(__name__)

class PLIPAgent:
    def __init__(self, api_url: str = "http://localhost:8000"):
        self.api_url = api_url
        self.output_dir = Path("plip_reports")
        self.output_dir.mkdir(exist_ok=True)

    async def analyze_structure(
        self, 
        input_structure: Union[str, Path], 
        pdb_id: Optional[str] = None,
        wait_for_result: bool = True,
        max_wait_time: int = 300  # 5 minutes timeout
    ) -> dict:
        """
        Analyze a protein structure using PLIP and save the report
        
        Args:
            input_structure: Path to PDB file or PDB ID
            pdb_id: Optional PDB ID for naming the output file
            wait_for_result: Whether to wait for analysis completion
            max_wait_time: Maximum time to wait for results in seconds
        
        Returns:
            dict: Analysis results
        """
        try:
            #print(f"[AGENT] Starting analysis for structure: {input_structure}")
            
            # Determine if input is a file path or PDB ID
            if isinstance(input_structure, (str, Path)) and os.path.exists(input_structure):
                #print(f"[AGENT] Processing local file")
                with open(input_structure, 'r') as f:
                    pdb_content = f.read()
                payload = {
                    "file_content": pdb_content,
                    "output_format": ["txt"],
                    "outpath": str(self.output_dir)
                }
                files = None
                body = json.dumps(payload)
                if not pdb_id:
                    pdb_id = Path(input_structure).stem
            else:
                #print(f"[AGENT] Processing PDB ID")
                pdb_id = input_structure
                payload = {
                    "pdb_id": pdb_id,
                    "output_format": ["txt"],
                    "outpath": str(self.output_dir)
                }
                files = None
                body = json.dumps(payload)

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

                if status['status'] == 'completed':
                    #print(f"[AGENT] Analysis completed successfully")
                    return status
                elif status['status'] == 'failed':
                    #print(f"[AGENT] Analysis failed with error: {status.get('error')}")
                    raise Exception(f"Analysis failed: {status.get('error')}")

                #print(f"[AGENT] Waiting 2 seconds before next check")
                time.sleep(2)


            raise TimeoutError(f"Analysis timed out after {max_wait_time} seconds")

        except Exception as e:

            logger.error(f"Error analyzing structure {pdb_id}: {str(e)}")
            raise



async def main():
    # Initialize the agent
    agent = PLIPAgent()

    # Analyze a local PDB file
    result = await agent.analyze_structure(
        input_structure="/Users/atabeyunlu/plip/4gv1.pdb",
        wait_for_result=True
    )
    #print(f"Analysis completed: {result}")

if __name__ == "__main__":
    asyncio.run(main())