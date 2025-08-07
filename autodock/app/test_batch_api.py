#!/usr/bin/env python3
# test_batch_api.py - Test client for the autodock API

import requests
import time
import json
import sys
import os
from pathlib import Path

# Default server URL
SERVER_URL = "http://localhost:8686"

def test_autodock_batch_api(protein_file, ligand_files, config=None):
    """Test the autodock API with batch processing"""
    print(f"Testing autodock batch API with {len(ligand_files)} ligands")
    print(f"Server URL: {SERVER_URL}")
    
    # Test ping endpoint
    print("\nTesting ping endpoint...")
    response = requests.get(f"{SERVER_URL}/ping")
    print(f"Response: {response.status_code}")
    if response.status_code == 200:
        print("Ping successful")
    else:
        print(f"Ping failed: {response.text}")
        return
    
    # Test inference endpoint with batch processing
    print("\nTesting inference endpoint with batch processing...")
    
    # Prepare files for upload using the new indexed format
    files = {
        'protein_file': open(protein_file, 'rb'),
    }

    # Add ligand files with proper indexed fields
    for i, ligand_file in enumerate(ligand_files):
        files[f'ligand_files_{i}'] = open(ligand_file, 'rb')
    
    # Prepare form data
    form_data = {}
    if config:
        form_data['config'] = json.dumps(config)
    
    # Print the request structure for debugging
    print("\nSending request with files:")
    for key in files:
        print(f"  {key}: {files[key].name}")
    
    # Send request
    response = requests.post(
        f"{SERVER_URL}/inference", 
        files=files,
        data=form_data
    )
    
    print(f"Response: {response.status_code}")
    if response.status_code == 202:
        print("Inference request accepted")
        result = response.json()
        task_id = result.get('task_id')
        print(f"Task ID: {task_id}")
        
        # Close all file handles
        for f in files.values():
            f.close()
        
        # Poll for results
        print("\nPolling for task status...")
        max_polls = 60  # Maximum number of times to poll
        poll_interval = 5  # Seconds between polls
        
        for i in range(max_polls):
            response = requests.get(f"{SERVER_URL}/task_status/{task_id}")
            
            if response.status_code == 200:
                status = response.json()
                print(f"Poll {i+1}: Status = {status.get('status')}")
                
                if status.get('status') == 'completed':
                    print("Task completed successfully!")
                    break
                elif status.get('status') == 'failed':
                    print(f"Task failed: {status.get('error', 'Unknown error')}")
                    break
            else:
                print(f"Error checking status: {response.status_code} - {response.text}")
                break
            
            time.sleep(poll_interval)
        
        # Download results if task completed
        if status.get('status') == 'completed':
            print("\nDownloading results...")
            response = requests.get(f"{SERVER_URL}/download_result/{task_id}")
            
            if response.status_code == 200:
                # Save the zip file
                output_file = f"autodock_batch_results_{task_id}.zip"
                with open(output_file, 'wb') as f:
                    f.write(response.content)
                print(f"Results downloaded to {output_file}")
            else:
                print(f"Error downloading results: {response.status_code} - {response.text}")
    else:
        print(f"Inference request failed: {response.text}")
        # Close all file handles
        for f in files.values():
            f.close()

if __name__ == "__main__":
    # Check for command line arguments
    if len(sys.argv) < 3:
        print("Usage: python test_batch_api.py [protein_file] [ligand_file1] [ligand_file2] ...")
        sys.exit(1)
    
    protein_file = sys.argv[1]
    ligand_files = sys.argv[2:]
    
    # Example configuration
    config = {
        "center_x": 0.0,
        "center_y": 0.0,
        "center_z": 0.0,
        "size_x": 20.0,
        "size_y": 20.0,
        "size_z": 20.0,
        "num_modes": 3,  # Reduced for testing
        "exhaustiveness": 2,  # Reduced for testing
        "num_runs": 3,  # Reduced for testing 
        "timeout": 60  # Short timeout for testing
    }
    
    test_autodock_batch_api(protein_file, ligand_files, config)