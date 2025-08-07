#!/usr/bin/env python
"""
Test script to verify AutoDock API is working locally.
This script tests the FastAPI endpoints directly.
"""
import os
import sys
import json
import requests
import time
from pathlib import Path

# Default settings
API_URL = "http://localhost:8686"
TEST_PROTEIN = Path("1ac8.pdb")  # Use relative path as test runs from /app directory
TEST_LIGAND = Path("ligand.sdf")  # Use relative path as test runs from /app directory

# ANSI colors for terminal output
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
BLUE = "\033[94m"
RESET = "\033[0m"

def print_success(msg):
    print(f"{GREEN}✓ {msg}{RESET}")

def print_error(msg):
    print(f"{RED}✗ {msg}{RESET}")

def print_warning(msg):
    print(f"{YELLOW}! {msg}{RESET}")

def print_info(msg):
    print(f"{BLUE}i {msg}{RESET}")

def download_pdb_if_needed(pdb_id, output_path):
    """Download a PDB file from RCSB PDB bank if it doesn't exist."""
    if output_path.exists():
        print_info(f"PDB file {output_path.name} already exists")
        return True
        
    pdb_id = pdb_id.lower()
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    print_info(f"Downloading {pdb_id.upper()} from RCSB PDB...")
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        with open(output_path, 'w') as f:
            f.write(response.text)
        
        print_success(f"Downloaded {pdb_id.upper()} to {output_path}")
        return True
        
    except Exception as e:
        print_error(f"Failed to download {pdb_id.upper()}: {e}")
        return False

def test_ping():
    """Test the /ping endpoint"""
    try:
        response = requests.get(f"{API_URL}/ping")
        if response.status_code == 200 and response.json()["message"] == "pong":
            print_success("Ping successful")
            return True
        else:
            print_error(f"Ping failed with status code {response.status_code}")
            return False
    except Exception as e:
        print_error(f"Ping failed with exception: {str(e)}")
        return False

def test_inference(protein_path, ligand_path):
    """Test the /inference endpoint with real files"""
    try:
        # Verify files exist
        if not protein_path.exists():
            print_error(f"Protein file not found: {protein_path}")
            return False
        if not ligand_path.exists():
            print_error(f"Ligand file not found: {ligand_path}")
            return False
            
        print_info(f"Testing inference with protein: {protein_path.name}, ligand: {ligand_path.name}")
        
        # Create config with very small box to test buffer overflow theory
        config = {
            "center_x": 0.0,
            "center_y": 0.0,
            "center_z": 0.0,
            "size_x": 10.0,          # Very small box to reduce grid size
            "size_y": 10.0,
            "size_z": 10.0,
            "num_runs": 2,           # Reduced from 5 to make testing faster
            "timeout": 120,
            "max_evaluations": 100000,  # Much reduced for testing
            "heuristics": 1,
            "exhaust": 4             # Lower exhaustiveness for testing
        }
        
        # Upload files to the inference endpoint
        with open(protein_path, "rb") as p_file, open(ligand_path, "rb") as l_file:
            files = {
                "protein_file": (protein_path.name, p_file),
                "ligand_files_0": (ligand_path.name, l_file)
            }
            data = {"config": json.dumps(config)}
            
            print_info("Submitting docking job...")
            response = requests.post(
                f"{API_URL}/inference",
                files=files,
                data=data
            )
            
            if response.status_code == 202:
                task_id = response.json()["task_id"]
                print_success(f"Job submitted successfully. Task ID: {task_id}")
                
                # Poll task status
                print_info("Polling for task completion...")
                max_polls = 180  # Maximum number of times to check (30 minutes at 10 second intervals)
                for i in range(max_polls):
                    time.sleep(10)  # Wait 10 seconds between polls
                    try:
                        status_response = requests.get(f"{API_URL}/task_status/{task_id}")
                        
                        if status_response.status_code == 200:
                            status_data = status_response.json()
                            status = status_data["status"]
                            print_info(f"Current status: {status} (poll {i+1}/{max_polls})")
                            
                            # Print more details if available
                            if "output_dir" in status_data:
                                print_info(f"Output directory: {status_data['output_dir']}")
                            
                            if status == "completed":
                                print_success("Task completed successfully!")
                                
                                # Try to download results
                                download_response = requests.get(f"{API_URL}/download_result/{task_id}")
                                if download_response.status_code == 200:
                                    # Save the zip file
                                    result_file = Path(f"autodock_result_{task_id}.zip")
                                    with open(result_file, "wb") as f:
                                        f.write(download_response.content)
                                    print_success(f"Results downloaded to {result_file}")
                                    return True
                                else:
                                    print_error(f"Failed to download results: {download_response.status_code}")
                                    print_error(f"Response: {download_response.text}")
                                    return False
                            
                            elif status == "failed":
                                error = status_data.get("error", "Unknown error")
                                print_error(f"Task failed: {error}")
                                return False
                        else:
                            print_warning(f"Failed to get status: {status_response.status_code}")
                            print_warning(f"Response: {status_response.text}")
                    except Exception as e:
                        print_warning(f"Error during polling: {str(e)}")
                
                print_error("Task did not complete within the timeout period (30 minutes)")
                return False
            else:
                print_error(f"Job submission failed with status code {response.status_code}")
                print_error(f"Response: {response.text}")
                return False
    except Exception as e:
        print_error(f"Inference test failed with exception: {str(e)}")
        return False

def main():
    """Main test function"""
    print_info("Starting AutoDock API tests")
    
    if not test_ping():
        print_error("Ping test failed, aborting further tests")
        return 1
    
    # Check for environment-specific test files
    protein_path = TEST_PROTEIN
    ligand_path = TEST_LIGAND
    
    # Download PDB file if needed
    if not download_pdb_if_needed("1ac8", protein_path):
        print_error("Failed to obtain protein file")
        return 1
    
    # Allow command-line overrides
    if len(sys.argv) > 1:
        protein_path = Path(sys.argv[1])
    if len(sys.argv) > 2:
        ligand_path = Path(sys.argv[2])
    
    if not test_inference(protein_path, ligand_path):
        print_error("Inference test failed")
        return 1
    
    print_success("All tests passed!")
    return 0

if __name__ == "__main__":
    sys.exit(main())