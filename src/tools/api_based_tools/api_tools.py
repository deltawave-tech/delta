import os
import time
import json
import logging
from pathlib import Path
from typing import List, Optional, Dict, Any, Union, Tuple
import requests
from dataclasses import dataclass
import zipfile
from pydantic import BaseModel, Field   
import numpy as np
logger = logging.getLogger(__name__)
from admet_ai import ADMETModel
from rdkit import Chem
import tempfile

def save_zipfile(result_path: Path, output_dir: Path, protein_name: str) -> None:
    output_dir = Path(output_dir) / protein_name
    print(f"Opening zip file from: {result_path}")
    print(f"Output directory set to: {output_dir}")
    
    with zipfile.ZipFile(result_path, 'r') as zip_ref:
        print(f"Files in zip: {zip_ref.namelist()}")
        for file_info in zip_ref.filelist:
            print(f"Processing file: {file_info.filename}")
            parts = file_info.filename.split('/')
            print(f"Split parts: {parts}")
            
            if len(parts) >= 3 and parts[-1].endswith('.sdf'):
                task_id = parts[0]
                rank_file = parts[-1]
                print(f"Task ID: {task_id}")
                print(f"Rank file: {rank_file}")
                
                # Replace forward slashes with a special token that won't be interpreted as a path separator
                smiles_parts = parts[1:-1]
                smiles = '[SLASH]'.join(smiles_parts)
                print(f"Processed SMILES: {smiles}")
                
                # Create path with encoded SMILES
                new_path = f"{task_id}/{smiles}/{rank_file}"
                print(f"New path: {new_path}")
                    
                target = Path(output_dir) / new_path
                print(f"Target path: {target}")
                target.parent.mkdir(parents=True, exist_ok=True)
                
                try:
                    with zip_ref.open(file_info.filename) as source, open(target, 'wb') as target_file:
                        content = source.read()
                        print(f"Read {len(content)} bytes from source file")
                        target_file.write(content)
                        print(f"Successfully wrote file to: {target}")
                except Exception as e:
                    print(f"Error writing file {target}: {str(e)}")

    print(f"Finished processing zip file. Output directory: {output_dir}")
    return output_dir


def get_original_smiles(path: str) -> str:
    """Convert back to original SMILES with proper slashes"""
    return path.replace('[SLASH]', '/')




def save_zipfile(result_path: Path, output_dir: Path, protein_name: str) -> None:
    output_dir = Path(output_dir) / protein_name
    with zipfile.ZipFile(result_path, 'r') as zip_ref:
        for file_info in zip_ref.filelist:
            parts = file_info.filename.split('/')
            if len(parts) >= 3 and parts[-1].endswith('.sdf'):
                task_id = parts[0]
                rank_file = parts[-1]
                
                # Replace forward slashes with a special token that won't be interpreted as a path separator
                smiles_parts = parts[1:-1]
                smiles = '[SLASH]'.join(smiles_parts)
                
                # Create path with encoded SMILES
                new_path = f"{task_id}/{smiles}/{rank_file}"
                    
                target = Path(output_dir) / new_path
                target.parent.mkdir(parents=True, exist_ok=True)
                with zip_ref.open(file_info.filename) as source, open(target, 'wb') as target_file:
                    target_file.write(source.read())

    return output_dir


def get_original_smiles(path: str) -> str:
    """Convert back to original SMILES with proper slashes"""
    return path.replace('[SLASH]', '/')

@dataclass
class DiffDockResult:
    """Structure for DiffDock results."""
    task_id: str
    output_dir: Path
    scores_path: Path
    pdb_files: List[Path]



class DockingResult(BaseModel):
    """Stores and analyzes results from a docking run"""
    protein_name: str
    output_dir: Path
    results_by_smiles: Dict[str, List[Path]] = Field(description="Dictionary of SMILES to their ranked files")

    class Config:
        arbitrary_types_allowed = True

    def _get_score_from_file(self, file_path: Path) -> float:
        """Extract confidence score from filename"""
        filename = file_path.stem
        confidence_part = filename.split('confidence')[1]
        # end = start.split('_')[0]
        return float(confidence_part)

    def get_smiles_scores(self, smiles: str) -> List[float]:
        """Get all scores for a specific SMILES"""
        return [self._get_score_from_file(f) for f in self.results_by_smiles[smiles]]

    @property
    def all_smiles(self) -> List[str]:
        """Get list of all SMILES"""
        return list(self.results_by_smiles.keys())

    def to_json(self) -> str:
        """Convert to JSON string with hierarchical structure"""
        result_dict = {
            "protein_name": str(self.protein_name),
            "results": {
                smiles: [str(f.name) for f in files]
                for smiles, files in self.results_by_smiles.items()
            }
        }
        return json.dumps(result_dict, indent=2)

    def to_best_values(self) -> str:
        """Get the best value from the docking results"""
        return self.get_smiles_scores(self.all_smiles[0])

    def to_markdown(self) -> str:
        """Create a markdown table with SMILES rows and rank columns"""
        # Header
        md = f"# Docking Results for {str(self.protein_name)}\n\n"
        
        # Fixed width for each column
        smiles_width = 60  # Adjust this value based on your longest SMILES
        score_width = 8    # Width for score columns
        
        # Create aligned header
        header = (f"| {'SMILES':<{smiles_width}} |" + 
                "".join(f" {'Rank '+str(i+1):<{score_width}} |" for i in range(10)))
        
        # Create separator with matching widths
        separator = (f"|{'-'*smiles_width}-|" + 
                    "".join(f"{'-'*score_width}-|" for i in range(10)))
        
        # Create rows with aligned columns
        rows = []
        for smiles in self.all_smiles:
            scores = self.get_smiles_scores(smiles)
            score_cells = [f"{scores[i]:.2f}" if i < len(scores) else "N/A" 
                        for i in range(10)]
            row = (f"| {smiles:<{smiles_width}} |" + 
                "".join(f" {score:<{score_width}} |" for score in score_cells))
            rows.append(row)
        
        # Combine all parts
        md += "\n".join([header, separator] + rows)
        return md
class DiffDockArgs(BaseModel):
    protein_file: str
    smiles: List[str]
    output_dir: str

class DiffDockAPI:
    """Client for DiffDock API."""

    def __init__(self, base_url: str = "http://34.170.69.34:8000") -> None:
        """Initialize DiffDock API client.

        :param base_url: Base URL for DiffDock API (default: production server)
        """
        self.base_url = base_url.rstrip('/')

    def ping(self) -> bool:
        """Check if API is alive.
        :return: True if API responds with pong
        """
        response = requests.get(f"{self.base_url}/ping")
        response.raise_for_status()
        return response.json().get('message') == 'pong'

    def parse_args(self, args: Dict[str, Any]) -> DiffDockArgs:
        return DiffDockArgs(**args)

    def get_task_status(self, task_id: str) -> dict:
        """Get status of a task.

        :param task_id: Task ID to check
        :return: Task status information
        """
        url = f"{self.base_url}/task_status/{task_id}"
        response = requests.get(url)
        response.raise_for_status()
        return response.json()

    def predict(
        self,
        protein_file: Union[str, Path],
        smiles: List[str],
        output_dir: Union[str, Path],
        timeout: int = 3600,  # 1 hour timeout
        polling_interval: int = 30,  # Check every 30 seconds
        to_markdown: bool = True,
        merge_pdb_sdf: bool = True
    ) -> DiffDockResult:
        """Run DiffDock prediction.

        :param protein_file: Path to protein PDB file
        :param smiles: List of SMILES strings for ligands
        :param output_dir: Directory to save results
        :param timeout: Maximum time to wait for results in seconds
        :param polling_interval: Time between download attempts in seconds
        :return: DiffDock results
        """
        # Ensure output directory exists
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Submit inference job
        print(f"Submitting inference for protein: {protein_file} with {len(smiles)} ligands")
        protein_file = clean_pdb_file(protein_file)
        task_id = self._submit_inference(protein_file, smiles)
        logger.info(f"Submitted DiffDock job: {task_id}")
        protein_name = Path(protein_file).name
        print(f"Processing protein: {protein_name}")

        # Poll for results
        start_time = time.time()
        while True:
            # Check task status
            status = self.get_task_status(task_id)
            logger.info(f"Task status: {status}")

            if status.get('status') == 'completed':
                try:
                    print(f"Task completed. Downloading results for task {task_id}")
                    result_path = self._download_results(task_id, output_dir)
                    print(f"Download completed. Result path: {result_path}")
                    if not result_path.exists():
                        print(f"Warning: Result path does not exist: {result_path}")
                        return None
                    
                    print(f"Saving results to {output_dir}")
                    saved_dir = save_zipfile(result_path, output_dir, protein_name)
                    print(f"Files saved to directory: {saved_dir}")
                    
                    print("Processing downloaded results...")
                    result = self._process_results(task_id, output_dir, protein_name)
                    if result is None:
                        print("Warning: _process_results returned None")
                        return None
                    if merge_pdb_sdf:
                        merged_pdb_paths = self.merge_pdb_sdf_rdkit(task_id, protein_file, saved_dir, output_dir)
                    if to_markdown:
                        print("Converting results to markdown")
                        result_md = result.to_markdown()
                        if merged_pdb_paths:
                            result_md += f"\n\n# Merged PDB files: {', '.join([str(p) for p in merged_pdb_paths])}"
                        return result_md
                    return result
                    
                except Exception as e:
                    print(f"Error processing results: {str(e)}")
                    import traceback
                    print(f"Traceback: {traceback.format_exc()}")
                    raise
            elif status.get('status') == 'failed':
                print(f"Task failed with status: {status}")
                raise RuntimeError(f"Task failed: {status.get('error', 'Unknown error')}")

            if time.time() - start_time > timeout:
                raise TimeoutError(
                    f"DiffDock job {task_id} timed out after {timeout} seconds. "
                    f"Last status: {status}"
                )

            time.sleep(polling_interval)

    def _submit_inference(self, protein_file: Union[str, Path], smiles: List[str]) -> str:
        """Submit inference job to DiffDock API.

        :param protein_file: Path to protein PDB file
        :param smiles: List of SMILES strings
        :return: Task ID
        """
        url = f"{self.base_url}/inference"
        # protein_name = Path(protein_file).name
        with open(protein_file, 'rb') as f:
            files = {
                'file': (Path(protein_file).name, f, 'application/octet-stream')
            }

            body = json.dumps({
                "smiles": smiles
            })

            data = {
                'body': body
            }

            response = requests.post(url, files=files, data=data)
            response.raise_for_status()

            result = response.json()
            return result['task_id']

    def _download_results(self, task_id: str, output_dir: Path) -> Optional[Path]:
        """Download results.
        
        :param task_id: Task ID
        :param output_dir: Output directory
        :return: Path to results if successful, None if not ready
        """
        url = f"{self.base_url}/download_result/{task_id}"
        print(f"[DEBUG-DOWNLOAD] Downloading results from: {url}")
        
        # Create a new session for downloading to avoid connection pool issues
        with requests.Session() as session:
            # Implement retry logic with exponential backoff
            max_retries = 3
            retry_delay = 5
            
            for attempt in range(max_retries):
                try:
                    # Stream the response with increased timeout
                    with session.get(url, stream=True, timeout=(30, 1800)) as response:
                        print(f"[DEBUG-DOWNLOAD] Response status code: {response.status_code}")
                        
                        # Check response status before raising for status
                        if response.status_code != 200:
                            print(f"[DEBUG-DOWNLOAD] Non-200 response: {response.text[:500]}")
                        
                        response.raise_for_status()
                        
                        # Save zip file
                        zip_path = output_dir / f"{task_id}_vina_output.zip"
                        print(f"[DEBUG-DOWNLOAD] Saving to: {zip_path}")
                        
                        total_size = 0
                        with open(zip_path, 'wb') as f:
                            for i, chunk in enumerate(response.iter_content(chunk_size=8192)):
                                if chunk:  # filter out keep-alive chunks
                                    chunk_size = len(chunk)
                                    total_size += chunk_size
                                    f.write(chunk)
                                    if i % 100 == 0:  # Log every ~800KB
                                        print(f"[DEBUG-DOWNLOAD] Downloaded {total_size/1024:.1f}KB so far...")
                        
                        print(f"[DEBUG-DOWNLOAD] Download complete, total size: {total_size/1024:.1f}KB")
                        
                        # Verify the file exists and has content
                        if zip_path.exists():
                            file_size = zip_path.stat().st_size
                            print(f"[DEBUG-DOWNLOAD] Zip file exists with size: {file_size/1024:.1f}KB")
                            if file_size == 0:
                                print(f"[DEBUG-DOWNLOAD] WARNING: Downloaded file is empty!")
                                if attempt < max_retries - 1:
                                    # If file is empty, try again
                                    raise ValueError("Downloaded zip file is empty")
                            else:
                                return zip_path
                        else:
                            print(f"[DEBUG-DOWNLOAD] ERROR: Zip file does not exist after download!")
                            raise FileNotFoundError("Zip file not created during download")
                
                except (requests.exceptions.Timeout, requests.exceptions.ConnectionError, 
                        ValueError, FileNotFoundError) as e:
                    if attempt < max_retries - 1:
                        wait_time = retry_delay * (2 ** attempt)  # Exponential backoff
                        print(f"[DEBUG-DOWNLOAD] Attempt {attempt+1} failed with {type(e).__name__}: {str(e)}")
                        print(f"[DEBUG-DOWNLOAD] Retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                    else:
                        print(f"[DEBUG-DOWNLOAD] All {max_retries} attempts failed. Last error: {str(e)}")
                        raise
                except Exception as e:
                    print(f"[DEBUG-DOWNLOAD] Error downloading results: {str(e)}")
                    import traceback
                    print(f"[DEBUG-DOWNLOAD] Traceback: {traceback.format_exc()}")
                    raise
            
            # This should not be reached if max_retries was exhausted
            raise RuntimeError("[DEBUG-DOWNLOAD] Failed to download results after exhausting all retries")
    
    def _process_results(self, task_id: str, 
                         output_dir: Path, 
                         protein_name: str,
                         top_k: int = 3) -> DockingResult:
        """Process the docking results from a directory"""
        task_dir = Path(output_dir) / protein_name / task_id
        print(f"Processing results from directory: {task_dir}")
        
        # Group files by SMILES
        results_by_smiles = {}
        
        # List all files first to debug
        all_files = list(task_dir.glob("**/rank*_confidence-*.sdf"))
        print(f"Found {len(all_files)} files matching pattern")
        
        for path in all_files:
            # Get SMILES from parent directory name
            smiles = path.parent.name
            # Replace any encoding used during saving (e.g., [SLASH] -> /)
            smiles = smiles.replace('[SLASH]', '/')
            print(f"Processing file: {path} with SMILES: {smiles}")
            
            if smiles not in results_by_smiles:
                results_by_smiles[smiles] = []
            results_by_smiles[smiles].append(path)
        
        print(f"Found results for {len(results_by_smiles)} SMILES strings")
        
        # Sort files by score for each SMILES - Fixed the sorting/slicing issue
        for smiles in results_by_smiles:
            paths = results_by_smiles[smiles]
            print(f"Sorting {len(paths)} results for SMILES: {smiles}")
            # Sort the paths list
            paths.sort(key=lambda p: float(p.stem.split('-')[1]))
            # Take top k results
            results_by_smiles[smiles] = paths[:top_k]
            print(f"Kept top {len(results_by_smiles[smiles])} results")
        
        return DockingResult(
            protein_name=protein_name, 
            output_dir=task_dir,
            results_by_smiles=results_by_smiles
        )
    def merge_pdb_sdf_rdkit(self, task_id: str, pdb_file: Union[str, Path], saved_dir: Path, output_dir: Path) -> List[Path]:
        """Merge PDB and SDF files using RDKit for multiple ligands.
        
        Args:
            task_id: Task ID
            pdb_file: Path to PDB file
            saved_dir: Directory containing SMILES subdirectories with SDF files
            output_dir: Output directory
            
        Returns:
            List of paths to merged PDB files
        """
        from rdkit import Chem
        
        # Use a temporary PDB file with hydrogens added
        temp_pdb = Path(tempfile.mktemp(suffix=".pdb"))
        merged_paths = []
        
        try:
            # Process protein once
            protein = Chem.MolFromPDBFile(str(pdb_file), removeHs=False)
            if not protein:
                # If protein loading fails, try cleaning it first
                cleaned_pdb = clean_pdb_file(pdb_file, temp_pdb)
                protein = Chem.MolFromPDBFile(str(cleaned_pdb), removeHs=False)
            
            if not protein:
                raise ValueError(f"Could not load protein from {pdb_file}")
            
            # Find all SMILES directories under the task_id directory
            task_dir = Path(saved_dir) / task_id
            smiles_dirs = [d for d in task_dir.iterdir() if d.is_dir()]
            
            if not smiles_dirs:
                logger.warning(f"No SMILES directories found in {task_dir}")
                return []
            
            # Process each SMILES directory
            for smiles_dir in smiles_dirs:
                # Get SMILES string from directory name
                smiles = smiles_dir.name
                
                # Find SDF file in the SMILES directory
                sdf_files = list(smiles_dir.glob("rank1_confidence-*.sdf"))
                if not sdf_files:
                    logger.warning(f"No SDF files found in {smiles_dir}")
                    continue
                    
                # Use the first SDF file (typically there's only one)
                sdf_file = sdf_files[0]
                
                # Load ligand
                supplier = Chem.SDMolSupplier(str(sdf_file), removeHs=False)
                if len(supplier) == 0:
                    logger.warning(f"No molecules found in {sdf_file}")
                    continue
                    
                ligand = supplier[0]
                
                # Combine the molecules
                combined = Chem.CombineMols(protein, ligand)
                
                # Prepare output path
                pdb_name = Path(pdb_file).stem
                smiles_short = smiles  # We might try to truncate it later
                plip_dir = Path(output_dir) / task_id
                plip_dir.mkdir(parents=True, exist_ok=True)
                
                output_path = plip_dir / f"{pdb_name}_{smiles_short}_merged.pdb"
                
                # Write output
                writer = Chem.PDBWriter(str(output_path))
                writer.write(combined)
                writer.close()
                
                merged_paths.append(output_path)
            
            return merged_paths
        finally:
            # Clean up temporary file
            if temp_pdb.exists():
                temp_pdb.unlink()

from ..tool_definitions import tool_registry    

@tool_registry.register("run_diffdock")
def run_diffdock(arguments: dict, output_dir_id: str) -> dict:
    print("DiffDock is being used")
    diffdock = DiffDockAPI()

    return diffdock.predict(
        protein_file=arguments["protein_path"],
        smiles=arguments["smiles"],
        output_dir=tool_registry.output_dir/output_dir_id/'docking'
    )

def get_admet_scores(smiles: Union[str, List[str]], check_validity: bool = True) -> Union[float, List[float]]:
    """Get ADMET scores for a SMILES string or list of SMILES strings
    
    Args:
        smiles: Single SMILES string or list of SMILES strings
        check_validity: If True, validates SMILES and removes duplicates
        
    Returns:
        Single score if input is string, list of scores if input is list
        
    Raises:
        ValueError: If any SMILES strings are invalid when check_validity=True
    """
    model = ADMETModel()
    
    if isinstance(smiles, str):
        if check_validity:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {smiles}")
        return model.predict([smiles])[0]
    
    elif isinstance(smiles, list):
        if check_validity:
            # Check validity and uniqueness
            valid_smiles = []
            seen = set()
            invalid = []
            
            for s in smiles:
                if s in seen:
                    continue
                mol = Chem.MolFromSmiles(s)
                if mol is None:
                    invalid.append(s)
                else:
                    valid_smiles.append(s)
                    seen.add(s)
            
            if invalid:
                raise ValueError(f"Invalid SMILES strings found: {invalid}")
            
            return model.predict(valid_smiles)
        return model.predict(smiles)
    

class VinaDockingResult(BaseModel):
    """Stores and analyzes results from a Vina docking run"""
    protein_name: str
    ligand_name: str
    output_dir: Path
    docking_files: Dict[str, Path] = Field(description="Dictionary of file types to their paths")
    
    class Config:
        arbitrary_types_allowed = True
    
    def _get_score_from_log(self) -> List[float]:
        """Extract binding energies from summary.json file"""
        scores = []
        import re
        
        # Get the ligand number from ligand_pdbqt filename - server uses 1-based indexing
        ligand_num = None
        ligand_pdbqt = self.docking_files.get('ligand_pdbqt')
        
        if ligand_pdbqt:
            if 'docking_result_best_pose' in str(ligand_pdbqt):
                # Special case for single molecule - use index 1
                ligand_num = "1"
                # Check for dlg file without index too
                dlg_file = self.output_dir / "docking_result.dlg"
                if dlg_file.exists() and 'docking_result_dlg' not in self.docking_files:
                    self.docking_files['docking_result_dlg'] = dlg_file
            else:
                # Match ligand_X pattern (server uses 1-based, so ligand_1, ligand_2, etc.)
                match = re.search(r'ligand_(\d+)', str(ligand_pdbqt))
                if match:
                    # The server uses 1-based indexing, so we can use this directly
                    ligand_num = match.group(1)
        
        # If we couldn't determine ligand number, we can't extract binding energy
        if not ligand_num:
            return scores
        
        # Look for summary.json in the output directory
        summary_path = self.output_dir / "summary.json"
        if not summary_path.exists():
            return scores
        
        # Load and parse summary.json
        try:
            with open(summary_path, 'r') as f:
                summary = json.load(f)
            
            # Get the specific ligand data
            ligand_key = f"ligand_{ligand_num}"
            if ligand_key in summary.get('ligands', {}):
                ligand_data = summary['ligands'][ligand_key]
                
                # Check if docking was successful
                if ligand_data.get('success', False):
                    # Extract binding energies from poses
                    for pose in ligand_data.get('poses', []):
                        if 'binding_energy' in pose:
                            scores.append(float(pose['binding_energy']))
                            
            # Special case for single molecule - also check the dlg file directly
            if ligand_num == "1" and not scores:
                # Try to find the binding energy directly in the dlg file
                dlg_file = self.docking_files.get('docking_result_dlg')
                if dlg_file and dlg_file.exists():
                    try:
                        with open(dlg_file, 'r') as f:
                            dlg_content = f.read()
                            # Look for binding energy patterns in the dlg file
                            energy_matches = re.findall(r'Estimated Free Energy of Binding\s*=\s*(-?\d+\.\d+)', dlg_content)
                            if energy_matches:
                                for energy in energy_matches:
                                    scores.append(float(energy))
                            # Also try to find alternative energy format
                            alt_matches = re.findall(r'binding energy\s*(-?\d+\.\d+)', dlg_content)
                            if alt_matches and not scores:
                                for energy in alt_matches:
                                    scores.append(float(energy))
                    except Exception as e:
                        print(f"Error parsing dlg file for energies: {e}")
                        
        except Exception as e:
            print(f"Error extracting binding energies from summary.json: {e}")
        
        return scores
    
    @property
    def binding_energies(self) -> List[float]:
        """Get list of binding energies"""
        return self._get_score_from_log()
    
    
    def to_json(self) -> str:
        """Convert to JSON string with hierarchical structure"""
        result_dict = {
            "protein_name": self.protein_name,
            "ligand_name": self.ligand_name,
            "binding_energies": self.binding_energies,
            "files": {k: str(v) for k, v in self.docking_files.items()},
            "output_dir": str(self.output_dir)
        }
        return json.dumps(result_dict, indent=2)
    
    def to_markdown(self) -> str:
        """Create a markdown table with results"""
        # Header
        md = f"# AutoDock Vina Results\n\n"
        md += f"## Protein: {self.protein_name}\n"
        md += f"## Ligand: {self.ligand_name}\n\n"
        
        # Binding energies section (show this first as it's most important)
        energies = self.binding_energies
        if energies:
            md += "### Binding Energies (kcal/mol)\n\n"
            md += "| Pose | Binding Energy | Interpretation |\n"
            md += "|------|---------------|----------------|\n"
            for i, energy in enumerate(energies):
                # Add interpretation (lower/more negative is better)
                interpretation = "Strong binding" if energy < -9.0 else \
                                "Good binding" if energy < -7.0 else \
                                "Moderate binding" if energy < -5.0 else "Weak binding"
                md += f"| {i+1} | {energy:.2f} | {interpretation} |\n"
            
            md += "\n**Note**: Lower (more negative) binding energies indicate stronger binding.\n"
        else:
            md += "### No binding energy data available\n\n"
        
        # Best pose file highlight
        best_pose = self.docking_files.get('best_pose')
        if best_pose:
            md += f"\n### Best Docking Pose\n\n"
            md += f"- **File**: {best_pose}\n"
        
        # Other files section
        md += "\n### All Files\n\n"
        for file_type, file_path in self.docking_files.items():
            if file_type != 'best_pose':  # Already showed best pose
                md += f"- **{file_type}**: {file_path}\n"
        
        return md

class VinaDockingResult(BaseModel):
    """Stores and analyzes results from a Vina docking run"""
    protein_name: str
    ligand_name: str
    output_dir: Path
    docking_files: Dict[str, Path] = Field(description="Dictionary of file types to their paths")
    
    class Config:
        arbitrary_types_allowed = True
    
    def _get_score_from_log(self) -> List[float]:
        """Extract binding energies from summary.json file"""
        scores = []
        import re
        
        # Get the ligand number from ligand_pdbqt filename - server uses 1-based indexing
        ligand_num = None
        ligand_pdbqt = self.docking_files.get('ligand_pdbqt')
        
        if ligand_pdbqt:
            if 'docking_result_best_pose' in str(ligand_pdbqt):
                # Special case for single molecule - use index 1
                ligand_num = "1"
                # Check for dlg file without index too
                dlg_file = self.output_dir / "docking_result.dlg"
                if dlg_file.exists() and 'docking_result_dlg' not in self.docking_files:
                    self.docking_files['docking_result_dlg'] = dlg_file
            else:
                # Match ligand_X pattern (server uses 1-based, so ligand_1, ligand_2, etc.)
                match = re.search(r'ligand_(\d+)', str(ligand_pdbqt))
                if match:
                    # The server uses 1-based indexing, so we can use this directly
                    ligand_num = match.group(1)
        
        # If we couldn't determine ligand number, we can't extract binding energy
        if not ligand_num:
            return scores
        
        # Look for summary.json in the output directory
        summary_path = self.output_dir / "summary.json"
        if not summary_path.exists():
            return scores
        
        # Load and parse summary.json
        try:
            with open(summary_path, 'r') as f:
                summary = json.load(f)
            
            # Get the specific ligand data
            ligand_key = f"ligand_{ligand_num}"
            if ligand_key in summary.get('ligands', {}):
                ligand_data = summary['ligands'][ligand_key]
                
                # Check if docking was successful
                if ligand_data.get('success', False):
                    # Extract binding energies from poses
                    for pose in ligand_data.get('poses', []):
                        if 'binding_energy' in pose:
                            scores.append(float(pose['binding_energy']))
                            
            # Special case for single molecule - also check the dlg file directly
            if ligand_num == "1" and not scores:
                # Try to find the binding energy directly in the dlg file
                dlg_file = self.docking_files.get('docking_result_dlg')
                if dlg_file and dlg_file.exists():
                    try:
                        with open(dlg_file, 'r') as f:
                            dlg_content = f.read()
                            # Look for binding energy patterns in the dlg file
                            energy_matches = re.findall(r'Estimated Free Energy of Binding\s*=\s*(-?\d+\.\d+)', dlg_content)
                            if energy_matches:
                                for energy in energy_matches:
                                    scores.append(float(energy))
                            # Also try to find alternative energy format
                            alt_matches = re.findall(r'binding energy\s*(-?\d+\.\d+)', dlg_content)
                            if alt_matches and not scores:
                                for energy in alt_matches:
                                    scores.append(float(energy))
                    except Exception as e:
                        print(f"Error parsing dlg file for energies: {e}")
                        
        except Exception as e:
            print(f"Error extracting binding energies from summary.json: {e}")
        
        return scores
    
    @property
    def binding_energies(self) -> List[float]:
        """Get list of binding energies"""
        return self._get_score_from_log()
    
    
    def to_json(self) -> str:
        """Convert to JSON string with hierarchical structure"""
        result_dict = {
            "protein_name": self.protein_name,
            "ligand_name": self.ligand_name,
            "binding_energies": self.binding_energies,
            "files": {k: str(v) for k, v in self.docking_files.items()},
            "output_dir": str(self.output_dir)
        }
        return json.dumps(result_dict, indent=2)
    
    def to_markdown(self) -> str:
        """Create a markdown table with results"""
        # Header
        md = f"# AutoDock Vina Results\n\n"
        md += f"## Protein: {self.protein_name}\n"
        md += f"## Ligand: {self.ligand_name}\n\n"
        
        # Binding energies section (show this first as it's most important)
        energies = self.binding_energies
        if energies:
            md += "### Binding Energies (kcal/mol)\n\n"
            md += "| Pose | Binding Energy | Interpretation |\n"
            md += "|------|---------------|----------------|\n"
            for i, energy in enumerate(energies):
                # Add interpretation (lower/more negative is better)
                interpretation = "Strong binding" if energy < -9.0 else \
                                "Good binding" if energy < -7.0 else \
                                "Moderate binding" if energy < -5.0 else "Weak binding"
                md += f"| {i+1} | {energy:.2f} | {interpretation} |\n"
            
            md += "\n**Note**: Lower (more negative) binding energies indicate stronger binding.\n"
        else:
            md += "### No binding energy data available\n\n"
        
        # Best pose file highlight
        best_pose = self.docking_files.get('best_pose')
        if best_pose:
            md += f"\n### Best Docking Pose\n\n"
            md += f"- **File**: {best_pose}\n"
        
        # Other files section
        md += "\n### All Files\n\n"
        for file_type, file_path in self.docking_files.items():
            if file_type != 'best_pose':  # Already showed best pose
                md += f"- **{file_type}**: {file_path}\n"
        
        return md

class VinaDockerAPI:
    """Client for Vina Docking API."""
    
    def __init__(self, base_url: str = "http://localhost:8686") -> None:
        """Initialize Vina Docking API client.
        
        :param base_url: Base URL for Vina API (default: local server)
        """
        self.base_url = base_url.rstrip('/')
        print(f"VinaDockerAPI initialized with base_url: {self.base_url}")
    
    def ping(self) -> bool:
        """Check if API is alive.
        :return: True if API responds with pong
        """
        ping_url = f"{self.base_url}/ping"
        print(f"Pinging Vina server at: {ping_url}")
        try:
            response = requests.get(ping_url)
            print(f"Ping response status code: {response.status_code}")
            if response.status_code == 200:
                json_resp = response.json()
                print(f"Ping response content: {json_resp}")
                return json_resp.get('message') == 'pong'
            else:
                print(f"Ping failed with status code: {response.status_code}")
                print(f"Response text: {response.text}")
                return False
        except Exception as e:
            print(f"Exception during ping: {str(e)}")
            return False
    
    def get_task_status(self, task_id: str) -> dict:
        """Get status of a task.
        
        :param task_id: Task ID to check
        :return: Task status information
        """
        url = f"{self.base_url}/task_status/{task_id}"
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    
    def predict(
        self,
        protein_file: Union[str, Path],
        ligand_file: Union[str, Path, List[Path]],
        output_dir: Union[str, Path],
        smiles_list: Optional[Union[str, List[str]]] = None,
        config: Optional[Dict[str, Any]] = None,
        timeout: int = 3600,  # Default timeout: 1 hour
        polling_interval: int = 30,
        to_markdown: bool = False,
        merge_pdbqt: bool = True,
        batch_id: Optional[int] = None
    ) -> Union[Dict[str, Any], List[Dict[str, Any]], str, List[str]]:
        """Run Vina docking prediction.
        
        :param protein_file: Path to protein PDB file
        :param ligand_file: Path to ligand SDF file or list of ligand SDF files
        :param output_dir: Directory to save results
        :param smiles_list: SMILES string or list of SMILES strings corresponding to ligand files
        :param config: Optional configuration dict for Vina (e.g., {"num_runs": 5, "local_search_method": "sw", "heuristics": 1})
        :param timeout: Maximum time to wait for results in seconds
        :param polling_interval: Time between download attempts in seconds
        :param to_markdown: If True, returns markdown representation of results
        :param merge_pdbqt: If True, creates merged complex files from protein.pdbqt and ligand.pdbqt
        :param batch_id: Optional batch ID for tracking multiple batches
        :return: Dictionary or list of dictionaries with smiles, binding_energy, and merged_path
        """
        # Start batch tracking
        batch_info = f"BATCH-{batch_id}" if batch_id is not None else "NO-BATCH"
        print(f"[DEBUG-{batch_info}] Starting Vina docking predict() method")
        
        # Ensure output directory exists
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"[DEBUG-{batch_info}] Output directory: {output_dir}")
        
        # Handle single vs multiple ligands
        is_batch = isinstance(ligand_file, list) and len(ligand_file) > 1
        batch_size = len(ligand_file) if isinstance(ligand_file, list) else 1
        print(f"[DEBUG-{batch_info}] Is batch mode: {is_batch}, Ligands count: {batch_size}")
        
        # Adjust timeout based on batch size (more ligands = more time needed)
        adjusted_timeout = timeout * max(1, batch_size // 10)
        print(f"[DEBUG-{batch_info}] Adjusted timeout: {adjusted_timeout} seconds")
        
        # Convert smiles to list if it's a single string
        if smiles_list is not None:
            if not isinstance(smiles_list, list):
                smiles_list = [smiles_list]
                print(f"[DEBUG-{batch_info}] Converted single SMILES to list")
            print(f"[DEBUG-{batch_info}] SMILES list length: {len(smiles_list)}")
        else:
            # If no SMILES provided, create placeholders
            if isinstance(ligand_file, list):
                smiles_list = [None] * len(ligand_file)
                print(f"[DEBUG-{batch_info}] Created placeholder SMILES list of length {len(ligand_file)}")
            else:
                smiles_list = [None]
                print(f"[DEBUG-{batch_info}] Created single placeholder SMILES")
        
        # Get ligand names for logging
        if isinstance(ligand_file, list):
            ligand_names = [Path(lf).name for lf in ligand_file]
            ligand_display = f"{len(ligand_file)} ligands: {', '.join(ligand_names[:3])}" + ("..." if len(ligand_file) > 3 else "")
        else:
            ligand_names = Path(ligand_file).name
            ligand_display = ligand_names
        
        # Create a specific output directory for this batch if batch_id is provided
        if batch_id is not None:
            batch_dir = output_dir / f"vina_docking_{batch_id}"
            batch_dir.mkdir(parents=True, exist_ok=True)
            output_dir = batch_dir
            print(f"[DEBUG-{batch_info}] Created batch directory: {batch_dir}")
        
        # Submit inference job - with retries
        max_submission_retries = 5
        submission_retry_delay = 2
        protein_file = clean_pdb_file(protein_file)
        print(f"[DEBUG-{batch_info}] Submitting Vina docking for protein: {protein_file} with {ligand_display}")
        
        # Use try/except to catch and handle submission failures
        task_id = None
        submission_error = None
        
        for attempt in range(max_submission_retries):
            try:
                # Use a small delay before each retry to avoid overwhelming the server
                if attempt > 0:
                    delay = submission_retry_delay * (2 ** attempt) # Exponential backoff
                    print(f"[DEBUG-{batch_info}] Retry submission attempt {attempt+1}/{max_submission_retries} after {delay}s delay")
                    time.sleep(delay)
                
                task_id = self._submit_inference(protein_file, ligand_file, config)
                logger.info(f"Submitted Vina docking job: {task_id}")
                protein_name = Path(protein_file).name
                print(f"[DEBUG-{batch_info}] Successfully submitted job, task_id: {task_id}")
                submission_error = None
                break
            except Exception as e:
                submission_error = str(e)
                logger.warning(f"Submission attempt {attempt+1} failed for batch {batch_id}: {str(e)}")
                print(f"[DEBUG-{batch_info}] Submission attempt {attempt+1} failed: {str(e)}")
                # Continue to next retry attempt
        
        # If all submission attempts failed, raise the final error
        if task_id is None:
            error_msg = f"Failed to submit Vina docking job for batch {batch_id} after {max_submission_retries} attempts: {submission_error}"
            logger.error(error_msg)
            print(f"[DEBUG-{batch_info}] {error_msg}")
            raise RuntimeError(error_msg)

        # Poll for results - with dedicated session for better connection management
        start_time = time.time()
        # Create a dedicated session for this batch's polling to avoid connection pool issues
        with requests.Session() as polling_session:
            while True:
                # Check task status with retries
                status_check_retries = 3
                status = None
                
                for status_attempt in range(status_check_retries):
                    try:
                        # Use polling_session for status checks
                        url = f"{self.base_url}/task_status/{task_id}"
                        response = polling_session.get(url, timeout=(10, 30))
                        response.raise_for_status()
                        status = response.json()
                        logger.info(f"[{batch_info}] Task status: {status}")
                        break
                    except (requests.exceptions.ConnectionError, requests.exceptions.Timeout) as e:
                        if status_attempt < status_check_retries - 1:
                            print(f"[DEBUG-{batch_info}] Status check attempt {status_attempt+1} failed: {str(e)}. Retrying...")
                            time.sleep(2)
                        else:
                            print(f"[DEBUG-{batch_info}] All status check attempts failed: {str(e)}")
                            # If all status checks fail, wait longer before trying again
                            time.sleep(polling_interval)
                            continue
                
                # If we couldn't get status after retries, continue to next polling cycle
                if status is None:
                    print(f"[DEBUG-{batch_info}] Could not get task status after {status_check_retries} attempts, continuing to wait...")
                    time.sleep(polling_interval)
                    continue
                
                if status.get('status') == 'completed':
                    try:
                        print(f"[DEBUG-{batch_info}] Task completed. Downloading results for task {task_id}")
                        result_path = self._download_results(task_id, output_dir)
                        print(f"[DEBUG-{batch_info}] Download completed. Result path: {result_path}")
                        if not result_path.exists():
                            print(f"[DEBUG-{batch_info}] Warning: Result path does not exist: {result_path}")
                            return None
                        
                        print(f"[DEBUG-{batch_info}] Processing downloaded results...")
                        results = self._process_results(task_id, output_dir, protein_name, ligand_names, result_path)
                        if results is None:
                            print(f"[DEBUG-{batch_info}] Warning: _process_results returned None")
                            return None
                        print(f"[DEBUG-{batch_info}] Results: {results}")
                        
                        # Merge PDBQT files if requested
                        merged_paths = []
                        if merge_pdbqt:
                            print(f"[DEBUG-{batch_info}] Merging protein and ligand PDBQT files into complex PDB files")
                            # Create unique subdirectory for merged results with batch ID if available
                            merge_output_dir = output_dir / f"vina_{task_id}_merged"
                            merged_paths = self.merge_vina_results(results, merge_output_dir)
                            print(f"[DEBUG-{batch_info}] Created {len(merged_paths)} merged complex files")
                        
                        # Create new format of results with the requested fields
                        if isinstance(results, list):
                            # Create list of dictionaries with the requested format
                            formatted_results = []
                            for i, result in enumerate(results):
                                # Get binding energy
                                binding_energies = result.binding_energies
                                binding_energy = min(binding_energies) if binding_energies else None
                                
                                # Create result dictionary
                                result_dict = {
                                    "smiles": smiles_list[i] if i < len(smiles_list) else None,
                                    "binding_energy": binding_energy,
                                    "merged_path": str(merged_paths[i]) if merge_pdbqt and i < len(merged_paths) else None,
                                    "batch_id": batch_id  # Include batch ID for better tracking
                                }
                                formatted_results.append(result_dict)
                            
                            # Return markdown or dictionary as requested
                            if to_markdown:
                                md_results = []
                                for result_dict in formatted_results:
                                    md = f"# Docking Result\n\n"
                                    md += f"- **SMILES**: {result_dict['smiles']}\n"
                                    md += f"- **Binding Energy**: {result_dict['binding_energy']} kcal/mol\n"
                                    md += f"- **Complex Path**: {result_dict['merged_path']}\n"
                                    md += f"- **Batch ID**: {result_dict['batch_id']}\n"
                                    md_results.append(md)
                                return md_results
                            else:
                                print(f"[DEBUG-{batch_info}] Successfully completed batch, returning {len(formatted_results)} results")
                                return formatted_results
                        else:
                            # Single result
                            binding_energies = results.binding_energies
                            binding_energy = min(binding_energies) if binding_energies else None
                            
                            result_dict = {
                                "smiles": smiles_list[0] if smiles_list else None,
                                "binding_energy": binding_energy,
                                "merged_path": str(merged_paths[0]) if merge_pdbqt and merged_paths else None,
                                "batch_id": batch_id  # Include batch ID for better tracking
                            }
                            
                            if to_markdown:
                                md = f"# Docking Result\n\n"
                                md += f"- **SMILES**: {result_dict['smiles']}\n"
                                md += f"- **Binding Energy**: {result_dict['binding_energy']} kcal/mol\n"
                                md += f"- **Complex Path**: {result_dict['merged_path']}\n"
                                md += f"- **Batch ID**: {result_dict['batch_id']}\n"
                                return md
                            else:
                                print(f"[DEBUG-{batch_info}] Successfully completed batch, returning single result")
                                return result_dict
                        
                    except Exception as e:
                        print(f"[DEBUG-{batch_info}] Error processing results: {str(e)}")
                        import traceback
                        print(f"[DEBUG-{batch_info}] Traceback: {traceback.format_exc()}")
                        raise
                elif status.get('status') == 'failed':
                    print(f"[DEBUG-{batch_info}] Task failed with status: {status}")
                    raise RuntimeError(f"Task failed: {status.get('error', 'Unknown error')}")
                
                if time.time() - start_time > adjusted_timeout:
                    raise TimeoutError(
                        f"Vina docking job {task_id} (batch {batch_info}) timed out after {adjusted_timeout} seconds. "
                        f"Last status: {status}"
                    )
                
                # Exponential backoff for polling interval
                current_wait = polling_interval
                # For larger batches or if waiting a long time, gradually increase polling interval
                elapsed_time = time.time() - start_time
                if elapsed_time > 300:  # After 5 minutes
                    current_wait = min(polling_interval * 2, 120)  # Max 2 minutes
                
                print(f"[DEBUG-{batch_info}] Waiting {current_wait} seconds before checking status again...")
                time.sleep(current_wait)
    
    def _submit_inference(self, protein_file: Union[str, Path], ligand_file: Union[str, Path, List[Path]], config: Optional[Dict[str, Any]] = None) -> str:
        """Submit inference job to Vina Docking API.
        
        :param protein_file: Path to protein PDB file
        :param ligand_file: Path to ligand SDF file or list of ligand SDF files
        :param config: Optional configuration dict for Vina
        :return: Task ID
        """
        url = f"{self.base_url}/inference"
        print(f"[DEBUG-SUBMIT] Starting _submit_inference to {url}")
        
        # Ensure ligand_file is a list
        if not isinstance(ligand_file, list):
            ligand_file = [ligand_file]
            print(f"[DEBUG-SUBMIT] Converted single ligand file to list")
            
        print(f"[DEBUG-SUBMIT] Processing {len(ligand_file)} ligand files")
        
        # Create a new session for each batch submission to avoid connection pool issues
        with requests.Session() as session:
            # Read all files into memory with proper context management
            files = {}
            
            try:
                # Read protein file into memory
                with open(protein_file, 'rb') as pf:
                    protein_content = pf.read()
                    files['protein_file'] = (Path(protein_file).name, protein_content, 'application/octet-stream')
                
                # Read all ligand files into memory
                for i, lf_path in enumerate(ligand_file):
                    print(f"[DEBUG-SUBMIT] Opening ligand file {i+1}/{len(ligand_file)}: {lf_path}")
                    with open(lf_path, 'rb') as lf:
                        ligand_content = lf.read()
                        # Use separate numbered keys for each file
                        files[f'ligand_files_{i}'] = (Path(lf_path).name, ligand_content, 'application/octet-stream')
                
                # Add config as form data, not as a file
                # Create separate data dictionary for non-file form fields
                data = {}
                if config:
                    data['config'] = json.dumps(config)
                
                print(f"[DEBUG-SUBMIT] Submitting inference to Vina server with {len(files)} files")
                print(f"[DEBUG-SUBMIT] File keys: {list(files.keys())}")
                print(f"[DEBUG-SUBMIT] Data: {data}")
                
                # First try to ping the server to check connectivity
                ping_result = self.ping()
                print(f"[DEBUG-SUBMIT] Ping result before submitting: {ping_result}")
                
                # Implement retry logic with exponential backoff
                max_retries = 3
                retry_delay = 5
                
                for attempt in range(max_retries):
                    try:
                        print(f"[DEBUG-SUBMIT] Sending POST request (attempt {attempt + 1}/{max_retries})...")
                        logger.info(f"with data: {data} and files names: {list(files.keys())}")
                        # Increase timeout: (connect_timeout, read_timeout)
                        response = session.post(url, files=files, data=data, timeout=(30, 1800))
                        print(f"[DEBUG-SUBMIT] Response status code: {response.status_code}")
                        
                        if response.status_code != 200 and response.status_code != 202:
                            print(f"[DEBUG-SUBMIT] Error response from server: {response.text}")
                        
                        response.raise_for_status()
                        
                        result = response.json()
                        task_id = result.get('task_id')
                        print(f"[DEBUG-SUBMIT] Successfully received task_id: {task_id}")
                        return task_id
                    except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as e:
                        if attempt < max_retries - 1:
                            wait_time = retry_delay * (2 ** attempt)  # Exponential backoff
                            print(f"[DEBUG-SUBMIT] Attempt {attempt+1} failed with {type(e).__name__}: {str(e)}")
                            print(f"[DEBUG-SUBMIT] Retrying in {wait_time} seconds...")
                            time.sleep(wait_time)
                        else:
                            print(f"[DEBUG-SUBMIT] All {max_retries} attempts failed. Last error: {str(e)}")
                            raise
                    except Exception as e:
                        print(f"[DEBUG-SUBMIT] Unexpected error: {type(e).__name__}: {str(e)}")
                        import traceback
                        print(f"[DEBUG-SUBMIT] Traceback: {traceback.format_exc()}")
                        raise
                
                # This should not be reached if max_retries was exhausted (it would raise in the loop)
                raise RuntimeError("[DEBUG-SUBMIT] Failed to submit inference after exhausting all retries")
            except Exception as e:
                print(f"[DEBUG-SUBMIT] Exception during inference submission: {str(e)}")
                import traceback
                print(f"[DEBUG-SUBMIT] Traceback: {traceback.format_exc()}")
                raise
    
    def _download_results(self, task_id: str, output_dir: Path) -> Optional[Path]:
        """Download results.
        
        :param task_id: Task ID
        :param output_dir: Output directory
        :return: Path to results if successful, None if not ready
        """
        url = f"{self.base_url}/download_result/{task_id}"
        print(f"[DEBUG-DOWNLOAD] Downloading results from: {url}")
        
        # Create a new session for downloading to avoid connection pool issues
        with requests.Session() as session:
            # Implement retry logic with exponential backoff
            max_retries = 3
            retry_delay = 5
            
            for attempt in range(max_retries):
                try:
                    # Stream the response with increased timeout
                    with session.get(url, stream=True, timeout=(30, 1800)) as response:
                        print(f"[DEBUG-DOWNLOAD] Response status code: {response.status_code}")
                        
                        # Check response status before raising for status
                        if response.status_code != 200:
                            print(f"[DEBUG-DOWNLOAD] Non-200 response: {response.text[:500]}")
                        
                        response.raise_for_status()
                        
                        # Save zip file
                        zip_path = output_dir / f"{task_id}_vina_output.zip"
                        print(f"[DEBUG-DOWNLOAD] Saving to: {zip_path}")
                        
                        total_size = 0
                        with open(zip_path, 'wb') as f:
                            for i, chunk in enumerate(response.iter_content(chunk_size=8192)):
                                if chunk:  # filter out keep-alive chunks
                                    chunk_size = len(chunk)
                                    total_size += chunk_size
                                    f.write(chunk)
                                    if i % 100 == 0:  # Log every ~800KB
                                        print(f"[DEBUG-DOWNLOAD] Downloaded {total_size/1024:.1f}KB so far...")
                        
                        print(f"[DEBUG-DOWNLOAD] Download complete, total size: {total_size/1024:.1f}KB")
                        
                        # Verify the file exists and has content
                        if zip_path.exists():
                            file_size = zip_path.stat().st_size
                            print(f"[DEBUG-DOWNLOAD] Zip file exists with size: {file_size/1024:.1f}KB")
                            if file_size == 0:
                                print(f"[DEBUG-DOWNLOAD] WARNING: Downloaded file is empty!")
                                if attempt < max_retries - 1:
                                    # If file is empty, try again
                                    raise ValueError("Downloaded zip file is empty")
                            else:
                                return zip_path
                        else:
                            print(f"[DEBUG-DOWNLOAD] ERROR: Zip file does not exist after download!")
                            raise FileNotFoundError("Zip file not created during download")
                
                except (requests.exceptions.Timeout, requests.exceptions.ConnectionError, 
                        ValueError, FileNotFoundError) as e:
                    if attempt < max_retries - 1:
                        wait_time = retry_delay * (2 ** attempt)  # Exponential backoff
                        print(f"[DEBUG-DOWNLOAD] Attempt {attempt+1} failed with {type(e).__name__}: {str(e)}")
                        print(f"[DEBUG-DOWNLOAD] Retrying in {wait_time} seconds...")
                        time.sleep(wait_time)
                    else:
                        print(f"[DEBUG-DOWNLOAD] All {max_retries} attempts failed. Last error: {str(e)}")
                        raise
                except Exception as e:
                    print(f"[DEBUG-DOWNLOAD] Error downloading results: {str(e)}")
                    import traceback
                    print(f"[DEBUG-DOWNLOAD] Traceback: {traceback.format_exc()}")
                    raise
            
            # This should not be reached if max_retries was exhausted
            raise RuntimeError("[DEBUG-DOWNLOAD] Failed to download results after exhausting all retries")
    
    def _process_results(self, task_id: str, output_dir: Path, protein_name: str, 
                         ligand_name: Union[str, List[str]], zip_path: Path) -> Union[VinaDockingResult, List[VinaDockingResult]]:
        """Process the Vina docking results from a zip file
        
        :param task_id: Task ID
        :param output_dir: Output directory
        :param protein_name: Name of the protein file
        :param ligand_name: Name of the ligand file or list of ligand names
        :param zip_path: Path to the downloaded zip file
        :return: Single VinaDockingResult or list of VinaDockingResult objects
        """
        print(f"[DEBUG-PROCESS] Starting _process_results for task_id: {task_id}")
        
        result_dir = output_dir / f"vina_{task_id}"
        result_dir.mkdir(parents=True, exist_ok=True)
        print(f"[DEBUG-PROCESS] Output directory: {result_dir}")
        
        # Extract files from zip
        try:
            print(f"[DEBUG-PROCESS] Extracting zip file: {zip_path}")
            with zipfile.ZipFile(zip_path, 'r') as zip_ref:
                files_in_zip = zip_ref.namelist()
                print(f"[DEBUG-PROCESS] Files in zip: {files_in_zip[:10]}..." if len(files_in_zip) > 10 else f"[DEBUG-PROCESS] Files in zip: {files_in_zip}")
                zip_ref.extractall(result_dir)
                print(f"[DEBUG-PROCESS] Extracted {len(files_in_zip)} files to {result_dir}")
        except Exception as e:
            print(f"[DEBUG-PROCESS] Error extracting zip file: {str(e)}")
            import traceback
            print(f"[DEBUG-PROCESS] Traceback: {traceback.format_exc()}")
            raise
        
        # Process the flat structure format with summary.json
        print(f"[DEBUG-PROCESS] Processing results with flat structure format from: {result_dir}")
        
        # Convert ligand_name to a list if it's not already
        if not isinstance(ligand_name, list):
            ligand_names = [ligand_name]
            print(f"[DEBUG-PROCESS] Converted single ligand name to list")
        else:
            ligand_names = ligand_name
            print(f"[DEBUG-PROCESS] Using list of {len(ligand_names)} ligand names")
            
        # Map ligand indices to original names
        ligand_name_map = {}
        for i, name in enumerate(ligand_names):
            ligand_name_map[i+1] = name
        print(f"[DEBUG-PROCESS] Created ligand name map: {ligand_name_map}")
        
        # Check for summary.json
        summary_path = result_dir / "summary.json"
        summary_data = {}
        if summary_path.exists():
            try:
                with open(summary_path, 'r') as f:
                    summary_data = json.load(f)
                print(f"[DEBUG-PROCESS] Found summary.json with {len(summary_data.get('ligands', {}))} ligands")
                print(f"[DEBUG-PROCESS] Summary keys: {list(summary_data.keys())}")
                if 'ligands' in summary_data:
                    print(f"[DEBUG-PROCESS] Ligand keys in summary: {list(summary_data['ligands'].keys())}")
            except Exception as e:
                print(f"[DEBUG-PROCESS] Error loading summary.json: {str(e)}")
                import traceback
                print(f"[DEBUG-PROCESS] Traceback: {traceback.format_exc()}")
        else:
            print(f"[DEBUG-PROCESS] summary.json not found in {result_dir}")
            # List all files in the directory to see what we have
            all_files = list(result_dir.glob("*"))
            print(f"[DEBUG-PROCESS] Files in result directory: {[f.name for f in all_files]}")
        
        # Find all ligand files (excluding best_pose files which we'll handle separately)
        ligand_files = [f for f in result_dir.glob("ligand_*.pdbqt") if "best_pose" not in str(f)]
        print(f"[DEBUG-PROCESS] Found {len(ligand_files)} ligand files: {[f.name for f in ligand_files]}")
        
        # Also check for best pose files from single molecule docking
        best_pose_files = list(result_dir.glob("docking_result_best_pose.pdbqt"))
        if best_pose_files and not ligand_files:
            print(f"[DEBUG-PROCESS] Found single molecule best pose file: {best_pose_files[0]}")
            # Create a virtual ligand file entry with index 1 for the best pose
            ligand_files = [best_pose_files[0]]
        
        # Process each ligand
        results = []
        print(f"[DEBUG-PROCESS] Starting processing of {len(ligand_files)} ligand files")
        
        # Process ligand files with special handling for the best_pose file case
        if ligand_files:
            # Sort the ligand files - handle special case for best_pose file 
            try:
                # First try standard sorting for ligand_1, ligand_2, etc.
                sorted_ligand_files = sorted(ligand_files, key=lambda x: int(x.stem.split('_')[1]) if x.stem.startswith('ligand_') else 1)
            except (ValueError, IndexError):
                # If that fails, handle best_pose file or other special cases
                print("Using alternative file sorting method for special filenames")
                sorted_ligand_files = ligand_files
                
            # Process each file
            for ligand_file in sorted_ligand_files:
                # Get 1-based index from filename (server uses 1-based) - handle special cases
                if ligand_file.stem.startswith('ligand_'):
                    # Regular case: ligand_1, ligand_2, etc.
                    try:
                        server_idx = int(ligand_file.stem.split('_')[1])
                    except (ValueError, IndexError):
                        print(f"Warning: Could not extract index from {ligand_file.stem}, using default index 1")
                        server_idx = 1
                elif 'docking_result_best_pose' in str(ligand_file):
                    # Special case for single molecule result
                    server_idx = 1
                else:
                    # Fallback for any other format - use index 1
                    server_idx = 1
                
                # Summary.json also uses 1-based (ligand_1, ligand_2, etc.)
                summary_idx = server_idx
                
                # Get ligand name from the map or use default
                # Our map is 1-based (position 1 = first SMILES), server is also 1-based
                current_ligand_name = ligand_name_map.get(server_idx, f"ligand_{server_idx}")
                
                # Collect related files
                docking_files = {
                    'ligand_pdbqt': ligand_file,
                }
                
                # Look for docking result files - server uses 1-based indexing
                dlg_file = result_dir / f"docking_result_{server_idx}.dlg"
                if dlg_file.exists():
                    docking_files['docking_result_dlg'] = dlg_file
                
                # Special case for best_pose file
                if 'best_pose' in str(ligand_file):
                    docking_files['best_pose'] = ligand_file
                    
                    # Check for the corresponding dlg file without index for single molecule case
                    single_dlg = result_dir / "docking_result.dlg"
                    if single_dlg.exists():
                        docking_files['docking_result_dlg'] = single_dlg
                        
                xml_file = result_dir / f"docking_result_{server_idx}.xml"
                if xml_file.exists():
                    docking_files['docking_result_xml'] = xml_file
                
                # Add protein file
                protein_pdbqt = result_dir / "protein.pdbqt"
                if protein_pdbqt.exists():
                    docking_files['protein_pdbqt'] = protein_pdbqt
                
                # Get status and best pose info from summary.json if available
                success = False
                if summary_data and 'ligands' in summary_data:
                    # Server uses 1-based indexing in summary.json (ligand_1, ligand_2, etc.)
                    ligand_key = f"ligand_{summary_idx}"
                    if ligand_key in summary_data['ligands']:
                        ligand_status = summary_data['ligands'][ligand_key]
                        success = ligand_status.get('success', False)
                        
                        # Add best pose file if available
                        if success and 'poses' in ligand_status:
                            for pose in ligand_status['poses']:
                                if 'filename' in pose:
                                    pose_file = result_dir / pose['filename']
                                    if pose_file.exists():
                                        docking_files['best_pose'] = pose_file
                                        break
                
                # Create result object
                result = VinaDockingResult(
                    protein_name=protein_name,
                    ligand_name=current_ligand_name,
                    output_dir=result_dir,
                    docking_files=docking_files
                )
                
                results.append(result)
        
        return results
        
    def merge_pdbqt_complex(self, protein_pdbqt: Path, ligand_pdbqt: Path, output_dir: Path, 
                            output_filename: Optional[str] = None) -> Path:
        """Merge protein.pdbqt and ligand.pdbqt files into a single complex file for PLIP analysis.
        
        Args:
            protein_pdbqt: Path to the protein PDBQT file
            ligand_pdbqt: Path to the ligand PDBQT file
            output_dir: Directory to save the merged file
            output_filename: Optional filename for the output file
            
        Returns:
            Path to the merged complex file
        """
        # Create output directory if it doesn't exist
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate output filename if not provided
        if output_filename is None:
            protein_name = protein_pdbqt.stem
            ligand_name = ligand_pdbqt.stem
            output_filename = f"{protein_name}_{ligand_name}_complex.pdb"
        
        output_path = output_dir / output_filename
        
        # Read protein file lines - keep ATOM records
        protein_lines = []
        with open(protein_pdbqt, 'r') as f:
            for line in f:
                if line.startswith('ATOM'):
                    # Convert PDBQT format to PDB format by removing charge info
                    # PDB format uses 80 character records; PDBQT adds charges at the end
                    # We'll keep just the standard PDB fields
                    protein_lines.append(line[:66] + '\n')
        
        # Find the last atom number from protein
        last_atom_number = 0
        if protein_lines:
            for line in reversed(protein_lines):
                if line.startswith('ATOM'):
                    try:
                        last_atom_number = int(line[6:11].strip())
                        break
                    except ValueError:
                        pass
        
        # Read ligand file lines
        ligand_lines = []
        raw_ligand_lines = []  # Store the original lines before renumbering
        with open(ligand_pdbqt, 'r') as f:
            ligand_content = f.read()
        
        # Check for the different file formats
        is_best_pose_format = "MODEL" in ligand_content and "USER" in ligand_content
        is_single_line_format = ligand_content.count('\n') < 10  # Very few line breaks indicates single-line format
        
        print(f"Processing {ligand_pdbqt}: is_best_pose_format={is_best_pose_format}, is_single_line_format={is_single_line_format}")
        
        if is_single_line_format:
            # Handle single-line format - split by USER tags and process
            print("Handling single-line best pose format")
            if "ATOM" in ligand_content:
                # Extract each ATOM entry from the single line
                # This uses regex to find all ATOM entries with their data
                import re
                atom_pattern = r'ATOM\s+\d+\s+\S+\s+\S+\s+\d+\s+[-]?\d+\.\d+\s+[-]?\d+\.\d+\s+[-]?\d+\.\d+'
                atom_matches = re.findall(atom_pattern, ligand_content)
                
                for atom_line in atom_matches:
                    # Convert ATOM to HETATM for ligands
                    hetatm_line = 'HETATM' + atom_line[4:66] + '\n'
                    raw_ligand_lines.append(hetatm_line)
        elif is_best_pose_format:
            # Split the content by lines for normal processing
            lines = ligand_content.split('\n')
            in_atom_section = False
            
            for line in lines:
                # Skip MODEL, USER, and REMARK lines until we find atom data
                if line.startswith(('MODEL', 'USER', 'REMARK')) and not in_atom_section:
                    # Look for the line that indicates start of coordinates
                    if 'x       y       z' in line:
                        in_atom_section = True
                    continue
                
                # Once we're in the atom section, process ATOM/HETATM lines
                if in_atom_section and line.startswith(('ATOM', 'HETATM')):
                    # Convert ATOM to HETATM for ligands
                    if line.startswith('ATOM'):
                        line = 'HETATM' + line[6:]
                    # Remove charge info
                    raw_ligand_lines.append(line[:66] + '\n')
                
                # Stop at TORSDOF or ENDMDL
                if line.startswith(('TORSDOF', 'ENDMDL')):
                    break
        else:
            # Regular PDBQT format
            for line in ligand_content.split('\n'):
                if line.startswith(('ATOM', 'HETATM')):
                    # Convert ATOM to HETATM for ligands
                    if line.startswith('ATOM'):
                        line = 'HETATM' + line[6:]
                    # Remove charge info
                    raw_ligand_lines.append(line[:66] + '\n')
        
        # If we didn't get any ligand lines, try a more aggressive extraction method
        if not raw_ligand_lines:
            print(f"Warning: No ligand atoms found using standard parsing. Trying advanced method for {ligand_pdbqt}")
            
            # Special handling for AutoDock Vina best_pose format
            if "_best_pose" in str(ligand_pdbqt) and "USER" in ligand_content and "ATOM" in ligand_content:
                print("Attempting specialized best_pose format extraction")
                
                try:
                    # First, see if we can find the USER line that describes the coordinate format,
                    # followed by ATOM records with proper coordinates
                    lines = ligand_content.split('\n')
                    if len(lines) <= 1 and "USER" in ligand_content:
                        # This is likely a single-line file with USER markers instead of newlines
                        # Try to split by USER markers
                        lines = ligand_content.split('USER')
                    
                    # Look for the coordinates section
                    in_atom_section = False
                    atom_lines = []
                    
                    for line in lines:
                        # Skip header until we find coordinate section
                        line = line.strip()
                        if not in_atom_section:
                            if "x       y       z" in line or "_______ _______ _______" in line:
                                in_atom_section = True
                                continue
                        elif line.startswith(('ATOM', 'HETATM')) and "UNL" in line:
                            # Extract atom information - this is a valid ATOM/HETATM line
                            atom_lines.append(line)
                        elif line.startswith(('TORSDOF', 'ENDMDL')):
                            # End of atom section
                            break
                    
                    if atom_lines:
                        print(f"Found {len(atom_lines)} atom lines in best_pose format")
                        for line in atom_lines:
                            try:
                                # Parse atom information - structure should match the standard PDBQT format
                                # Example: ATOM      1  C   UNL     1      -2.434   2.237  -1.145 -0.29 +0.00    -0.046 C
                                # Convert to HETATM for ligand
                                parts = line.split()
                                if len(parts) >= 11:
                                    atom_num = int(parts[1])
                                    atom_type = parts[2]
                                    res_name = parts[3]
                                    res_num = parts[4]
                                    x = float(parts[5])
                                    y = float(parts[6])
                                    z = float(parts[7])
                                    
                                    # Create proper PDB format HETATM record
                                    hetatm_line = f"HETATM{atom_num:5d} {atom_type:^4s} {res_name}  {int(res_num):4d}    {x:8.3f}{y:8.3f}{z:8.3f}\n"
                                    raw_ligand_lines.append(hetatm_line)
                            except (ValueError, IndexError) as e:
                                print(f"Error parsing atom line: {line}. Error: {e}")
                except Exception as e:
                    print(f"Error during specialized best_pose parsing: {e}")
            
            # If still no ligand lines, try regex approaches
            if not raw_ligand_lines:
                # This regex looks for any pattern that resembles atomic coordinates
                import re
                
                # Try a more specific pattern that matches the format from the best_pose file
                # Example: ATOM      1  C   UNL     1      -2.434   2.237  -1.145 -0.29 +0.00    -0.046 C  
                best_pose_pattern = r'ATOM\s+(\d+)\s+(\S+)\s+UNL\s+(\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)'
                best_pose_matches = re.findall(best_pose_pattern, ligand_content)
                
                if best_pose_matches:
                    print(f"Found {len(best_pose_matches)} atoms using best_pose pattern")
                    for match in best_pose_matches:
                        atom_num, atom_type, res_num, x, y, z = match
                        try:
                            # Create a standard PDB HETATM record
                            hetatm_line = f"HETATM{int(atom_num):5d} {atom_type:^4s} UNL  {int(res_num):4d}    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}\n"
                            raw_ligand_lines.append(hetatm_line)
                        except ValueError as e:
                            print(f"Error converting coordinate match: {match}. Error: {e}")
                else:
                    # Try a more general pattern as fallback
                    print("No matches with best_pose pattern, trying general coordinate pattern")
                    coord_pattern = r'(?:ATOM|HETATM)?\s*(\d+)\s+(\S+)\s+UNL\s+\d+\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)'
                    coord_matches = re.findall(coord_pattern, ligand_content)
                    
                    if coord_matches:
                        print(f"Found {len(coord_matches)} atoms using general pattern")
                        for match in coord_matches:
                            atom_num, atom_type, x, y, z = match
                            try:
                                # Create a standard PDB HETATM record
                                hetatm_line = f"HETATM{int(atom_num):5d} {atom_type:^4s} UNL     1    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}\n"
                                raw_ligand_lines.append(hetatm_line)
                            except ValueError as e:
                                print(f"Error converting coordinate match: {match}. Error: {e}")
                    else:
                        # Last resort: try to find any 3 consecutive floating point numbers that could be coordinates
                        print("No matches with general pattern, trying last resort pattern")
                        # Find sets of 3 floats that might be x,y,z coordinates
                        float_pattern = r'(\S+)\s+(\S+)\s+UNL\s+\d+\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)'
                        float_matches = re.findall(float_pattern, ligand_content)
                        
                        if float_matches:
                            print(f"Found {len(float_matches)} potential coordinate sets")
                            for i, match in enumerate(float_matches):
                                try:
                                    atom_type, atom_num, x, y, z = match
                                    if not atom_num.isdigit():
                                        atom_num = i + 1  # Use index if not a valid number
                                    else:
                                        atom_num = int(atom_num)
                                    # Create a standard PDB HETATM record
                                    hetatm_line = f"HETATM{int(atom_num):5d} {atom_type:^4s} UNL     1    {float(x):8.3f}{float(y):8.3f}{float(z):8.3f}\n"
                                    raw_ligand_lines.append(hetatm_line)
                                except (ValueError, IndexError) as e:
                                    print(f"Error converting float coordinates: {match}. Error: {e}")
                        else:
                            print("Failed to extract any atom coordinates from file")
        
        # Now renumber the ligand atoms to continue from where protein atoms left off
        for line in raw_ligand_lines:
            if line.startswith('HETATM'):
                try:
                    # Split line into parts and extract the necessary components
                    parts = line.split()
                    if len(parts) >= 7:  # HETATM, atom_num, atom_type, residue_name, chain_id, residue_num, x, y, z
                        # Get all components
                        record_type = "HETATM"
                        atom_type = parts[2]
                        res_name = parts[3]
                        chain_id = ""  # Empty chain ID for ligands
                        res_num = parts[4] if len(parts) > 4 else "1"
                        x = float(parts[5]) if len(parts) > 5 else 0.0
                        y = float(parts[6]) if len(parts) > 6 else 0.0
                        z = float(parts[7])
                        
                        # Increment the atom number
                        last_atom_number += 1
                        
                        # Create a new line with proper PDB formatting
                        # PDB format: columns 1-6: record name, 7-11: atom serial number, 13-16: atom name,
                        # 18-20: residue name, 22: chain ID, 23-26: residue sequence number,
                        # 31-38: x coordinate, 39-46: y coordinate, 47-54: z coordinate
                        new_line = f"{record_type}{last_atom_number:5d} {atom_type:^4s} {res_name}  {int(res_num):4d}    {x:8.3f}{y:8.3f}{z:8.3f}\n"
                        ligand_lines.append(new_line)
                    else:
                        # If not enough parts, try a more permissive approach with regex
                        import re
                        match = re.match(r'^HETATM\s+\d+\s+(\S+)\s+(\S+)\s+(\S+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)', line)
                        if match:
                            atom_type, res_name, res_num, x, y, z = match.groups()
                            last_atom_number += 1
                            try:
                                # Convert coordinates to float
                                x_float = float(x)
                                y_float = float(y)
                                z_float = float(z)
                                new_line = f"HETATM{last_atom_number:5d} {atom_type:^4s} {res_name}  {int(res_num):4d}    {x_float:8.3f}{y_float:8.3f}{z_float:8.3f}\n"
                                ligand_lines.append(new_line)
                            except (ValueError, TypeError):
                                # Just use original line with new atom number if conversion fails
                                new_line = re.sub(r'^HETATM\s+\d+', f"HETATM{last_atom_number:5d}", line)
                                ligand_lines.append(new_line)
                        else:
                            # Last resort: just replace the atom number
                            match = re.match(r'^HETATM\s+\d+', line)
                            if match:
                                last_atom_number += 1
                                new_line = re.sub(r'^HETATM\s+\d+', f"HETATM{last_atom_number:5d}", line)
                                ligand_lines.append(new_line)
                            else:
                                print(f"Could not parse HETATM line format: {line}")
                                ligand_lines.append(line)
                except Exception as e:
                    print(f"Error renumbering atom: {line}. Error: {e}")
                    # If we can't parse it, just add the original line
                    ligand_lines.append(line)
        
        # Write combined file - protein followed by ligand
        with open(output_path, 'w') as f:
            f.writelines(protein_lines)
            f.writelines(ligand_lines)
            f.write('END\n')  # Add standard PDB terminator
        
        # Check if we successfully merged atoms
        if not protein_lines:
            print(f"Warning: No protein atoms found in {protein_pdbqt}")
        if not ligand_lines:
            print(f"Warning: No ligand atoms found in {ligand_pdbqt}")
        else:
            print(f"Successfully extracted {len(ligand_lines)} ligand atoms from {ligand_pdbqt}")
        
        return output_path
    
    def merge_vina_results(self, result: Union[VinaDockingResult, List[VinaDockingResult]], 
                          output_dir: Path) -> List[Path]:
        """Merge protein and ligand files from Vina docking results into complexes for PLIP analysis.
        
        Args:
            result: VinaDockingResult or list of VinaDockingResult objects
            output_dir: Directory to save merged complex files
            
        Returns:
            List of paths to merged complex files
        """
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        merged_paths = []
        
        # Handle single result or list of results
        results = [result] if not isinstance(result, list) else result
        
        for res in results:
            # Get protein and ligand files
            protein_pdbqt = res.docking_files.get('protein_pdbqt')
            
            # Try best pose first, then default ligand
            ligand_pdbqt = res.docking_files.get('best_pose', res.docking_files.get('ligand_pdbqt'))
            
            if protein_pdbqt and ligand_pdbqt and protein_pdbqt.exists() and ligand_pdbqt.exists():
                # Create output filename
                protein_name = res.protein_name.split('.')[0]
                ligand_name = res.ligand_name.split('.')[0]
                output_filename = f"{protein_name}_{ligand_name}_complex.pdb"
                
                # Merge and add to results
                merged_path = self.merge_pdbqt_complex(
                    protein_pdbqt=protein_pdbqt,
                    ligand_pdbqt=ligand_pdbqt,
                    output_dir=output_dir,
                    output_filename=output_filename
                )
                
                merged_paths.append(merged_path)
        
        return merged_paths

def clean_pdb_file(input_pdb: Union[str, Path], output_pdb: Union[str, Path] = None) -> Path:
    """Clean PDB file to contain only protein ATOM records.
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output PDB file (if None, creates a cleaned file in the same directory)
        
    Returns:
        Path to cleaned PDB file
    """
    if output_pdb is None:
        input_path = Path(input_pdb)
        output_pdb = input_path.parent / f"cleaned_{input_path.name}"
    else:
        output_pdb = Path(output_pdb)
        
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            # Only write ATOM records (protein atoms) and exclude HETATM and ANISOU
            if line.startswith('ATOM  ') and not line.startswith('ANISOU'):
                f_out.write(line)
        # Write END to properly terminate the PDB file
        f_out.write('END\n')
    
    return output_pdb
def clean_pdb_file(input_pdb: Union[str, Path], output_pdb: Union[str, Path] = None) -> Path:
    """Clean PDB file to contain only protein ATOM records.
    
    Args:
        input_pdb: Path to input PDB file
        output_pdb: Path to output PDB file (if None, creates a cleaned file in the same directory)
        
    Returns:
        Path to cleaned PDB file
    """
    if output_pdb is None:
        input_path = Path(input_pdb)
        output_pdb = input_path.parent / f"cleaned_{input_path.name}"
    else:
        output_pdb = Path(output_pdb)
        
    with open(input_pdb, 'r') as f_in, open(output_pdb, 'w') as f_out:
        for line in f_in:
            # Only write ATOM records (protein atoms) and exclude HETATM and ANISOU
            if line.startswith('ATOM  ') and not line.startswith('ANISOU'):
                f_out.write(line)
        # Write END to properly terminate the PDB file
        f_out.write('END\n')
    
    return output_pdb

    
class NvidiaGenMolAPI:
    """Client for NVIDIA NIM GenMol API."""

    def __init__(self, api_key: str = None) -> None:
        """Initialize NVIDIA Cloud GenMol API client.

        :param api_key: API key for NVIDIA Cloud API (can also be set via environment variable)
        """
        self.base_url = "https://health.api.nvidia.com/v1/biology/nvidia/genmol/generate"
        self.api_key = api_key or os.environ.get("NVIDIA_API_KEY")
        if not self.api_key:
            logger.warning("NVIDIA API key not provided. Will attempt to use environment variable at runtime.")
            
        # Create a session for connection reuse
        self.session = requests.Session()

    def _get_auth_header(self) -> dict:
        """Get the authorization headers for NVIDIA Cloud API."""
        api_key = self.api_key or os.environ.get("NVIDIA_API_KEY")
        if not api_key:
            raise ValueError("NVIDIA API key must be provided either as a parameter or via NVIDIA_API_KEY environment variable")
        return {
            "Authorization": f"Bearer {api_key}",
            "Accept": "application/json",
            "Content-Type": "application/json"
        }

    def generate(
        self,
        input_smiles: str = None,
        num_molecules: int = 30,
        temperature: float = 1.0,
        noise: float = 0.0,
        step_size: float = 1.0,
        scoring: str = "QED",
        unique_only: bool = True,
        output_dir: Union[str, Path] = None
    ) -> List[str]:
        """Generate molecules using NVIDIA NIM GenMol API.

        Args:
            input_smiles: Input SMILES string for guided generation (required by API)
            num_molecules: Number of molecules to generate
            temperature: Temperature for generation (higher = more diverse)
            noise: Randomness factor (higher = more diverse)
            step_size: Diffusion step size
            scoring: Scoring method to use ("QED" or "LogP")
            unique_only: Return only unique molecules
            output_dir: Optional directory to save generated molecules

        Returns:
            List of generated molecules as SMILES strings
        """
        # If no input SMILES is provided, use a default scaffold
        if not input_smiles:
            # Default scaffold if no SMILES provided
            input_smiles = "[*{50-100}]"
            logger.info(f"Using default scaffold SMILES: {input_smiles}")

        # Prepare payload for the API
        payload = {
            "smiles": input_smiles,
            "num_molecules": num_molecules,
            "temperature": temperature,
            "noise": noise,
            "step_size": step_size,
            "scoring": scoring,
            "unique_only": unique_only
        }

        logger.info(f"Calling NVIDIA GenMol API with parameters: {payload}")

        try:
            # Make API request
            response = self.session.post(
                self.base_url,
                headers=self._get_auth_header(),
                json=payload
            )
            response.raise_for_status()
            
            # Process response
            response_data = response.json()
            
            # Log the first part of the response for debugging
            logger.debug(f"GenMol API response: {str(response_data)[:500]}...")
            
            # Extract generated molecules
            generated_molecules = []
            
            if "molecules" in response_data:
                molecules = response_data["molecules"]
                for mol in molecules:
                    if isinstance(mol, dict) and "smiles" in mol:
                        generated_molecules.append(mol["smiles"])
                    elif isinstance(mol, str):
                        generated_molecules.append(mol)
            
            # If scores are provided, log them
            if "scores" in response_data:
                scores = response_data["scores"]
                logger.info(f"Generated {len(generated_molecules)} molecules with scores: {scores[:5]}...")
            else:
                logger.info(f"Generated {len(generated_molecules)} molecules")
            
            # Save to output directory if provided
            if output_dir and generated_molecules:
                output_dir = Path(output_dir)
                output_dir.mkdir(parents=True, exist_ok=True)
                
                # Save all molecules to a text file
                with open(output_dir / "generated_molecules.txt", "w") as f:
                    for i, smiles in enumerate(generated_molecules):
                        f.write(f"{smiles}\n")
                
                # Save individual molecules as SDF files using RDKit
                try:
                    from rdkit import Chem
                    from rdkit.Chem import AllChem
                    
                    for i, smiles in enumerate(generated_molecules):
                        try:
                            mol = Chem.MolFromSmiles(smiles)
                            if mol:
                                # Add hydrogens and generate 3D coordinates
                                mol = Chem.AddHs(mol)
                                AllChem.EmbedMolecule(mol, randomSeed=42)
                                AllChem.MMFFOptimizeMolecule(mol)
                                
                                # Write to SDF file
                                sdf_path = output_dir / f"molecule_{i}.sdf"
                                writer = Chem.SDWriter(str(sdf_path))
                                writer.write(mol)
                                writer.close()
                        except Exception as e:
                            logger.warning(f"Error creating SDF for molecule {i}: {e}")
                except ImportError:
                    logger.warning("RDKit not available, skipping SDF file creation")
            
            return generated_molecules
            
        except requests.exceptions.RequestException as e:
            logger.error(f"Error calling NVIDIA GenMol API: {e}")
            if hasattr(e, 'response') and e.response:
                logger.error(f"Response status: {e.response.status_code}")
                logger.error(f"Response body: {e.response.text}")
            raise