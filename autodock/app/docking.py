# docking.py
import subprocess
from pathlib import Path
import logging
import shutil
import os
import traceback
import re
import time
import json
import tempfile
from typing import Union, List

# Import models
from models import AutoDockConfig

logger = logging.getLogger(__name__)

def reset_gpu_state():
    """
    Minimally reset the GPU state to clear memory and resources between runs.
    
    This function performs basic GPU monitoring and stats reset without
    attempting to create or destroy CUDA contexts, which could conflict
    with AutoDock-GPU's own CUDA management.
    """
    try:
        logger.info("Performing minimal GPU state reset")
        
        # Just run nvidia-smi with stats reset
        # This doesn't interfere with existing CUDA contexts but helps
        # keep GPU monitoring accurate
        try:
            subprocess.run(
                ["nvidia-smi", "--gpu-reset-stats"],
                stdout=subprocess.DEVNULL, 
                stderr=subprocess.DEVNULL, 
                check=False,
                timeout=5
            )
            logger.info("GPU monitoring stats reset")
        except Exception as e:
            logger.debug(f"nvidia-smi reset failed: {str(e)}")
        
        # Deliberately NOT using cudaDeviceReset() or other aggressive methods
        # that could interfere with AutoDock-GPU's CUDA context management
        
        # Adding a small delay to allow any pending GPU operations to complete
        time.sleep(0.5)
        
        return True
            
    except Exception as e:
        logger.warning(f"Failed to reset GPU state: {str(e)}")
        logger.debug(traceback.format_exc())
        return False

def run_autodock_gpu(protein_pdbqt: Path, ligand_pdbqts: Union[Path, list[Path]], maps_fld: Path, output_dir: Path, config: AutoDockConfig) -> Path:
    """Run docking with AutoDock-GPU v1.5.3.
    
    Handles both single and batch ligand docking.
    
    Args:
        protein_pdbqt: Path to the prepared protein PDBQT file
        ligand_pdbqts: Path or list of paths to prepared ligand PDBQT files
        maps_fld: Path to the grid maps fld file
        output_dir: Directory to save results
        config: Configuration parameters
        
    Returns:
        Path to the output directory
    """
    try:
        # Ensure we work with absolute paths
        output_dir = output_dir.absolute()
        work_dir = output_dir / "work"
        work_dir.mkdir(exist_ok=True)
        maps_fld = maps_fld.absolute()
        protein_pdbqt = protein_pdbqt.absolute()
        
        # Handle both single ligand and list of ligands
        if not isinstance(ligand_pdbqts, list):
            # If single Path object, convert to a list
            ligand_pdbqts = [ligand_pdbqts]
        
        # Now we can safely convert all paths to absolute
        ligand_pdbqts = [p.absolute() for p in ligand_pdbqts]
        
        output_dlg = output_dir / "docking_result.dlg"
        best_pdbqt = output_dir / "best_pose.pdbqt"

        # Get docking parameters from config
        num_runs = getattr(config, 'num_runs', 10)
        seed = getattr(config, 'seed', 42)
        grid_volume = getattr(config, 'size_x', 20.0) * getattr(config, 'size_y', 20.0) * getattr(config, 'size_z', 20.0)
        max_evaluations = getattr(config, 'max_evaluations', max(2500000, int(2500000 * (grid_volume / 8000))))
        
        # Make sure max_evaluations is not 1 (a bug was occurring where it was defaulting to 1)
        if max_evaluations <= 1:
            max_evaluations = 2500000
            logger.warning(f"Invalid max_evaluations value detected, using default of {max_evaluations}")
        
        logger.info(f"Docking parameters: num_runs={num_runs}, seed={seed}, max_evals={max_evaluations}")
        logger.info(f"Config details: {config.dict() if hasattr(config, 'dict') else str(config)}")

        # Copy grid maps to work directory
        work_maps_fld = work_dir / maps_fld.name
        if maps_fld.resolve() != work_maps_fld.resolve():
            logger.info(f"Copying maps.fld from {maps_fld} to {work_maps_fld}")
            shutil.copy2(maps_fld, work_maps_fld)
        
        logger.info(f"Copying map files from {maps_fld.parent} to {work_dir}")
        map_count = 0
        for map_file in maps_fld.parent.glob("protein.*.map"):
            target_file = work_dir / map_file.name
            if map_file.resolve() != target_file.resolve():
                shutil.copy2(map_file, target_file)
                map_count += 1
        logger.info(f"Copied {map_count} map files to work directory")
        
        # Copy protein to work directory
        work_protein = work_dir / "protein.pdbqt"
        if protein_pdbqt.resolve() != work_protein.resolve():
            logger.info(f"Copying protein from {protein_pdbqt} to {work_protein}")
            shutil.copy2(protein_pdbqt, work_protein)
            
        # Copy ligands to work directory with unique names
        work_ligands = []
        for i, ligand_pdbqt in enumerate(ligand_pdbqts):
            # Use original filename if single ligand, otherwise add index to prevent conflicts
            if len(ligand_pdbqts) == 1:
                work_ligand = work_dir / ligand_pdbqt.name
            else:
                basename = ligand_pdbqt.stem
                work_ligand = work_dir / f"{basename}_{i+1}.pdbqt"
                
            if ligand_pdbqt.resolve() != work_ligand.resolve():
                logger.info(f"Copying ligand from {ligand_pdbqt} to {work_ligand}")
                shutil.copy2(ligand_pdbqt, work_ligand)
            work_ligands.append(work_ligand)

        # Locate the autodock executable
        autodock_path = "/opt/AutoDock-GPU/bin/autodock_gpu_128wi"
        if not os.path.exists(autodock_path):
            # Look for the executable in PATH
            autodock_path = shutil.which("autodock_gpu_128wi") or autodock_path
            logger.info(f"Using AutoDock-GPU executable: {autodock_path}")
            if not os.path.exists(autodock_path):
                raise FileNotFoundError(f"AutoDock-GPU executable not found at {autodock_path}")
                
        # Create base command with file arguments first - use filenames only, not full paths
        # Some command-line tools are sensitive to argument order
        if len(work_ligands) == 1:
            # Single ligand mode - add both ffile and lfile parameters at the beginning
            cmd = [autodock_path, 
                   "-ffile", work_maps_fld.name, 
                   "-lfile", work_ligands[0].name,
                   "-resnam", "docking_result", 
                   "-nrun", str(num_runs), "-seed", str(seed), "-xmloutput", "1",
                   "-dlgoutput", "1", "-dlg2stdout", "0", "-lsmet", getattr(config, 'local_search_method', 'ad'),
                   "-autostop", "1" if getattr(config, 'autostop', True) else "0", "-heuristics", str(getattr(config, 'heuristics', 1)),
                   "-nev", str(max_evaluations), "-D", "1"]
            logger.info(f"Running single ligand docking with -ffile {work_maps_fld.name} -lfile {work_ligands[0].name}")
        else:
            # Batch mode - create file list in the correct format according to AutoDock-GPU documentation
            # Format: Grid map file followed by pairs of ligand file and result name
            ligand_list_file = work_dir / "ligand_list.txt"
            with open(ligand_list_file, "w") as f:
                # First line: Grid map file with relative path (use ./ prefix)
                f.write(f"./{work_maps_fld.name}\n")
                
                # Then for each ligand: ligand file followed by result name
                for i, ligand_path in enumerate(work_ligands):
                    # Ligand file with relative path (use ./ prefix)
                    f.write(f"./{ligand_path.name}\n")
                    # Result name (must be on its own line)
                    f.write(f"Ligand{i+1}\n")
            
            # Log the contents of ligand_list.txt to verify what's being passed to AutoDock-GPU
            with open(ligand_list_file, "r") as f:
                logger.info(f"ligand_list.txt contents:\n{f.read()}")
            
            # Verify that all files in the list exist in the work directory
            if not work_maps_fld.exists():
                logger.error(f"Grid map file not found: {work_maps_fld}")
                raise FileNotFoundError(f"Grid map file not found: {work_maps_fld}")
                
            for ligand_path in work_ligands:
                if not ligand_path.exists():
                    logger.error(f"Ligand file not found: {ligand_path}")
                    raise FileNotFoundError(f"Ligand file not found: {ligand_path}")
            
            # For batch mode, create the command with file arguments first (using filenames only)
            cmd = [autodock_path,
                   "-B", ligand_list_file.name,
                   "-resnam", "docking_result", 
                   "-nrun", str(num_runs), "-seed", str(seed), "-xmloutput", "1",
                   "-dlgoutput", "1", "-dlg2stdout", "0", "-lsmet", getattr(config, 'local_search_method', 'ad'),
                   "-autostop", "1" if getattr(config, 'autostop', True) else "0", "-heuristics", str(getattr(config, 'heuristics', 1)),
                   "-nev", str(max_evaluations), "-D", "1"]
            logger.info(f"Running batch docking with -B {ligand_list_file.name} containing {len(work_ligands)} ligands")
        
        # Verify all required map files exist and are accessible
        logger.info("Verifying map files in work directory...")
        try:
            # Read maps.fld to understand required map types
            with open(work_maps_fld, 'r') as f:
                fld_content = f.read()
                logger.debug(f"Maps FLD content (first 200 chars): {fld_content[:200]}...")
                
                # Log all referenced map files
                map_refs = []
                for line in fld_content.split('\n'):
                    if ".map" in line:
                        map_refs.append(line.strip())
                logger.debug(f"Map files referenced in FLD: {map_refs}")
        except Exception as e:
            logger.warning(f"Could not verify maps.fld content: {str(e)}")
        
        # Verify actual map files on disk
        map_files = list(work_dir.glob("protein.*.map"))
        logger.info(f"Found {len(map_files)} map files in work directory:")
        for map_file in map_files:
            map_size = os.path.getsize(map_file)
            logger.info(f"  {map_file.name} ({map_size} bytes)")
            if map_size == 0:
                logger.warning(f"  WARNING: {map_file.name} has zero size!")

        # Extract atom types from all ligand PDBQT files
        ligand_atom_types = set()
        minimal_pdbqt_detected = False
        
        # Process all ligands to collect all atom types that need to be mapped
        logger.info(f"Extracting atom types from {len(work_ligands)} ligand(s)")
        for i, ligand_path in enumerate(work_ligands):
            try:
                with open(ligand_path, 'r') as f:
                    content = f.read()
                    if i == 0:  # Log only first ligand content to avoid verbose output
                        logger.debug(f"First ligand PDBQT content (first 200 chars): {content[:200]}...")
                    
                    # Check for minimal PDBQT file
                    if "Minimal PDBQT file" in content:
                        # For minimal files, explicitly add a standard type
                        ligand_atom_types.add("A")
                        minimal_pdbqt_detected = True
                        logger.warning(f"Detected minimal PDBQT file in {ligand_path.name}, using atom type 'A'")
                    
                    # Extract atom types normally
                    for line in content.splitlines():
                        if line.startswith(("ATOM", "HETATM")):
                            atom_type = line[77:79].strip()
                            if atom_type:
                                ligand_atom_types.add(atom_type)
                                
                logger.info(f"Ligand {i+1}: {ligand_path.name} - atom types found: {sorted(ligand_atom_types)}")
            except Exception as e:
                logger.warning(f"Error extracting atom types from {ligand_path.name}: {str(e)}")
        
        # Add "ATO" to the list of atom types to handle - AutoDock-GPU sometimes interprets
        # atom types this way from minimal PDBQT files
        if minimal_pdbqt_detected:
            ligand_atom_types.add("ATO")
            logger.warning("Added 'ATO' atom type for minimal PDBQT handling")
            
        # Read maps.fld to see which atom types have maps available (more reliable than directory listing)
        available_map_types = set()
        try:
            with open(work_maps_fld, 'r') as f:
                fld_content = f.read()
                for line in fld_content.split('\n'):
                    if "protein." in line and ".map" in line:
                        map_name = line.split("protein.")[1].split(".map")[0]
                        available_map_types.add(map_name)
        except Exception as e:
            logger.warning(f"Error reading maps.fld, falling back to directory listing: {str(e)}")
            available_map_types = {m.stem.split('.')[1] for m in work_dir.glob("protein.*.map")}
        
        logger.info(f"Ligand atom types: {', '.join(sorted(ligand_atom_types))}")
        logger.info(f"Available map types: {', '.join(sorted(available_map_types))}")
        
        # Check if all required atom types have map files
        missing_types = ligand_atom_types - available_map_types - {'ATO', 'H'}  # ATO and H are handled specially
        if missing_types:
            logger.warning(f"Missing map files for atom types: {', '.join(sorted(missing_types))}")
        
        # Create derivtype mappings for atom types not in available maps
        deriv_mappings = []
        for at in ligand_atom_types:
            # Handle H atom type as a special case - always map to HD
            if at == 'H':
                mapping = 'HD'
                logger.info(f"Mapping H to HD as standard practice")
                deriv_mappings.append(f"{at}={mapping}")
                continue
                
            # Skip atom types that already have maps (except special cases)
            if at in available_map_types and at != 'ATO':
                continue
                
            if at in ['MG', 'ZN', 'MN', 'CA']:
                # Metal ions map to Fe
                mapping = 'Fe'
                logger.warning(f"Mapping metal {at} to Fe")
            elif at == 'ATO':
                # Handle ATO explicitly for minimal PDBQT files
                mapping = 'A'
                logger.warning(f"Mapping ATO to generic A type")
            elif at.startswith('O') and at not in ['OA', 'O']:
                # Oxygen variants
                mapping = 'OA'
                logger.warning(f"Mapping oxygen variant {at} to OA")
            elif at.startswith('N') and at not in ['NA', 'N']:
                # Nitrogen variants
                mapping = 'N'
                logger.warning(f"Mapping nitrogen variant {at} to N")
            elif at in ['CG', 'CG0', 'G0', 'G1'] or at.startswith('CG') or at.startswith('G0') or at.startswith('G1'):
                # Special case for graphene-related atom types
                mapping = 'C'  # Map to C (carbon) as suggested by AutoDock-GPU code
                logger.warning(f"Mapping graphene-related atom type {at} to C")
            elif at.startswith('C') and at != 'C':
                # Carbon variants
                mapping = 'C'
                logger.warning(f"Mapping carbon variant {at} to C")
            elif at.startswith('S') and at not in ['SA', 'S']:
                # Sulfur variants
                mapping = 'SA'
                logger.warning(f"Mapping sulfur variant {at} to SA")
            elif at.startswith('H') and at not in ['HD', 'H']:
                # Hydrogen variants
                mapping = 'HD'
                logger.warning(f"Mapping hydrogen variant {at} to HD")
            else:
                # Skip if the type is already available
                if at in available_map_types:
                    continue
                # Generic fallback for unknown types
                mapping = 'A'
                logger.warning(f"Mapping unusual atom type {at} to generic A")
                
            deriv_mappings.append(f"{at}={mapping}")
            
        if deriv_mappings:
            # Try using individual -T flags for each mapping (instead of comma-separated list)
            # This may fix parsing issues in AutoDock-GPU
            for mapping in deriv_mappings:
                cmd.extend(["-T", mapping])
            logger.warning(f"Using derivtype mappings with individual -T flags: {deriv_mappings}")
        
        
        # Setup environment variables for CUDA
        env = os.environ.copy()
        
        # Ensure proper CUDA environment setup
        env.update({
            "CUDA_DEVICE_ORDER": "PCI_BUS_ID", 
            "CUDA_VISIBLE_DEVICES": "0",
            "NVIDIA_VISIBLE_DEVICES": "0"
        })
        
        # Set up LD_LIBRARY_PATH for CUDA libs
        cuda_lib_paths = [
            "/usr/local/cuda-11.5-compat/lib64",
            "/usr/local/cuda-11.5/lib64",
            "/usr/local/cuda/lib64"
        ]
        
        # Add all existing CUDA paths to LD_LIBRARY_PATH
        cuda_paths = []
        for path in cuda_lib_paths:
            if os.path.exists(path):
                cuda_paths.append(path)
                
        if cuda_paths:
            current_ld_path = env.get('LD_LIBRARY_PATH', '')
            new_ld_path = ':'.join(cuda_paths)
            if current_ld_path:
                new_ld_path = f"{new_ld_path}:{current_ld_path}"
            env['LD_LIBRARY_PATH'] = new_ld_path
            logger.info(f"Set LD_LIBRARY_PATH to include CUDA: {new_ld_path}")
        else:
            logger.warning(f"No CUDA library paths found, using system default")
        
        # Enhanced logging to debug H atom type issues
        logger.info(f"Ligand atom types detected: {sorted(ligand_atom_types)}")
        logger.info(f"Available map types: {sorted(available_map_types)}")
        logger.info(f"Using derivtype mappings: {deriv_mappings}")
        
        # Log the final command with extra detail
        logger.info(f"Running AutoDock-GPU from {work_dir} with command:")
        logger.info(f"  {' '.join(cmd)}")
        
        # Run the command
        with open(output_dir / "autodock_stdout.log", "w") as stdout_file, \
             open(output_dir / "autodock_stderr.log", "w") as stderr_file:
            
            # Use Popen to capture output in real-time
            process = subprocess.Popen(
                cmd, 
                stdout=subprocess.PIPE, 
                stderr=subprocess.PIPE, 
                text=True, 
                cwd=str(work_dir), 
                env=env
            )
            
            # Set up timeout
            timeout = getattr(config, 'timeout', 300)
            start_time = time.time()
            
            # Process output in real-time
            while process.poll() is None:
                # Check for timeout
                if time.time() - start_time > timeout:
                    process.terminate()
                    raise TimeoutError(f"AutoDock-GPU timed out after {timeout} seconds")
                
                # Read stdout
                stdout_line = process.stdout.readline()
                if stdout_line:
                    stdout_file.write(stdout_line)
                    stdout_file.flush()
                    if "ERROR" in stdout_line or "error" in stdout_line.lower():
                        logger.error(f"AutoDock-GPU output: {stdout_line.strip()}")
                    elif "WARNING" in stdout_line:
                        logger.warning(f"AutoDock-GPU output: {stdout_line.strip()}")
                
                # Read stderr
                stderr_line = process.stderr.readline()
                if stderr_line:
                    stderr_file.write(stderr_line)
                    stderr_file.flush()
                    logger.error(f"AutoDock-GPU error: {stderr_line.strip()}")
                
                # Avoid CPU spinning
                if not stdout_line and not stderr_line:
                    time.sleep(0.1)
            
            # Get final output
            stdout, stderr = process.communicate()
            stdout_file.write(stdout)
            stderr_file.write(stderr)
        
        # Check return code
        if process.returncode != 0:
            logger.error(f"AutoDock-GPU failed with exit code {process.returncode}")
            with open(output_dir / "autodock_stderr.log", "r") as f:
                stderr = f.read()
            
            # Save additional debug info
            with open(output_dir / "gpu_debug.log", "w") as f:
                f.write(f"Environment variables:\n")
                for k, v in env.items():
                    if "CUDA" in k or "NVIDIA" in k or "LD_LIBRARY" in k:
                        f.write(f"{k}={v}\n")
                f.write(f"\nCUDA device check:\n")
                try:
                    cuda_check = subprocess.run("nvidia-smi", capture_output=True, text=True)
                    f.write(f"nvidia-smi output:\n{cuda_check.stdout}\n")
                except:
                    f.write("Failed to run nvidia-smi\n")
                
            raise Exception(f"AutoDock-GPU failed with exit code {process.returncode}: {stderr}")

        # Copy all result files to output directory
        # In batch mode, result files will have unique names based on result names (Ligand1, Ligand2, etc.)
        if len(work_ligands) == 1:
            result_files = list(work_dir.glob("docking_result*"))
        else:
            # For batch mode with the new format, collect all Ligand* files
            result_files = []
            for i in range(1, len(work_ligands) + 1):
                ligand_results = list(work_dir.glob(f"Ligand{i}*"))
                if ligand_results:
                    result_files.extend(ligand_results)
            
            # Also add any other standard result files
            result_files.extend(list(work_dir.glob("docking_result*")))
            
            # If no Ligand* files found, try the old pattern as fallback (ligand_*)
            if not result_files:
                logger.warning("No result files found with Ligand* pattern, trying other patterns")
                result_files = list(work_dir.glob("docking_result*")) + list(work_dir.glob("ligand_*"))
        
        logger.info(f"Copying {len(result_files)} result files to output directory")
        for output_file in result_files:
            shutil.copy2(output_file, output_dir / output_file.name)

        # Process results for each ligand
        # For batch mode, we'll have multiple DLG files, one per ligand
        dlg_files = []
        
        # Look for all potential dlg file patterns
        if len(work_ligands) == 1:
            # For single ligand, look for docking_result.dlg
            ligand_dlg = list(output_dir.glob("docking_result.dlg"))
            if ligand_dlg:
                dlg_files.extend(ligand_dlg)
                logger.info(f"Found single ligand result file: {ligand_dlg[0].name}")
        else:
            # For batch mode, look for numbered files
            for i in range(1, len(work_ligands) + 1):
                ligand_dlg = list(output_dir.glob(f"docking_result_{i}.dlg"))
                if ligand_dlg:
                    dlg_files.extend(ligand_dlg)
        
        # If we still haven't found any DLG files, look for any files matching docking_result*.dlg
        if not dlg_files:
            all_dlg = list(output_dir.glob("docking_result*.dlg"))
            if all_dlg:
                dlg_files.extend(all_dlg)
                logger.info(f"Found {len(all_dlg)} DLG files with generic pattern")
        
        logger.info(f"Found total of {len(dlg_files)} DLG result files")
        # Track results for each ligand for summary reporting
        ligand_results = {}
        
        # Process each DLG file
        for dlg_file in dlg_files:
            logger.info(f"Processing results from DLG file: {dlg_file.name}")
            
            # Determine ligand name and output file for best pose
            # Extract ligand number from filename (e.g., "docking_result_1.dlg" -> "1")
            # Parse filename to extract number
            file_parts = dlg_file.stem.split('_')
            if len(file_parts) >= 3 and file_parts[-1].isdigit():
                ligand_number = int(file_parts[-1])
                ligand_name = f"ligand_{ligand_number}"
            else:
                # Fallback if format is unexpected
                ligand_name = dlg_file.stem
            
            best_pose_file = output_dir / f"{ligand_name}_best_pose.pdbqt"
            logger.info(f"Processing ligand {ligand_name} from file {dlg_file.name}")
            
            try:
                best_energy = float('inf')
                best_pose = None
                
                # Try to extract from DLG file
                with open(dlg_file, 'r') as f:
                    content = f.read()
                    # Find all models and their energies
                    models = []
                    current_model = None
                    model_lines = []
                    energy = None
                    
                    for line in content.splitlines():
                        if line.startswith("DOCKED: MODEL"):
                            if current_model is not None and energy is not None:
                                models.append((current_model, energy, model_lines))
                            current_model = line.split()[2]
                            model_lines = [line]
                            energy = None
                        elif line.startswith("DOCKED: ENDMDL"):
                            model_lines.append(line)
                            if current_model is not None and energy is not None:
                                models.append((current_model, energy, model_lines))
                                current_model = None
                                model_lines = []
                                energy = None
                        elif line.startswith("DOCKED: USER    Estimated Free Energy of Binding"):
                            energy = float(line.split("=")[1].split()[0])
                            model_lines.append(line)
                        elif current_model is not None:
                            model_lines.append(line)
                
                # Find model with lowest energy
                if models:
                    logger.info(f"Found {len(models)} docked models in DLG file for {ligand_name}")
                    best_model = min(models, key=lambda x: x[1])
                    best_energy = best_model[1]
                    best_pose = [line[8:] for line in best_model[2] if line.startswith("DOCKED:")]
                    
                    # Save best pose to PDBQT file
                    if best_pose:
                        with open(best_pose_file, 'w') as f:
                            f.writelines(best_pose)
                        logger.info(f"Best pose for {ligand_name} saved with energy {best_energy} kcal/mol")
                        
                        # Store result for summary
                        ligand_results[ligand_name] = {
                            "success": True,
                            "binding_energy": best_energy,
                            "pose_file": best_pose_file.name
                        }
                    else:
                        logger.warning(f"No pose content extracted for {ligand_name}")
                        ligand_results[ligand_name] = {
                            "success": False,
                            "error": "No pose content extracted"
                        }
                else:
                    logger.warning(f"No models found in DLG file for {ligand_name}")
                    ligand_results[ligand_name] = {
                        "success": False,
                        "error": "No docking models found"
                    }
                
                # Try to extract from XML file if DLG parsing failed
                if ligand_name not in ligand_results or not ligand_results[ligand_name].get("success", False):
                    logger.warning(f"No poses extracted from DLG for {ligand_name}, trying XML file")
                    
                    # Determine XML file name based on ligand
                    xml_file = None
                    if len(work_ligands) == 1:
                        xml_files = list(output_dir.glob("docking_result*.xml"))
                    else:
                        xml_files = list(output_dir.glob(f"{ligand_name}*.xml"))
                    
                    if xml_files:
                        xml_file = xml_files[0]
                        logger.info(f"Found XML file for {ligand_name}: {xml_file.name}")
                        try:
                            import xml.etree.ElementTree as ET
                            tree = ET.parse(xml_file)
                            root = tree.getroot()
                            lowest_energy = float('inf')
                            best_run = None
                            
                            for run in root.findall(".//run"):
                                energy_elem = run.find("free_NRG_binding")
                                if energy_elem is not None:
                                    energy = float(energy_elem.text)
                                    if energy < lowest_energy:
                                        lowest_energy = energy
                                        best_run = run
                                        
                            if best_run is not None:
                                logger.info(f"Best pose found in XML for {ligand_name} with energy {lowest_energy} kcal/mol")
                                
                                # Store result for summary
                                ligand_results[ligand_name] = {
                                    "success": True,
                                    "binding_energy": lowest_energy,
                                    "pose_file": f"{ligand_name}_from_xml.pdbqt"  # Placeholder
                                }
                        except Exception as e:
                            logger.error(f"Error parsing XML file for {ligand_name}: {str(e)}")
                            logger.error(traceback.format_exc())
                
            except Exception as e:
                logger.error(f"Error extracting best pose for {ligand_name}: {str(e)}")
                logger.error(traceback.format_exc())
                ligand_results[ligand_name] = {
                    "success": False,
                    "error": str(e)
                }
        
        if not dlg_files:
            logger.warning("No DLG files found after docking")
        
        # Create summary files with results from all ligands
        
        # JSON summary
        json_summary = {
            "method": "autodock-gpu",
            "num_ligands_processed": len(dlg_files),
            "num_ligands_submitted": len(work_ligands),
            "success": any(result.get("success", False) for result in ligand_results.values()),
            "ligands": {}
        }
        
        # Even if we didn't find DLG files, we should report what we attempted
        if len(dlg_files) == 0 and len(work_ligands) > 0:
            logger.warning(f"No DLG files found but {len(work_ligands)} ligands were submitted")
            # Add all submitted ligands to the summary with error status
            for i, ligand_path in enumerate(work_ligands):
                ligand_name = f"ligand_{i+1}"
                json_summary["ligands"][ligand_name] = {
                    "success": False,
                    "error": "No docking results found"
                }
        
        # Add detailed info for each ligand
        for ligand_name, result in ligand_results.items():
            if result.get("success", False):
                json_summary["ligands"][ligand_name] = {
                    "success": True,
                    "poses": [{
                        "rank": 1,
                        "binding_energy": round(result["binding_energy"], 2),
                        "filename": result["pose_file"]
                    }]
                }
            else:
                json_summary["ligands"][ligand_name] = {
                    "success": False,
                    "error": result.get("error", "Unknown error")
                }
        
        with open(output_dir / "summary.json", "w") as f:
            json.dump(json_summary, f, indent=2)
        
        # Text summary
        with open(output_dir / "summary.txt", "w") as f:
            f.write(f"Docking completed with AutoDock-GPU v1.5.3\n")
            f.write(f"Number of ligands: {len(work_ligands)}\n")
            f.write(f"Number of runs per ligand: {num_runs}\n\n")
            
            for ligand_name, result in ligand_results.items():
                if result.get("success", False):
                    f.write(f"{ligand_name}: Best binding energy = {result['binding_energy']:.2f} kcal/mol\n")
                else:
                    f.write(f"{ligand_name}: Failed - {result.get('error', 'Unknown error')}\n")
        
        # Verify results for each ligand
        for ligand_name in ligand_results:
            if ligand_results[ligand_name].get("success", False):
                pose_file = ligand_results[ligand_name].get("pose_file")
                if pose_file and not (output_dir / pose_file).exists():
                    logger.warning(f"Missing expected pose file for {ligand_name}: {pose_file}")
        
        # Also add the method.txt file
        with open(output_dir / "method.txt", "w") as f:
            f.write(f"autodock-gpu\n")
            f.write(f"Configuration: {{\n")
            f.write(f"  \"num_runs\": {num_runs},\n")
            f.write(f"  \"local_search_method\": \"{getattr(config, 'local_search_method', 'ad')}\",\n")
            f.write(f"  \"heuristics\": {getattr(config, 'heuristics', 1)},\n")
            f.write(f"  \"max_evaluations\": {max_evaluations},\n")
            f.write(f"  \"batch_mode\": {len(work_ligands) > 1}\n")
            f.write(f"}}")
        
        # Clean up GPU resources after docking to prevent memory leaks
        logger.info("Cleaning up GPU resources after docking")
        reset_gpu_state()
        
        # Note: We don't delete the work directory here anymore
        # It's needed by the API for post-processing and downloading results
        # The work directory will be cleaned up by the API after the results have been
        # downloaded or when the task is complete
        # 
        # try:
        #     logger.info(f"Removing work directory to free resources: {work_dir}")
        #     shutil.rmtree(work_dir)
        #     logger.info("Work directory successfully removed")
        # except Exception as cleanup_error:
        #     # Non-fatal error - log it but don't fail the run
        #     logger.warning(f"Failed to remove work directory: {str(cleanup_error)}")
        #     logger.debug(traceback.format_exc())
            
        logger.info("Docking completed successfully")
        return output_dir

    except Exception as e:
        logger.error(f"Docking failed: {str(e)}")
        logger.error(traceback.format_exc())
        
        # Try to clean up GPU resources even after an error
        try:
            logger.info("Attempting to clean up GPU resources after docking failure")
            reset_gpu_state()
            
            # We don't delete the work directory on failure either
            # It might be needed for debugging
            # logger.info(f"Removing work directory after failure to free resources: {work_dir}")
            # shutil.rmtree(work_dir)
            # logger.info("Work directory successfully removed")
        except Exception as cleanup_error:
            # Non-fatal error - log it but continue with the exception
            logger.warning(f"Failed to clean up resources after error: {str(cleanup_error)}")
            
        raise