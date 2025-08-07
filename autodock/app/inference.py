# inference.py
from pathlib import Path
import logging
import shutil
import os
from typing import Dict, Any, Optional, List
import json
import traceback
import sys

# Import project modules
from models import AutoDockConfig
from preprocessing import prepare_protein, prepare_ligand
from grid import prepare_grid_maps
from docking import run_autodock_gpu, reset_gpu_state

logger = logging.getLogger(__name__)

def run_autodock_docking(protein_content: bytes, ligand_content: bytes, config: AutoDockConfig, output_dir: Path) -> Path:
    """
    Orchestrate AutoDock-GPU docking process.
    
    Args:
        protein_content: The raw protein PDB file content as bytes
        ligand_content: The raw ligand SDF file content as bytes
        config: Configuration parameters for the docking
        output_dir: Directory where results will be stored
        
    Returns:
        Path to the output directory containing all results
    """
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        work_dir = output_dir / "work"
        work_dir.mkdir(exist_ok=True)
        
        # Save debug info
        if getattr(config, 'verbose', False):
            with open(output_dir / "debug_info.txt", "w") as f:
                f.write(f"Configuration: {config.dict()}\n")
                f.write(f"Output directory: {output_dir}\n")
                f.write(f"Work directory: {work_dir}\n")

        logger.info("Preparing protein and ligand")
        protein_pdbqt, _ = prepare_protein(protein_content, work_dir)
        ligand_pdbqt = prepare_ligand(ligand_content, work_dir)

        logger.info("Generating grid maps")
        maps_fld = prepare_grid_maps(protein_pdbqt, work_dir, config)
        
        # Reset GPU state before running docking to prevent memory leaks
        logger.info("Resetting GPU state before docking")
        reset_gpu_state()

        logger.info("Running AutoDock-GPU v1.5.3")
        result_dir = run_autodock_gpu(protein_pdbqt, ligand_pdbqt, maps_fld, output_dir, config)

        # Copy input files to output
        for file, src in [("protein.pdbqt", protein_pdbqt), ("ligand.pdbqt", ligand_pdbqt)]:
            target = output_dir / file
            if not target.exists() and src.resolve() != target.resolve():
                shutil.copy2(src, target)

        # Write method information
        with open(output_dir / "method.txt", "w") as f:
            f.write("autodock-gpu-v1.5.3\n")
            f.write(f"Configuration: {config.dict()}")
        
        # Create summary.json with standardized output format
        try:
            summary = create_result_summary(output_dir)
            with open(output_dir / "summary.json", "w") as f:
                json.dump(summary, f, indent=2)
        except Exception as e:
            logger.error(f"Failed to create summary.json: {str(e)}")
        
        # Work directory is now automatically cleaned up in the docking.py function
        # This helps prevent resource accumulation and memory issues
        logger.info("Docking completed successfully")
            
        return output_dir

    except Exception as e:
        logger.error(f"Docking process failed: {str(e)}")
        logger.error(traceback.format_exc())
        # Save error information to output directory
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            with open(output_dir / "error.txt", "w") as f:
                f.write(f"Error: {str(e)}\n")
                f.write(traceback.format_exc())
        except:
            pass
        raise

def run_batch_autodock_docking(protein_content: bytes, ligand_contents: List[bytes], config: AutoDockConfig, output_dir: Path) -> Path:
    """
    Orchestrate batch AutoDock-GPU docking process for multiple ligands.
    
    Args:
        protein_content: The raw protein PDB file content as bytes
        ligand_contents: List of raw ligand SDF file contents as bytes
        config: Configuration parameters for the docking
        output_dir: Directory where results will be stored
        
    Returns:
        Path to the output directory containing all results
    """
    try:
        # Create output directories
        output_dir.mkdir(parents=True, exist_ok=True)
        work_dir = output_dir / "work"
        work_dir.mkdir(exist_ok=True)
        
        # Save debug info if verbose mode is enabled
        if getattr(config, 'verbose', False):
            with open(output_dir / "debug_info.txt", "w") as f:
                f.write(f"Batch Docking Configuration: {config.dict()}\n")
                f.write(f"Number of ligands: {len(ligand_contents)}\n")
                f.write(f"Output directory: {output_dir}\n")
                f.write(f"Work directory: {work_dir}\n")

        # Prepare protein once for all ligands
        logger.info("Preparing protein")
        protein_pdbqt, _ = prepare_protein(protein_content, work_dir)
        
        # Generate grid maps once for all ligands
        logger.info("Generating grid maps")
        maps_fld = prepare_grid_maps(protein_pdbqt, work_dir, config)
        
        # Reset GPU state before preparing ligands to ensure a clean state
        logger.info("Resetting GPU state before batch docking")
        reset_gpu_state()
        
        # Prepare each ligand with unique names
        logger.info(f"Preparing {len(ligand_contents)} ligands")
        ligand_pdbqts = []
        for i, ligand_content in enumerate(ligand_contents):
            # Use ligand_{index} naming convention for multiple ligands
            ligand_name = f"ligand_{i+1}"
            ligand_pdbqt = prepare_ligand(ligand_content, work_dir, ligand_name)
            ligand_pdbqts.append(ligand_pdbqt)
            logger.info(f"Prepared ligand {i+1}/{len(ligand_contents)}: {ligand_pdbqt}")
        
        # Run batch docking
        logger.info(f"Running AutoDock-GPU batch docking with {len(ligand_pdbqts)} ligands")
        result_dir = run_autodock_gpu(protein_pdbqt, ligand_pdbqts, maps_fld, output_dir, config)
        
        # Copy essential input files to output directory for reference
        shutil.copy2(protein_pdbqt, output_dir / "protein.pdbqt")
        for i, ligand_pdbqt in enumerate(ligand_pdbqts):
            filename = ligand_pdbqt.name
            if not (output_dir / filename).exists():
                shutil.copy2(ligand_pdbqt, output_dir / filename)
        
        # Note: The docking module now handles result processing and summary creation
        
        # Work directory is now automatically cleaned up in the docking.py function
        # This helps prevent resource accumulation and memory issues
        logger.info("Batch docking completed successfully")
            
        return output_dir
        
    except Exception as e:
        logger.error(f"Batch docking process failed: {str(e)}")
        logger.error(traceback.format_exc())
        # Save error information to output directory
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            with open(output_dir / "error.txt", "w") as f:
                f.write(f"Error: {str(e)}\n")
                f.write(traceback.format_exc())
        except:
            pass
        raise

def run_sequential_autodock_docking(protein_content: bytes, ligand_contents: List[bytes], config: AutoDockConfig, output_dir: Path) -> Path:
    """
    Process multiple ligands sequentially using single ligand mode.
    
    Args:
        protein_content: The raw protein PDB file content as bytes
        ligand_contents: List of raw ligand SDF file contents as bytes
        config: Configuration parameters for the docking
        output_dir: Directory where results will be stored
        
    Returns:
        Path to the output directory containing all results
    """
    try:
        # Create output directories
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Save debug info if verbose mode is enabled
        if getattr(config, 'verbose', False):
            with open(output_dir / "debug_info.txt", "w") as f:
                f.write(f"Sequential Docking Configuration: {config.dict()}\n")
                f.write(f"Number of ligands: {len(ligand_contents)}\n")
                f.write(f"Output directory: {output_dir}\n")
        
        # Process each ligand individually
        ligand_results = {}
        for i, ligand_content in enumerate(ligand_contents):
            # Create a subfolder for each ligand
            ligand_name = f"ligand_{i+1}"
            ligand_dir = output_dir / ligand_name
            
            logger.info(f"Processing ligand {i+1}/{len(ligand_contents)}: {ligand_name}")
            
            # Reset GPU state before each sequential ligand processing
            if i > 0:  # No need to reset before the first one (already reset in run_autodock_docking)
                logger.info(f"Resetting GPU state before processing ligand {i+1}")
                reset_gpu_state()
            
            try:
                # Process this individual ligand
                result_dir = run_autodock_docking(
                    protein_content,
                    ligand_content,
                    config,
                    ligand_dir
                )
                
                # Copy ligand results to main output directory with unique names
                if result_dir and result_dir.exists():
                    # Extract binding energy from the result summary
                    ligand_summary = create_result_summary(result_dir)
                    best_energy = None
                    if ligand_summary["poses"]:
                        best_energy = ligand_summary["poses"][0].get("binding_energy")
                    
                    # Copy best pose with renamed file
                    best_pose_file = result_dir / "best_pose.pdbqt"
                    if best_pose_file.exists():
                        ligand_pose_file = output_dir / f"{ligand_name}_best_pose.pdbqt"
                        shutil.copy2(best_pose_file, ligand_pose_file)
                    
                    # Copy docking result dlg with renamed file
                    dlg_file = result_dir / "docking_result.dlg"
                    if dlg_file.exists():
                        ligand_dlg_file = output_dir / f"{ligand_name}_docking_result.dlg"
                        shutil.copy2(dlg_file, ligand_dlg_file)
                    
                    # Copy any XML results
                    xml_files = list(result_dir.glob("*.xml"))
                    for xml_file in xml_files:
                        shutil.copy2(xml_file, output_dir / f"{ligand_name}_{xml_file.name}")
                    
                    # Store ligand result info
                    ligand_results[ligand_name] = {
                        "success": True,
                        "binding_energy": best_energy,
                        "output_dir": str(result_dir),
                        "pose_file": f"{ligand_name}_best_pose.pdbqt" if best_pose_file.exists() else None
                    }
                    logger.info(f"Successfully processed {ligand_name} with energy {best_energy}")
                else:
                    ligand_results[ligand_name] = {
                        "success": False,
                        "error": "Docking failed to produce results"
                    }
                    logger.warning(f"Failed to get results for {ligand_name}")
            
            except Exception as e:
                ligand_results[ligand_name] = {
                    "success": False,
                    "error": str(e)
                }
                logger.error(f"Error processing {ligand_name}: {str(e)}")
                logger.error(traceback.format_exc())
        
        # Create a ligand list file
        with open(output_dir / "ligand_list.txt", "w") as f:
            for ligand_name in ligand_results.keys():
                f.write(f"{ligand_name}\n")
        
        # Create combined JSON summary file
        json_summary = {
            "method": "autodock-gpu-sequential",
            "num_ligands_processed": len(ligand_contents),
            "success": any(result.get("success", False) for result in ligand_results.values()),
            "ligands": {}
        }
        
        # Add detailed info for each ligand
        for ligand_name, result in ligand_results.items():
            if result.get("success", False):
                json_summary["ligands"][ligand_name] = {
                    "success": True,
                    "poses": [{
                        "rank": 1,
                        "binding_energy": result.get("binding_energy"),
                        "filename": result.get("pose_file")
                    }]
                }
            else:
                json_summary["ligands"][ligand_name] = {
                    "success": False,
                    "error": result.get("error", "Unknown error")
                }
        
        with open(output_dir / "summary.json", "w") as f:
            json.dump(json_summary, f, indent=2)
        
        # Create combined text summary
        with open(output_dir / "summary.txt", "w") as f:
            f.write(f"Sequential docking completed with AutoDock-GPU v1.5.3\n")
            f.write(f"Number of ligands: {len(ligand_contents)}\n")
            f.write(f"Number of runs per ligand: {getattr(config, 'num_runs', 10)}\n\n")
            
            for ligand_name, result in ligand_results.items():
                if result.get("success", False):
                    energy = result.get("binding_energy")
                    if energy is not None:
                        f.write(f"{ligand_name}: Best binding energy = {energy:.2f} kcal/mol\n")
                    else:
                        f.write(f"{ligand_name}: Completed, but energy information not available\n")
                else:
                    f.write(f"{ligand_name}: Failed - {result.get('error', 'Unknown error')}\n")
        
        # Create method.txt file
        with open(output_dir / "method.txt", "w") as f:
            f.write(f"autodock-gpu-sequential\n")
            f.write(f"Multiple ligands processed individually in single mode\n")
            f.write(f"Configuration: {config.dict()}")
        
        # Work directories are now automatically cleaned up in the docking.py function
        # This helps prevent resource accumulation and memory issues
        logger.info(f"Sequential docking completed for {len(ligand_contents)} ligands")
        
        return output_dir
    
    except Exception as e:
        logger.error(f"Sequential docking process failed: {str(e)}")
        logger.error(traceback.format_exc())
        # Save error information to output directory
        try:
            output_dir.mkdir(parents=True, exist_ok=True)
            with open(output_dir / "error.txt", "w") as f:
                f.write(f"Error: {str(e)}\n")
                f.write(traceback.format_exc())
        except:
            pass
        raise

def create_result_summary(output_dir: Path) -> Dict[str, Any]:
    """Create a standardized summary of docking results."""
    summary = {
        "method": "autodock-gpu-v1.5.3",
        "success": True,
        "poses": []
    }
    
    # Try to extract binding energy information
    best_energy = None
    
    # Check for summary.txt file
    summary_file = output_dir / "summary.txt"
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            for line in f:
                if "Best binding energy:" in line:
                    try:
                        best_energy = float(line.split(":")[1].strip().split()[0])
                        break
                    except:
                        pass
    
    # If no summary.txt, try to extract from dlg file
    if best_energy is None:
        dlg_file = output_dir / "docking_result.dlg"
        if dlg_file.exists():
            try:
                with open(dlg_file, 'r') as f:
                    for line in f:
                        if "DOCKED: USER    Estimated Free Energy of Binding" in line:
                            try:
                                best_energy = float(line.split("=")[1].split()[0])
                                break
                            except:
                                pass
            except:
                pass
    
    # Check for pose files in the following order:
    # 1. Legacy format: best_pose.pdbqt
    # 2. Base name with best_pose: docking_result_best_pose.pdbqt
    # 3. Ligand-specific format: ligand_X_best_pose.pdbqt
    # 4. Main docking result file: docking_result.pdbqt
    
    # 1. Check legacy format first
    best_pose_file = output_dir / "best_pose.pdbqt"
    if best_pose_file.exists():
        summary["poses"].append({
            "rank": 1,
            "binding_energy": best_energy,
            "filename": "best_pose.pdbqt"
        })
    else:
        # 2. Check for docking_result_best_pose.pdbqt
        best_pose_file = output_dir / "docking_result_best_pose.pdbqt"
        if best_pose_file.exists():
            summary["poses"].append({
                "rank": 1,
                "binding_energy": best_energy,
                "filename": "docking_result_best_pose.pdbqt"
            })
        else:
            # 3. Check for ligand specific formats (ligand_X_best_pose.pdbqt)
            ligand_pose_files = list(output_dir.glob("ligand_*_best_pose.pdbqt"))
            if ligand_pose_files:
                # Use the first one found (typically ligand_1_best_pose.pdbqt for single ligand)
                pose_file = ligand_pose_files[0]
                summary["poses"].append({
                    "rank": 1,
                    "binding_energy": best_energy,
                    "filename": pose_file.name
                })
            else:
                # 4. Finally check for the main docking result file
                result_file = output_dir / "docking_result.pdbqt"
                if result_file.exists():
                    summary["poses"].append({
                        "rank": 1,
                        "binding_energy": best_energy,
                        "filename": "docking_result.pdbqt"
                    })
    
    return summary

def main():
    """
    Main function for testing purposes.
    Usage example:
    ```
    from inference import main
    main()
    ```
    """
    import argparse
    parser = argparse.ArgumentParser(description='Run AutoDock docking')
    parser.add_argument('--protein', type=str, required=True, help='Path to protein PDB file')
    parser.add_argument('--ligands', type=str, nargs='+', required=True, help='Path(s) to ligand SDF file(s)')
    parser.add_argument('--output', type=str, required=True, help='Path to output directory')
    parser.add_argument('--num-runs', type=int, default=10, help='Number of docking runs')
    parser.add_argument('--center-x', type=float, default=0.0, help='X coordinate of search box center')
    parser.add_argument('--center-y', type=float, default=0.0, help='Y coordinate of search box center')
    parser.add_argument('--center-z', type=float, default=0.0, help='Z coordinate of search box center')
    parser.add_argument('--size-x', type=float, default=20.0, help='Width of search box in X dimension')
    parser.add_argument('--size-y', type=float, default=20.0, help='Height of search box in Y dimension')
    parser.add_argument('--size-z', type=float, default=20.0, help='Depth of search box in Z dimension')
    parser.add_argument('--timeout', type=int, default=300, help='Timeout in seconds')
    parser.add_argument('--verbose', action='store_true', help='Enable verbose logging')
    
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Load input files
    with open(args.protein, 'rb') as f:
        protein_content = f.read()
    
    # Load all ligand files
    ligand_contents = []
    for ligand_path in args.ligands:
        with open(ligand_path, 'rb') as f:
            ligand_contents.append(f.read())
            
    # Create config from args
    config = AutoDockConfig(
        num_runs=args.num_runs,
        center_x=args.center_x,
        center_y=args.center_y,
        center_z=args.center_z,
        size_x=args.size_x,
        size_y=args.size_y,
        size_z=args.size_z,
        timeout=args.timeout,
        verbose=args.verbose
    )
    
    # Run docking
    output_dir = Path(args.output)
    
    # Choose between single and batch mode based on number of ligands
    if len(ligand_contents) == 1:
        print(f"Running single ligand docking")
        result_dir = run_autodock_docking(protein_content, ligand_contents[0], config, output_dir)
    else:
        # Using native batch mode for multiple ligands
        print(f"Running batch docking with {len(ligand_contents)} ligands")
        result_dir = run_batch_autodock_docking(protein_content, ligand_contents, config, output_dir)
        
        # Sequential processing alternative - commented out but preserved for fallback
        # print(f"Running sequential docking with {len(ligand_contents)} ligands")
        # result_dir = run_sequential_autodock_docking(protein_content, ligand_contents, config, output_dir)
    
    print(f"Docking completed. Results saved to: {result_dir}")
    
    # Print summary
    summary_file = result_dir / "summary.json"
    if summary_file.exists():
        with open(summary_file, 'r') as f:
            summary = json.load(f)
            print("\nDocking Summary:")
            print(f"Method: {summary['method']}")
            print(f"Success: {summary['success']}")
            
            # Handle both single and batch mode summary formats
            if 'ligands' in summary:
                # Batch mode summary
                print(f"Number of ligands: {summary.get('num_ligands_processed', 'N/A')}")
                print("\nLigand results:")
                for ligand_name, ligand_data in summary['ligands'].items():
                    if ligand_data.get('success', False):
                        pose = ligand_data['poses'][0]
                        print(f"  {ligand_name}: {pose.get('binding_energy', 'N/A')} kcal/mol ({pose.get('filename', 'N/A')})")
                    else:
                        print(f"  {ligand_name}: Failed - {ligand_data.get('error', 'Unknown error')}")
            elif 'poses' in summary:
                # Single ligand summary
                if summary['poses']:
                    print(f"Best binding energy: {summary['poses'][0].get('binding_energy', 'N/A')} kcal/mol")
                    print(f"Best pose file: {summary['poses'][0].get('filename', 'N/A')}")

if __name__ == "__main__":
    main()