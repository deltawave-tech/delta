# grid.py
import subprocess
from pathlib import Path
import logging
import shutil
import os
import traceback
from typing import Tuple

# Import models
from models import AutoDockConfig

logger = logging.getLogger(__name__)

def calculate_protein_center_and_size(protein_pdbqt: Path) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Calculate protein geometric center and size."""
    try:
        coords = {'x': [], 'y': [], 'z': []}
        with open(protein_pdbqt, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    coords['x'].append(float(line[30:38].strip()))
                    coords['y'].append(float(line[38:46].strip()))
                    coords['z'].append(float(line[46:54].strip()))
        if not coords['x']:
            logger.warning("No coordinates found, using defaults")
            return (0.0, 0.0, 0.0), (40.0, 40.0, 40.0)
        center = tuple((min(coords[d]) + max(coords[d])) / 2 for d in ['x', 'y', 'z'])
        buffer = 10.0
        size = tuple(max(coords[d]) - min(coords[d]) + buffer for d in ['x', 'y', 'z'])
        logger.info(f"Center: {center}, Size: {size}")
        return center, size
    except Exception as e:
        logger.error(f"Center/size calculation failed: {str(e)}")
        return (0.0, 0.0, 0.0), (40.0, 40.0, 40.0)

def run_autogrid(protein_pdbqt: Path, work_dir: Path, config: AutoDockConfig) -> Path:
    """Generate grid maps using AutoGrid."""
    try:
        # Ensure we have absolute paths for reliability
        work_dir = work_dir.absolute()
        protein_pdbqt = protein_pdbqt.absolute()
        
        # Calculate grid parameters
        default_center = all(getattr(config, f'center_{d}', 0.0) == 0.0 for d in ['x', 'y', 'z'])
        default_size = all(getattr(config, f'size_{d}', 40.0) == 40.0 for d in ['x', 'y', 'z'])
        center, size = calculate_protein_center_and_size(protein_pdbqt) if default_center or default_size else (
            (config.center_x, config.center_y, config.center_z),
            (config.size_x, config.size_y, config.size_z)
        )
        spacing = 0.375
        npts = tuple(int((s / spacing) // 2 * 2) for s in size)
        logger.info(f"Grid center: {center}, Points: {npts}")

        # Determine atom types
        receptor_atom_types = set()
        with open(protein_pdbqt, 'r') as f:
            for line in f:
                if line.startswith(("ATOM", "HETATM")):
                    atom_type = line[77:79].strip()
                    if atom_type:
                        receptor_atom_types.add(atom_type)
        logger.info(f"Receptor atom types: {', '.join(sorted(receptor_atom_types))}")
        
        # Standard ligand types that should always be included
        standard_ligand_types = {'A', 'C', 'N', 'HD', 'OA', 'SA', 'NA', 'HS', 'F', 'Cl', 'Br', 'I', 'P', 'S'}
        logger.info(f"Standard ligand types: {', '.join(sorted(standard_ligand_types))}")

        # Copy protein to work directory if needed
        work_protein = work_dir / "protein.pdbqt"
        if protein_pdbqt.resolve() != work_protein.resolve():
            logger.info(f"Copying protein from {protein_pdbqt} to {work_protein}")
            shutil.copy2(protein_pdbqt, work_protein)
        
        # Create GPF file
        gpf_path = work_dir / "protein.gpf"
        maps_fld_path = work_dir / "protein.maps.fld"
        logger.info(f"Creating GPF file at {gpf_path}")
        with open(gpf_path, "w") as f:
            f.write(f"npts {npts[0]} {npts[1]} {npts[2]}\n")
            f.write(f"gridfld {maps_fld_path}\n")
            f.write(f"spacing {spacing}\n")
            f.write(f"receptor_types {' '.join(sorted(receptor_atom_types))}\n")
            f.write(f"ligand_types {' '.join(sorted(standard_ligand_types))}\n")
            f.write(f"receptor {work_protein}\n")
            f.write(f"gridcenter {center[0]} {center[1]} {center[2]}\n")
            f.write(f"smooth 0.5\n")
            for atom_type in sorted(standard_ligand_types):
                f.write(f"map {work_dir / f'protein.{atom_type}.map'}\n")
            f.write(f"elecmap {work_dir / 'protein.e.map'}\n")
            f.write(f"dsolvmap {work_dir / 'protein.d.map'}\n")
            f.write(f"dielectric -0.1465\n")

        # Run AutoGrid
        autogrid_path = "/usr/local/bin/autogrid4"
        if not os.path.exists(autogrid_path):
            autogrid_path = shutil.which("autogrid4") or "/opt/mgltools/bin/autogrid4"
            if not os.path.exists(autogrid_path):
                raise FileNotFoundError(f"AutoGrid executable not found at {autogrid_path}")
                
        cmd = [autogrid_path, "-p", str(gpf_path), "-l", str(work_dir / "autogrid.log")]
        logger.info(f"Running AutoGrid: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, cwd=str(work_dir))
        logger.debug(f"AutoGrid stdout: {result.stdout}")
        
        # Verify map files were created
        maps_fld = work_dir / "protein.maps.fld"
        if not maps_fld.exists():
            raise FileNotFoundError(f"Grid map generation failed at {maps_fld}")
            
        # Log each generated map file for debugging
        map_files = list(work_dir.glob("protein.*.map"))
        logger.info(f"Generated {len(map_files)} map files:")
        for map_file in map_files:
            logger.info(f"  {map_file.name} ({os.path.getsize(map_file)} bytes)")
            
        return maps_fld

    except Exception as e:
        logger.error(f"AutoGrid failed: {str(e)}")
        logger.error(traceback.format_exc())
        raise

def prepare_grid_maps(protein_pdbqt: Path, work_dir: Path, config: AutoDockConfig) -> Path:
    """Prepare grid maps for docking.
    
    Maps are generated in the work directory and used directly from there.
    """
    try:
        logger.info(f"Generating grid maps for {protein_pdbqt}")
        maps_fld = run_autogrid(protein_pdbqt, work_dir, config)
        
        # Verify map files were created successfully
        map_files = list(work_dir.glob("protein.*.map"))
        logger.info(f"Generated {len(map_files)} map files in work directory")
        
        return maps_fld
    except Exception as e:
        logger.error(f"Grid map preparation failed: {str(e)}")
        raise