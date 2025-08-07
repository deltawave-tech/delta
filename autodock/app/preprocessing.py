#!/usr/bin/env python3
# preprocessing.py
import subprocess
from pathlib import Path
import logging
import shutil
from typing import Tuple
import os
import traceback
import re

# Handle external dependencies
from meeko.preparation import PDBQTWriterLegacy, MoleculePreparation
from rdkit import Chem
from rdkit.Chem import AllChem

logger = logging.getLogger(__name__)

def validate_protein_size(pdbqt_path: Path) -> bool:
    """Check if protein size is within safe limits."""
    atom_count = 0
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_count += 1
    return atom_count < 50000  # Increased limit for large proteins

def validate_ligand_size(pdbqt_path: Path) -> bool:
    """Check if ligand size and rotatable bonds are within limits."""
    atom_count = 0
    torsion_count = 0
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom_count += 1
            elif line.startswith('TORSDOF'):
                torsion_count = int(line.split()[1])
    return atom_count < 1000 and torsion_count < 32

def prepare_protein(pdb_content: bytes, work_dir: Path) -> Tuple[Path, Path]:
    """Convert PDB to PDBQT using Meeko and MGLTools."""
    try:
        work_dir = work_dir.absolute()
        pdb_path = work_dir / "protein.pdb"
        pdbqt_path = work_dir / "protein.pdbqt"
        
        logger.info(f"Saving protein PDB to {pdb_path}")
        with open(pdb_path, "wb") as f:
            f.write(pdb_content)

        # Clean PDB file to contain only relevant ATOM/HETATM records
        clean_pdb_path = work_dir / "clean_protein.pdb"
        logger.info(f"Cleaning PDB file to {clean_pdb_path}")
        with open(pdb_path, 'r') as f, open(clean_pdb_path, 'w') as out:
            for line in f:
                if line.startswith(('ATOM', 'HETATM', 'TER', 'END')):
                    out.write(line)
        
        # Use the cleaned PDB if it's not empty, otherwise use original
        if clean_pdb_path.stat().st_size > 0:
            input_pdb = clean_pdb_path
            logger.info(f"Using cleaned PDB file ({clean_pdb_path.stat().st_size} bytes)")
        else:
            input_pdb = pdb_path
            logger.warning(f"Cleaned PDB file is empty, using original ({pdb_path.stat().st_size} bytes)")

        # First try MGLTools for preparation (the preferred method)
        try:
            logger.info("Preparing protein with MGLTools prepare_receptor4.py")
            
            # Try multiple possible locations for MGLTools
            mgltools_paths = [
                ("/opt/mgltools/bin/pythonsh", "/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"),
                ("/usr/local/mgltools/bin/pythonsh", "/usr/local/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py"),
                (shutil.which("pythonsh"), shutil.which("prepare_receptor4.py"))
            ]
            
            # Find the first valid path
            pythonsh_path = None
            prepare_receptor_path = None
            for p_path, pr_path in mgltools_paths:
                if p_path and pr_path and os.path.exists(p_path) and os.path.exists(pr_path):
                    pythonsh_path = p_path
                    prepare_receptor_path = pr_path
                    break
            
            if not pythonsh_path or not prepare_receptor_path:
                logger.warning("Could not find MGLTools executables, will try alternative methods")
                raise FileNotFoundError("MGLTools not found")
                
            # Run MGLTools
            cmd = [
                pythonsh_path,
                prepare_receptor_path,
                "-r", str(input_pdb), 
                "-o", str(pdbqt_path), 
                "-A", "hydrogens"
            ]
            logger.info(f"Running MGLTools: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            logger.info(f"MGLTools output: {result.stdout}")
            
            if not pdbqt_path.exists() or os.path.getsize(pdbqt_path) == 0:
                raise Exception("MGLTools completed but produced no output")
                
        except Exception as e:
            logger.warning(f"MGLTools preparation failed: {str(e)}")
            logger.warning("Falling back to Meeko for protein preparation")
            
            try:
                # Use Meeko/RDKit as fallback
                logger.info(f"Preparing protein with Meeko/RDKit from {input_pdb}")
                mol = Chem.MolFromPDBFile(str(input_pdb), removeHs=False)
                if mol is None:
                    raise Exception(f"RDKit could not read {input_pdb}")
                
                logger.info("Adding hydrogens with RDKit")
                mol = Chem.AddHs(mol, addCoords=True)
                Chem.SanitizeMol(mol)
                
                logger.info("Preparing molecule with Meeko")
                preparator = MoleculePreparation()
                preparator.prepare(mol)
                # Use the write_pdbqt_string method directly from the preparator
                # instead of trying to access preparator.mol_info which might not exist in newer versions
                logger.info(f"Writing Meeko output to {pdbqt_path}")
                with open(pdbqt_path, 'w') as f:
                    f.write(preparator.write_pdbqt_string())
                    
                if not pdbqt_path.exists() or os.path.getsize(pdbqt_path) == 0:
                    raise Exception("Meeko preparation completed but produced no output")
            
            except Exception as e:
                logger.error(f"Meeko fallback also failed: {str(e)}")
                logger.error(traceback.format_exc())
                raise Exception(f"All protein preparation methods failed: {str(e)}")

        # Fix unknown atom types in the PDBQT file
        logger.info("Fixing unknown atom types in PDBQT file")
        with open(pdbqt_path, 'r') as f:
            lines = f.readlines()
            
        fixed_lines = []
        unknown_types_count = 0
        
        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and "____" in line:
                unknown_types_count += 1
                atom_name = line[12:16].strip()
                element = atom_name[:2].upper() if atom_name[:2] in ["CL", "BR"] else atom_name[0]
                
                # Map element to appropriate atom type
                atom_type = {'C': 'C', 'N': 'N', 'O': 'OA', 'S': 'SA', 'H': 'HD', 'P': 'P',
                            'F': 'F', 'CL': 'CL', 'BR': 'BR', 'I': 'I'}.get(element, 'A')
                            
                # Extract and preserve charge
                try:
                    charge = float(line[70:76].strip() or 0.0)
                except ValueError:
                    charge = 0.0
                    
                # Create fixed line with proper atom type
                fixed_lines.append(f"{line[:70]}{charge:6.3f} {atom_type:<2s}\n")
            else:
                fixed_lines.append(line)
                
        if unknown_types_count > 0:
            logger.info(f"Fixed {unknown_types_count} unknown atom types")
            
            # Write fixed file
            with open(pdbqt_path, 'w') as f:
                f.writelines(fixed_lines)

        # Validate protein size
        if not validate_protein_size(pdbqt_path):
            logger.error("Protein exceeds safe size limits")
            raise Exception("Protein exceeds safe size limits (>50,000 atoms)")
            
        # Count atoms for logging
        atom_count = 0
        with open(pdbqt_path, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_count += 1
                    
        logger.info(f"Protein prepared at {pdbqt_path} with {atom_count} atoms")
        return pdbqt_path, work_dir / "dummy.maps.fld"  # Dummy maps_fld

    except Exception as e:
        logger.error(f"Protein preparation failed: {str(e)}")
        logger.error(traceback.format_exc())
        raise

def prepare_ligand(sdf_content: bytes, work_dir: Path, ligand_name: str = "ligand") -> Path:
    """Convert SDF to PDBQT using RDKit and Meeko's Python API.

    Args:
        sdf_content: Raw bytes of the SDF file
        work_dir: Directory to save files to
        ligand_name: Base name for the ligand files (default: "ligand")

    Returns:
        Path to the prepared PDBQT file

    Raises:
        Exception: If preparation fails at any step
    """
    try:
        work_dir = work_dir.absolute()
        sdf_path = work_dir / f"{ligand_name}.sdf"
        pdbqt_path = work_dir / f"{ligand_name}.pdbqt"

        # Save original SDF
        logger.info(f"Saving ligand SDF to {sdf_path}")
        with open(sdf_path, "wb") as f:
            f.write(sdf_content)

        # Validate and preprocess with RDKit
        logger.info(f"Reading SDF with RDKit from {sdf_path}")
        supplier = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        mol = next(supplier, None)
        if mol is None:
            logger.error(f"RDKit failed to parse SDF file {sdf_path}")
            with open(sdf_path, 'r') as f:
                logger.debug(f"SDF content (first 500 chars): {f.read(500)}")
            raise ValueError("Invalid SDF file: RDKit could not parse")

        logger.info("Adding hydrogens with RDKit")
        mol = Chem.AddHs(mol, addCoords=True)

        if mol.GetNumConformers() == 0:
            logger.info("Generating 3D coordinates with RDKit")
            AllChem.EmbedMolecule(mol, useRandomCoords=True)
            AllChem.UFFOptimizeMolecule(mol)

        # Convert to PDBQT with Meeko
        logger.info("Preparing ligand with Meeko")
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        # Use the write_pdbqt_string method directly from the preparator
        # instead of trying to access mol_info which might not exist in newer versions
        
        logger.info(f"Writing PDBQT to {pdbqt_path}")
        with open(pdbqt_path, 'w') as f:
            f.write(preparator.write_pdbqt_string())

        if not pdbqt_path.exists() or pdbqt_path.stat().st_size == 0:
            raise Exception("Meeko produced no output")

        # Validate size and log atom count
        atom_count = 0
        with open(pdbqt_path, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_count += 1
        if not validate_ligand_size(pdbqt_path):
            logger.warning("Ligand may be too large or have too many rotatable bonds")
        logger.info(f"Ligand prepared at {pdbqt_path} with {atom_count} atoms")

        return pdbqt_path

    except Exception as e:
        logger.error(f"Ligand preparation failed: {str(e)}")
        logger.error(traceback.format_exc())
        raise