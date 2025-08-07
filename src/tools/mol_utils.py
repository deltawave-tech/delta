import rdkit
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
from pathlib import Path
import base64
import io
from typing import Dict, List, Optional, Union, Tuple
import re

def merge_ligand_protein(ligand_input, protein_pdb_path, output_path):
    """
    Merge a ligand (from SDF or SMILES) with a protein PDB file.
    
    Args:
        ligand_input (str): Either path to SDF file or SMILES string
        protein_pdb_path (str): Path to protein PDB file
        output_path (str): Path for output merged PDB file
    """
    # Step 1: Read the ligand (detect if it's a SMILES or SDF)
    ligand = None
    
    # Check if input is a file path that exists
    if Path(ligand_input).exists():
        # Try reading as SDF
        try:
            ligand = Chem.SDMolSupplier(ligand_input)[0]
        except:
            raise ValueError("Could not read ligand SDF file")
    else:
        # Try reading as SMILES
        try:
            ligand = Chem.MolFromSmiles(ligand_input)
            # Generate 3D coordinates for SMILES input
            ligand = Chem.AddHs(ligand)  # Add hydrogens
            success = Chem.EmbedMolecule(ligand, randomSeed=42)  # Generate 3D coordinates
            if success == -1:
                raise ValueError("Could not generate 3D coordinates for SMILES")
            Chem.MMFFOptimizeMolecule(ligand)  # Energy minimize
        except:
            raise ValueError("Could not process SMILES input")
    
    if ligand is None:
        raise ValueError("Could not process ligand input")
    
    # Step 2: Convert ligand to PDB format
    ligand_pdb = Chem.MolToPDBBlock(ligand)
    
    # Step 3: Read the protein PDB file
    try:
        with open(protein_pdb_path, 'r') as f:
            protein_pdb = f.read()
    except:
        raise ValueError("Could not read protein PDB file")
    
    # Step 4: Merge the two PDB contents
    protein_pdb = protein_pdb.replace('END\n', '')
    merged_pdb = protein_pdb + ligand_pdb + "END\n"
    
    # Step 5: Write the merged PDB file
    try:
        with open(output_path, 'w') as f:
            f.write(merged_pdb)
        print(f"Successfully created merged PDB file at: {output_path}")
    except:
        raise ValueError(f"Could not write to output path: {output_path}")

def merge_pdb_sdf_rdkit(sdf_file, pdb_file, output_file):
    # Read PDB file
    print(f"PDB file: {str(pdb_file)}")
    protein = Chem.MolFromPDBFile(str(pdb_file))
    print(f"Protein: {protein}")
    # Read SDF file
    print(f"SDF file: {str(sdf_file)}")
    ligand = Chem.SDMolSupplier(str(sdf_file))[0]
    print(f"Ligand: {ligand}")
    # Combine the molecules
    combined = Chem.CombineMols(protein, ligand)
    print(f"Combined: {combined}")
    # Write output as PDB
    writer = Chem.PDBWriter(str(output_file))
    writer.write(combined)
    writer.close()

def get_largest_fragment(mol: Union[Chem.Mol, str, List[Chem.Mol], List[str]]) -> Union[str, List[str]]:
    """
    Get the largest fragment from a molecule and return the SMILES of the largest fragment.
    
    Parameters:
    -----------
    mol : Union[Chem.Mol, str, List[Chem.Mol], List[str]]
        Input can be a single SMILES string, single RDKit mol object, 
        list of SMILES strings, or list of RDKit mol objects.
        
    Returns:
    --------
    Union[str, List[str]]
        SMILES representation of the largest fragment(s).
    """
    # Handle list input
    if isinstance(mol, list):
        return [get_largest_fragment(m) for m in mol]
    
    # Convert SMILES to mol if needed
    if isinstance(mol, str):
        mol = Chem.MolFromSmiles(mol)
    
    # Check for valid molecule
    if mol is None:
        raise ValueError("Invalid molecule input")
    
    # Get fragments and find the largest
    frags = Chem.GetMolFrags(mol, asMols=True)
    if not frags:
        return Chem.MolToSmiles(mol)  # Return original if no fragments
    
    largest_frag = max(frags, key=lambda x: x.GetNumAtoms())
    return Chem.MolToSmiles(largest_frag)

def generate_mol_image(smiles: str, size: Tuple[int, int] = (300, 300)) -> Dict[str, str]:
    """Generate image data for a molecule in Claude's vision API format.
    
    Args:
        smiles: SMILES string of the molecule
        size: Tuple of (width, height) for the image. Default 1000x1000 for optimal token usage.
        
    Returns:
        Dictionary with image data in Claude's vision API format
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
        
    img = Draw.MolToImage(mol, size=size)
    
    # Convert PIL image to base64 string
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    img_str = base64.b64encode(buffered.getvalue()).decode()
    
    return {
        "type": "image",
        "source": {
            "type": "base64",
            "media_type": "image/png",
            "data": img_str
        }
    }

def format_molecules_output(molecules: Dict[str, List[str]], include_images: bool = True) -> List[Dict[str, str]]:
    """Format molecules dictionary into a list of content blocks for Claude's vision API.
    
    Args:
        molecules: Dictionary with molecule categories as keys and lists of SMILES as values
        include_images: Whether to include images
        
    Returns:
        List of content blocks (text and images) for Claude's vision API
    """
    content_blocks = []
    
    # Start with the opening molecules tag
    content_blocks.append({
        "type": "text",
        "text": "<molecules>"
    })
    
    for category, smiles_list in molecules.items():
        if not smiles_list:  # Skip empty categories
            continue
            
        # Add category header with the standardized format
        content_blocks.append({
            "type": "text",
            "text": f"{category}_molecules:"
        })
        
        for idx, smiles in enumerate(smiles_list, 1):
            # Add SMILES text with numbered format
            content_blocks.append({
                "type": "text",
                "text": f"{idx}. {smiles}"
            })
            
            # Add molecule image if requested
            if include_images:
                try:
                    # Validate SMILES before generating image
                    if is_valid_smiles(smiles):
                        img_block = generate_mol_image(smiles)
                        if img_block:
                            content_blocks.append(img_block)
                    else:
                        print(f"Skipping image generation for invalid SMILES: {smiles}")
                except Exception as e:
                    print(f"Error generating image for SMILES {smiles}: {str(e)}")
            
        # Add empty line between categories
        content_blocks.append({
            "type": "text",
            "text": ""
        })
    
    # End with the closing molecules tag
    content_blocks.append({
        "type": "text",
        "text": "</molecules>"
    })
        
    return content_blocks

def is_valid_smiles(smiles: str) -> bool:
    """Check if a string is a valid SMILES representation.
    
    Args:
        smiles: String to check
        
    Returns:
        True if the string is a valid SMILES, False otherwise
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except:
        return False

def get_mol(smiles_or_mol):
    '''
    Loads SMILES/molecule into RDKit's object
    '''
    if isinstance(smiles_or_mol, str):
        if len(smiles_or_mol) == 0:
            return None
        mol = Chem.MolFromSmiles(smiles_or_mol)
        if mol is None:
            return None
        try:
            Chem.SanitizeMol(mol)
        except ValueError:
            return None
        return mol
    return smiles_or_mol

from src.tools.tool_definitions import tool_registry
from rdkit.Chem import Descriptors
from rdkit.Chem import RDConfig
import os
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

@tool_registry.register("validity")
def validity(arguments: dict, output_dir_id: str = None) -> dict:
    print("Validity is being used")
    """
    Get the validity of one or more molecules and return valid/invalid status
    """
    smiles_input = arguments["smiles"]

    # Handle single SMILES string
    if isinstance(smiles_input, str):
        gen = get_mol(smiles_input)
        return {"validity": True if gen is not None else False}

    # Handle list of SMILES strings
    elif isinstance(smiles_input, list):
        results = []
        for smiles in smiles_input:
            gen = get_mol(smiles)
            results.append(True if gen is not None else False)
        return {"validity": results}

    return {"error": "Invalid input format for smiles"}

@tool_registry.register("QED")
def QED(arguments: dict, output_dir_id: str = None) -> dict:
    print("QED is being used")
    smiles_input = arguments["smiles"]

    # Handle single SMILES string
    if isinstance(smiles_input, str):
        mol = get_mol(smiles_input)
        qed = Descriptors.qed(mol) if mol is not None else None
        return {"QED": qed}

    # Handle list of SMILES strings
    elif isinstance(smiles_input, list):
        results = []
        for smiles in smiles_input:
            mol = get_mol(smiles)
            qed = Descriptors.qed(mol) if mol is not None else None
            results.append(qed)
        return {"QED": results}

    return {"error": "Invalid input format for smiles"}

@tool_registry.register("SA")
def SA(arguments: dict, output_dir_id: str = None) -> dict:
    print("SA is being used")
    smiles_input = arguments["smiles"]

    # Handle single SMILES string
    if isinstance(smiles_input, str):
        mol = get_mol(smiles_input)
        sa_score = sascorer.calculateScore(mol) if mol is not None else None
        return {"SA": sa_score}

    # Handle list of SMILES strings
    elif isinstance(smiles_input, list):
        results = []
        for smiles in smiles_input:
            mol = get_mol(smiles)
            sa_score = sascorer.calculateScore(mol) if mol is not None else None
            results.append(sa_score)
        return {"SA": results}

    return {"error": "Invalid input format for smiles"}

@tool_registry.register("logP")
def logP(arguments: dict, output_dir_id: str = None) -> dict:
    print("logP is being used")
    smiles_input = arguments["smiles"]

    # Handle single SMILES string
    if isinstance(smiles_input, str):
        mol = get_mol(smiles_input)
        logp = Descriptors.MolLogP(mol) if mol is not None else None
        return {"logP": logp}

    # Handle list of SMILES strings
    elif isinstance(smiles_input, list):
        results = []
        for smiles in smiles_input:
            mol = get_mol(smiles)
            logp = Descriptors.MolLogP(mol) if mol is not None else None
            results.append(logp)
        return {"logP": results}

    return {"error": "Invalid input format for smiles"}

@tool_registry.register("molecular_weight")
def molecular_weight(arguments: dict, output_dir_id: str = None) -> dict:
    print("molecular_weight is being used")
    smiles_input = arguments["smiles"]

    # Handle single SMILES string
    if isinstance(smiles_input, str):
        mol = get_mol(smiles_input)
        mw = Descriptors.MolWt(mol) if mol is not None else None
        return {"molecular_weight": mw}

    # Handle list of SMILES strings
    elif isinstance(smiles_input, list):
        results = []
        for smiles in smiles_input:
            mol = get_mol(smiles)
            mw = Descriptors.MolWt(mol) if mol is not None else None
            results.append(mw)
        return {"molecular_weight": results}

    return {"error": "Invalid input format for smiles"}

@tool_registry.register("molecular_similarity")
def molecular_similarity(arguments: dict, output_dir_id: str = None) -> dict:
    print("molecular_similarity is being used")

    from rdkit.Chem import AllChem
    from rdkit import DataStructs

    # Extract SMILES arguments
    ref_smiles = arguments.get("reference_smiles")
    query_smiles = arguments.get("query_smiles")

    if ref_smiles is None or query_smiles is None:
        return {"error": "Both reference_smiles and query_smiles must be provided"}

    # Convert to lists if single SMILES provided
    if isinstance(ref_smiles, str):
        ref_smiles = [ref_smiles]
    if isinstance(query_smiles, str):
        query_smiles = [query_smiles]

    # Calculate similarity for each pair
    results = []

    for r_smiles in ref_smiles:
        r_mol = get_mol(r_smiles)
        if r_mol is None:
            results.append({
                "reference": r_smiles,
                "error": "Invalid reference SMILES"
            })
            continue

        # Generate Morgan fingerprint for reference
        r_fp = AllChem.GetMorganFingerprintAsBitVect(r_mol, 2, 2048)

        for q_smiles in query_smiles:
            q_mol = get_mol(q_smiles)
            if q_mol is None:
                results.append({
                    "reference": r_smiles,
                    "query": q_smiles,
                    "error": "Invalid query SMILES"
                })
                continue

            # Generate Morgan fingerprint for query
            q_fp = AllChem.GetMorganFingerprintAsBitVect(q_mol, 2, 2048)

            # Calculate Tanimoto similarity
            similarity = DataStructs.TanimotoSimilarity(r_fp, q_fp)

            results.append({
                "reference": r_smiles,
                "query": q_smiles,
                "similarity": similarity
            })

    return {
        "results": results,
        "count": len(results)
    }

@tool_registry.register("visualize_molecules")
def visualize_molecules(arguments: dict, output_dir_id: str) -> dict:
    """
    Visualize molecules from SMILES strings and return content blocks with images.
    """
    from mol_utils import generate_mol_image, is_valid_smiles

    smiles_list = arguments["smiles"]
    labels = arguments.get("labels", [])

    # Ensure labels list matches smiles list length
    if labels and len(labels) != len(smiles_list):
        labels = [f"Molecule {i+1}" for i in range(len(smiles_list))]
    elif not labels:
        labels = [f"Molecule {i+1}" for i in range(len(smiles_list))]

    results = []

    for i, (smiles, label) in enumerate(zip(smiles_list, labels)):
        # Check SMILES validity
        if not is_valid_smiles(smiles):
            results.append({
                "smiles": smiles,
                "label": label,
                "valid": False,
                "message": "Invalid SMILES string"
            })
            continue

        # Generate molecule image
        try:
            img_block = generate_mol_image(smiles)
            if img_block:
                results.append({
                    "smiles": smiles,
                    "label": label,
                    "valid": True,
                    "image": img_block
                })
            else:
                results.append({
                    "smiles": smiles,
                    "label": label,
                    "valid": False,
                    "message": "Could not generate image"
                })
        except Exception as e:
            results.append({
                "smiles": smiles,
                "label": label,
                "valid": False,
                "message": f"Error generating image: {str(e)}"
            })

    return {
        "results": results,
        "message": f"Visualized {sum(1 for r in results if r.get('valid', False))}/{len(smiles_list)} molecules"
    }

def _parse_plip_text(plip_text: str) -> List[Dict[str, str]]:
    """Convert raw PLIP CSV text into a list of interaction dictionaries."""
    if not plip_text:
        return []
    lines = [l.strip() for l in plip_text.splitlines() if l.strip()]
    if not lines:
        return []
    header = [h.strip() for h in lines[0].split(',')]
    interactions = []
    for line in lines[1:]:
        values = [v.strip() for v in line.split(',')]
        if len(values) < len(header):
            values += ["" for _ in range(len(header) - len(values))]
        interaction = {header[i]: values[i] for i in range(len(header))}
        interactions.append(interaction)
    return interactions

def _transform_ligand_results_to_chembl_format(ligands: List[Dict], agent_name=None, iteration: int =0) -> List[Dict]:
    """
    Transform ligand results from get_mol_gen_report format to ChEMBL-like format.
    
    Args:
        ligands: List of ligand dictionaries from get_mol_gen_report
        
    Returns:
        List of transformed ligand dictionaries in ChEMBL format
    """
    transformed_results = []
    direct_fields = ['smiles', 'qed', 'sa', 'logp', 'mw', 'docking', 'ligand_path']
    
    for idx, ligand in enumerate(ligands):
        # Create new result dictionary with transformed structure
        result = {field: ligand.get(field) for field in direct_fields}
        if agent_name:
            result["friendly_id"] = tool_registry._id_generator.generate_id(agent_name, idx, iteration, parent_id = None)
        # Transform PLIP data to interactions
        plip_data = ligand.get('plip')
        if plip_data and isinstance(plip_data, dict):
            # Extract text field and convert to interactions
            plip_text = plip_data.get('text', '')
            if plip_text:
                # Use existing _parse_plip_text function to convert CSV to list of dicts
                interactions = _parse_plip_text(plip_text)
                result['plip_interactions'] = interactions
            else:
                result['plip_interactions'] = []
        else:
            result['plip_interactions'] = []
            
        transformed_results.append(result)
    
    return transformed_results
