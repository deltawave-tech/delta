from src.tools.api_based_tools.api_tools import VinaDockerAPI
from src.tools.tool_definitions import tool_registry
from src.tools.mol_utils import get_largest_fragment
from rdkit import Chem

@tool_registry.register("run_vina_docker")
def run_vina_docker(arguments: dict, output_dir_id: str) -> dict:
    print("Vina Docker is being used")
    vina_docker = VinaDockerAPI()

    # Create output directory with optional batch_id for parallel processing
    batch_id = arguments.get("batch_id")
    if batch_id is not None:
        output_dir = tool_registry.output_dir/output_dir_id/f'vina_docking_{batch_id}'
        print(f"Processing batch {batch_id} in directory: {output_dir}")
    else:
        output_dir = tool_registry.output_dir/output_dir_id/'vina_docking'
    
    output_dir.mkdir(parents=True, exist_ok=True)

    # Prepare the protein file
    protein_file = arguments["protein_path"]

    # Determine the source of the ligand(s)
    ligand_files = []
    smiles_list = []

    # Case 1: Ligand file path is provided
    if "ligand_path" in arguments and arguments["ligand_path"]:
        if isinstance(arguments["ligand_path"], list):
            ligand_files = arguments["ligand_path"]
        else:
            ligand_files = [arguments["ligand_path"]]

    # Case 2: Single SMILES string is provided
    elif "ligand_smiles" in arguments and arguments["ligand_smiles"]:
        smiles = arguments["ligand_smiles"]
        if isinstance(smiles, list):
            # Handle list of SMILES strings
            smiles = get_largest_fragment(smiles)
            smiles_list = smiles
            for i, s in enumerate(smiles):
                print(f"Converting SMILES {i+1}/{len(smiles)}: {s}")
                mol = Chem.MolFromSmiles(s)
                if mol is None:
                    print(f"Warning: Invalid SMILES string: {s} - skipping")
                    continue
                
                # Add hydrogens (optional but helpful for docking)
                mol = Chem.AddHs(mol)
                
                # Generate 3D coordinates
                try:
                    from rdkit.Chem import AllChem
                    AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coordinates
                    AllChem.MMFFOptimizeMolecule(mol)  # Energy minimize with MMFF
                    print(f"  Generated 3D coordinates for molecule {i+1}")
                except Exception as e:
                    print(f"  Warning: Could not generate 3D coordinates for {s}: {e}")
                    # If 3D generation fails, skip this molecule
                    continue
                
                ligand_file = output_dir / f"temp_ligand_{i}.sdf"
                writer = Chem.SDWriter(str(ligand_file))
                writer.write(mol)
                writer.close()
                ligand_files.append(ligand_file)
                print(f"Created temporary SDF file: {ligand_file}")
        else:
            # Handle single SMILES string
            smiles_list = [smiles]
            print(f"Converting SMILES to SDF: {smiles}")
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"error": f"Invalid SMILES string: {smiles}"}
            
            # Add hydrogens (optional but helpful for docking)
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            try:
                from rdkit.Chem import AllChem
                AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coordinates
                AllChem.MMFFOptimizeMolecule(mol)  # Energy minimize with MMFF
                print(f"  Generated 3D coordinates")
            except Exception as e:
                print(f"  Warning: Could not generate 3D coordinates for {smiles}: {e}")
                # If 3D generation fails, return an error
                return {"error": f"Could not generate 3D coordinates for: {smiles}"}
            
            ligand_file = output_dir / "temp_ligand.sdf"
            writer = Chem.SDWriter(str(ligand_file))
            writer.write(mol)
            writer.close()
            ligand_files.append(ligand_file)
            print(f"Created temporary SDF file: {ligand_file}")

    # Case 3: List of SMILES is provided
    elif "smiles" in arguments and arguments["smiles"] and len(arguments["smiles"]) > 0:
        smiles_list = arguments["smiles"]
        smiles_list = get_largest_fragment(smiles_list)
        for i, s in enumerate(smiles_list):
            print(f"Converting SMILES {i+1}/{len(smiles_list)}: {s}")
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                print(f"Warning: Invalid SMILES string: {s} - skipping")
                continue
            
            # Add hydrogens (optional but helpful for docking)
            mol = Chem.AddHs(mol)
            
            # Generate 3D coordinates
            try:
                from rdkit.Chem import AllChem
                AllChem.EmbedMolecule(mol, randomSeed=42)  # Generate 3D coordinates
                AllChem.MMFFOptimizeMolecule(mol)  # Energy minimize with MMFF
                print(f"  Generated 3D coordinates for molecule {i+1}")
            except Exception as e:
                print(f"  Warning: Could not generate 3D coordinates for {s}: {e}")
                # If 3D generation fails, skip this molecule
                continue
                
            ligand_file = output_dir / f"temp_ligand_{i}.sdf"
            writer = Chem.SDWriter(str(ligand_file))
            writer.write(mol)
            writer.close()
            ligand_files.append(ligand_file)
            print(f"Created temporary SDF file: {ligand_file}")

    # Error if no ligand source is provided
    else:
        return {"error": "Either ligand_path, ligand_smiles, or smiles list must be provided"}

    # Error if no valid ligands were processed
    if not ligand_files:
        return {"error": "No valid ligands were processed"}

    # Run Vina docking with the prepared ligand files
    # Always pass as a list for consistency
    print(f"Calling vina_docker.predict with {len(ligand_files)} ligand files")
    result = vina_docker.predict(
        protein_file=protein_file,
        ligand_file=ligand_files,  # Always pass as a list, even for a single file
        smiles_list=smiles_list,   # Pass the SMILES list
        output_dir=output_dir
    )
    print("Vina Docker result: ", result)

    
    return result