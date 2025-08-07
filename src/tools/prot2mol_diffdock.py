from tool_definitions import tool_registry
from rdkit import Chem
import json
from pathlib import Path
from mol_utils import SA, QED, logP, molecular_weight, merge_pdb_sdf_rdkit
import glob

@tool_registry.register("mol_gen")
def get_mol_gen_report(arguments: dict, output_dir_id: str) -> dict:
    # Extract parameters from arguments dictionary
    protein_path = arguments.get("protein_path")
    molecules = arguments.get("molecules")
    protein_sequence = arguments.get("protein_sequence")
    thresholds = arguments.get("thresholds", {})  # Default to empty dict if None
    run_id = output_dir_id

    # Use output_dir directly from tool_registry
    output_dir = tool_registry.output_dir / run_id
    print(f"Output directory: {output_dir}")
    protein_path = Path(protein_path)
    protein_name = protein_path.name

    # Generate molecules if not provided
    if not molecules:
        # Use prot2mol to generate molecules based on protein sequence
        num_molecules = arguments.get("num_molecules", 10)
        smiles_list = tool_registry.execute('generate_molecules', {'protein_sequence': protein_sequence, 'num_molecules': num_molecules}, run_id = run_id)
        if not smiles_list or len(smiles_list) == 0:
            return {'error': 'No valid molecules generated'}
        molecules = smiles_list

    print(f"Molecules: {molecules}")
    
    # Check molecule validity
    molecule_validity = {}
    for smiles in molecules:
        try:
            mol = Chem.MolFromSmiles(smiles)
            molecule_validity[smiles] = mol is not None
        except:
            molecule_validity[smiles] = False
    
    # Get molecular properties for all molecules - do this regardless of thresholds
    molecular_properties = {}
    # Force calculation of properties even if thresholds aren't defined
    try:
        sa_values = SA({"smiles": molecules})["SA"] if SA else None
        qed_values = QED({"smiles": molecules})["QED"] if QED else None
        logp_values = logP({"smiles": molecules})["logP"] if logP else None
        mw_values = molecular_weight({"smiles": molecules})["molecular_weight"] if molecular_weight else None

        for i, mol_smiles in enumerate(molecules):
            molecular_properties[mol_smiles] = {
                "sa": sa_values[i] if sa_values and i < len(sa_values) else None,
                "qed": qed_values[i] if qed_values and i < len(qed_values) else None,
                "logp": logp_values[i] if logp_values and i < len(logp_values) else None,
                "molecular_weight": mw_values[i] if mw_values and i < len(mw_values) else None
            }
        print(f"Calculated properties for {len(molecular_properties)} molecules")
    except Exception as e:
        print(f"Error calculating molecular properties: {e}")
        # Continue with empty properties rather than failing

    # Filter molecules based on properties only if thresholds are specified
    filtered_molecules = []
    molecules_from_args = arguments.get("molecules") is not None
    
    if thresholds and any([SA, QED, logP, molecular_weight]) and molecular_properties and not molecules_from_args:
        max_attempts = 3  # Limit the number of generation attempts
        attempts = 0

        while attempts < max_attempts:
            filtered_molecules = []

            for i, mol_smiles in enumerate(molecules):
                # Check if molecule passes all threshold criteria
                passes_filters = True
                props = molecular_properties.get(mol_smiles, {})

                # Synthetic accessibility
                if SA and 'sa' in thresholds:
                    sa_value = props.get('sa')
                    sa_min = thresholds['sa'].get('min', float('-inf'))
                    sa_max = thresholds['sa'].get('max', float('inf'))
                    if sa_value is None or not (sa_min <= sa_value <= sa_max):
                        passes_filters = False

                # Drug-likeness (QED)
                if passes_filters and QED and 'qed' in thresholds:
                    qed_value = props.get('qed')
                    qed_min = thresholds['qed'].get('min', float('-inf'))
                    qed_max = thresholds['qed'].get('max', float('inf'))
                    if qed_value is None or not (qed_min <= qed_value <= qed_max):
                        passes_filters = False

                # LogP (lipophilicity)
                if passes_filters and logP and 'logp' in thresholds:
                    logp_value = props.get('logp')
                    logp_min = thresholds['logp'].get('min', float('-inf'))
                    logp_max = thresholds['logp'].get('max', float('inf'))
                    if logp_value is None or not (logp_min <= logp_value <= logp_max):
                        passes_filters = False

                # Molecular weight
                if passes_filters and molecular_weight and 'molecular_weight' in thresholds:
                    mw_value = props.get('molecular_weight')
                    mw_min = thresholds['molecular_weight'].get('min', float('-inf'))
                    mw_max = thresholds['molecular_weight'].get('max', float('inf'))
                    if mw_value is None or not (mw_min <= mw_value <= mw_max):
                        passes_filters = False

                if passes_filters:
                    filtered_molecules.append(mol_smiles)

            if filtered_molecules:
                molecules = filtered_molecules
                print(f"After filtering: {len(molecules)} molecules remain")
                break
            else:
                print(f"No molecules left after filtering (attempt {attempts+1}/{max_attempts}). Regenerating...")
                attempts += 1
                # Generate new molecules
                num_molecules = arguments.get("num_molecules", 10)
                smiles_list = tool_registry.execute('generate_molecules', {'protein_sequence': protein_sequence, 'num_molecules': num_molecules}, run_id = run_id)
                if not smiles_list or len(smiles_list) == 0:
                    return {'error': 'No valid molecules generated after filtering attempts'}
                molecules = smiles_list

                # Update property values for new molecules
                try:
                    sa_values = SA({"smiles": molecules})["SA"] if SA else None
                    qed_values = QED({"smiles": molecules})["QED"] if QED else None
                    logp_values = logP({"smiles": molecules})["logP"] if logP else None
                    mw_values = molecular_weight({"smiles": molecules})["molecular_weight"] if molecular_weight else None

                    for i, mol_smiles in enumerate(molecules):
                        molecular_properties[mol_smiles] = {
                            "sa": sa_values[i] if sa_values and i < len(sa_values) else None,
                            "qed": qed_values[i] if qed_values and i < len(qed_values) else None,
                            "logp": logp_values[i] if logp_values and i < len(logp_values) else None,
                            "molecular_weight": mw_values[i] if mw_values and i < len(mw_values) else None
                        }
                except Exception as e:
                    print(f"Error calculating properties for regenerated molecules: {e}")

        if not molecules:
            return {'error': 'No molecules passed the filtering criteria after multiple attempts'}
    elif thresholds:
        # If thresholds are specified but property calculation failed
        print("Warning: Thresholds were specified but property calculation failed. Using all molecules.")
    else:
        print("No thresholds specified. Using all molecules without filtering.")
        
    # Add log message if filtering was skipped because molecules were supplied as arguments
    if molecules_from_args and thresholds and any([SA, QED, logP, molecular_weight]):
        print("Filtering bypassed because molecules were supplied as arguments, not generated de novo.")

    # Check if docking results already exist
    docking_dir = output_dir / 'docking' / f"cleaned_{protein_name}"
    rank1_files = glob.glob(str(docking_dir / '**/rank1_*.sdf'), recursive=True)

    # Check if all current molecules have been docked
    existing_smiles = set()
    need_docking = False

    # Only process if there are existing docking results
    if rank1_files:
        # Extract SMILES from directory names
        for rank1_file in rank1_files:
            ligand_dir = Path(rank1_file).parent
            # The parent directory name is the SMILES string (potentially encoded)
            encoded_smiles = ligand_dir.name
            # Replace any encoding used during saving (e.g., [SLASH] -> /)
            smiles = encoded_smiles.replace('[SLASH]', '/')
            existing_smiles.add(smiles)

        # Check if all current molecules exist in previous docking results
        for mol in molecules:
            if mol not in existing_smiles:
                need_docking = True
                break
    else:
        # No existing results at all
        need_docking = True

    # If any molecules need docking or no previous results exist, run DiffDock
    if need_docking:
        print(f"Running DiffDock for {len(molecules)} molecules ({len(existing_smiles.intersection(molecules))} already docked)...")
        docking_results = tool_registry.execute('run_diffdock', {
            'protein_path': protein_path,
            'smiles': molecules
        }, run_id = run_id)

        # Check for docking results again
        rank1_files = glob.glob(str(docking_dir / '**/rank1_*.sdf'), recursive=True)
        if not rank1_files:
            return {'error': f"No rank1*.sdf files found in {docking_dir} after docking"}
    else:
        print(f"All {len(molecules)} molecules have already been docked. Using existing results.")
        docking_results = {"note": "Using existing docking results"}

    # Process all ligands
    standardized_results = []
    all_ligand_data = []
    all_complex_paths = []

    # Create a mapping from directory names to original SMILES
    # This improves the matching between docking results and original molecules
    dir_to_smiles_map = {}
    for smiles in molecules:
        # Create the same encoding for the directory name that would be used during docking
        sanitized_smiles = smiles.replace('/', '[SLASH]')
        dir_to_smiles_map[sanitized_smiles] = smiles

    # First phase: Create all complexes and collect data
    for rank1_file in rank1_files:
        ligand_path = Path(rank1_file).resolve()
        # Extract ligand name from directory path
        ligand_dir = ligand_path.parent
        ligand_name = ligand_dir.name

        # Extract confidence score from filename
        confidence = None
        filename = ligand_path.name
        try:
            if 'confidence-' in filename:
                # Handle negative confidence (confidence-0.48.sdf)
                confidence_part = filename.split('confidence-')[-1].split('.sdf')[0]
                try:
                    confidence = -float(confidence_part)
                except ValueError:
                    print(f"Could not parse negative confidence from {filename}")
            elif 'confidence' in filename:
                # Handle positive confidence (confidence0.48.sdf)
                confidence_part = filename.split('confidence')[-1].split('.sdf')[0]
                try:
                    confidence = float(confidence_part)
                except ValueError:
                    print(f"Could not parse positive confidence from {filename}")
        except Exception as e:
            print(f"Error extracting confidence from filename {filename}: {e}")

        # Get SMILES for this ligand using improved matching
        smiles = None
        # First try direct match from directory name map
        if ligand_name in dir_to_smiles_map:
            smiles = dir_to_smiles_map[ligand_name]
        else:
            # Fall back to the previous matching method
            for mol_smiles in molecules:
                sanitized_smiles = mol_smiles.replace('/', '[SLASH]')
                if ligand_name == sanitized_smiles or ligand_name in sanitized_smiles:
                    smiles = mol_smiles
                    break

            # If still no match, try a more flexible approach
            if not smiles:
                for mol_smiles in molecules:
                    if ligand_name in mol_smiles or any(part in ligand_name for part in mol_smiles.split('=')):
                        smiles = mol_smiles
                        break

        # Create complex path for this ligand
        complex_path = protein_path.parent / f"{protein_path.stem}_{ligand_name}_merged{protein_path.suffix}"

        # Check if complex already exists
        try:
            if not complex_path.exists():
                print(f"Creating complex for ligand: {ligand_name}")
                cleaned_protein_path = protein_path.parent / f"cleaned_{protein_path.stem}.pdb"
                # Merge the protein and ligand into a complex for interaction analysis
                merge_pdb_sdf_rdkit(ligand_path, cleaned_protein_path, complex_path)
            else:
                print(f"Complex already exists for ligand: {ligand_name}")
        except Exception as e:
            print(f"Error creating complex for {ligand_name}: {e}")
            # Use the protein path as a fallback if complex creation fails
            complex_path = protein_path

        # Get property values for this molecule
        qed_value = None
        sa_value = None
        logp_value = None
        mw_value = None

        if smiles and smiles in molecular_properties:
            props = molecular_properties[smiles]
            qed_value = props.get('qed')
            sa_value = props.get('sa')
            logp_value = props.get('logp')
            mw_value = props.get('molecular_weight')
        else:
            print(f"No properties found for smiles: {smiles}")

        # Store ligand data and complex path for later
        ligand_data = {
            'ligand_name': ligand_name,
            'ligand_path': str(ligand_path),
            'complex_path': str(complex_path),
            'smiles': smiles,
            'qed': qed_value,
            'sa': sa_value,
            'logp': logp_value,
            'mw': mw_value,
            'docking': confidence  # Rename to match the final output format
        }
        all_ligand_data.append(ligand_data)
        all_complex_paths.append(str(complex_path))

    # Second phase: Run PLIP on all complexes only if not done previously
    print(f"Checking for existing PLIP results for {len(all_complex_paths)} complexes")

    # Create a dictionary to track which complexes need analysis
    complexes_to_analyze = []
    existing_plip_results = {}

    # Check for existing PLIP result files
    for complex_path in all_complex_paths:
        complex_file = Path(complex_path)
        # Construct expected PLIP results file path
        plip_result_file = output_dir / 'plip_results' / f"{complex_file.stem}_plip_report.json"

        if plip_result_file.exists():
            print(f"Found existing PLIP results for {complex_file.name}")
            try:
                # Load existing results from file
                with open(plip_result_file, 'r') as f:
                    existing_plip_results[complex_path] = json.load(f)
            except (json.JSONDecodeError, IOError) as e:
                print(f"Error loading existing PLIP results for {complex_file.name}: {e}")
                # If there's an error loading the file, add to list for re-analysis
                complexes_to_analyze.append(complex_path)
        else:
            # No existing results, add to list for analysis
            complexes_to_analyze.append(complex_path)

    # Only run PLIP on complexes that need analysis
    if complexes_to_analyze:
        print(f"Running PLIP on {len(complexes_to_analyze)} new complexes")
        try:
            new_plip_results = tool_registry.execute('get_plip_report', {
                'input_structure': complexes_to_analyze,
                'wait_for_result': True,
                'max_wait_time': 600  # Increase timeout for multiple complexes
            }, run_id = run_id)

            # Save new results to files for future use
            if new_plip_results and 'results' in new_plip_results:
                plip_results_dir = output_dir / 'plip_results'
                plip_results_dir.mkdir(parents=True, exist_ok=True)

                if isinstance(new_plip_results['results'], dict):
                    # If results is a dictionary with complex paths as keys
                    for path, result in new_plip_results['results'].items():
                        existing_plip_results[path] = result
                        # Save individual result to file
                        complex_file = Path(path)
                        result_file = plip_results_dir / f"{complex_file.stem}_plip_report.json"
                        with open(result_file, 'w') as f:
                            json.dump(result, f)
                elif isinstance(new_plip_results['results'], list) and len(new_plip_results['results']) == len(complexes_to_analyze):
                    # If results is a list in the same order as input complexes
                    for i, path in enumerate(complexes_to_analyze):
                        if i < len(new_plip_results['results']):
                            existing_plip_results[path] = new_plip_results['results'][i]
                            # Save individual result to file
                            complex_file = Path(path)
                            result_file = plip_results_dir / f"{complex_file.stem}_plip_report.json"
                            with open(result_file, 'w') as f:
                                json.dump(new_plip_results['results'][i], f)
        except Exception as e:
            print(f"Error running PLIP analysis: {e}")
            # Continue without PLIP results rather than failing completely
    else:
        print("All complexes already have PLIP results. Using cached data.")

    # Now we have combined results in existing_plip_results
    plip_results = {"results": existing_plip_results}

    # Create a flag to indicate if any new PLIP analysis was performed
    new_analysis_performed = len(complexes_to_analyze) > 0

    # Third phase: Associate PLIP results with ligands and create standardized output
    for ligand_data in all_ligand_data:
        complex_path = ligand_data['complex_path']
        plip_output = None

        # Find matching PLIP result for this complex
        if complex_path in existing_plip_results:
            # For previously analyzed complexes, include full results only if new analysis was performed
            if new_analysis_performed:
                plip_output = existing_plip_results[complex_path]
            else:
                # Create a summary instead of full report for previously analyzed complexes
                # when no new analysis was performed
                if isinstance(existing_plip_results[complex_path], dict):
                    # Extract key metrics from the full report
                    interaction_count = 0
                    interaction_types = set()

                    # Process the PLIP report structure to extract summary information
                    if 'interactions' in existing_plip_results[complex_path]:
                        for interaction_type, interactions in existing_plip_results[complex_path]['interactions'].items():
                            if isinstance(interactions, list):
                                interaction_count += len(interactions)
                                interaction_types.add(interaction_type)

                    plip_output = {
                        "summary": f"Complex previously analyzed: {interaction_count} interactions of types {', '.join(interaction_types)}",
                        "full_report_available": True,
                        "complex_path": complex_path
                    }
                else:
                    plip_output = {"summary": "Complex previously analyzed (see existing report)", "full_report_available": True}

        # Get validity status
        smiles = ligand_data['smiles']
        is_valid = molecule_validity.get(smiles, False) if smiles else False

        # Create standardized result entry with requested format
        result = {
            'smiles': ligand_data['smiles'] or None,
            'qed': ligand_data['qed'] or None,
            'sa': ligand_data['sa'] or None,
            'logp': ligand_data['logp'] or None,
            'mw': ligand_data['mw'] or None,
            'docking': ligand_data['docking'] or None,
            'plip': plip_output or None,
            'ligand_path': ligand_data['ligand_path'],
            'complex_path': ligand_data['complex_path'],
            'is_valid': is_valid
        }

        standardized_results.append(result)

    return {
        'ligands': standardized_results,
        'docking_summary': docking_results,
        'plip_summary': plip_results
    }


