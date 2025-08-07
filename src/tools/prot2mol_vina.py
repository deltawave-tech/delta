from src.tools.tool_definitions import tool_registry
from rdkit import Chem
import json
from pathlib import Path
from src.tools.mol_utils import SA, QED, logP, molecular_weight
from src.tools.prot2mol import generate_molecules
from src.tools.vina import run_vina_docker
from src.tools.api_based_tools.plip import get_plip_report
from src.tools.mol_utils import _transform_ligand_results_to_chembl_format

@tool_registry.register("vina_mol_gen")
def get_vina_mol_gen_report(arguments: dict, output_dir_id: str) -> dict:
    from pathlib import Path
    # Extract parameters from arguments dictionary
    protein_path = arguments.get("protein_path")
    molecules = arguments.get("molecules", [])
    protein_sequence = arguments.get("protein_sequence")
    generate_de_novo = arguments.get("generate_de_novo", False)
    thresholds = arguments.get("thresholds", {})  # Default to empty dict if None
    run_id = output_dir_id

    # Validate input parameters
    if not molecules and not generate_de_novo:
        return {'error': 'Either molecules must be provided or generate_de_novo must be set to True'}

    # Use output_dir directly from tool_registry
    output_dir = tool_registry.output_dir / run_id
    print(f"Output directory: {output_dir}")
    protein_path = Path(protein_path)

    # Generate molecules if requested and not provided
    if generate_de_novo is True:
        print("Generating de novo molecules")
        # Use prot2mol to generate molecules based on protein sequence
        num_molecules = arguments.get("num_molecules", 10)
        smiles_list = tool_registry.execute('generate_molecules', {'protein_sequence': protein_sequence, 'num_molecules': num_molecules}, run_id = run_id)
        if not smiles_list or len(smiles_list) == 0:
            return {'error': 'No valid molecules generated'}
        molecules.extend(smiles_list)

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
    # AND molecules were generated de novo
    filtered_molecules = []
    
    if thresholds and any([SA, QED, logP, molecular_weight]) and molecular_properties and generate_de_novo == True:
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
        
    if not generate_de_novo and thresholds and any([SA, QED, logP, molecular_weight]):
        print("Filtering bypassed because molecules were supplied directly, not generated de novo.")

    # Process all molecules with vina_docker in batch
    docking_results_list = []
    complex_paths = []

    # Create a directory for vina docking results
    vina_dir = output_dir / 'vina_docking'
    vina_dir.mkdir(parents=True, exist_ok=True)

    # Handle multi-fragment molecules using RDKit's canonicalization
    sanitized_molecules = []
    smiles_map = {}  # Map processed SMILES to original SMILES
    for smiles in molecules:
        try:
            smi_to_process = smiles
            if '.' in smiles:
                smi_to_process = max(smiles.split('.'), key=len)
            
            mol = Chem.MolFromSmiles(smi_to_process)
            if mol:
                canonical_smiles = Chem.MolToSmiles(mol)
                sanitized_molecules.append(canonical_smiles)
                if canonical_smiles not in smiles_map:
                    smiles_map[canonical_smiles] = smiles # Map canonical of largest fragment to original
            else:
                print(f"Warning: Could not parse SMILES: {smiles} - skipping")
        except Exception as e:
            print(f"Error processing SMILES {smiles}: {e}")
    
    # Run Vina docking on all molecules at once
    print(f"Running Vina docking on {len(sanitized_molecules)} valid molecules")
    if not sanitized_molecules:
        return {
            'error': 'No valid molecules to dock after sanitization',
            'status': 'failed',
            'details': f"Original molecules count: {len(molecules)}, None were valid after sanitization.",
            'result': []
        }
    
    import concurrent.futures
    import threading
    import time

    docking_result_list = []
    batch_size = 10
    
    if len(sanitized_molecules) > batch_size:
        batches = [sanitized_molecules[i:i + batch_size] for i in range(0, len(sanitized_molecules), batch_size)]
        print(f"Processing {len(sanitized_molecules)} molecules in {len(batches)} batches of up to {batch_size}")

        # Set max concurrent batches.
        # The server has 12 vCPUs, so we can run several batches in parallel.
        # This value can be tuned based on server performance.
        max_concurrent = 6
        
        batch_semaphore = threading.Semaphore(max_concurrent)

        def process_batch_with_semaphore(batch_smiles, batch_idx):
            try:
                batch_semaphore.acquire()
                # print(f"Starting batch {batch_idx} processing (acquired semaphore)")
                if batch_idx > 0: # Stagger batch submissions
                    delay = 1 * batch_idx # Reduced delay slightly
                    # print(f"Delaying batch {batch_idx} submission by {delay}s")
                    time.sleep(delay)
                
                batch_docking_result = tool_registry.execute(
                    'run_vina_docker',
                    {'protein_path': str(protein_path),
                        'smiles': batch_smiles,
                        'batch_id': batch_idx}, # Pass batch_id for tracking if run_vina_docker uses it
                    run_id=run_id
                )
                # print(f"Batch {batch_idx} completed with {len(batch_smiles)} molecules.")
                return batch_docking_result, None
            except Exception as e:
                print(f"Error in batch {batch_idx}: {str(e)}")
                return None, str(e) # Return error to handle
            finally:
                batch_semaphore.release()
                # print(f"Released semaphore for batch {batch_idx}")

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            futures = [executor.submit(process_batch_with_semaphore, batch, i) for i, batch in enumerate(batches)]
            for future in concurrent.futures.as_completed(futures):
                result, error = future.result()
                if error:
                    print(f"A batch failed: {error}") # Log or handle as needed
                if result:
                    if isinstance(result, list):
                        docking_result_list.extend(result)
                    elif isinstance(result, dict): # if a single batch returns a dict (e.g. combined result)
                        docking_result_list.append(result) 
                        # This might need adjustment based on run_vina_docker's actual return for single batches
                        # The downstream logic expects a flat list of dicts, one per molecule.
    else: # Process as a single batch if not exceeding batch_size
        # print(f"Processing {len(sanitized_molecules)} molecules in a single batch.")
        direct_result = tool_registry.execute(
            'run_vina_docker',
            {'protein_path': str(protein_path),
                'smiles': sanitized_molecules},
            run_id=run_id
        )
        if isinstance(direct_result, list):
            docking_result_list.extend(direct_result)
        elif isinstance(direct_result, dict):
                # This assumes if vina returns a single dict for few molecules, it's a combined result
                # that needs to be deconstructed or the downstream logic needs to handle it.
                # For simplicity, let's assume it returns a list even for small inputs,
                # or the processing logic below can handle a list containing one complex dict.
                # The current downstream logic (isinstance(docking_result_list, list): for res in docking_result_list)
                # should work if it's a list of dicts.
            docking_result_list.append(direct_result)

    # helper → return best (lowest) energy from list or scalar
    def _best_energy(val):
        if isinstance(val, list) and val:
            try:
                return min(float(v) for v in val if v is not None)
            except (ValueError, TypeError):
                return None
        elif val is not None:
            try:
                return float(val)
            except (ValueError, TypeError):
                return None
        return None

    best_energy_map   = {}               # original_smiles → float
    complex_path_map  = {}               # original_smiles → str

    # normalise the result into the two maps above
    for res_item in docking_result_list: # res_item is a dict for one molecule or a combined dict
        if not isinstance(res_item, dict):
            # print(f"Warning: Unexpected item in docking_result_list: {type(res_item)}")
            continue

        # Scenario A: res_item is a dictionary for a single SMILES (typical from batch item)
        if 'smiles' in res_item and 'binding_energy' in res_item:
            res_smi   = res_item.get('smiles')
            orig_smi  = smiles_map.get(res_smi, res_smi) # Map canonical back to original
            if orig_smi:
                best_energy_map[orig_smi]  = _best_energy(res_item.get('binding_energy'))
                complex_path_map[orig_smi] = res_item.get('merged_path')
        
        # Scenario B: res_item is a combined dictionary (e.g., if a batch returns one dict for multiple SMILES)
        elif 'binding_energy' in res_item and isinstance(res_item['binding_energy'], dict):
            for res_smi_key, energies in res_item['binding_energy'].items():
                orig_smi  = smiles_map.get(res_smi_key, res_smi_key)
                if orig_smi:
                    best_energy_map[orig_smi]  = _best_energy(energies)
                    if 'merged_path' in res_item and isinstance(res_item['merged_path'], dict):
                        complex_path_map[orig_smi] = res_item['merged_path'].get(res_smi_key)
                    elif isinstance(res_item.get('merged_path'), str): # if merged_path is single string for this combined dict
                            complex_path_map[orig_smi] = res_item.get('merged_path')

    # Run PLIP on all complex files
    plip_paths   = [p for p in complex_path_map.values() if p and Path(p).exists()]
    plip_results_map = {}
    if plip_paths:
        try:
            plip_out = tool_registry.execute('get_plip_report', {
                'input_structure': plip_paths,
                'wait_for_result': True,
                'max_wait_time': 600
            }, run_id=run_id)
            if plip_out and 'results' in plip_out and isinstance(plip_out['results'], dict):
                plip_results_map = plip_out['results']
        except Exception as e:
            print(f"Error running PLIP analysis: {e}")

    # Compile all results based on successfully docked molecules
    standardized_results = []
    
    for original_smiles in molecules:
        binding_energy = best_energy_map.get(original_smiles)
        complex_path = complex_path_map.get(original_smiles)

        plip_output = None
        if complex_path and complex_path in plip_results_map:
            plip_output = plip_results_map[complex_path]
        
        props = molecular_properties.get(original_smiles, {})
        is_valid = molecule_validity.get(original_smiles, False)

        result = {
            'smiles': original_smiles,
            'qed': props.get('qed'),
            'sa': props.get('sa'),
            'logp': props.get('logp'),
            'mw': props.get('molecular_weight'),
            'docking': binding_energy,
            'plip': plip_output,
            'is_valid': is_valid
        }
        
        standardized_results.append(result)
    
    # Also check if we have summary.json data as a fallback
    if not standardized_results:
        # Instead of using summary.json as fallback, return molecules with error message
        return {
            'ligands': [],
            'error': 'Docking process failed to generate valid results',
            'molecules': molecules,
            'docking_summary': {'method': 'Autodock Vina', 'status': 'failed'},
            'plip_summary': {'results': {}}
        }
    transformed_ligands = _transform_ligand_results_to_chembl_format(standardized_results)
    
    return {
        'ligands': transformed_ligands,
        # 'docking_summary': docking_results,
        # 'plip_summary': plip_results
    }    