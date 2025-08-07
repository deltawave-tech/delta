from src.tools.tool_definitions import tool_registry
from rdkit import Chem
import json
from pathlib import Path
from src.tools.mol_utils import SA, QED, logP, molecular_weight
from src.tools.vina import run_vina_docker
from src.tools.api_based_tools.plip import get_plip_report
from src.tools.mol_utils import _transform_ligand_results_to_chembl_format


@tool_registry.register("vina_report")
def get_vina_report(arguments: dict, output_dir_id: str) -> dict:
    # Extract parameters from arguments dictionary
    from pathlib import Path
    import concurrent.futures
    import threading
    import time
    # Extract parame
    protein_path = arguments.get("protein_path")
    molecules = arguments.get("molecules", [])
    run_id = output_dir_id
    if not molecules:
        return {'error': 'No molecules provided'}
    # Use output_dir directly from tool_registry
    output_dir = tool_registry.output_dir / run_id
    print(f"Output directory: {output_dir}")
    protein_path = Path(protein_path)
    
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

    # Process all molecules with vina_docker in batch
    docking_result_list = []

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
                if batch_idx > 0: # Stagger batch submissions
                    delay = 1 * batch_idx
                    time.sleep(delay)
                
                batch_docking_result = tool_registry.execute(
                    'run_vina_docker',
                    {'protein_path': str(protein_path), 'smiles': batch_smiles, 'batch_id': batch_idx},
                    run_id=run_id
                )
                return batch_docking_result, None
            except Exception as e:
                print(f"Error in batch {batch_idx}: {str(e)}")
                return None, str(e)
            finally:
                batch_semaphore.release()

        with concurrent.futures.ThreadPoolExecutor(max_workers=max_concurrent) as executor:
            futures = [executor.submit(process_batch_with_semaphore, batch, i) for i, batch in enumerate(batches)]
            for future in concurrent.futures.as_completed(futures):
                result, error = future.result()
                if error:
                    print(f"A batch failed: {error}")
                if result:
                    if isinstance(result, list):
                        docking_result_list.extend(result)
                    elif isinstance(result, dict):
                        docking_result_list.append(result)
    else: # Process as a single batch
        direct_result = tool_registry.execute(
            'run_vina_docker',
            {'protein_path': str(protein_path), 'smiles': sanitized_molecules},
            run_id=run_id
        )
        if isinstance(direct_result, list):
            docking_result_list.extend(direct_result)
        elif isinstance(direct_result, dict):
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

    best_energy_map   = {}
    complex_path_map  = {}

    # normalise the result into the two maps above
    for res_item in docking_result_list:
        if not isinstance(res_item, dict):
            continue

        if 'smiles' in res_item and 'binding_energy' in res_item:
            res_smi   = res_item.get('smiles')
            orig_smi  = smiles_map.get(res_smi, res_smi)
            if orig_smi:
                best_energy_map[orig_smi]  = _best_energy(res_item.get('binding_energy'))
                complex_path_map[orig_smi] = res_item.get('merged_path')
        
        elif 'binding_energy' in res_item and isinstance(res_item['binding_energy'], dict):
            for res_smi_key, energies in res_item['binding_energy'].items():
                orig_smi  = smiles_map.get(res_smi_key, res_smi_key)
                if orig_smi:
                    best_energy_map[orig_smi]  = _best_energy(energies)
                    if 'merged_path' in res_item and isinstance(res_item['merged_path'], dict):
                        complex_path_map[orig_smi] = res_item['merged_path'].get(res_smi_key)
                    elif isinstance(res_item.get('merged_path'), str):
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

