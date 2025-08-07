import requests
import pandas as pd
import os
import time
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import concurrent.futures


def create_session():
    """Create a requests session with retries and connection pooling"""
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.5, status_forcelist=[500, 502, 503, 504])
    session.mount('https://', HTTPAdapter(max_retries=retries, pool_connections=100, pool_maxsize=100))
    return session

def fetch_page(session, url, params):
    """Helper function to fetch a single page of results"""
    try:
        response = session.get(url, params=params)
        response.raise_for_status()
        return response.json()
    except Exception as e:
        print(f"Error fetching page: {e}")
        return None

def fetch_activities(target_chembl_ids, assay_types, pchembl_threshold_for_download, active_inactive_compounds, sa, qed, logP, molecular_weight):
    print("Starting to fetch activities...")
    base_url = "https://www.ebi.ac.uk/chembl/api/data/activity.json"
    
    params = {
        'target_chembl_id__in': ','.join(target_chembl_ids),
        'assay_type__in': ','.join(assay_types),
        'pchembl_value__isnull': 'false',
        'only': 'molecule_chembl_id,pchembl_value,target_chembl_id',
        'limit': 1000
    }
    
    session = create_session()
    
    # Create dictionary to store all activities first
    all_compounds_dict = {}    # molecule_id -> list of activities
    
    offset = 0
    
    while True:
        params['offset'] = offset
        result = fetch_page(session, base_url, params)
        
        if not result or 'activities' not in result:
            break
            
        activities_found = len(result['activities'])
        if activities_found == 0:
            break
            
        for activity in result['activities']:
            try:
                molecule_id = activity['molecule_chembl_id']
                pchembl_value = float(activity['pchembl_value'])
                
                # Store all activities for this molecule
                if molecule_id not in all_compounds_dict:
                    all_compounds_dict[molecule_id] = []
                all_compounds_dict[molecule_id].append(activity)
                
            except (ValueError, TypeError):
                continue
        
        # Break if we've processed all available results
        if activities_found < params['limit']:
            break
            
        offset += params['limit']
    
    # Now classify molecules based on their average pchembl_value
    active_compounds_dict = {}
    inactive_compounds_dict = {}
    
    for molecule_id, activities in all_compounds_dict.items():
        avg_pchembl = sum(float(a['pchembl_value']) for a in activities) / len(activities)
        
        if avg_pchembl >= pchembl_threshold_for_download:
            if len(active_compounds_dict) < active_inactive_compounds:
                active_compounds_dict[molecule_id] = activities
        else:
            if len(inactive_compounds_dict) < active_inactive_compounds:
                inactive_compounds_dict[molecule_id] = activities
                
        # Check if we have enough compounds
        if len(active_compounds_dict) >= active_inactive_compounds and len(inactive_compounds_dict) >= active_inactive_compounds:
            break

    # Aggregate the activities
    aggregated_actives = []
    aggregated_inactives = []

    for molecule_id, activities in active_compounds_dict.items():
        avg_pchembl = sum(float(a['pchembl_value']) for a in activities) / len(activities)
        base_activity = activities[0].copy()
        base_activity['pchembl_value'] = avg_pchembl
        aggregated_actives.append(base_activity)

    for molecule_id, activities in inactive_compounds_dict.items():
        avg_pchembl = sum(float(a['pchembl_value']) for a in activities) / len(activities)
        base_activity = activities[0].copy()
        base_activity['pchembl_value'] = avg_pchembl
        aggregated_inactives.append(base_activity)

    # Create DataFrame with aggregated results
    df = pd.DataFrame(aggregated_actives + aggregated_inactives)
    if not df.empty:
        df['activity_class'] = df['pchembl_value'].apply(
            lambda x: 'active' if float(x) >= pchembl_threshold_for_download else 'inactive'
        )

    print(f"\nFound {len(aggregated_actives)} unique active and {len(aggregated_inactives)} unique inactive compounds")
    if len(aggregated_actives) < active_inactive_compounds or len(aggregated_inactives) < active_inactive_compounds:
        print(f"Notice: Could not find {active_inactive_compounds} unique compounds for both active and inactive classes")
    
    return df

def fetch_molecule_data_batch(compound_ids, session):
    """Fetch SMILES and SDF for a batch of compounds"""
    results = []
    for compound_id in compound_ids:
        url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{compound_id}.json"
        try:
            response = session.get(url)
            response.raise_for_status()
            data = response.json()
            
            molecule_structures = data.get('molecule_structures', {})
            smiles = molecule_structures.get('canonical_smiles')
            molfile = molecule_structures.get('molfile')  # This is the SDF/MOL format
            
            results.append((compound_id, smiles, molfile))
            
        except Exception as e:
            print(f"Error fetching {compound_id}: {e}")
            results.append((compound_id, None, None))
            
    return results

def check_and_download_molecule_data(compound_ids):
    print("Starting to download molecular data (SMILES and SDF)...")
    
    session = create_session()
    
    # Split compounds into batches of 50
    batch_size = 50
    batches = [compound_ids[i:i + batch_size] for i in range(0, len(compound_ids), batch_size)]
    
    molecule_data = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        future_to_batch = {
            executor.submit(fetch_molecule_data_batch, batch, session): batch 
            for batch in batches
        }
        
        for future in concurrent.futures.as_completed(future_to_batch):
            batch_results = future.result()
            molecule_data.extend([r for r in batch_results if r[1] is not None])
    
    print("Finished downloading molecular data.")
    return molecule_data

def read_chembl_ids_from_file(file_path):
    if os.path.exists(file_path):
        print(f"Reading ChEMBL IDs from {file_path}...")  # File read message
        with open(file_path, 'r') as file:
            chembl_ids = [line.strip() for line in file.readlines() if line.strip()]
            return chembl_ids
    else:
        print(f"File {file_path} does not exist.")
        return []

def fetch_all_protein_targets():
    """
    Fetches all protein targets from ChEMBL database.
    Returns a list of ChEMBL IDs for proteins.
    First checks if cached file exists, if not downloads and saves for future use.
    """
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    cache_file = os.path.join(base_dir, 'training_files', 'all_protein_targets.txt')
    
    # Check if cached file exists
    if os.path.exists(cache_file):
        print(f"Loading protein targets from cache: {cache_file}")
        with open(cache_file, 'r') as f:
            targets = [line.strip() for line in f if line.strip()]
        print(f"Loaded {len(targets)} protein targets from cache")
        return targets
    
    print("Fetching all protein targets from ChEMBL...")
    base_url = "https://www.ebi.ac.uk/chembl/api/data/target.json"
    params = {
        'target_type': 'SINGLE PROTEIN',
        'only': 'target_chembl_id'
    }
    
    targets = []
    while True:
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()
            
            data = response.json()
            if 'targets' in data:
                targets.extend([t['target_chembl_id'] for t in data['targets']])
                
            if 'page_meta' in data and data['page_meta']['next']:
                params['offset'] = data['page_meta']['offset'] + data['page_meta']['limit']
                time.sleep(0.1)  # Add delay to avoid overwhelming the API
            else:
                break
                
        except requests.exceptions.RequestException as e:
            print(f"Error fetching targets: {e}")
            break
    
    print(f"Found {len(targets)} protein targets")
    
    # Save targets to cache file
    os.makedirs(os.path.dirname(cache_file), exist_ok=True)
    with open(cache_file, 'w') as f:
        f.write('\n'.join(targets))
    print(f"Saved protein targets to: {cache_file}")
    
    return targets

def download_activity_data(target_chembl_id, assay_types, pchembl_threshold_for_download, active_inactive_compounds, output_dir: str,
                           sa ,
                           qed ,
                           logP ,
                           molecular_weight ):
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    target_chembl_ids = []
    

    if target_chembl_id:
        target_chembl_ids.extend(target_chembl_id.split(','))
    
    assay_types = assay_types.split(',')
    
    for chembl_id in target_chembl_ids:
        output_dir = output_dir / "active_inactive_compounds"
        filename = f"{chembl_id}_activity_data.csv"
        output_path = os.path.join(output_dir, filename)
        
        if os.path.exists(output_path):
            print(f"File {output_path} already exists. Skipping download.")
            continue

        data = fetch_activities([chembl_id], assay_types, pchembl_threshold_for_download, active_inactive_compounds, sa, qed, logP, molecular_weight)
        
        if not data.empty:
            compound_ids = data['molecule_chembl_id'].unique().tolist()
            molecule_data = check_and_download_molecule_data(compound_ids)
            
            if molecule_data:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                
                # Create DataFrame with both SMILES and SDF data
                molecule_df = pd.DataFrame(molecule_data, 
                                         columns=["molecule_chembl_id", "canonical_smiles", "molfile"])
                data = data.merge(molecule_df, on='molecule_chembl_id')
                data['SA_Score'] = data['canonical_smiles'].apply(lambda smiles: sa({"smiles": smiles})["SA"])
                data['QED'] = data['canonical_smiles'].apply(lambda smiles: qed({"smiles": smiles})["QED"])
                data['logP'] = data['canonical_smiles'].apply(lambda smiles: logP({"smiles": smiles})["logP"])
                data['molecular_weight'] = data['canonical_smiles'].apply(lambda smiles: molecular_weight({"smiles": smiles})["molecular_weight"])
                            # Create SDF directory
                sdf_dir = output_dir / "sdf_files" / chembl_id
                os.makedirs(sdf_dir, exist_ok=True)
                
                # Add column for SDF file paths
                data['sdf_path'] = ''
                
                # Save each molecule's SDF and store the path
                for idx, row in data.iterrows():
                    if row['molfile']:  # Check if SDF data exists
                        molecule_id = row['molecule_chembl_id']
                        activity_class = row['activity_class']
                        sdf_filename = f"{molecule_id}_{activity_class}.sdf"
                        sdf_path = sdf_dir / sdf_filename
                        with open(sdf_path, 'w') as f:
                            f.write(row['molfile'])
                        # Store relative path to make it more portable
                        data.at[idx, 'sdf_path'] = os.path.join("sdf_files", chembl_id, sdf_filename)
                
                # Save the CSV with SMILES and SDF paths
                csv_data = data[['canonical_smiles', 'activity_class', 'pchembl_value', 'SA_Score', 'QED', 'logP', 'molecular_weight', 'sdf_path']]
                csv_path = output_path
                csv_data.to_csv(csv_path, index=False)
                
                print(f"Activity data and molecular files for {chembl_id} saved to {output_dir}")
                return csv_data
            else:
                return_msg = f"No molecular data found for {chembl_id}."
                print(return_msg)
                return return_msg
        else:
            return_msg = f"No activity data found for {chembl_id}."
            print(return_msg)
            return return_msg


from src.tools.tool_definitions import tool_registry

from src.tools.mol_utils import SA, QED, logP, molecular_weight, _parse_plip_text
from rdkit import Chem

@tool_registry.register("search_chembl_activity")
def search_chembl_activity(arguments: dict, output_dir_id: str) -> dict:
    # print("ChEMBL API is being used")
    output_dir = tool_registry.output_dir / output_dir_id
    
    # Get ChEMBL data
    chembl_df = download_activity_data(
        arguments["target_chembl_id"],
        assay_types='B',
        pchembl_threshold_for_download=6,
        active_inactive_compounds=5,
        output_dir=output_dir,
        sa=SA,
        qed=QED,
        logP=logP,
        molecular_weight=molecular_weight,
    )

    # Round numeric columns
    numeric_cols = [
        "pchembl_value",
        "SA_Score",
        "QED",
        "logP",
        "molecular_weight",
    ]
    chembl_df[numeric_cols] = chembl_df[numeric_cols].round(2)

    # Rename columns
    column_mapping = {
        'SA_Score': 'sa',
        'QED': 'qed',
        'logP': 'logp',
        'molecular_weight': 'mw',
        'canonical_smiles': 'smiles'
    }
    chembl_df = chembl_df.rename(columns=column_mapping)
    
    protein_path = arguments.get("protein_path")
    if protein_path:
        try:
            # ----------  VINA-BASED DOCKING (replaces DiffDock block) ----------
            molecules = chembl_df["smiles"].tolist()

            # 1) sanitise & canonicalise SMILES so Vina (and RDKit) are happy
            sanitized_molecules = []
            smiles_map = {}                          # canonical → original
            for smi in molecules:
                try:
                    smi_to_process = smi
                    if '.' in smi:
                        smi_to_process = max(smi.split('.'), key=len)
                    
                    mol = Chem.MolFromSmiles(smi_to_process) # Ensure Chem is imported
                    if mol:
                        can = Chem.MolToSmiles(mol)
                        sanitized_molecules.append(can)
                        # Map canonical from largest fragment back to original full SMILES
                        if can not in smiles_map:
                           smiles_map[can] = smi
                except Exception:
                    # print(f"Could not sanitize SMILES: {smi}") # Optional: log failed sanitizations
                    pass                             # ignore unparsable SMILES
            
            if not sanitized_molecules:
                print("No valid molecules to dock after sanitization.")
                chembl_df["docking"] = None
                chembl_df["plip_interactions"] = [[] for _ in range(len(chembl_df))]
            else:
                # 2) run AutoDock-Vina (with batching)
                import concurrent.futures
                import threading
                import time

                docking_result_list = []
                batch_size = 10 # As in get_vina_mol_gen_report
                
                if len(sanitized_molecules) > batch_size:
                    batches = [sanitized_molecules[i:i + batch_size] for i in range(0, len(sanitized_molecules), batch_size)]
                    print(f"Processing {len(sanitized_molecules)} molecules in {len(batches)} batches of up to {batch_size}")

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
                                run_id=output_dir_id
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
                        run_id=output_dir_id
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

                # 3) normalise the result into the two maps above
                # This part processes the aggregated docking_result_list
                processed_smiles_from_docking = set()

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
                            processed_smiles_from_docking.add(orig_smi)
                    
                    # Scenario B: res_item is a combined dictionary (e.g., if a batch returns one dict for multiple SMILES)
                    # This was the structure assumed in the original non-batched version for 'docking_result'
                    elif 'binding_energy' in res_item and isinstance(res_item['binding_energy'], dict):
                        for res_smi_key, energies in res_item['binding_energy'].items():
                            orig_smi  = smiles_map.get(res_smi_key, res_smi_key)
                            if orig_smi:
                                best_energy_map[orig_smi]  = _best_energy(energies)
                                if 'merged_path' in res_item and isinstance(res_item['merged_path'], dict):
                                    complex_path_map[orig_smi] = res_item['merged_path'].get(res_smi_key)
                                elif isinstance(res_item.get('merged_path'), str): # if merged_path is single string for this combined dict
                                     complex_path_map[orig_smi] = res_item.get('merged_path')
                                processed_smiles_from_docking.add(orig_smi)


                # 4) add the Vina score column in the same order as chembl_df
                docking_vals      = []
                plip_input_paths  = []
                for smi_original in chembl_df["smiles"]: # Iterate over original SMILES from DataFrame
                    docking_vals.append(best_energy_map.get(smi_original))
                    plip_input_paths.append(complex_path_map.get(smi_original))
                
                chembl_df["docking"] = docking_vals

                # 5) run PLIP on the merged complexes produced by Vina
                # Ensure Path is imported from pathlib
                from pathlib import Path
                plip_paths   = [p for p in plip_input_paths if p and Path(p).exists()]
                plip_results_map = {} # Changed from plip_results to avoid confusion
                if plip_paths:
                    plip_out = tool_registry.execute(
                        'get_plip_report',
                        {'input_structure': plip_paths, # This should be a list of paths
                         'wait_for_result': True,
                         'max_wait_time': 600}, # Increased timeout
                        run_id=output_dir_id
                    )
                    if plip_out and 'results' in plip_out and isinstance(plip_out['results'], dict):
                        plip_results_map = plip_out['results'] # results is a dict path_str -> plip_data

                # 6) attach PLIP interactions to the dataframe
                plip_vals = []
                for complex_p_str in plip_input_paths: # These are paths corresponding to chembl_df order
                    if complex_p_str and complex_p_str in plip_results_map:
                        data = plip_results_map[complex_p_str]
                        # Assuming _parse_plip_text is defined elsewhere or needs to be included
                        if isinstance(data, dict) and 'text' in data:
                            plip_vals.append(_parse_plip_text(data['text']))
                        elif isinstance(data, list): # If PLIP data is already parsed list
                            plip_vals.append(data)
                        else: # Fallback if structure is unexpected
                            plip_vals.append([])
                    else:
                        plip_vals.append([])
                chembl_df["plip_interactions"] = plip_vals
            # ----------  END VINA-BASED DOCKING  -------------------------------

        except Exception as e:
            import traceback
            print(f"Error in docking/PLIP analysis: {str(e)}")
            print(traceback.format_exc())
            # Ensure columns are added even on failure, matching original length
            num_rows = len(chembl_df)
            chembl_df["docking"] = [None] * num_rows
            chembl_df["plip_interactions"] = [[] for _ in range(num_rows)] # list of empty lists

    # Attach friendly IDs for downstream reference
    agent_name = arguments.get("agent_name",'Database Agent')
    iteration = arguments.get("run_iteration",1) 

    friendly_ids = []
    for idx in range(len(chembl_df)):
        friendly_id = tool_registry._id_generator.generate_id(
            agent_name=agent_name, 
            molecule_number = idx,
            iteration=iteration, 
            parent_id=None  # These are original molecules from database
        )
        friendly_ids.append(friendly_id)
    
    chembl_df["friendly_id"] = friendly_ids
    # chembl_df.rename(columns={"canonical_smiles": "smiles"}, inplace=True) # Already done
    # Convert to records and return
    records = chembl_df.to_dict(orient="records")
    return {"ligands": records}

# Helper function, ensure it's defined or imported if not part of this file's global scope
# def _parse_plip_text(plip_report_text: str) -> list:
#     # Placeholder: Actual implementation would parse the PLIP text report
#     # For example, split by lines and extract relevant interaction types
#     interactions = []
#     if isinstance(plip_report_text, str):
#         for line in plip_report_text.splitlines():
#             if "interaction" in line.lower(): # Example simple parsing
#                 interactions.append(line.strip())
#     return interactions


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Download ChEMBL activity data and SMILES")
#     parser.add_argument('--target_chembl_id', type=str, help="Target ChEMBL ID(s) to search for, comma-separated")
#     parser.add_argument('--assay_type', type=str, default='B', help="Assay type(s) to search for, comma-separated")
#     parser.add_argument('--pchembl_threshold_for_download', type=float, default=6, help="Threshold for pChembl value to determine active/inactive")
#     parser.add_argument('--active_inactive_compounds', type=int, default=50, help="Number of active and inactive compounds to download")
#     args = parser.parse_args()
#     download_activity_data(args.target_chembl_id, args.assay_type, args.pchembl_threshold_for_download, args.active_inactive_compounds)
