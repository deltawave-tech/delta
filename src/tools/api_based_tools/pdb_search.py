import requests
import json
import os
from pathlib import Path

def search_pdb(query_text, search_type="gene_name"):
    """
    Search PDB database using text search with configurable search type
    
    Parameters:
    - query_text: Text to search for
    - search_type: Type of search to perform. Options include:
        - "gene_name": Search by gene name
        - "text": General text search
        - "structure_id": Search by PDB ID
    """
    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    
    # Define different search parameters based on search type
    search_params = {
        "gene_name": {
            "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
            "operator": "exact_match"
        },
        "text": {
            "attribute": "text",
            "operator": "contains_words"
        },
        "structure_id": {
            "attribute": "rcsb_id",
            "operator": "exact_match"
        }
    }
    
    # Get search parameters for the requested search type
    params = search_params.get(search_type, search_params["text"])
    
    # Construct the query
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": params["attribute"],
                "operator": params["operator"],
                "value": query_text
            }
        },
        "return_type": "entry"
    }

    # Make the POST request
    response = requests.post(url, json=query)
    
    if response.status_code != 200:
        print(f"Error: {response.status_code}")
        return None
    
    return response.json()

def get_structure_details(pdb_id):
    """
    Get additional details about a specific PDB structure and save selected fields to JSON
    """
    url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    response = requests.get(url)
    
    if response.status_code != 200:
        print(f"Error getting details for {pdb_id}: {response.status_code}")
        return None
    
    full_details = response.json()
    
    # Get uniprot_id using existing function
    uniprot_id = get_uniprot_id(pdb_id)
    
    # Extract only the desired fields
    filtered_details = {
        'rcsb_id': full_details.get('rcsb_id'),
        'uniprot_id': uniprot_id,
        'struct': full_details.get('struct', {}),
        'citation': {
            'pdbx_database_id_doi': full_details.get('citation', [{}])[0].get('pdbx_database_id_doi'),
            'title': full_details.get('citation', [{}])[0].get('title'),
            'year': full_details.get('citation', [{}])[0].get('year')
        },
        'rcsb_entry_container_identifiers': {
            'pubmed_id': full_details.get('rcsb_entry_container_identifiers', {}).get('pubmed_id')
        }
    }
    
    # Create pdb_files directory if it doesn't exist
    os.makedirs('pdb_files', exist_ok=True)
    
    # Save filtered details to JSON file in pdb_files folder
    filename = os.path.join('pdb_files', f"{pdb_id}.json")
    with open(filename, 'w') as f:
        json.dump(filtered_details, f, indent=4)
    
    return filtered_details

def get_structure_file(pdb_id, file_format="pdb", output_dir: str = ''):
    """
    Download structure file in specified format and save to pdb_files directory
    
    Parameters:
    - pdb_id: PDB ID of the structure
    - file_format: Format to download ('pdb', 'cif', or 'xml')
    
    Returns:
    - Path to saved file or None if download failed
    """
    # Format-specific file extensions
    format_extensions = {
        "pdb": "pdb",
        "cif": "cif",
        "xml": "xml"
    }
    
    if file_format not in format_extensions:
        print(f"Unsupported file format: {file_format}")
        return None
    
    # Construct download URL
    url = f"https://files.rcsb.org/download/{pdb_id}.{format_extensions[file_format]}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise exception for bad status codes
        
        # Create pdb_files directory if it doesn't exist
        output_dir = output_dir/'pdb_files'
        os.makedirs(output_dir, exist_ok=True)
        # Save file
        filename = os.path.join(output_dir, f"{pdb_id}.{format_extensions[file_format]}")
        with open(filename, 'wb') as f:
            f.write(response.content)
        
        return filename
        
    except requests.exceptions.RequestException as e:
        print(f"Error downloading structure file for {pdb_id}: {e}")
        return None

def get_uniprot_id(pdb_id):
    """
    Get UniProt ID for a given PDB structure
    
    Parameters:
    - pdb_id: PDB ID of the structure
    
    Returns:
    - UniProt ID or None if not found
    """
    # First get entity information to find entity_id
    entity_url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
    entity_response = requests.get(entity_url)
    
    if entity_response.status_code != 200:
        print(f"Error getting entity details for {pdb_id}: {entity_response.status_code}")
        return None
    
    # Get the first entity ID (usually 1 for single-protein structures)
    entity_id = "1"  # Default to 1 if not found
    
    # Get UniProt information
    uniprot_url = f"https://data.rcsb.org/rest/v1/core/uniprot/{pdb_id}/{entity_id}"
    uniprot_response = requests.get(uniprot_url)
    
    if uniprot_response.status_code != 200:
        print(f"Error getting UniProt details for {pdb_id}: {uniprot_response.status_code}")
        return None
    
    uniprot_data = uniprot_response.json()
    
    # Extract UniProt ID from the response
    if uniprot_data and len(uniprot_data) > 0:
        uniprot_id = uniprot_data[0].get('rcsb_uniprot_container_identifiers', {}).get('uniprot_id')
        return uniprot_id
    
    return None

def get_alphafold_models(uniprot_accession, sequence_checksum=None):
    """
    Get AlphaFold models for a UniProt accession
    
    Parameters:
    - uniprot_accession: UniProt accession (e.g., "Q8NFU0")
    - sequence_checksum: Optional CRC64 checksum of the UniProt sequence
    
    Returns:
    - JSON response with model information or None if failed
    """
    url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}"
    
    params = {}
    if sequence_checksum:
        params['sequence_checksum'] = sequence_checksum
    
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        print(f"Error getting AlphaFold models for {uniprot_accession}: {e}")
        return None

def download_alphafold_structure(uniprot_accession, file_format="pdb", output_dir="", version=None):
    """
    Download AlphaFold structure file in specified format
    
    Parameters:
    - uniprot_accession: UniProt accession
    - file_format: Format to download ('pdb', 'cif', or 'bcif')
    - output_dir: Directory to save the file
    - version: Specific version to download (uses latest if None)
    
    Returns:
    - Path to saved file or None if download failed
    """
    # Get model information first
    models = get_alphafold_models(uniprot_accession)
    
    if not models or len(models) == 0:
        print(f"No AlphaFold models found for {uniprot_accession}")
        return None
    
    # Use the first model (should be the main one)
    model = models[0]
    
    # Map file formats to URLs
    format_urls = {
        "pdb": model.get("pdbUrl"),
        "cif": model.get("cifUrl"), 
        "bcif": model.get("bcifUrl")
    }
    
    if file_format not in format_urls:
        print(f"Unsupported file format: {file_format}")
        return None
    
    download_url = format_urls[file_format]
    if not download_url:
        print(f"No {file_format} URL found for {uniprot_accession}")
        return None
    
    try:
        response = requests.get(download_url)
        response.raise_for_status()
        
        # Create output directory if it doesn't exist
        if output_dir:
            output_dir = Path(output_dir) / "af2_db_entries"
            output_dir.mkdir(parents=True, exist_ok=True)
        else:
            output_dir = Path("af2_db_entries")
            output_dir.mkdir(exist_ok=True)
        
        # Extract filename from URL or create one
        entry_id = model.get("entryId", f"AF-{uniprot_accession}-F1")
        version_str = f"_v{model.get('latestVersion', 4)}"
        filename = output_dir / f"{entry_id}-model{version_str}.{file_format}"
        
        with open(filename, 'wb') as f:
            f.write(response.content)
        
        return str(filename)
        
    except requests.exceptions.RequestException as e:
        print(f"Error downloading AlphaFold structure for {uniprot_accession}: {e}")
        return None

def save_alphafold_metadata(uniprot_accession, output_dir=""):
    """
    Get and save AlphaFold model metadata to JSON file
    
    Parameters:
    - uniprot_accession: UniProt accession
    - output_dir: Directory to save the metadata file
    
    Returns:
    - Path to saved JSON file or None if failed
    """
    models = get_alphafold_models(uniprot_accession)
    
    if not models:
        return None
    
    # Create output directory if it doesn't exist
    if output_dir:
        output_dir = Path(output_dir) / "af2_db_entries"
        output_dir.mkdir(parents=True, exist_ok=True)
    else:
        output_dir = Path("af2_db_entries")
        output_dir.mkdir(exist_ok=True)
    
    # Save metadata to JSON file
    filename = output_dir / f"{uniprot_accession}_alphafold_metadata.json"
    with open(filename, 'w') as f:
        json.dump(models, f, indent=4)
    
    return str(filename)

def main():
    # Search for AKT1 structures
    results = search_pdb("AKT1")
    
    if results and 'result_set' in results:
        print(f"Found {len(results['result_set'])} structures")
        
        # Print details for first 5 structures and save files
        for item in results['result_set'][:5]:
            pdb_id = item['identifier']
            
            # Get and save structure details JSON
            details = get_structure_details(pdb_id)
            
            # Get UniProt ID
            uniprot_id = get_uniprot_id(pdb_id)
            
            # Download PDB file
            pdb_file = get_structure_file(pdb_id, "pdb")
            
            if details:
                print(f"\nPDB ID: {pdb_id}")
                print(f"UniProt ID: {uniprot_id if uniprot_id else 'N/A'}")
                print(f"Structure file saved to: {pdb_file}")
                print(f"JSON details saved to: {pdb_id}.json")
                print("Title:", details.get('struct', {}).get('title', 'N/A'))
                print("Method:", details.get('exptl', [{}])[0].get('method', 'N/A'))
                print("Resolution:", details.get('refine', [{}])[0].get('ls_d_res_high', 'N/A'))

from src.tools.tool_definitions import tool_registry

@tool_registry.register("get_pdb_file")
def get_pdb_file(arguments: dict, output_dir_id: str) -> dict:
    print("PDB file download is being used")
    os.makedirs(tool_registry.output_dir/output_dir_id, exist_ok=True)

    temp_file = get_structure_file(
        pdb_id=arguments["pdb_id"],
        file_format=arguments.get("file_format", "pdb"),
        output_dir=tool_registry.output_dir/output_dir_id
    )

    if not temp_file:
        return {
            "success": False,
            "error": f"Failed to download PDB file for {arguments['pdb_id']}"
        }

    # target_file = tool_registry.pdb_dir / Path(temp_file).name
    # shutil.move(temp_file, target_file)
    return {
        "success": True,
        "file_path": str(temp_file)
    }

@tool_registry.register("af_pdb_file")
def get_alphafold_file(arguments: dict, output_dir_id: str) -> dict:
    """
    Download AlphaFold structure file for a given UniProt accession
    
    Arguments:
    - accession: UniProt accession (required)
    - file_format: Format to download ('pdb', 'cif', 'bcif') (optional, defaults to 'pdb')
    - include_metadata: Whether to also save metadata JSON (optional, defaults to True)
    """
    print("AlphaFold file download is being used")
    os.makedirs(tool_registry.output_dir/output_dir_id, exist_ok=True)
    
    uniprot_accession = arguments.get("accession")
    if not uniprot_accession:
        return {
            "success": False,
            "error": "uniprot_accession is required"
        }
    
    file_format = arguments.get("file_format", "pdb")
    include_metadata = arguments.get("include_metadata", True)
    
    output_dir = tool_registry.output_dir / output_dir_id
    
    # Download structure file
    structure_file = download_alphafold_structure(
        uniprot_accession=uniprot_accession,
        file_format=file_format,
        output_dir=output_dir
    )
    
    if not structure_file:
        return {
            "success": False,
            "error": f"Failed to download AlphaFold structure for {uniprot_accession}"
        }
    
    result = {
        "success": True,
        "structure_file": structure_file,
        "accession": uniprot_accession,
        "file_format": file_format
    }
    
    # Optionally save metadata
    if include_metadata:
        metadata_file = save_alphafold_metadata(uniprot_accession, output_dir)
        if metadata_file:
            result["metadata_file"] = metadata_file
    
    return result