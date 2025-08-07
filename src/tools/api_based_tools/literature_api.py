import json
import os
import urllib.parse
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List

import requests
from bs4 import BeautifulSoup

class UnpaywallAPI:
    def __init__(self, email: str):
        self.base_url = "https://api.unpaywall.org/v2"
        self.email = email
    
    def get_article_by_doi(self, doi: str) -> Optional[Dict]:
        """
        Get article information by DOI
        
        Args:
            doi (str): The DOI of the article
            
        Returns:
            dict: Article information or None if not found
        """
        url = f"{self.base_url}/{doi}?email={self.email}"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error fetching DOI {doi}: {str(e)}")
            return None
    
    def search_articles(self, query: str, is_oa: Optional[bool] = None, page: int = 1) -> List[Dict]:
        """
        Search for articles using keywords
        
        Args:
            query (str): Search query
            is_oa (bool, optional): Filter for open access articles
            page (int): Page number for pagination (50 results per page)
            
        Returns:
            list: List of matching articles
        """
        params = {
            'query': query,
            'email': self.email,
            'page': page
        }
        
        if is_oa is not None:
            params['is_oa'] = str(is_oa).lower()
            
        url = f"{self.base_url}/search"
        
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error searching articles: {str(e)}")
            return []

class UniProtAPI:
    def __init__(self, output_dir: Path):
        self.base_url = "https://www.ebi.ac.uk/proteins/api"
        self.headers = {"Accept": "application/json"}
        # Create uniprot_entries directory if it doesn't exist
        self.entries_dir = output_dir / "uniprot_entries"
        os.makedirs(self.entries_dir, exist_ok=True)

    def search_proteins(self, query: str, filters=None, size=3, offset=0):
        """
        Search for proteins using UniProt's REST API

        Args:
            query (str): Text search query (e.g., gene name like 'AKT1' or accession like 'P31750')
            filters (dict): Optional filters
            size (int): Number of results per page
            offset (int): Starting point for pagination

        Returns:
            dict: First matching protein entry with extracted features, or None if not found
        """
        endpoint = "https://rest.uniprot.org/uniprotkb/search"

        # Build the query string
        query_parts = []
        if query:
            # Check if query looks like an accession ID
            if any(query.startswith(prefix) and query[1:].isdigit() for prefix in ('P', 'Q', 'O')):
                query_parts.append(f"accession:{query}")
            else:
                query_parts.append(f"gene:{query}")

        # Add filters if provided
        if filters:
            if 'protein' in filters:
                query_parts.append(f"protein_name:\"{filters['protein']}\"")
            if 'gene' in filters:
                query_parts.append(f"gene:\"{filters['gene']}\"")
            if 'length' in filters:
                query_parts.append(f"length:{filters['length']}")
            if 'mass' in filters:
                query_parts.append(f"mass:{filters['mass']}")
            if 'reviewed' in filters:
                query_parts.append("reviewed:true" if filters['reviewed'] else "reviewed:false")

        # Combine all query parts
        final_query = " AND ".join(query_parts)

        params = {
            'query': final_query,
            'size': size,
            'offset': offset,
            'format': 'json'
        }

        response = requests.get(endpoint, headers=self.headers, params=params)

        if response.status_code == 200:
            results = response.json()
            if results.get('results'):
                # Get the first matching entry
                entry = results['results'][0]
                accession = entry.get('primaryAccession')

                # Extract features as before
                extracted_info = {
                    'function': None,
                    'domain': [],
                    'sequence': None,
                    'activity_regulation': None,
                    'domains_and_sites': [],
                    'interpro_ids': [],
                    'pfam_ids': [],
                    'bindingdb_ids': [],
                    'chembl_ids': [],
                    'drugbank_ids': [],
                    'pdb_entries': []
                }

                # Extract function, domain, and activity regulation from comments
                for comment in entry.get('comments', []):
                    comment_type = comment.get('commentType')
                    if comment_type == 'FUNCTION':
                        texts = comment.get('texts', [])
                        if texts:
                            extracted_info['function'] = texts[0].get('value')
                    elif comment_type == 'DOMAIN':
                        texts = comment.get('texts', [])
                        for text in texts:
                            extracted_info['domain'].append(text.get('value'))
                    elif comment_type == 'ACTIVITY REGULATION':
                        texts = comment.get('texts', [])
                        if texts:
                            extracted_info['activity_regulation'] = texts[0].get('value')

                # Extract detailed domain information from features
                for feature in entry.get('features', []):
                    feature_type = feature.get('type')
                    if feature_type in ['DOMAIN', 'REGION', 'MOTIF', 'SITE']:
                        domain_info = {
                            'type': feature_type,
                            'description': feature.get('description'),
                            'location': {
                                'start': feature.get('location', {}).get('start', {}).get('value'),
                                'end': feature.get('location', {}).get('end', {}).get('value')
                            }
                        }
                        extracted_info['domains_and_sites'].append(domain_info)
                    elif feature_type == 'BINDING':
                        # Extract binding site information
                        site_info = {
                            'description': feature.get('description', ''),
                            'ligand': feature.get('ligand', {}).get('name', ''),
                            'ligand_id': feature.get('ligand', {}).get('id', ''),
                            'ligand_label': feature.get('ligandPart', {}).get('name', ''),
                            'position': {
                                'start': feature.get('location', {}).get('start', {}).get('value'),
                                'end': feature.get('location', {}).get('end', {}).get('value'),
                            },
                            'evidence': [ev.get('evidenceCode', '') for ev in feature.get('evidences', [])]
                        }
                        # Initialize binding_sites list if it doesn't exist
                        if 'binding_sites' not in extracted_info:
                            extracted_info['binding_sites'] = []
                        extracted_info['binding_sites'].append(site_info)

                # Extract sequence
                sequence = entry.get('sequence', {})
                if sequence:
                    extracted_info['sequence'] = sequence.get('value', sequence.get('sequence'))

                # Extract database references
                for ref in entry.get('uniProtKBCrossReferences', []):
                    db = ref.get('database')
                    if db == 'InterPro':
                        extracted_info['interpro_ids'].append(ref.get('id'))
                    elif db == 'Pfam':
                        extracted_info['pfam_ids'].append(ref.get('id'))
                    elif db == 'BindingDB':
                        extracted_info['bindingdb_ids'].append(ref.get('id'))
                    elif db == 'ChEMBL':
                        extracted_info['chembl_ids'].append(ref.get('id'))
                    elif db == 'DrugBank':
                        extracted_info['drugbank_ids'].append(ref.get('id'))
                    elif db == 'PDB':
                        pdb_entry = {
                            'id': ref.get('id'),
                            'method': None,
                            'resolution': None,
                            'chains': None
                        }

                        # Extract PDB properties
                        for prop in ref.get('properties', []):
                            prop_key = prop.get('key')
                            if prop_key == 'Method':
                                pdb_entry['method'] = prop.get('value')
                            elif prop_key == 'Resolution':
                                pdb_entry['resolution'] = prop.get('value')
                            elif prop_key == 'Chains':
                                pdb_entry['chains'] = prop.get('value')
                        extracted_info['pdb_entries'].append(pdb_entry)

                # If no PDB entries found, try AlphaFold DB
                if not extracted_info['pdb_entries']:
                    alphafold_entry = self._get_alphafold_structure(accession)
                    if alphafold_entry:
                        extracted_info['pdb_entries'].append(alphafold_entry)

                # Save to cache
                file_path = os.path.join(self.entries_dir, f"{accession}.json")
                simplified_entry = {
                    'accession': accession,
                    'protein_name': entry.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
                    'extracted_features': extracted_info
                }

                # Sort and filter PDB entries before saving
                simplified_entry['extracted_features']['pdb_entries'] = self._sort_and_filter_pdb_entries(
                    simplified_entry['extracted_features']['pdb_entries']
                )

                with open(file_path, 'w') as f:
                    json.dump(simplified_entry, f, indent=2)

                return simplified_entry

        return None

    def get_protein_by_accession(self, accession):
        """
        Get protein information by UniProt accession ID

        Args:
            accession (str): UniProt accession ID (e.g., 'P31749')

        Returns:
            dict: Protein information with extracted features
        """
        endpoint = f"https://rest.uniprot.org/uniprotkb/{accession}"

        response = requests.get(endpoint, headers=self.headers)

        if response.status_code == 200:
            hit = response.json()

            # Use the same extraction logic as search_proteins
            extracted_info = {
                'function': None,
                'domain': [],
                'sequence': None,
                'activity_regulation': None,
                'domains_and_sites': [],
                'interpro_ids': [],
                'pfam_ids': [],
                'bindingdb_ids': [],
                'chembl_ids': [],
                'drugbank_ids': [],
                'pdb_entries': []
            }

            # Extract function, domain, and activity regulation from comments
            for comment in hit.get('comments', []):
                comment_type = comment.get('commentType')
                if comment_type == 'FUNCTION':
                    texts = comment.get('texts', [])
                    if texts:
                        extracted_info['function'] = texts[0].get('value')
                elif comment_type == 'DOMAIN':
                    texts = comment.get('texts', [])
                    for text in texts:
                        extracted_info['domain'].append(text.get('value'))
                elif comment_type == 'ACTIVITY REGULATION':
                    texts = comment.get('texts', [])
                    if texts:
                        extracted_info['activity_regulation'] = texts[0].get('value')

            # Extract detailed domain information from features
            for feature in hit.get('features', []):
                feature_type = feature.get('type')
                if feature_type in ['DOMAIN', 'REGION', 'MOTIF', 'SITE']:
                    domain_info = {
                        'type': feature_type,
                        'description': feature.get('description'),
                        'location': {
                            'start': feature.get('location', {}).get('start', {}).get('value'),
                            'end': feature.get('location', {}).get('end', {}).get('value')
                        }
                    }
                    extracted_info['domains_and_sites'].append(domain_info)
                elif feature_type == 'BINDING':
                    # Extract binding site information
                    site_info = {
                        'description': feature.get('description', ''),
                        'ligand': feature.get('ligand', {}).get('name', ''),
                        'ligand_id': feature.get('ligand', {}).get('id', ''),
                        'ligand_label': feature.get('ligandPart', {}).get('name', ''),
                        'position': {
                            'start': feature.get('location', {}).get('start', {}).get('value'),
                            'end': feature.get('location', {}).get('end', {}).get('value'),
                        },
                        'evidence': [ev.get('evidenceCode', '') for ev in feature.get('evidences', [])]
                    }
                    # Initialize binding_sites list if it doesn't exist
                    if 'binding_sites' not in extracted_info:
                        extracted_info['binding_sites'] = []
                    extracted_info['binding_sites'].append(site_info)

            # Extract sequence
            sequence = hit.get('sequence', {})
            if sequence:
                extracted_info['sequence'] = sequence.get('value', sequence.get('sequence'))

            # Extract database references
            for ref in hit.get('uniProtKBCrossReferences', []):
                db = ref.get('database')
                if db == 'InterPro':
                    extracted_info['interpro_ids'].append(ref.get('id'))
                elif db == 'Pfam':
                    extracted_info['pfam_ids'].append(ref.get('id'))
                elif db == 'BindingDB':
                    extracted_info['bindingdb_ids'].append(ref.get('id'))
                elif db == 'ChEMBL':
                    extracted_info['chembl_ids'].append(ref.get('id'))
                elif db == 'DrugBank':
                    extracted_info['drugbank_ids'].append(ref.get('id'))
                elif db == 'PDB':
                    pdb_entry = {
                        'id': ref.get('id'),
                        'method': None,
                        'resolution': None,
                        'chains': None
                    }
                    for prop in ref.get('properties', []):
                        prop_key = prop.get('key')
                        if prop_key == 'Method':
                            pdb_entry['method'] = prop.get('value')
                        elif prop_key == 'Resolution':
                            pdb_entry['resolution'] = prop.get('value')
                        elif prop_key == 'Chains':
                            pdb_entry['chains'] = prop.get('value')
                    extracted_info['pdb_entries'].append(pdb_entry)

            # If no PDB entries found, try AlphaFold DB
            if not extracted_info['pdb_entries']:
                alphafold_entry = self._get_alphafold_structure(accession)
                if alphafold_entry:
                    extracted_info['pdb_entries'].append(alphafold_entry)

            # Save to JSON file
            file_path = os.path.join(self.entries_dir, f"{accession}.json")
            simplified_entry = {
                'accession': accession,
                'protein_name': hit.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value'),
                'extracted_features': extracted_info
            }
            with open(file_path, 'w') as f:
                json.dump(simplified_entry, f, indent=2)

            return simplified_entry
        else:
            raise Exception(f"Error fetching protein: {response.status_code} - {response.text}")

    def get_protein_field(self, accession: str, field: str):
        """
        Get specific field from protein information by UniProt accession ID

        Args:
            accession (str): UniProt accession ID (e.g., 'P31749')
            field (str): Field to extract (e.g., 'sequence', 'function', etc.)

        Returns:
            Any: Value of the requested field
        """
        # First check if we have a cached file
        file_path = os.path.join(self.entries_dir, f"{accession}.json")

        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                data = json.load(f)
                if field in data.get('extracted_features', {}):
                    return data['extracted_features'][field]

        # If not cached or field not found, fetch from API
        protein_data = self.get_protein_by_accession(accession)
        if protein_data and 'extracted_features' in protein_data:
            return protein_data['extracted_features'].get(field)

        return None

    def get_binding_sites(self, accession: str):
        """
        Get binding site information for a protein from UniProtKB
        Returns a simple list of amino acid positions
        
        Args:
            accession (str): UniProt accession ID (e.g., 'P05067')
            
        Returns:
            list: List of binding site positions (e.g. ["A174", "Y176"])
        """
        # First, get the protein sequence
        sequence = self._get_protein_sequence(accession)
        if not sequence:
            raise Exception(f"Could not retrieve sequence for {accession}")
        
        # Check if we have cached binding site info
        file_path = os.path.join(self.entries_dir, f"{accession}.json")
        
        # If cached data exists and contains binding sites, return it
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                data = json.load(f)
                if 'binding_sites' in data:
                    return data['binding_sites']
        
        # If not cached, fetch from UniProt API
        endpoint = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        
        response = requests.get(endpoint, headers=self.headers)
        
        if response.status_code == 200:
            data = response.json()
            binding_positions = []
            
            # Extract binding site information from features
            for feature in data.get('features', []):
                if feature.get('type') == 'Binding site':
                    # Get position information
                    position = feature.get('location', {}).get('position', {})
                    begin = feature.get('location', {}).get('start', {}).get('value')
                    end = feature.get('location', {}).get('end', {}).get('value')
                    
                    # Determine position value
                    position_value = position.get('value') if position else None
                    
                    # If position_value exists, use it for both begin and end
                    if position_value:
                        begin = end = position_value
                    
                    # Get simplified position info with amino acid
                    if begin and isinstance(begin, int) and 0 < begin <= len(sequence):
                        # Get amino acid at this position (subtract 1 for 0-indexing)
                        amino_acid = sequence[begin-1]
                        binding_positions.append(f"{amino_acid}{begin}")
            
            # Save the binding sites to the cache file
            with open(file_path, 'w') as f:
                json.dump({"binding_sites": binding_positions}, f, indent=2)
            
            # Return the list directly
            return binding_positions
        else:
            raise Exception(f"Error fetching binding sites: {response.status_code} - {response.text}")
            
    def _get_protein_sequence(self, accession: str) -> str:
        """Get the amino acid sequence for a protein"""
        # Check if we have cached data
        file_path = os.path.join(self.entries_dir, f"{accession}_sequence.json")
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                cached_data = json.load(f)
                if 'sequence' in cached_data:
                    return cached_data['sequence']
        
        # If not cached, fetch from UniProt API
        api_url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
        response = requests.get(api_url, headers=self.headers)
        
        if response.status_code == 200:
            data = response.json()
            sequence = data.get('sequence', {}).get('value', '')
            
            # Save sequence to a separate cache file
            with open(file_path, 'w') as f:
                json.dump({'sequence': sequence}, f, indent=2)
            
            return sequence
        else:
            return ""

    def _sort_and_filter_pdb_entries(self, pdb_entries, max_entries=5):
        """
        Sort and filter PDB entries based on preferred criteria
        """
        def get_resolution(entry):
            # Convert resolution string to float, handling missing or invalid values
            res_str = entry.get('resolution', '')
            try:
                return float(res_str.split()[0]) if res_str else float('inf')
            except (ValueError, IndexError):
                return float('inf')

        def get_sequence_length(entry):
            # Extract sequence length from chains field
            chains = entry.get('chains', '')
            try:
                # Get the sequence length (should be same for all chains)
                range_part = chains.split('=')[1]
                start, end = map(int, range_part.split('-'))
                return end - start + 1
            except (ValueError, IndexError):
                return 0

        def get_chain_count(entry):
            # Count number of chains (fewer is better)
            chains = entry.get('chains', '')
            try:
                # Split by '=' and take the left part to get chain identifiers
                chain_ids = chains.split('=')[0]
                # Count chains that are separated by '/'
                return len(chain_ids.split('/'))
            except (ValueError, IndexError):
                return float('inf')

        # Sort entries based on multiple criteria
        sorted_entries = sorted(
            pdb_entries,
            key=lambda x: (
                0 if x.get('method') == 'X-ray' else 1,  # Prefer X-ray

                get_chain_count(x),                      # Fewer chains is better
                -get_sequence_length(x),                # Longer sequence is better
                get_resolution(x),                        # Lower resolution is better
            )
        )

        return sorted_entries[:max_entries]

    def _get_alphafold_models(self, uniprot_accession):
        """
        Get AlphaFold models for a UniProt accession
        
        Parameters:
        - uniprot_accession: UniProt accession (e.g., "Q8NFU0")
        
        Returns:
        - JSON response with model information or None if failed
        """
        url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_accession}"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            return response.json()
        except requests.exceptions.RequestException as e:
            print(f"Error getting AlphaFold models for {uniprot_accession}: {e}")
            return None

    def _download_alphafold_structure(self, uniprot_accession, file_format="pdb"):
        """
        Download AlphaFold structure file in specified format
        
        Parameters:
        - uniprot_accession: UniProt accession
        - file_format: Format to download ('pdb', 'cif', or 'bcif')
        
        Returns:
        - Path to saved file or None if download failed
        """
        # Get model information first
        models = self._get_alphafold_models(uniprot_accession)
        
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
            
            # Create pdb_files directory for AlphaFold models (same as regular PDB files)
            pdb_dir = Path(self.entries_dir).parent / "pdb_files"
            pdb_dir.mkdir(exist_ok=True)
            
            # Extract filename from URL or create one
            entry_id = model.get("entryId", f"AF-{uniprot_accession}-F1")
            version_str = f"_v{model.get('latestVersion', 4)}"
            filename = pdb_dir / f"{entry_id}-model{version_str}.{file_format}"
            
            with open(filename, 'wb') as f:
                f.write(response.content)
            
            return str(filename)
            
        except requests.exceptions.RequestException as e:
            print(f"Error downloading AlphaFold structure for {uniprot_accession}: {e}")
            return None

    def _get_alphafold_structure(self, uniprot_accession):
        """
        Get AlphaFold structure information and download PDB file
        
        Parameters:
        - uniprot_accession: UniProt accession
        
        Returns:
        - Dict with AlphaFold structure info in PDB entry format, or None if failed
        """
        # Get model information
        models = self._get_alphafold_models(uniprot_accession)
        
        if not models or len(models) == 0:
            return None
        
        # Use the first model
        model = models[0]
        
        # Download the PDB file
        pdb_file_path = self._download_alphafold_structure(uniprot_accession, "pdb")
        
        if not pdb_file_path:
            return None
        
        # Create PDB entry-like structure for consistency
        alphafold_entry = {
            'id': model.get("entryId", f"AF-{uniprot_accession}-F1"),
            'method': 'AlphaFold',
            'resolution': None,  # AlphaFold doesn't have resolution
            'chains': f"A={model.get('uniprotStart', 1)}-{model.get('uniprotEnd', 333)}",
            'source': 'alphafold',
            'pdb_file_path': pdb_file_path,
            'model_version': model.get('latestVersion', 4),
            'confidence_score': None,  # Could add confidence metrics if needed
            'organism': model.get('organismScientificName', ''),
            'gene': model.get('gene', ''),
            'uniprot_accession': uniprot_accession
        }
        
        print(f"Found AlphaFold model for {uniprot_accession}: {alphafold_entry['id']}")
        return alphafold_entry

from src.tools.tool_definitions import tool_registry

@tool_registry.register("search_uniprot")
def search_uniprot(arguments: dict, output_dir_id: str) -> dict:
    print("UniProt API is being used")
    uniprot = UniProtAPI(tool_registry.output_dir / output_dir_id)

    protein_data = uniprot.search_proteins(query=arguments["query"])
    if not protein_data or 'extracted_features' not in protein_data:
        return {"error": f"No results found for query: {arguments['query']}"}

    features = protein_data['extracted_features']
    tool_registry.previous_sequence = features.get('sequence')  # Store for later use

    # Include accession and protein_name in extracted_information
    extracted_info = {
        "accession": protein_data.get('accession'),
        "protein_name": protein_data.get('protein_name'),
        **features  # Include all the extracted features
    }

    return {
        "available_fields": list(extracted_info.keys()),
        "extracted_information": extracted_info
    }

def fetch_pubmed_data(search_term,
                      max_results=100, 
                      max_pdf_count = 10,
                      save_dir="pubmed_articles", 
                      sort_by="relevance"):
    """
    Fetch PubMed articles with abstracts and PDFs when available
    
    Optimized for performance with:
    1. Batched API requests (iCite)
    2. Parallel processing with ThreadPoolExecutor
    3. Smarter search strategy with fewer queries
    4. Time limits to prevent excessive runtime
    5. Performance metrics tracking
    
    Args:
        search_term (str): PubMed search query with AND/OR operators
        max_results (int): Maximum number of article results to fetch
        max_pdf_count (int): Maximum number of PDFs to download
        save_dir (str): Directory to save PDFs and metadata
        sort_by (str): Sort results by "relevance" or "pub_date"
        
    Returns:
        dict: Information about downloaded articles, including:
            - articles: List of article metadata
            - num_pdfs: Number of PDFs downloaded
            - num_articles: Number of articles fetched
            - performance: Performance metrics
    """
    import time
    from concurrent.futures import ThreadPoolExecutor
    # Record start time for performance monitoring
    overall_start_time = time.time()
    
    # Create separate directories for PDFs, full text, and metadata
    save_path = Path(save_dir)
    pdf_path = save_path / "pdfs"
    text_path = save_path / "texts"
    metadata_path = save_path / "metadata"
    save_path.mkdir(parents=True, exist_ok=True)
    pdf_path.mkdir(parents=True, exist_ok=True)
    text_path.mkdir(parents=True, exist_ok=True)
    metadata_path.mkdir(parents=True, exist_ok=True)

    # Initialize articles list at the function level
    final_output = {}
    articles = []
    downloaded_pmids = set()
    pdf_count = 0
    text_count = 0

    unpaywall = UnpaywallAPI("atababeyunlu36@gmail.com")  # Replace with your email
    # Setup base URLs
    search_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    fetch_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    pmc_text_base_url = "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_JSON/PMC{}/unicode"

    def construct_search_url(term, max_res, date_filter=None, prioritize_free=True):
        # Add free full text filter if requested - this prioritizes PMC articles and other free PDFs
        search_term = term
        if prioritize_free:
            search_term = f"({term}) AND free full text[filter]"
            
        url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
            f"db=pubmed"
            f"&term={urllib.parse.quote_plus(search_term)}"
            f"&retmax={max_res}"
            f"&retmode=json"
            f"&sort={sort_by}"
            f"&tool=your_tool_name"
            f"&email=your_email@example.com"
        )
        if date_filter:
            url += f"&datetype=pdat&mindate={date_filter[0]}&maxdate={date_filter[1]}"
        return url

    def get_citation_metrics(pmids):
        """Fetch citation metrics from iCite using individual requests as fallback"""
        # First try the batch endpoint, which might not be working correctly
        batch_url = "https://icite.od.nih.gov/api/pubs"
        individual_url = "https://icite.od.nih.gov/api/pubs"
        
        # Limit the number of PMIDs to a manageable number to reduce load
        pmids_to_fetch = pmids[:50]  # Just get metrics for the first 50
        pmids_list = [str(pmid) for pmid in pmids_to_fetch]
        
        try:
            # First try batch endpoint
            payload = {"pmids": pmids_list}
            response = requests.post(batch_url, json=payload)
            
            if response.status_code == 200:
                data = response.json()
                if "data" in data:
                    return data["data"]
            
            # If batch fails, fall back to individual requests with a limit
            print("Batch citation retrieval failed, falling back to individual requests...")
            citation_data = []
            
            # Only process a limited number to avoid too many requests
            for pmid in pmids_list[:10]:  # Limit to 10 requests as fallback
                try:
                    current_url = f"{individual_url}/{str(pmid)}"
                    response = requests.get(current_url)
                    if response.status_code == 200:
                        data = response.json()
                        citation_data.append(data)
                except Exception:
                    pass  # Skip individual failures
                    
                # Be respectful of the API
                time.sleep(0.1)
                
            return citation_data
            
        except Exception as e:
            print(f"Warning: Unable to fetch citation metrics: {e}")
            return []

    # Helper function for downloading PDFs (to be used with ThreadPoolExecutor)
    def download_pdf(article_info):
        nonlocal pdf_count
        
        pmid = article_info['pmid']
        title = article_info['title']
        pdf_url = article_info['pdf_url']
        
        # Skip if we already have too many PDFs
        if pdf_count >= max_pdf_count:
            return None
            
        # Skip known problematic domains that often block automated downloads
        if any(blocked in pdf_url.lower() for blocked in [
            'sciencedirect.com', 'elsevier.com', 'tandfonline.com', 
            'wiley.com/doi/pdfdirect', 'springer.com', 'aacrjournals.org',
            'thelancet.com'
        ]):
            print(f"Skipping likely-blocked publisher URL for PMID {pmid}: {pdf_url}")
            return article_info
            
        try:
            # Custom handling based on URL patterns
            if 'ncbi.nlm.nih.gov/pmc' in pdf_url:
                # PMC articles - use specific headers and handling
                headers = {
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
                    'Accept': 'application/pdf',
                    'Referer': 'https://www.ncbi.nlm.nih.gov/'
                }
                # Try alternate URL format if needed
                if 'nihms' in pdf_url:
                    pmc_id = article_info.get('pmc_id')
                    alternate_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/"
                    pdf_response = requests.get(alternate_url, headers=headers, timeout=20)
                else:
                    pdf_response = requests.get(pdf_url, headers=headers, timeout=20)
                
            elif 'mdpi.com' in pdf_url:
                # MDPI specific handling
                headers = {
                    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/15.1 Safari/605.1.15',
                    'Accept': 'application/pdf,application/x-pdf',
                    'Referer': 'https://www.mdpi.com/'
                }
                pdf_response = requests.get(pdf_url, headers=headers, timeout=20)
                
            elif 'frontiersin.org' in pdf_url:
                # Frontiers journals
                headers = {
                    'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:95.0) Gecko/20100101 Firefox/95.0',
                    'Accept': 'application/pdf,application/x-pdf',
                    'Referer': 'https://www.frontiersin.org/'
                }
                pdf_response = requests.get(pdf_url, headers=headers, timeout=20)
                
            elif 'plos' in pdf_url:
                # PLOS journals - very accessible
                headers = {
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
                    'Accept': 'application/pdf,application/x-pdf',
                    'Referer': 'https://journals.plos.org/'
                }
                pdf_response = requests.get(pdf_url, headers=headers, timeout=20)
                
            else:
                # Default handling for other publishers
                headers = {
                    'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/91.0.4472.124 Safari/537.36',
                    'Accept': 'application/pdf,application/x-pdf,text/html,application/xhtml+xml,application/xml',
                    'Accept-Language': 'en-US,en;q=0.9',
                    'Referer': 'https://pubmed.ncbi.nlm.nih.gov/'
                }
                pdf_response = requests.get(pdf_url, headers=headers, timeout=20)
            
            pdf_response.raise_for_status()
            
            # Check if the response is actually a PDF (by content type or first few bytes)
            content_type = pdf_response.headers.get('Content-Type', '').lower()
            if not ('pdf' in content_type or pdf_response.content[:4] == b'%PDF'):
                print(f"Response for PMID {pmid} is not a PDF (content type: {content_type})")
                article_info['has_pdf'] = False
                article_info['pdf_url'] = None
                return article_info

            # Save PDF file in pdf subdirectory
            safe_title = title[:50] if title else str(pmid)
            pdf_filename = f"{pmid}_{sanitize_filename(safe_title)}.pdf"
            pdf_file_path = pdf_path / pdf_filename
            
            with open(pdf_file_path, 'wb') as f:
                f.write(pdf_response.content)
                
            print(f"Successfully downloaded PDF for PMID: {pmid}")
            
            # Update article info with successful download
            article_info['has_pdf'] = True
            return article_info
            
        except Exception as e:
            print(f"Failed to download PDF for PMID {pmid}: {str(e)}")
            article_info['has_pdf'] = False
            article_info['pdf_url'] = None
            return article_info

    # Helper function to get PDF URL from DOI using Unpaywall (for ThreadPoolExecutor)
    def get_pdf_url(article_info):
        doi = article_info.get('doi')
        pmc_id = article_info.get('pmc_id')
        
        # First try direct PMC URL if available - these have highest success rate
        if pmc_id:
            # Construct direct PMC PDF URL
            pmc_pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/PMC{pmc_id}/pdf/nihms-{pmc_id}.pdf"
            article_info['pdf_url'] = pmc_pdf_url
            article_info['pdf_source'] = 'pmc_direct'
            article_info['has_pdf'] = True
            return article_info
        
        # If no PMC ID or DOI, can't get PDF
        if not doi:
            return article_info
            
        # Skip DOIs with potential formatting issues
        if ' ' in doi or '/' not in doi or len(doi) > 50:
            return article_info
            
        try:
            unpaywall_data = unpaywall.get_article_by_doi(doi)
            if unpaywall_data:
                # First try to get the direct PDF URL
                pdf_url = (unpaywall_data.get('best_oa_location') or {}).get('url_for_pdf')
                
                # If no direct PDF URL, try the landing page URL as fallback
                if not pdf_url:
                    pdf_url = (unpaywall_data.get('best_oa_location') or {}).get('url')
                    
                article_info['pdf_url'] = pdf_url
                article_info['pdf_source'] = 'unpaywall'
                article_info['has_pdf'] = bool(pdf_url)
            return article_info
        except Exception as e:
            # Don't log 404 errors which are expected for many DOIs
            if "404" not in str(e):
                print(f"Error fetching Unpaywall data for DOI {doi}: {str(e)}")
            return article_info

    def get_pmc_full_text(pmc_id: str) -> str | None:
        """Fetch full text of an article from PubMed Central."""
        try:
            resp = requests.get(pmc_text_base_url.format(pmc_id))
            resp.raise_for_status()
            article = resp.json()
            if len(article) != 1 or len(article[0]["documents"]) != 1:
                return None
            document = article[0]["documents"][0]
            full_text = "\n".join(
                passage["text"]
                for passage in document["passages"]
                if passage["infons"].get("section_type") != "REF"
            )
            return full_text
        except Exception as e:  # pragma: no cover - network failure
            print(f"Error fetching PMC text for {pmc_id}: {e}")
            return None

    def download_full_text(article_info):
        nonlocal text_count

        pmc_id = article_info.get("pmc_id")
        if not pmc_id or text_count >= max_pdf_count:
            return article_info

        text = get_pmc_full_text(pmc_id)
        if not text:
            return article_info

        safe_title = article_info.get("title", "")[:50] if article_info.get("title") else str(pmc_id)
        txt_filename = f"{pmc_id}_{sanitize_filename(safe_title)}.txt"
        txt_file_path = text_path / txt_filename
        try:
            with open(txt_file_path, "w", encoding="utf-8") as f:
                f.write(text)
            article_info["has_full_text"] = True
            article_info["text_path"] = str(txt_file_path)
            text_count += 1
        except Exception as e:  # pragma: no cover - disk failure
            print(f"Failed to save full text for PMC {pmc_id}: {e}")
        return article_info
            
    def process_batch(pmids):
        nonlocal pdf_count, downloaded_pmids, articles
        
        if pdf_count >= max_pdf_count:
            return pdf_count
            
        start_batch_time = time.time()
        print(f"Processing batch of {len(pmids)} PMIDs...")
        
        # Process in smaller chunks to avoid "Request-URI Too Long" errors
        # PubMed recommends no more than ~200 IDs per request
        chunk_size = 100
        all_articles_data = []
        
        for i in range(0, len(pmids), chunk_size):
            chunk_pmids = pmids[i:i+chunk_size]
            print(f"Processing chunk {i//chunk_size + 1}/{(len(pmids) + chunk_size - 1)//chunk_size} ({len(chunk_pmids)} PMIDs)...")
            
            # Use a single request to fetch metadata for each chunk of PMIDs
            fetch_params = {
                "db": "pubmed",
                "id": ",".join(chunk_pmids),
                "retmode": "xml",
                "tool": "your_tool_name",
                "email": "your_email@example.com"
            }
            
            fetch_response = requests.get(fetch_base_url, params=fetch_params)
            fetch_response.raise_for_status()
            
            chunk_soup = BeautifulSoup(fetch_response.content, 'xml')
            all_articles_data.append(chunk_soup)
            
            # Don't overwhelm the server
            if i + chunk_size < len(pmids):
                time.sleep(0.5)  
        
        # Process each chunk instead of trying to combine them
        all_articles = []
        for chunk_soup in all_articles_data:
            all_articles.extend(chunk_soup.find_all('PubmedArticle'))
            
        print(f"Total articles found across all chunks: {len(all_articles)}")
        
        # Fetch citation metrics in a single batch
        citation_data = {str(item['pmid']): item for item in get_citation_metrics(pmids)}
        
        # Extract basic metadata from all articles first
        candidate_articles = []
        
        for article in all_articles:
            if pdf_count >= max_pdf_count:
                break
                
            pmid = article.find('PMID').text if article.find('PMID') else None
            
            # Skip if already processed
            if pmid in downloaded_pmids:
                continue
                
            downloaded_pmids.add(pmid)
            
            # Extract DOI and PMC ID
            doi = None
            pmc_id = None
            article_ids = article.find_all('ArticleId')
            
            for id_elem in article_ids:
                if id_elem.get('IdType') == 'doi':
                    doi = id_elem.text
                elif id_elem.get('IdType') == 'pmc':
                    pmc_id = id_elem.text.replace('PMC', '')
            
            # Extract title
            title = article.find('ArticleTitle').text if article.find('ArticleTitle') else "Untitled"
            
            # Extract journal info
            journal = article.find('Journal')
            journal_title = journal.find('Title').text if journal and journal.find('Title') else None
            
            # Extract publication date
            pub_date = article.find('PubDate')
            year = pub_date.find('Year').text if pub_date and pub_date.find('Year') else None
            
            # Extract abstract
            abstract = None
            abstract_elem = article.find('Abstract')
            if abstract_elem:
                abstract_parts = []
                abstract_texts = abstract_elem.find_all('AbstractText')
                for abstract_text in abstract_texts:
                    text = abstract_text.text.strip()
                    abstract_parts.append(text)
                abstract = ' '.join(abstract_parts)
            
            # Skip full text retrieval from PMC since we only care about PDFs/abstracts
            content = abstract
            
            # Initialize article data (we'll add PDF info later)
            article_data = {
                'pmid': pmid,
                'pmc_id': pmc_id,
                'doi': doi,
                'title': title,
                'journal': journal_title,
                'year': year,
                'content': content,
                'has_full_text': False,  # We're skipping full text
                'pdf_url': None,
                'has_pdf': False
            }
            
            # Add citation metrics if available
            if pmid in citation_data:
                metrics = citation_data[pmid]
                article_data.update({
                    'citation_count': metrics.get('citation_count', 0),
                    'relative_citation_ratio': metrics.get('relative_citation_ratio', None),
                    'citations_per_year': metrics.get('citations_per_year', 0)
                })
            else:
                article_data.update({
                    'citation_count': 0,
                    'relative_citation_ratio': None,
                    'citations_per_year': 0
                })
                
            # Only add to candidates if there's DOI (we need it for PDF)
            if doi:
                candidate_articles.append(article_data)
                
        # Use ThreadPoolExecutor to parallelize Unpaywall requests
        with ThreadPoolExecutor(max_workers=10) as executor:
            # First, get all PDF URLs in parallel
            articles_with_urls = list(executor.map(get_pdf_url, candidate_articles))
            
            # Filter to only articles with PDF URLs
            articles_with_pdf_urls = [a for a in articles_with_urls if a.get('pdf_url')]
            print(f"Found {len(articles_with_pdf_urls)} articles with PDF URLs")
            
            # Download PDFs in parallel (limited to what we need)
            needed_pdfs = max_pdf_count - pdf_count
            download_candidates = articles_with_pdf_urls[:needed_pdfs*2]  # Get 2x what we need in case some fail
            
            if download_candidates:
                futures = []
                for article in download_candidates:
                    futures.append(executor.submit(download_pdf, article))
                
                # Process results as they complete
                successful_downloads = []
                for future in futures:
                    try:
                        result = future.result()
                        if result and result.get('has_pdf'):
                            successful_downloads.append(result)
                            pdf_count += 1
                            if pdf_count >= max_pdf_count:
                                break
                    except Exception as e:
                        print(f"Error in PDF download thread: {str(e)}")
                
                # Update articles list with successful downloads first
                articles.extend(successful_downloads)

            # Attempt to fetch full text for articles without PDFs
            remaining_articles = [a for a in articles_with_urls if not a.get('has_pdf')]
            if remaining_articles:
                futures = [executor.submit(download_full_text, art) for art in remaining_articles]
                updated = [f.result() for f in futures]
            else:
                updated = remaining_articles

            for art in updated:
                articles.append(art)
                if len(articles) >= max_results:
                    break
        
        # Save metadata for all articles
        for article_data in articles:
            if article_data.get('content'):
                safe_title = article_data.get('title', '')[:50] if article_data.get('title') else str(article_data.get('pmid'))
                filename = f"{article_data.get('pmid')}_{sanitize_filename(safe_title)}.json"
                file_path = metadata_path / filename
                with open(file_path, 'w', encoding='utf-8') as f:
                    json.dump(article_data, f, indent=2, ensure_ascii=False)
        
        # Filter and sort articles by citation count
        articles.sort(key=lambda x: x.get('citation_count', 0), reverse=True)
        
        batch_time = time.time() - start_batch_time
        print(f"Batch processing completed in {batch_time:.1f} seconds, found {pdf_count} PDFs")
        
        return pdf_count

    try:
        import time
        start_time = time.time()
        TIME_LIMIT = 180  # 3-minute time limit
        
        print(f"Searching for up to {max_pdf_count} PDFs related to: {search_term}")
        
        # Use a larger retmax to get more results in a single query
        initial_retmax = 500  # Get a larger batch in one go
        current_year = datetime.now().year
        
        # First try recent papers (2020-current)
        search_url = construct_search_url(search_term, initial_retmax, (2020, current_year))
        search_response = requests.get(search_url)
        search_response.raise_for_status()
        search_data = search_response.json()
        pmids = search_data.get('esearchresult', {}).get('idlist', [])
        
        total_results = int(search_data.get('esearchresult', {}).get('count', 0))
        print(f"Found {total_results} results from 2020-{current_year}")
        
        # If we got very few results, try searching all years
        if len(pmids) < 50:  # Arbitrary threshold to decide if we need more results
            print(f"Finding more results by searching all years...")
            search_url = construct_search_url(search_term, initial_retmax)
            search_response = requests.get(search_url)
            search_response.raise_for_status()
            search_data = search_response.json()
            pmids = search_data.get('esearchresult', {}).get('idlist', [])
            total_results = int(search_data.get('esearchresult', {}).get('count', 0))
            print(f"Found {total_results} results from all years")
        
        if not pmids:
            print("No results found")
            return {"articles": [], "num_pdfs": 0, "num_articles": 0}
            
        # Process all PMIDs at once
        process_batch(pmids)
        
        # Check if we've reached the time limit
        if time.time() - start_time > TIME_LIMIT:
            print(f"Exceeded {TIME_LIMIT} second time limit, stopping")
        elif pdf_count < max_pdf_count and total_results > len(pmids):
            # We have more results available and haven't reached our target PDF count
            print(f"Fetched {pdf_count} PDFs but need {max_pdf_count}. Fetching more results...")
            # Get the next batch of results
            next_retmax = 500
            next_retstart = len(pmids)
            
            # Use date filter if it worked well initially
            date_filter = (2020, current_year) if len(pmids) >= 50 else None
            search_url = construct_search_url(search_term, next_retmax, date_filter)
            search_url += f"&retstart={next_retstart}"
            
            search_response = requests.get(search_url)
            search_response.raise_for_status()
            search_data = search_response.json()
            next_pmids = search_data.get('esearchresult', {}).get('idlist', [])
            
            if next_pmids:
                process_batch(next_pmids)
        # Calculate total execution time
        total_time = time.time() - overall_start_time
        
        # Create performance summary
        final_output = {
            "articles": articles,
            "num_pdfs": pdf_count,
            "num_full_texts": text_count,
            "num_articles": len(articles),
            "performance": {
                "total_time_seconds": total_time,
                "avg_time_per_pdf": total_time / pdf_count if pdf_count > 0 else 0,
                "pdfs_per_minute": (pdf_count / total_time) * 60,
            }
        }
        
        # Print performance summary
        print(f"\n{'='*50}")
        print(f"PERFORMANCE SUMMARY:")
        print(f"{'='*50}")
        print(f"Total execution time: {total_time:.2f} seconds")
        print(f"PDFs downloaded: {pdf_count}")
        print(f"Full texts downloaded: {text_count}")
        print(f"Articles processed: {len(articles)}")
        print(f"Average time per PDF: {final_output['performance']['avg_time_per_pdf']:.2f} seconds")
        print(f"PDFs per minute: {final_output['performance']['pdfs_per_minute']:.2f}")
        print(f"{'='*50}")
        
        return final_output

    except requests.exceptions.RequestException as e:
        print(f"Error fetching data: {str(e)}")
        total_time = time.time() - overall_start_time
        print(f"Execution failed after {total_time:.2f} seconds")
        return {
            "articles": articles,
            "num_pdfs": pdf_count,
            "num_full_texts": text_count,
            "num_articles": len(articles),
            "error": str(e)
        }

def sanitize_filename(filename):
    """Remove invalid characters from filename."""
    return "".join(c for c in filename if c.isalnum() or c in (' ', '-', '_')).strip()


from paperqa import Settings, ask
from paperqa.agents.main import agent_query
from paperqa.agents.models import QueryRequest
import asyncio
import nest_asyncio
nest_asyncio.apply()
import logging

@tool_registry.register("search_pubmed_agent")
def search_pubmed_agent(arguments: dict, output_dir_id: str) -> dict:
    print("PubMed API is being used")
    save_dir = tool_registry.output_dir / output_dir_id / "pubmed_articles"
    print(f"Save dir: {save_dir}")
    content = fetch_pubmed_data(
        arguments["search_term"],
        max_pdf_count=arguments.get("max_pdf_count", 20),
        sort_by="relevance",
        save_dir=save_dir,
    )

    total_docs = content.get("num_pdfs", 0) + content.get("num_full_texts", 0)
    if total_docs > 0:
        print(f"Documents found: {total_docs}")
        
        # Instead of creating a new event loop, use asyncio.run which properly manages the event loop
        try:
            # asyncio.run handles creating and closing the event loop properly
            answer = asyncio.run(agent_query(
                QueryRequest(   
                    query=arguments['query'],
                    settings=Settings(temperature=0.5, paper_directory=save_dir, concurrency=1, llm = 'claude-3-7-sonnet-20250219', summary_llm='claude-3-7-sonnet-20250219'),
                )
            ))
            
            print(f'-'*100)
            logging.info('Answer received')
            print(answer.session.answer)
            print(f'-'*100)
            answer2return = f'Here are the summary of the relevants chunks: {answer.session.answer}'
            return {
                "answer": answer2return,
            }
        except Exception as e:
            print(f"Error during agent query: {str(e)}")
            return {
                "error": f"Error processing query: {str(e)}"
            }
    else:
        print("No documents found")
        return {
            "error": "No documents found"
        }

