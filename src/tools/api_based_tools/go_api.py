import requests
import json
import os
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List, Union, Literal
from tqdm import tqdm
import logging
logger = logging.getLogger(__name__)
class GOTermAPI:
    """
    API client for retrieving Gene Ontology (GO) term information
    """
    
    def __init__(self, cache_dir: str = "go_cache"):
        """
        Initialize the GO API client
        
        Args:
            cache_dir: Directory to cache GO term data
        """
        self.base_url = "https://api.geneontology.org/api"
        self.quick_go_url = "https://www.ebi.ac.uk/QuickGO/services/ontology/go"
        self.cache_dir = Path(cache_dir)
        os.makedirs(self.cache_dir, exist_ok=True)
        
    def get_term(self, go_id: str) -> Optional[Dict]:
        """
        Get information about a specific GO term
        
        Args:
            go_id: GO ID (e.g., 'GO:0006915')
            
        Returns:
            dict: Information about the GO term or None if not found
        """
        # Check cache first
        cache_file = self.cache_dir / f"{go_id.replace(':', '_')}.json"
        if cache_file.exists():
            with open(cache_file, 'r') as f:
                return json.load(f)
        
        # Format the GO ID properly (removing 'GO:' prefix if present)
        formatted_id = go_id
        if go_id.startswith("GO:"):
            formatted_id = go_id[3:]
        
        # Make API request
        url = f"{self.base_url}/ontology/term/GO:{formatted_id}"
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            # Cache the result
            with open(cache_file, 'w') as f:
                json.dump(data, f, indent=2)
                
            return data
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error fetching GO term {go_id}: {str(e)}")
            
            # Try QuickGO as fallback
            try:
                quick_go_url = f"{self.quick_go_url}/terms/GO:{formatted_id}"
                fallback_response = requests.get(quick_go_url)
                fallback_response.raise_for_status()
                fallback_data = fallback_response.json()
                
                if "results" in fallback_data and len(fallback_data["results"]) > 0:
                    # Cache the result
                    with open(cache_file, 'w') as f:
                        json.dump(fallback_data["results"][0], f, indent=2)
                    
                    return fallback_data["results"][0]
            except:
                pass
                
            return None
    
    def search_terms(
        self, 
        query: str, 
        ontology_type: Optional[Literal["biological_process", "molecular_function", "cellular_component"]] = None,
        limit: int = 20
    ) -> List[Dict]:
        """
        Search for GO terms by keyword
        
        Args:
            query: The search query
            ontology_type: Filter by specific GO aspect (biological_process, molecular_function, cellular_component)
            limit: Maximum number of results to return
            
        Returns:
            list: List of matching GO terms
        """
        # Map user-friendly names to GO aspect codes
        aspect_map = {
            "biological_process": "P",
            "molecular_function": "F", 
            "cellular_component": "C"
        }
        
        # Build the URL with parameters
        params = {
            "query": query,
            "limit": limit
        }
        
        # Try QuickGO API first as it has better search functionality
        try:
            search_url = f"{self.quick_go_url}/search"
            response = requests.get(search_url, params=params)
            response.raise_for_status()
            data = response.json()
            
            results = []
            for hit in data.get("results", []):
                # Filter by aspect if requested
                if ontology_type and hit.get("aspect") != aspect_map.get(ontology_type):
                    continue
                    
                results.append(hit)
                
                # Cache individual terms
                cache_file = self.cache_dir / f"{hit['id'].replace(':', '_')}.json"
                with open(cache_file, 'w') as f:
                    json.dump(hit, f, indent=2)
            
            return results[:limit]
        
        except requests.exceptions.RequestException as e:
            logger.warning(f"QuickGO search failed: {str(e)}")
            
            # Fall back to Gene Ontology API
            try:
                search_url = f"{self.base_url}/search/entity/autocomplete/go"
                response = requests.get(search_url, params={"q": query, "rows": limit})
                response.raise_for_status()
                data = response.json()
                
                results = []
                for doc in data.get("docs", []):
                    # Filter by aspect if requested
                    if ontology_type:
                        go_id = doc.get("id")
                        term_data = self.get_term(go_id)
                        if term_data and term_data.get("aspect") != aspect_map.get(ontology_type):
                            continue
                            
                    results.append({
                        "id": doc.get("id"),
                        "name": doc.get("annotation_class_label"),
                        "aspect": doc.get("aspect"),
                        "definition": doc.get("description")
                    })
                
                return results[:limit]
            except requests.exceptions.RequestException as e:
                logger.warning(f"Gene Ontology API search failed: {str(e)}")
                return []
    
    def get_parents(self, go_id: str) -> List[Dict]:
        """
        Get parent terms of a GO term
        
        Args:
            go_id: GO ID (e.g., 'GO:0006915')
            
        Returns:
            list: List of parent GO terms
        """
        term_data = self.get_term(go_id)
        if not term_data:
            return []
            
        # Try to extract parents from the response
        # Different APIs return different structures
        parents = []
        
        # Check if we're using QuickGO or GO API format
        if "parentIds" in term_data:
            # QuickGO format
            parent_ids = term_data.get("parentIds", [])
            for parent_id in parent_ids:
                parent_data = self.get_term(parent_id)
                if parent_data:
                    parents.append({
                        "id": parent_id,
                        "name": parent_data.get("name"),
                        "relation": "is_a"  # Default relation
                    })
        elif "parents" in term_data:
            # GO API format
            for parent in term_data.get("parents", []):
                parents.append({
                    "id": parent.get("id"),
                    "name": parent.get("label"),
                    "relation": parent.get("relation", "is_a")
                })
        
        return parents
    
    def get_children(self, go_id: str) -> List[Dict]:
        """
        Get child terms of a GO term
        
        Args:
            go_id: GO ID (e.g., 'GO:0006915')
            
        Returns:
            list: List of child GO terms
        """
        # Format the GO ID properly (removing 'GO:' prefix if present)
        formatted_id = go_id
        if go_id.startswith("GO:"):
            formatted_id = go_id[3:]
            
        url = f"{self.quick_go_url}/terms/GO:{formatted_id}/children"
        
        try:
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            results = []
            for child in data.get("results", []):
                results.append({
                    "id": child.get("id"),
                    "name": child.get("name"),
                    "relation": child.get("relation", "is_a")
                })
            
            return results
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error fetching children for {go_id}: {str(e)}")
            return []
    
    def get_annotations(self, go_id: str, taxon_id: Optional[int] = None, limit: int = 100) -> List[Dict]:
        """
        Get protein annotations for a GO term
        
        Args:
            go_id: GO ID (e.g., 'GO:0006915')
            taxon_id: Taxonomy ID to filter by organism (e.g., 9606 for human)
            limit: Maximum number of results to return
            
        Returns:
            list: List of proteins annotated with this GO term
        """
        # Format the GO ID properly
        formatted_id = go_id
        if not go_id.startswith("GO:"):
            formatted_id = f"GO:{go_id}"
            
        params = {
            "goId": formatted_id,
            "limit": limit
        }
        
        if taxon_id:
            params["taxonId"] = taxon_id
            
        url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
        
        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            data = response.json()
            
            results = []
            for annotation in data.get("results", []):
                results.append({
                    "protein_id": annotation.get("geneProductId"),
                    "protein_name": annotation.get("geneProductSimpleName"),
                    "taxon_id": annotation.get("taxonId"),
                    "taxon_name": annotation.get("taxonName"),
                    "evidence_code": annotation.get("goEvidence"),
                    "evidence_type": annotation.get("evidenceType"),
                    "reference": annotation.get("reference")
                })
            
            return results
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error fetching annotations for {go_id}: {str(e)}")
            return []
    
    def export_as_csv(self, go_ids: List[str], output_dir: str = "go_exports"):
        """
        Export GO term information and annotations to CSV files
        
        Args:
            go_ids: List of GO IDs to export
            output_dir: Directory to save CSV files
        """
        os.makedirs(output_dir, exist_ok=True)
        
        terms_data = []
        annotations_data = []
        
        for go_id in tqdm(go_ids, desc="Exporting GO terms"):
            term_data = self.get_term(go_id)
            if term_data:
                # Add term info
                terms_data.append({
                    "id": term_data.get("id"),
                    "name": term_data.get("name"),
                    "aspect": term_data.get("aspect"),
                    "definition": term_data.get("definition", {}).get("text") if "definition" in term_data else term_data.get("definition")
                })
                
                # Get annotations
                annotations = self.get_annotations(go_id)
                for annotation in annotations:
                    annotation["go_id"] = go_id
                    annotation["go_term"] = term_data.get("name")
                    annotations_data.append(annotation)
        
        # Write to CSV
        if terms_data:
            terms_df = pd.DataFrame(terms_data)
            terms_file = os.path.join(output_dir, "go_terms.csv")
            terms_df.to_csv(terms_file, index=False)
            logger.info(f"GO terms data is written: {terms_file}")
            
        if annotations_data:
            annotations_df = pd.DataFrame(annotations_data)
            annotations_file = os.path.join(output_dir, "go_annotations.csv")
            annotations_df.to_csv(annotations_file, index=False)
            logger.info(f"GO annotations data is written: {annotations_file}")


class GOEnrichmentAPI:
    """
    API client for performing GO term enrichment analysis
    """
    
    def __init__(self, cache_dir: str = "go_enrichment_cache"):
        """
        Initialize the GO Enrichment API client
        
        Args:
            cache_dir: Directory to cache enrichment results
        """
        self.base_url = "http://pantherdb.org/services/oai/pantherdb/enrich/overrep"
        self.cache_dir = Path(cache_dir)
        os.makedirs(self.cache_dir, exist_ok=True)
    
    def perform_enrichment(
        self,
        gene_list: List[str],
        organism: Literal["human", "mouse", "rat", "zebrafish", "fly", "worm", "yeast"] = "human",
        ontology: Literal["biological_process", "molecular_function", "cellular_component"] = "biological_process",
        correction: Literal["fdr", "bonferroni", "none"] = "fdr",
        reference_list: Optional[List[str]] = None
    ) -> Dict:
        """
        Perform GO term enrichment analysis on a list of genes
        
        Args:
            gene_list: List of gene symbols or IDs
            organism: Organism to use for enrichment
            ontology: GO aspect to analyze
            correction: Multiple testing correction method
            reference_list: Custom background gene list (if None, use all genes in organism)
            
        Returns:
            dict: Enrichment analysis results
        """
        # Map user-friendly names to API parameters
        organism_map = {
            "human": "9606",
            "mouse": "10090",
            "rat": "10116",
            "zebrafish": "7955",
            "fly": "7227",
            "worm": "6239",
            "yeast": "4932"
        }
        
        ontology_map = {
            "biological_process": "GO:0008150",
            "molecular_function": "GO:0003674",
            "cellular_component": "GO:0005575"
        }
        
        # Create cache key
        genes_str = ",".join(sorted(gene_list))
        ref_str = ",".join(sorted(reference_list)) if reference_list else "default"
        cache_key = f"{organism}_{ontology}_{correction}_{hash(genes_str)}_{hash(ref_str)}"
        cache_file = self.cache_dir / f"{cache_key}.json"
        
        # Check cache
        if cache_file.exists():
            with open(cache_file, 'r') as f:
                return json.load(f)
        
        # Prepare request data
        params = {
            "geneInputList": ",".join(gene_list),
            "organism": organism_map.get(organism, "9606"),
            "annotDataSet": ontology_map.get(ontology, "GO:0008150"),
            "enrichmentTestType": "FISHER",
            "correction": correction.upper()
        }
        
        if reference_list:
            params["refInputList"] = ",".join(reference_list)
        
        try:
            response = requests.post(self.base_url, data=params)
            response.raise_for_status()
            data = response.json()
            
            # Process and structure the response
            result = {
                "query": {
                    "genes": gene_list,
                    "organism": organism,
                    "ontology": ontology,
                    "correction": correction
                },
                "enrichment": []
            }
            
            if "results" in data and "result" in data["results"]:
                for item in data["results"]["result"]:
                    term_data = {
                        "go_id": item.get("term", {}).get("id"),
                        "name": item.get("term", {}).get("label"),
                        "p_value": float(item.get("pValue", 0)),
                        "fold_enrichment": float(item.get("foldEnrichment", 0)),
                        "expected": float(item.get("expectedMatchCount", 0)),
                        "observed": int(item.get("numberInTerm", 0)),
                        "genes": item.get("matchingTerms", [])
                    }
                    result["enrichment"].append(term_data)
            
            # Sort by p-value
            result["enrichment"].sort(key=lambda x: x["p_value"])
            
            # Cache the result
            with open(cache_file, 'w') as f:
                json.dump(result, f, indent=2)
            
            return result
        
        except requests.exceptions.RequestException as e:
            logger.warning(f"Error performing enrichment: {str(e)}")
            return {
                "query": {
                    "genes": gene_list,
                    "organism": organism,
                    "ontology": ontology,
                    "correction": correction
                },
                "enrichment": [],
                "error": str(e)
            }
    
    def export_as_csv(self, enrichment_result: Dict, output_file: str = "go_enrichment.csv"):
        """
        Export enrichment results to CSV
        
        Args:
            enrichment_result: Result from perform_enrichment
            output_file: Path to output CSV file
        """
        if not enrichment_result or "enrichment" not in enrichment_result:
            logger.warning("No enrichment data to export")
            return
        
        enrichment_data = enrichment_result["enrichment"]
        if not enrichment_data:
            logger.warning("Enrichment data is empty")
            return
        
        df = pd.DataFrame(enrichment_data)
        df.to_csv(output_file, index=False)
        logger.info(f"Enrichment data is written: {output_file}") 