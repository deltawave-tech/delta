import requests
import json
from typing import Optional, Dict, List

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
