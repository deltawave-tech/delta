from socket import socket
from typing import Dict, Any
from ..models.molecule import MoleculeResponse
from ..utils.json_utils import convert_keys
from ..utils.socket_utils import send_request

def validate_molecule(sock: socket, session_id: str, molecule_data: dict) -> MoleculeResponse:
    """
    Validate a molecule against the schema.
    
    Args:
        sock: Socket connection
        session_id: ID of the session
        molecule_data: Molecule data to validate
        
    Returns:
        MoleculeResponse or raises error if validation fails
    """
    # Use the provided molecule data
    request = {
        "tag": "ValidateMolecule",
        "contents": [session_id, molecule_data]
    }
    response = send_request(sock, request)
    if isinstance(response, dict) and response.get("tag") == "Success":
        converted_contents = convert_keys(response["contents"])
        return MoleculeResponse(**converted_contents)
    raise RuntimeError(str(response))

def create_empty_molecule() -> Dict[str, Any]:
    """Helper function to create empty molecule structure"""
    return {
        "structure": {},
        "properties": {},
        "computedProperties": {}
    }
