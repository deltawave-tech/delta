from socket import socket
from typing import Dict, Any, Optional, List, Union
from ..models.molecule import MoleculeResponse
from ..models.transform import TransformRequest, MethodDetails
from ..utils.json_utils import convert_keys, convert_keys_to_camel
from ..utils.socket_utils import send_request

def apply_transform(
    sock: socket,
    session_id: str,
    transform_type: str,
    molecule_data: Union[dict, List[dict]],
    agent: str,
    rationale: str,
    parent_transformation_id: Optional[str] = None,
    user_message: Optional[str] = None,
    method_details: Optional[Dict[str, Any]] = None,
    protein_sequence: Optional[str] = None,
    protein_path: Optional[str] = None,
    iteration: Optional[int] = None,
    user_id: Optional[str] = None
) -> Dict[str, Any]:
    """
    Apply a transformation to one or more molecules in a session.

    Note: Each molecule must have its own rationale field to capture the reasoning
    for that specific molecule.

    Args:
        sock: Socket connection
        session_id: ID of the session
        transform_type: Type of transformation to apply
        molecule_data: Single molecule data dict or list of molecule data dicts (each should have a rationale)
        agent: Name of the agent performing the transformation
        rationale: Rationale for the transformation
        parent_transformation_id: Optional ID of the parent transformation
        user_message: Optional user message
        method_details: Optional method details
        protein_sequence: Optional protein sequence to store with the transformation
        protein_path: Optional protein path to store with the transformation
        iteration: Optional iteration number for the transformation
        user_id: Optional user ID to associate with the transformation (helps with retrieval)

    Returns:
        Dictionary containing transformationId and moleculeIds
    """
    # Handle single molecule or list of molecules
    if not isinstance(molecule_data, list):
        molecules = [molecule_data or {
            "structure": {},
            "properties": {},
            "computed_properties": {}
        }]
    else:
        molecules = molecule_data

    # Ensure each molecule has the required structure and verify rationale exists
    for i, mol in enumerate(molecules):
        if not mol:
            molecules[i] = {
                "structure": {},
                "properties": {},
                "computed_properties": {}
            }

        # Verify each molecule has a rationale
        if not molecules[i].get("rationale"):
            raise ValueError(f"Molecule at index {i} is missing a rationale field. Each molecule must have its own rationale.")

    # Create the transform request using snake_case for consistency
    transform_request = {
        "type": transform_type,
        "agent": agent,
        "user_message": user_message,
        "method_details": method_details,
        "rationale": rationale,
        "protein_sequence": protein_sequence,
        "protein_path": protein_path,
        "iteration": iteration
    }
    
    # Convert to camelCase before sending to the API
    transform_request = convert_keys_to_camel(transform_request)

    parent_id = parent_transformation_id

    # Convert molecule data to camelCase as well
    camel_case_molecules = [convert_keys_to_camel(mol) for mol in molecules]
    
    # Add user_id to the request if provided
    transform_request["userId"] = user_id
    
    # Construct the request array format that the server expects
    request = {
        "tag": "ApplyTransform",
        "contents": [session_id, parent_id, transform_request] + camel_case_molecules
    }
    response = send_request(sock, request)
    if isinstance(response, dict) and response.get("tag") == "Success":
        # New format: Response is a dict with transformation (including ID) and molecules (with IDs)
        contents = response["contents"]
        if isinstance(contents, dict):
            # Just return the new format directly - clients will access transformation.transformationId
            # and the molecule IDs directly from the molecules array
            return contents
        # # Old format compatibility (string)
        # elif isinstance(contents, str):
        #     return {"transformation": {"transformationId": contents}, "molecules": []}
        # else:
        #     # For backward compatibility with older API versions
        #     converted_contents = convert_keys(contents)
        #     try:
        #         molecule = MoleculeResponse(**converted_contents)
        #         # Convert to new format
        #         return {
        #             "transformation": {"transformationId": "unknown"},
        #             "molecules": [{"moleculeId": getattr(molecule, "moleculeId", "unknown")}]
        #         }
        #     except Exception as e:
        #         raise RuntimeError(f"Failed to parse response: {e}")
    raise RuntimeError(str(response))

def rollback_transformation(sock: socket, session_id: str, transform_id: str) -> Dict[str, Any]:
    request = {
        "tag": "RollbackTransformation",
        "contents": [session_id, transform_id]
    }
    response = send_request(sock, request)
    if isinstance(response, dict) and response.get("tag") == "Success":
        return convert_keys(response["contents"])
    raise RuntimeError(str(response))

def reapply_transformation(
    sock: socket,
    session_id: str,
    transform_id: str,
    molecule_data: Union[dict, List[dict]],
    agent: str,
    user_message: Optional[str] = None,
    method_details: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    """
    Reapply a transformation with potentially updated parameters.

    Args:
        sock: Socket connection
        session_id: ID of the session
        transform_id: ID of the transformation to reapply
        molecule_data: Single molecule data dict or list of molecule data dicts (each should have a rationale)
        agent: Name of the agent performing the transformation
        user_message: Optional user message
        method_details: Optional method details

    Returns:
        Dictionary containing transformationId and moleculeIds
    """
    # Handle single molecule or list of molecules
    if not isinstance(molecule_data, list):
        molecules = [molecule_data or {
            "structure": {},
            "properties": {},
            "computedProperties": {}
        }]
    else:
        molecules = molecule_data

    # Ensure each molecule has the required structure and verify rationale exists
    for i, mol in enumerate(molecules):
        if not mol:
            molecules[i] = {
                "structure": {},
                "properties": {},
                "computedProperties": {}
            }

        # Verify each molecule has a rationale
        if not molecules[i].get("rationale"):
            raise ValueError(f"Molecule at index {i} is missing a rationale field. Each molecule must have its own rationale.")

    transform_request = {
        "transformationId": transform_id,
        "agent": agent,
        "userMessage": user_message,
        "methodDetails": method_details,
        "rationale": None  # No transformation-level rationale, only per-molecule
    }

    # Construct the request array format that the server expects
    request = {
        "tag": "ReapplyTransformation",
        "contents": [session_id, transform_id, transform_request] + molecules
    }
    response = send_request(sock, request)
    if isinstance(response, dict) and response.get("tag") == "Success":
        # New format: Response is a dict with transformation (including ID) and molecules (with IDs)
        contents = response["contents"]
        if isinstance(contents, dict):
            # Just return the new format directly
            return contents
        # For backward compatibility
        converted_contents = convert_keys(contents)
        try:
            return {
                "transformation": {"transformationId": transform_id},
                "molecules": [{"moleculeId": converted_contents.get("moleculeId", "unknown")}]
            }
        except Exception as e:
            raise RuntimeError(f"Failed to parse response: {e}")
    raise RuntimeError(str(response))

def create_interaction(
    sock: socket,
    session_id: str,
    user_message: str,
    agent: str,
    parent_transformation_id: Optional[str] = None,
    molecule_data: Optional[Union[dict, List[dict]]] = None
) -> Dict[str, Any]:
    """
    Create a user interaction transformation.

    Args:
        sock: Socket connection
        session_id: ID of the session
        user_message: User message
        agent: Name of the agent performing the transformation
        parent_transformation_id: Optional ID of the parent transformation
        molecule_data: Optional molecule data (single dict or list of dicts, each should have a rationale)

    Returns:
        Dictionary containing transformationId and moleculeIds
    """
    # Use default empty molecule if none provided
    if molecule_data is None:
        raise ValueError("molecule_data is required and each molecule must include a rationale field")

    # Handle single molecule or list of molecules
    if not isinstance(molecule_data, list):
        molecules = [molecule_data]
    else:
        molecules = molecule_data

    # Ensure each molecule has the required structure and verify rationale exists
    for i, mol in enumerate(molecules):
        if not mol:
            raise ValueError(f"Molecule at index {i} cannot be None or empty")

        # Verify each molecule has a rationale
        if not molecules[i].get("rationale"):
            raise ValueError(f"Molecule at index {i} is missing a rationale field. Each molecule must have its own rationale.")

    transform_request = {
        "type": "user-interaction",
        "agent": agent,
        "userMessage": user_message,
        "methodDetails": None,
        "rationale": None  # No transformation-level rationale, only per-molecule
    }

    # Use a default parent transformation ID if none is provided
    parent_id = parent_transformation_id

    # Construct the request array format that the server expects
    request = {
        "tag": "ApplyTransform",
        "contents": [session_id, parent_id, transform_request] + molecules
    }
    response = send_request(sock, request)
    if isinstance(response, dict) and response.get("tag") == "Success":
        # New format: Response is a dict with transformation (including ID) and molecules (with IDs)
        contents = response["contents"]
        if isinstance(contents, dict):
            # Just return the new format directly
            return contents
        # For backward compatibility
        converted_contents = convert_keys(contents)
        try:
            return {
                "transformation": {"transformationId": "unknown"},
                "molecules": [{"moleculeId": converted_contents.get("moleculeId", "unknown")}]
            }
        except Exception as e:
            raise RuntimeError(f"Failed to parse response: {e}")
    raise RuntimeError(str(response))

# The update_transformation function has been removed as part of the refactoring.
# All transformation interactions should now be done through apply_transform,
# which includes the ability to specify a parent_transformation_id to create
# chains of transformations.
