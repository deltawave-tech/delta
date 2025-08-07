from socket import socket
from typing import Optional, Dict, Any, Union, List
from ..models.transform import TransformationChain, TransformationResponse, TransformRequest
from ..utils.json_utils import convert_keys
from ..utils.socket_utils import send_request

def create_session(
    sock: socket, 
    user_id: str, 
    initial_transform: Union[TransformRequest, Dict[str, Any]], 
    initial_molecules: Optional[Union[dict, List[dict]]] = None
) -> Dict[str, Any]:
    """
    Create a new session with the specified user ID, initial transformation, and optional molecules.
    
    Args:
        sock: Socket connection
        user_id: ID of the user creating the session
        initial_transform: Initial transformation to start the session with (TransformRequest or dict)
        initial_molecules: Optional initial molecules (single dict or list of dicts)
        
    Returns:
        Dict with sessionId, transformationId and optional moleculeIds
    """
    # Handle both TransformRequest and dict inputs for backward compatibility
    if isinstance(initial_transform, TransformRequest):
        transform_data = {
            "type": initial_transform.type,
            "agent": initial_transform.agent,
            "userMessage": initial_transform.user_message,
            "methodDetails": initial_transform.method_details,
            "rationale": initial_transform.rationale,
            "userId": user_id  # Include user ID in the transformation data
        }
    else:
        # If a dict is provided, use it directly but ensure it has userId
        transform_data = initial_transform.copy()  # Make a copy to avoid modifying the original
        transform_data["userId"] = user_id  # Add userId to the transformation data
    
    # Build the contents object
    contents = {
        "user_id": user_id,  # Keep this for backward compatibility
        "transformation": transform_data
    }
    
    # Handle initial molecules if provided
    if initial_molecules is not None:
        # Convert single molecule to list
        if not isinstance(initial_molecules, list):
            molecules = [initial_molecules]
        else:
            molecules = initial_molecules
        
        # Add molecules to contents
        contents["molecules"] = molecules
    
    # Create the request
    request = {
        "tag": "CreateSession",
        "contents": contents
    }
    
    response = send_request(sock, request)
    if isinstance(response, dict) and response.get("tag") == "Success":
        # New format: Response with sessionId, transformation (with ID) and molecules (with IDs)
        contents = response["contents"]
        if isinstance(contents, dict):
            # Just return the contents directly
            return contents
        # Fallback for unexpected formats
        return {"sessionId": "unknown", "transformation": {}, "molecules": []}
    raise RuntimeError(str(response))

def get_chain(
    sock: socket, 
    session_id: str = None, 
    molecule_id: str = None,
    chain_type: str = "transformation",
    include_molecules: bool = False
) -> TransformationChain:
    """
    Get a chain of transformations or molecules.
    
    Args:
        sock: Socket connection
        session_id: ID of the session (optional)
        molecule_id: ID of the molecule (optional)
        chain_type: Type of chain to retrieve ("transformation" or "molecule")
        include_molecules: Whether to include full molecule data in the response (default: False)
        
    Returns:
        TransformationChain object with full molecule data included in results.molecules if requested
    """
    # Build the request contents
    contents = {}
    
    # Add session ID if provided
    if session_id is not None:
        contents["sessionId"] = session_id
        
    # Add molecule ID if provided
    if molecule_id is not None:
        contents["moleculeId"] = molecule_id
        
    # Add chain type if not the default
    if chain_type and chain_type != "transformation":
        contents["chainType"] = chain_type
    
    # Add includeMolecules parameter if set to True
    if include_molecules:
        contents["includeMolecules"] = True
    
    if not contents:
        raise ValueError("Either session_id or molecule_id must be provided")
    
    # Create the request
    request = {
        "tag": "GetChain",
        "contents": contents
    }
    
    response = send_request(sock, request)
    if isinstance(response, dict) and response.get("tag") == "Success":
        print(f"Chain response type: {type(response['contents'])}")
        
        # Handle string response (often just an ID in test mode)
        if isinstance(response["contents"], str):
            print(f"Received string response: {response['contents']}")
            # Create an empty chain
            return TransformationChain(
                steps=[],
                results={},
                status="ok"
            )
        
        # Handle dictionary response (normal mode)
        elif isinstance(response["contents"], dict):
            try:
                converted_contents = convert_keys(response["contents"])
                
                # Initialize results dictionary to store molecule data if provided
                results = converted_contents.get("results", {})
                
                # Check if molecules are included in the response outside of results
                if "molecules" in converted_contents:
                    # Add molecules to results if not already present
                    if "molecules" not in results:
                        results["molecules"] = converted_contents["molecules"]
                
                # Check if chain field exists and is a list
                if "chain" in converted_contents and isinstance(converted_contents["chain"], list):
                    # Parse each step in the chain
                    steps = []
                    
                    for step in converted_contents["chain"]:
                        if not isinstance(step, dict):
                            print(f"Unexpected step format: {type(step)}, {step}")
                            continue
                        
                        converted_step = convert_keys(step)
                        
                        # Handle molecule chain format (has moleculeId, smiles, transformationId)
                        if "molecule_id" in converted_step:
                            try:
                                # For molecule chains, create a simplified TransformationResponse
                                steps.append(TransformationResponse(
                                    transformation_id=converted_step.get("transformation_id", ""),
                                    transformation_type="molecule",
                                    user_message=f"Molecule: {converted_step.get('smiles', '')}",
                                    agent_response=None,
                                    method_details=None,
                                    input_molecule_ids=[],
                                    output_molecule_ids=[converted_step.get("molecule_id", "")],
                                    agent="system",
                                    rationale=None,
                                    parent_transformation_id=None
                                ))
                            except Exception as e:
                                print(f"Error creating molecule step: {e}")
                        else:
                            # For transformation chains, parse as normal
                            try:
                                steps.append(TransformationResponse(**converted_step))
                            except Exception as e:
                                print(f"Error parsing transformation step: {e}")
                    
                    # Create a TransformationChain with the steps and molecule data in results
                    return TransformationChain(
                        steps=steps,
                        results=results,
                        status=converted_contents.get("status", "ok")
                    )
                else:
                    print(f"Missing or invalid chain field: {converted_contents}")
                    return TransformationChain(steps=[], results={}, status="error")
            except Exception as e:
                print(f"Error parsing chain: {e}")
                return TransformationChain(steps=[], results={}, status="error")
        
        # Handle other response types
        else:
            print(f"Unexpected response format: {type(response['contents'])}")
            print(f"Response: {response['contents']}")
            return TransformationChain(steps=[], results={}, status="error")
    
    raise RuntimeError(str(response))

# The get_transformation_chain function has been removed as part of the refactoring.
# Use get_chain with molecule_id parameter instead to retrieve transformation chains
# related to a specific molecule.

def get_user_transformation_chains(
    sock: socket, 
    user_id: str, 
    include_molecules: bool = True
) -> List[Dict[str, Any]]:
    """
    Retrieve all transformation chains for a given user.
    
    This function finds all "leaf" transformations (those with no children)
    for the specified user, and then returns the complete transformation chain
    for each one, optionally including all associated molecule data.
    
    Args:
        sock: Socket connection
        user_id: ID of the user
        include_molecules: Whether to include full molecule data (default: True)
        
    Returns:
        List of transformation chains, each including the full chain of transformations
        and optionally the associated molecule data
    """
    # Build the request
    request = {
        "tag": "GetUserTransformationChains",
        "contents": {
            "userId": user_id,
            "includeMolecules": include_molecules
        }
    }
    
    # Send the request
    response = send_request(sock, request)
    
    # Process the response
    if isinstance(response, dict) and response.get("tag") == "Success":
        return response["contents"]
    
    # Handle error cases
    raise RuntimeError(str(response))

def log_transformation(sock: socket, transformation: TransformRequest) -> str:
    """
    Log a transformation to the system without applying it.
    
    Args:
        sock: Socket connection
        transformation: Transformation to log
        
    Returns:
        Transformation ID
    """
    transform_data = {
        "type": transformation.type,
        "agent": transformation.agent,
        "userMessage": transformation.user_message,
        "methodDetails": transformation.method_details,
        "rationale": transformation.rationale
    }
    
    request = {
        "tag": "LogTransformation",
        "contents": transform_data
    }
    response = send_request(sock, request)
    if isinstance(response, dict) and response.get("tag") == "Success":
        return response["contents"]
    raise RuntimeError(str(response))
