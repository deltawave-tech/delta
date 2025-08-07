from dataclasses import dataclass
from typing import Optional, Dict, Any, List

@dataclass
class MethodDetails:
    """Details about a transformation method"""
    method_id: str
    parameters: Optional[Dict[str, Any]] = None

@dataclass
class TransformRequest:
    """Request for a molecule transformation"""
    type: str
    agent: str
    user_message: Optional[str] = None
    method_details: Optional[Dict[str, Any]] = None
    rationale: Optional[str] = None
    protein_sequence: Optional[str] = None
    protein_path: Optional[str] = None
    iteration: Optional[int] = None

@dataclass
class TransformationResponse:
    """Response from a transformation operation"""
    transformation_id: str
    transformation_type: str
    user_message: Optional[str]
    agent_response: Optional[str]
    method_details: Optional[Dict[str, Any]]
    input_molecule_ids: List[str]  # Updated: Now a list of molecule IDs
    output_molecule_ids: List[str]
    agent: str                    # Identifies the agent (e.g., "Medicinal Chemist", "AI_expert")
    rationale: Optional[str] = None # Captures the reason for the change
    parent_transformation_id: Optional[str] = None  # New: ID of the parent transformation
    protein_sequence: Optional[str] = None # NEW: Stores a protein sequence associated with the transformation
    protein_path: Optional[str] = None # NEW: Stores a file path to the protein structure file
    iteration: Optional[int] = None # NEW: Stores the iteration number for this transformation

@dataclass
class TransformationChain:
    """Chain of transformations applied to molecules"""
    steps: List[TransformationResponse]
    results: Dict[str, Any]
    status: str

def method_details_to_json(details: MethodDetails) -> Dict[str, Any]:
    return {
        "methodId": details.method_id,
        "parameters": details.parameters or {}
    }

def method_details_from_json(data: Dict[str, Any]) -> MethodDetails:
    return MethodDetails(
        method_id=data.get("methodId", ""),
        parameters=data.get("parameters", {})
    )

def transform_request_to_json(request: TransformRequest) -> Dict[str, Any]:
    result = {
        "type": request.type,
        "agent": request.agent,
        "userMessage": request.user_message,
        "methodDetails": request.method_details,
        "rationale": request.rationale
    }
    
    if request.protein_sequence is not None:
        result["proteinSequence"] = request.protein_sequence
    
    if request.protein_path is not None:
        result["proteinPath"] = request.protein_path
    
    if request.iteration is not None:
        result["iteration"] = request.iteration
    
    return result

def transform_request_from_json(data: Dict[str, Any]) -> TransformRequest:
    return TransformRequest(
        type=data.get("type", ""),
        agent=data.get("agent", ""),
        user_message=data.get("userMessage"),
        method_details=data.get("methodDetails"),
        rationale=data.get("rationale"),
        protein_sequence=data.get("proteinSequence"),
        protein_path=data.get("proteinPath"),
        iteration=data.get("iteration")
    )

def transformation_response_to_json(response: TransformationResponse) -> Dict[str, Any]:
    result = {
        "transformationId": response.transformation_id,
        "transformationType": response.transformation_type,
        "userMessage": response.user_message,
        "agentResponse": response.agent_response,
        "methodDetails": response.method_details,
        "inputMoleculeIds": response.input_molecule_ids,
        "outputMoleculeIds": response.output_molecule_ids,
        "agent": response.agent,
        "rationale": response.rationale
    }
    
    if response.parent_transformation_id is not None:
        result["parentTransformationId"] = response.parent_transformation_id
    
    if response.protein_sequence is not None:
        result["proteinSequence"] = response.protein_sequence
    
    if response.protein_path is not None:
        result["proteinPath"] = response.protein_path
    
    if response.iteration is not None:
        result["iteration"] = response.iteration
    
    return result

def transformation_response_from_json(data: Dict[str, Any]) -> TransformationResponse:
    # Handle backward compatibility with old input_molecule_id
    if "inputMoleculeId" in data and "inputMoleculeIds" not in data:
        input_molecule_id = data.get("inputMoleculeId")
        input_molecule_ids = [input_molecule_id] if input_molecule_id else []
    else:
        input_molecule_ids = data.get("inputMoleculeIds", [])
    
    return TransformationResponse(
        transformation_id=data.get("transformationId", ""),
        transformation_type=data.get("transformationType", ""),
        user_message=data.get("userMessage"),
        agent_response=data.get("agentResponse"),
        method_details=data.get("methodDetails"),
        input_molecule_ids=input_molecule_ids,
        output_molecule_ids=data.get("outputMoleculeIds", []),
        agent=data.get("agent", ""),
        rationale=data.get("rationale"),
        parent_transformation_id=data.get("parentTransformationId"),
        protein_sequence=data.get("proteinSequence"),
        protein_path=data.get("proteinPath"),
        iteration=data.get("iteration")
    )

def transformation_chain_to_json(chain: TransformationChain) -> Dict[str, Any]:
    return {
        "steps": [transformation_response_to_json(step) for step in chain.steps],
        "results": chain.results,
        "status": chain.status
    }

def transformation_chain_from_json(data: Dict[str, Any]) -> TransformationChain:
    return TransformationChain(
        steps=[transformation_response_from_json(step) for step in data.get("steps", [])],
        results=data.get("results", {}),
        status=data.get("status", "")
    )
