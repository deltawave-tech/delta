"""
GURNEMANZ utilities for IPC integration.

This module provides functions for interacting with the GURNEMANZ transformation
chain logging system, which helps reduce verbosity in shared agent history by
replacing lengthy tool outputs with clean JSON references.
"""

import json
import logging
from typing import Dict, List, Any, Tuple, Optional

# Import GURNEMANZ IPC client libraries
from ipc.src.client.session_ops import create_session, get_chain
from ipc.src.client.transform_ops import apply_transform
from ipc.src.utils.socket_utils import connect_socket, send_request

# Configure logging
logging.basicConfig(level=logging.INFO,
                   format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def format_gurnemanz_for_llm(gurnemanz_data: Dict[str, Any]) -> str:
    """
    Convert GURNEMANZ JSON data to a structured text format that's easier for LLMs to read.

    Args:
        gurnemanz_data: Dictionary containing GURNEMANZ data with transformation and molecules

    Returns:
        Formatted string representation of the GURNEMANZ data
    """
    output = []

    if "transformation" in gurnemanz_data:
        trans = gurnemanz_data["transformation"]
        output.append("## Transformation")
        output.append(f"- Type: {trans.get('type', 'unknown')}")
        output.append(f"- Agent: {trans.get('agent', 'unknown')}")
        output.append(f"- Iteration: {trans.get('iteration', 'unknown')}")

        if "transformationId" in trans:
            output.append(f"- ID: {trans['transformationId']}")

        if "userMessage" in trans:
            output.append(f"- Message: {trans['userMessage']}")

        if "rationale" in trans:
            rationale = trans["rationale"].strip()
            if rationale:
                output.append(f"- Rationale: {rationale}")
            else:
                output.append("- Rationale: No rationale provided")

        # Add protein path information if present
        if "proteinPath" in trans:
            output.append(f"- Protein Path: {trans['proteinPath']}")

        # Add protein sequence information if present
        if "proteinSequence" in trans:
            output.append(f"- Protein Sequence: {trans['proteinSequence']}")

        if "methodDetails" in trans and trans["methodDetails"]:
            method = trans["methodDetails"]
            output.append("- Method Details:")
            output.append(f"  * Method ID: {method.get('methodId', 'unknown')}")
            if "parameters" in method and method["parameters"]:
                for param, value in method["parameters"].items():
                    output.append(f"  * Parameter - {param}: {value}")

    if "molecules" in gurnemanz_data and gurnemanz_data["molecules"]:
        logger.info(f"GURNEMANZ: gurnemanz_data {gurnemanz_data}")
        print(f"GURNEMANZ: gurnemanz_data {gurnemanz_data}")
        output.append("\n## Molecules")
        for i, mol in enumerate(gurnemanz_data["molecules"]):
            output.append(f"\n### Molecule {i+1}")
            logger.info(f"GURNEMANZ: Molecule {i+1}: {mol}")
            print(f"GURNEMANZ: Molecule {i+1}: {mol}")
            if "moleculeId" in mol:
                output.append(f"- ID: {mol['moleculeId']}")
            # Get the nested originalData
            original_data = mol.get("originalData", {})
            print('GURNEMANZ: original_data', original_data)
            # Structure and SMILES
            if "structure" in original_data and "smiles" in original_data["structure"]:
                output.append(f"- SMILES: {original_data['structure']['smiles']}")
            # Properties
            if "properties" in original_data and original_data["properties"]:
                props = original_data["properties"]
                output.append("- Properties:")
                if "status" in props:
                    output.append(f"  * Status: {props['status']}")
                if "parentMoleculeId" in props:
                    output.append(f"  * Parent ID: {props['parentMoleculeId']}")
                for prop, value in props.items():
                    if prop not in ["status", "parentMoleculeId"]:
                        output.append(f"  * {prop}: {value}")
            # Computed Properties
            if "computedProperties" in original_data and original_data["computedProperties"]:
                comp_props = original_data["computedProperties"]
                output.append("- Computed Properties:")
                for prop, value in comp_props.items():
                    output.append(f"  * {prop}: {value}")
            # Rationale
            if "rationale" in original_data:
                output.append(f"- Rationale: {original_data['rationale']}")
            # Activity (if it exists in your data)
            if "activity" in original_data:
                output.append(f"- Activity: {original_data['activity']}")
            if "protein_data" in original_data:
                output.append(f"- Protein Data: {original_data['protein_data']}")
            if "friendly_id" in original_data or "friendlyId" in original_data:
                output.append(f"- friendly_id: {original_data['friendlyId'] if 'friendlyId' in original_data else original_data['friendly_id']}")
            if "parent_friendly_id" in original_data or "parentFriendlyId" in original_data:
                output.append(f"- parent_friendly_id: {original_data['parentFriendlyId'] if 'parentFriendlyId' in original_data else original_data['parent_friendly_id']}")
    return "\n".join(output)

def initialize_session(user_id: str) -> Tuple[str, str]:
    """
    Initialize a new GURNEMANZ session.

    Args:
        user_id: Unique identifier for the user or run

    Returns:
        Tuple containing (session_id, initial_transform_id)
    """
    # Use the GURNEMANZ IPC API
    with connect_socket() as sock:
        initial_transform = {
            "type": "session-init",
            "agent": "session-initializer",
            "user_message": "Starting a new molecule design session",
            "rationale": "Initializing lab workflow"
        }
        response = create_session(sock, user_id, initial_transform)

        # Extract data from the response in the new format
        session_id = response["sessionId"]
        initial_transform_id = response["transformation"]["transformationId"]

        logger.info(f"Initialized GURNEMANZ session {session_id} for user {user_id}")
        return session_id, initial_transform_id

def submit_to_gurnemanz(data: Dict[str, Any], session_id: str, parent_transformation_id: Optional[str] = None) -> Dict[str, Any]:
    """
    Submit transformation and molecule data to GURNEMANZ.

    Args:
        data: Dictionary containing transformation and molecules data
        session_id: Session ID to use (from initialize_session)
        parent_transformation_id: ID of the parent transformation

    Returns:
        Dictionary with the GURNEMANZ response including assigned IDs
    """
    # Extract necessary data
    transformation = data.get("transformation", {})
    molecule_data = data.get("molecules", [])

    logger.info(f"Starting GURNEMANZ submission with transformation type: {transformation.get('type')}")
    logger.info(f"Session ID: {session_id}")
    logger.info(f"Parent transformation ID: {parent_transformation_id}")
    # logger.info(f"Number of molecules: {len(molecules)}")

    # Ensure each molecule has a rationale (required by GURNEMANZ)
    # If a molecule doesn't have a rationale, use the transformation rationale
    transformation_rationale = transformation.get("rationale", "No explicit rationale provided")
    for i, molecule in enumerate(molecule_data):
        if not molecule.get("rationale"):
            logger.warning(f"Molecule at index {i} is missing a rationale. Using transformation rationale.")
            molecule_data[i]["rationale"] = transformation_rationale
        else:
            # Sanitize rationale to avoid JSON parsing issues with special characters
            # Just use a simple string method (not HTML escape) to ensure the string is valid JSON
            rationale = molecule_data[i]["rationale"]
            if isinstance(rationale, str):
                # Simple sanitization to prevent JSON parsing errors
                # Replace characters that might cause issues in JSON
                rationale = rationale.replace('\\', '\\\\')  # Escape backslashes first
                rationale = rationale.replace('"', '\\"')    # Escape double quotes
                rationale = rationale.replace("'", "\\'")    # Escape single quotes
                rationale = rationale.replace('\n', ' ')     # Replace newlines with spaces
                rationale = rationale.replace('\r', ' ')     # Replace carriage returns
                rationale = rationale.replace('\t', ' ')     # Replace tabs
                molecule_data[i]["rationale"] = rationale

    # Use the real GURNEMANZ IPC API
    with connect_socket() as sock:
        try:
            # Extract the required fields from the transformation data
            transform_type = transformation["type"]
            agent = transformation["agent"]
            user_message = transformation.get("user_message")
            method_details = transformation.get("method_details")
            rationale = transformation.get("rationale")

            # Safely extract optional fields
            iteration = transformation.get("iteration")
            protein_sequence = transformation.get("protein_sequence")
            protein_path = transformation.get("protein_path")

            # Log what we're about to send
            logger.info(f"Using transform_type: {transform_type}")
            logger.info(f"Using agent: {agent}")
            if iteration is not None:
                logger.info(f"Using iteration: {iteration}")
            if protein_sequence is not None:
                logger.info(f"Using protein_sequence of length: {len(protein_sequence)}")

            # Use the apply_transform function directly with the extracted fields
            # This ensures we're using the exact API format expected by GURNEMANZ
            # Note: Following the exact parameter order from transform_ops.py
            response = apply_transform(
                sock=sock,
                session_id=session_id,
                transform_type=transform_type,
                molecule_data=molecule_data,
                agent=agent,
                rationale=rationale or "",
                parent_transformation_id=parent_transformation_id,
                user_message=user_message,
                method_details=method_details,
                protein_sequence=protein_sequence,
                protein_path=protein_path,
                iteration=iteration
            )

            # Log the response
            logger.debug(f"Response from apply_transform: {response}")

            # The response from apply_transform is already in the format we need
            # Extract transformation ID for logging
            if isinstance(response, dict) and "transformation" in response and "transformationId" in response["transformation"]:
                transform_id = response["transformation"]["transformationId"]
                logger.info(f"Received transformation ID: {transform_id}")

                # Log molecule IDs if present
                if "molecules" in response:
                    logger.info(f"Received {len(response['molecules'])} molecules in response")
                    for i, mol in enumerate(response["molecules"]):
                        logger.debug(f"Molecule {i} assigned ID: {mol.get('moleculeId')}")

                logger.info(f"Successfully submitted transformation to GURNEMANZ")
                return response
            else:
                logger.error(f"Unexpected response format from apply_transform: {response}")
                # Still return the response as-is to allow for debugging
                return response
        except Exception as e:
            logger.error(f"Error in submit_to_gurnemanz: {str(e)}")
            logger.error(f"Request data: {json.dumps(data, indent=2)}")
            raise

def get_transformation_chain_with_molecules(session_id: str, transformation_id: str) -> Dict[str, Any]:
    """
    Retrieve a specific transformation chain with full molecule data.

    Args:
        session_id: Current session identifier
        transformation_id: ID of the specific transformation to retrieve

    Returns:
        Dictionary containing the transformation chain with full molecule data
    """
    # Use the GURNEMANZ IPC API with the new include_molecules parameter
    with connect_socket() as sock:
        # Get the chain with molecules data for this transformation
        logger.info(f"Retrieving transformation chain for ID {transformation_id} with molecules")

        chain = get_chain(
            sock,
            session_id=session_id,
            include_molecules=True
        )

        # Find the specific transformation in the chain
        target_step_index = None
        for i, step in enumerate(chain.steps):
            if step.transformation_id == transformation_id:
                target_step_index = i
                break

        if target_step_index is not None:
            # Extract the steps up to and including the target transformation
            relevant_steps = chain.steps[:target_step_index + 1]

            # Collect all molecule IDs referenced in these steps
            molecule_ids = set()
            for step in relevant_steps:
                molecule_ids.update(step.input_molecule_ids)
                molecule_ids.update(step.output_molecule_ids)

            # Extract the molecule data for all referenced molecules
            molecules = {}
            if "molecules" in chain.results:
                all_molecules = chain.results["molecules"]
                for mol_id in molecule_ids:
                    if mol_id in all_molecules:
                        molecules[mol_id] = all_molecules[mol_id]

            # Create the result with the transformation chain and molecules
            result = {
                "steps": [step.__dict__ for step in relevant_steps],
                "molecules": molecules
            }

            logger.info(f"Found transformation chain with {len(relevant_steps)} steps and {len(molecules)} molecules")
            return result

        logger.warning(f"Transformation {transformation_id} not found in session {session_id}")
        return {"steps": [], "molecules": {}}

class MoleculeIDGenerator:
    """Simple molecule ID generator that tracks generation based on parent IDs"""
    
    def __init__(self):
        # Agent name to code mapping
        self.agent_codes = {
            "Database Agent": "DA",
            "AI Expert": "AI",
            "Medicinal Chemist": "MC",
            "Literature Agent": "LA",
            "Ranking Agent": "RA",
            "Principal Researcher": "PR",
            "Scientific Critic": "SC",
            "Summary Parser": "SP"
        }
        
        # Track molecule counts per agent per iteration
        self.counters = {}
    
    def generate_id(self, agent_name: str, molecule_number: int, iteration: int, parent_id: str = None) -> str:
        """
        Generate a molecule ID. Generation is automatically determined from parent.
        
        Args:
            agent_name: Full agent name (e.g., "Database Agent")
            iteration: Current iteration number
            parent_id: Parent molecule's ID (optional). If provided, generation = parent_generation + 1
        
        Returns:
            Friendly ID string (e.g., "DA:I1:N3:G0")
        """
        # Get agent code
        agent_code = self.agent_codes.get(agent_name, agent_name[:2].upper())
        
        # Create counter key
        # counter_key = f"{agent_code}:I{iteration}"
        
        # # Initialize counter if needed
        # if counter_key not in self.counters:
        #     self.counters[counter_key] = 0
        
        # # Increment molecule number
        # self.counters[counter_key] += 1
        # molecule_number = self.counters[counter_key]
        molecule_number = molecule_number
        # Determine generation from parent
        generation = self._get_generation(parent_id)
        
        return f"{agent_code}:I{iteration}:N{molecule_number}:G{generation}"
    
    def _get_generation(self, parent_id: str = None) -> int:
        """
        Determine generation based on parent ID.
        If no parent, generation = 0
        If parent exists, generation = parent_generation + 1
        """
        if not parent_id:
            return 0
        
        # Parse parent ID to extract generation
        try:
            parts = parent_id.split(':')
            if len(parts) == 4 and parts[3].startswith('G'):
                parent_generation = int(parts[3][1:])
                return parent_generation + 1
        except (ValueError, IndexError):
            # If parsing fails, assume it's a new molecule
            return 0
        
        return 0
    
    def parse_id(self, molecule_id: str) -> dict:
        """Parse a molecule ID into its components"""
        try:
            parts = molecule_id.split(':')
            if len(parts) != 4:
                return None
            
            return {
                "agent_code": parts[0],
                "iteration": int(parts[1][1:]),
                "molecule_number": int(parts[2][1:]),
                "generation": int(parts[3][1:])
            }
        except (ValueError, IndexError):
            return None
    
    def reset_counters(self):
        """Reset all molecule counters"""
        self.counters = {}
    
    def get_lineage_info(self, molecule_id: str, parent_id: str = None) -> dict:
        """Get information about a molecule's lineage"""
        parsed = self.parse_id(molecule_id)
        if not parsed:
            return None
        
        info = {
            "id": molecule_id,
            "generation": parsed["generation"],
            "is_original": parsed["generation"] == 0,
            "parent_id": parent_id
        }
        
        # Add parent info if available
        if parent_id:
            parent_parsed = self.parse_id(parent_id)
            if parent_parsed:
                info["parent_agent"] = parent_parsed["agent_code"]
                info["parent_generation"] = parent_parsed["generation"]
        
        return info