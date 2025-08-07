from src.tools.tool_definitions import tool_registry
from src.tools.gurnemanz_utils import format_gurnemanz_for_llm

@tool_registry.register("gurnemanz_apply")
def gurnemanz_apply(arguments: dict, output_dir_id: str) -> dict:
    """Submit transformation and molecule data to GURNEMANZ, format the response, and store it in session data."""
    from src.tools.gurnemanz_utils import submit_to_gurnemanz
    import json
    import os

    print("GURNEMANZ Apply tool is being used")
    client_name = arguments.get("client_name", "gemini")
    # Load session information
    session_file = os.path.join(tool_registry.output_dir, output_dir_id, "gurnemanz_session.json")
    print(f"Session file: {session_file}")
    try:
        with open(session_file, 'r') as f:
            session_data = json.load(f)
        session_id = session_data.get("session_id")
        parent_transformation_id = session_data.get("last_transform_id")
        if not session_id:
            raise ValueError("No session_id found in session file")
        if not parent_transformation_id:
            raise ValueError("No last_transform_id found in session file")
        print(f"Using session_id: {session_id} from session file")
        print(f"Using parent_transformation_id: {parent_transformation_id} from session file")
    except Exception as e:
        raise ValueError(f"Error loading session file: {str(e)}")

    # Prepare data
    data = arguments

    # Ensure molecules is at least an empty list if not present
    if "molecules" not in data:
        data["molecules"] = []

    for i, mol in enumerate(data.get("molecules", [])):
        if not mol.get("rationale"):
            trans_rationale = data["transformation"].get("rationale", "No explicit rationale")
            print(f"Warning: Molecule at index {i} is missing a rationale. Using transformation rationale.")
            data["molecules"][i]["rationale"] = trans_rationale

    try:
        data_copy = data.copy()
        print(f"Parser agent JSON: {json.dumps(data_copy)}")
        print(f"Using session_id: {session_id}")
        print(f"Using parent_transformation_id: {parent_transformation_id}")

        # Submit to GURNEMANZ and get the raw response
        response = submit_to_gurnemanz(data_copy, session_id, parent_transformation_id)

        # Update the session file with the new transformation ID (existing logic)
        if "transformation" in response and "transformationId" in response["transformation"]:
            session_data["last_transform_id"] = response["transformation"]["transformationId"]
            print(f"Updated last_transform_id to {session_data['last_transform_id']}")

        # Make a copy of the response and format it
        formatted_content = format_gurnemanz_for_llm(response)
        if client_name == "gemini":
            formatted_data = {
                "role": "model",
                "parts": [{"text": formatted_content}],
            }
        else:
            formatted_data = {
                "role": "assistant",
                "content": formatted_content,
            }

        # Store the formatted data and raw JSON separately in session_data
        session_data["latest_response"] = formatted_data
        session_data["latest_raw_json"] = response  # Store raw JSON separately

        # Save the updated session data to the file
        with open(session_file, 'w') as f:
            json.dump(session_data, f)
        print("Stored formatted Gurnemanz response in session data")

        # Return the original, unformatted response to the Parser Agent
        return response

    except Exception as e:
        error_msg = str(e)
        print(f"Error in GURNEMANZ Apply tool: {error_msg}")
        return {"error": error_msg}
