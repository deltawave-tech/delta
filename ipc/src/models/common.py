from typing import Dict, Any, TypeVar, Optional, Union
from datetime import datetime
import uuid
from uuid import UUID

# Type variable for generic return types
T = TypeVar('T')

# Custom types for clarity
JsonDict = Dict[str, Any]
ErrorResponse = Dict[str, Union[str, Dict[str, str]]]

def ensure_uuid(value: Optional[str]) -> Optional[UUID]:
    """
    Convert a string to UUID if not None

    Args:
        value: String representation of UUID or None

    Returns:
        UUID object if value is valid UUID string, None if value is None

    Raises:
        ValueError: If string is not a valid UUID
    """
    if value is None:
        return None
    return uuid.UUID(value) if isinstance(value, str) else value

def ensure_datetime(value: Optional[str]) -> Optional[datetime]:
    """
    Convert an ISO timestamp string to datetime if not None

    Args:
        value: ISO format timestamp string or None

    Returns:
        datetime object if value is valid timestamp, None if value is None

    Raises:
        ValueError: If string is not a valid ISO format
    """
    if value is None:
        return None
    return datetime.fromisoformat(value) if isinstance(value, str) else value

def create_empty_molecule() -> JsonDict:
    """
    Create an empty molecule structure

    Returns:
        Dictionary with empty structure, properties, and computedProperties
    """
    return {
        "structure": {},
        "properties": {},
        "computedProperties": {}
    }

def format_response_error(response: ErrorResponse) -> str:
    """
    Format error response for exception messages

    Args:
        response: Error response dictionary from server

    Returns:
        Formatted error message string

    Example:
        >>> format_response_error({"tag": "Failure", "contents": {"type": "ValidationError", "message": "Invalid input"}})
        'Error: ValidationError - Invalid input'
    """
    if isinstance(response, dict):
        if response.get("tag") == "Failure":
            contents = response.get("contents", {})
            return f"Error: {contents.get('type', 'unknown')} - {contents.get('message', 'No message')}"
    return str(response)
