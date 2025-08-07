from typing import Any
import re

def camel_to_snake(name: str) -> str:
    """Convert camelCase to snake_case."""
    name = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', name).lower()

def snake_to_camel(name: str) -> str:
    """Convert snake_case to camelCase."""
    components = name.split('_')
    return components[0] + ''.join(x.title() for x in components[1:])

def convert_keys_to_snake(obj: Any) -> Any:
    """Convert all dictionary keys from camelCase to snake_case recursively."""
    if isinstance(obj, dict):
        return {camel_to_snake(key): convert_keys_to_snake(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_keys_to_snake(element) for element in obj]
    return obj

def convert_keys_to_camel(obj: Any) -> Any:
    """Convert all dictionary keys from snake_case to camelCase recursively."""
    if isinstance(obj, dict):
        return {snake_to_camel(key): convert_keys_to_camel(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_keys_to_camel(element) for element in obj]
    return obj

# Alias for backward compatibility
convert_keys = convert_keys_to_snake
