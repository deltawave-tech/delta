# Gurnemanz IPC Client

## Overview

- **Transport**: Unix domain socket (`/tmp/gurnemanz.sock`)
- **Core Operations**: Session management, transformation logging, molecule tracking
- **Use Case**: Multi-agent molecular design workflows with full provenance

## Usage Example

```python
from ipc.src.utils.socket_utils import connect_socket
from ipc.src.client.session_ops import create_session
from ipc.src.client.transform_ops import apply_transform
from ipc.src.models.transform import TransformRequest

# Connect and create session
with connect_socket() as sock:
    # 1. Create session
    initial_transform = TransformRequest(
        type="session-init",
        agent="user",
        user_message="Start molecular design",
        rationale="Initialize design session"
    )

    session_response = create_session(sock, "user-123", initial_transform)
    session_id = session_response["sessionId"]
    transform_id = session_response["transformationId"]

    # 2. Apply transformation with molecule
    result = apply_transform(
        sock,
        session_id=session_id,
        transform_type="molecule-generation",
        molecule_data={
            "structure": {"smiles": "CCO"},
            "properties": {"name": "Ethanol"},
            "computedProperties": {"logP": 0.25},
            "rationale": "Generated ethanol as starting molecule"
        },
        agent="gen-agent",
        parent_transformation_id=transform_id,
        user_message="Generated initial molecule"
    )
```

## JSON Template Specs

### TransformRequest
```json
{
    "type": "string",              // Required: transformation type
    "agent": "string",             // Required: agent name
    "user_message": "string",      // Required: human-readable message
    "rationale": "string"          // Required: reasoning for transformation
}
```

### Molecule Data
```json
{
    "structure": {
        "smiles": "string"         // Required: SMILES notation
    },
    "properties": {
        "name": "string"           // Required: molecule name
    },
    "computedProperties": {        // Optional: calculated properties
        "logP": "number",
        "molecular_weight": "number"
    },
    "rationale": "string"          // Required: why this molecule was created
}
```

### Session Response
```json
{
    "sessionId": "uuid",           // Session identifier
    "transformationId": "uuid",    // Initial transformation ID
    "originalRequest": "object"    // Echo of request with IDs
}
```

### Transform Response
```json
{
    "transformationId": "uuid",    // New transformation ID
    "moleculeIds": ["uuid"],       // Generated molecule IDs
    "success": "boolean"
}
```
