# Gurnemanz Documentation

> GURNEMANZ:
> Du siehst, mein Sohn,
> zum Raum wird hier die Zeit.
>
> (You see, son, here time transforms into space.)
>
> -- Parsifal, Act 1

## Table of Contents

1. [Overview](#overview)
2. [Architecture](#architecture)
3. [Core Domain Model](#core-domain-model)
4. [Transformation System](#transformation-system)
5. [Session Management](#session-management)
6. [Socket API](#socket-api)
7. [Multi-Agent Support](#multi-agent-support)
8. [Error Handling](#error-handling)
9. [Configuration and Environment](#configuration-and-environment)

## Overview

Gurnemanz is a lean, in-memory molecular transformation system designed to:

- Validate and manage molecular structures
- Apply transformations to molecules in a type-safe manner
- Track transformation chains and history
- Support multiple agents working on the same molecules
- Provide fast, in-memory operations without external dependencies
- Communicate with external components through a socket API

The system uses a layered architecture with strong separation of concerns, phantom types for compile-time safety, and efficient in-memory storage for all molecular data and transformations.

## Architecture

Gurnemanz follows a clean architecture pattern with three main layers:

### Domain Layer

The domain layer contains the core business logic and domain-specific types:

- `Domain.Molecule.Types`: Core molecule representation
- `Domain.Transform.Types`: Transformation definitions and operations
- `Domain.Session.Types`: Session management
- `Domain.Base.Types`: Basic domain types and validation markers

This layer has no dependencies on infrastructure components.

### Infrastructure Layer

The infrastructure layer provides concrete implementations of interfaces defined in the domain layer:

- `Infrastructure.API`: Socket-based API implementation
- `Infrastructure.Config`: Configuration handling
- `Infrastructure.Persistence`: Logging mechanisms

This layer depends on the domain layer but not on the application layer.

### Application Layer

The application layer glues everything together and handles the app lifecycle:

- `App.Types`: Application monad and environment
- `App.Env`: Environment initialization and shutdown
- `App.Main`: Application entry point

This layer depends on both domain and infrastructure layers.

## Core Domain Model

### Molecule Representation

Molecules are represented by the `ValidatedMolecule a` type, which uses a phantom type parameter to track validation status:

```haskell
data ValidatedMolecule a = ValidatedMolecule
    { moleculeId :: Text                -- Unique identifier for the molecule
    , parentMoleculeId :: Maybe Text    -- The parent molecule this was derived from
    , transformationId :: Maybe Text    -- The transformation that created this molecule
    , rationale :: Maybe Text           -- Captures the reason for the molecule creation/modification
    , validationInfo :: Maybe Text
    , smiles :: Maybe Text
    , inchiKey :: Maybe Text
    , molecularFormula :: Maybe Text
    , molecularWeight :: Maybe Double
    , xlogp :: Maybe Double
    , hbd :: Maybe Int
    , hba :: Maybe Int
    , tpsa :: Maybe Double
    , rb :: Maybe Int
    , qed :: Maybe Double
    , sas :: Maybe Double
    , lipinskiCount :: Maybe Double
    , s3Path :: Maybe Text
    , dockingScore :: Maybe Double
    , favorite :: Maybe Bool
    , plipInteractions :: Maybe Text
    , status :: Maybe Text
    , properties :: Map Text Value
    }
```

The phantom type `a` can be either `Valid` or `Draft`, allowing compile-time safety checks on operations that require validated molecules.

### Type Safety

Type safety is enforced using phantom types that prevent invalid operations:

```haskell
-- Abstract marker types
data Valid
data Draft

-- Only molecules marked as Valid can be used in sensitive operations
applyTransform :: Transformation -> ValidatedMolecule Valid -> m (ValidatedMolecule Valid)
```

This approach ensures that operations can't be performed on invalid molecular structures.

## Transformation System

The transformation system is the core of Gurnemanz, managing how molecules are modified.

### Transformation

A single transformation step is represented by the `Transformation` type:

```haskell
data Transformation = Transformation
    { transformationId :: TransformationId
    , transformationType :: Text
    , userMessage :: Maybe Text
    , agentResponse :: Maybe Text
    , methodDetails :: Maybe MethodDetails
    , inputMoleculeIds :: [MoleculeId]  -- Support for multiple input molecules
    , outputMoleculeIds :: [MoleculeId]
    , agent :: Text                     -- Identifies the agent (e.g., "Medicinal Chemist")
    }
```

Key features:
- Each transformation has a unique ID
- Transformations can be associated with user messages and agent responses
- Method details capture the specific operation performed
- Multi-agent support through agent attribution and parent-child relationships
- Input/output molecule IDs track relationships between molecules and support multiple inputs

### Transformation Chain

Transformations are organized into chains using the `TransformationChain` type, with explicit parent-child relationships:

```haskell
data TransformationChain = TransformationChain
    { _steps :: NonEmpty Transformation
    , _status :: ChainStatus
    }

data ChainStatus
    = RolledBack
    | Active
```

Parent-child relationships between transformations are tracked in the API layer through explicit parent transformation IDs. When a new transformation is applied:

```haskell
ApplyTransform sid parentTransformId transformReq molecules
```

This creates a directed graph where each transformation can reference its parent, enabling:

Key operations on chains:
- Adding transformations to the chain with explicit parent links
- Rolling back to a previous state
- Tracking chain status (active or rolled back)
- Retrieving the complete history of a specific molecule
- Tracing molecule lineage through multiple transformations

### Transformation Operations

Core transformation operations include:

- `createTransformation`: Create a new transformation
- `addTransformation`: Add a transformation to a chain
- `rollbackTo`: Roll back a chain to a previous state
- `createInteraction`: Create a user interaction transformation
- `updateTransformation`: Update a transformation with agent response
- `getChainForMolecule`: Get the chain of transformations leading to a molecule

## Session Management

Sessions manage the state of interactions between users, agents, and molecules.

### Session State

The session state is represented by the `SessionState` type:

```haskell
data SessionState = SessionState
    { sessionId :: SessionId
    , userId :: UserId
    , transformationChain :: TransformationChain
    }
```

Each session has:
- A unique session ID
- An associated user ID
- A transformation chain tracking all operations

### Session Operations

Sessions support the following operations:

- `startSession`: Start a new session with an initial transformation
- `logTransformation`: Log a transformation to the chain
- `getTransformationChain`: Get the transformations for a molecule
- `rollbackTo`: Roll back to a previous state

### Free Monad Pattern

Session operations use a free monad pattern through the `SessionF` type:

```haskell
data SessionF next where
    StartSession :: UserId -> Transformation -> (SessionId -> next) -> SessionF next
    LogTransformation :: Transformation -> (TransformationId -> next) -> SessionF next
    GetTransformationChain :: MoleculeId -> ([Transformation] -> next) -> SessionF next
    RollbackTo :: TransformationId -> (TransformationChain -> next) -> SessionF next
```

This pattern allows for flexible interpretation of session operations in in-memory context.

## Socket API

Gurnemanz communicates with external components through a Unix Domain Socket API. This is the primary way to integrate with the system for multi-agent workflows.

### Socket Server

The socket server is implemented in `Infrastructure.API.Socket`:

```haskell
runSocketServer :: AppEnv -> APIConfig -> IO ()
```

This function starts a socket server that listens for client connections and handles requests on the specified socket path.

### Request/Response Protocol

Communication uses a JSON-based request/response protocol:

```haskell
data Request
    = CreateSession (Maybe TransformRequest) (Maybe [Value]) -- (transformation data, initial molecules)
    | ValidateMolecule SessionId Value
    | ApplyTransform SessionId UUID TransformRequest [Value] -- SessionId, parentTransformationId, transformation, molecules
    | RollbackTransformation SessionId Value
    | GetChain SessionId (Maybe UUID) (Maybe UUID) -- SessionId, optional transformationId, optional moleculeId
    | ...

data Response
    = Success Value
    | Failure Error
```

Requests are tagged to indicate the operation type, and responses include either a success value or an error message with details.

### Integration

The socket API provides a simple JSON-based protocol for external integration, supporting all core molecular transformation and multi-agent operations through Unix domain sockets.

## Multi-Agent Support

Gurnemanz supports multiple agents working on the same molecules through several mechanisms.

### Agent Information

Transformations include agent-specific fields:

```haskell
data Transformation = Transformation
    { ...
    , agent :: Text          -- Identifies the agent performing the transformation
    , ...
    }
```

Each molecule can include its own rationale:

```haskell
data ValidatedMolecule a = ValidatedMolecule
    { ...
    , rationale :: Maybe Text    -- The reasoning for this specific molecule
    , ...
    }
```

These fields identify:
- Which agent performed each transformation
- The reasoning behind each molecule creation or modification

### Asynchronous Updates

Transformations can be updated asynchronously:

```haskell
updateTransformation
    :: (MonadError Error m, MonadLogger m)
    => Transformation
    -> Text  -- ^ Agent response
    -> Maybe MethodDetails
    -> [MoleculeId]  -- ^ Output molecules
    -> m Transformation
```

This allows agents to:
1. Start a transformation (e.g., request a computation)
2. Later update it with results (e.g., when computation completes)

### Collaborative Workflows

The system supports complex multi-agent collaborative workflows through:
- Agent identification in transformations
- Parent-child relationship tracking
- Asynchronous transformation updates
- Session-based collaboration context

## Error Handling

Gurnemanz uses a structured error handling system:

```haskell
data Error
    = ValidationError Text
    | ComputationError Text
    | TransformationError Text
    | DataError Text
```

Each error type corresponds to a specific category of failures, and errors can be caught and handled at appropriate levels.

### Error Propagation

Errors are propagated through the `MonadError` typeclass:

```haskell
rollbackTo
    :: (MonadError Error m, MonadLogger m)
    => TransformationId
    -> TransformationChain
    -> m TransformationChain
```

This allows for consistent error handling across different contexts.

### Error Logging

Errors are logged through the `MonadLogger` typeclass:

```haskell
logError :: Text -> m ()
```

This ensures that errors are properly recorded for debugging and monitoring.

## Configuration and Environment

Gurnemanz uses a YAML-based configuration system:

```haskell
data AppConfig = AppConfig
    { chemConfig :: !ChemConfig
    , apiConfig :: !APIConfig
    , logConfig :: !LogConfig
    , environment :: !Environment
    }
```

### Environment Types

Three environment types are supported:

```haskell
data Environment
    = Development
    | Staging
    | Production
```

Configuration values can vary based on the environment.

### Resource Management

Application resources are managed through the `Resources` type:

```haskell
data Resources = Resources
    { resourceLogHandle :: !LogHandle
    }
```

Resources are initialized at startup and properly shut down on application termination:

```haskell
initializeResources :: AppConfig -> m Resources
shutdownResources :: AppEnv -> m ()
```

This ensures proper cleanup of logging handles and other resources.

---
