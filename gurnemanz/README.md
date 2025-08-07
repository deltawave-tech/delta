# GURNEMANZ

<img src="gurne_mane.png" width="50%">

> Du siehst, mein Sohn, zum Raum wird hier die Zeit.
>
> *(You see, son, here time transforms into space.)*
>
> -- Gralsritter Gurnemanz: Parsifal, Act 1

A lean, in-memory molecular transformation system with type-safe operations and multi-agent collaboration support.

## Description

Gurnemanz is an in-memory molecular provenance service that logs the evolutionary history of all molecular candidates in multi-agent design workflows. It maintains two distinct ordering structures:

1. **Transformation Chains**: A totally ordered temporal chain that records the immutable, linear sequence of all transformations as they occur within a session. Each transformation (agent-driven modification) is attributed to the specific agent that originated it, providing a complete chronological account of every action.

2. **Molecule Chains**: A direct lineage forest maintained across all molecules, where each molecule explicitly references its single parent. This enables rapid traversal of molecular ancestry and captures the parent-child relationships as agents modify and derive new molecules.

The system models relationships using a directed hypergraph where transformations act as hyperedges connecting input molecules (domain) to output molecules (codomain), precisely capturing complex events like a single parent molecule yielding multiple derivatives.

## Quick Start

```bash
# Build and run
cabal build
cabal run gurnemanz

# Or with custom config
GURNEMANZ_CONFIG=./config.yaml cabal run gurnemanz
```

## Configuration

Example `config.yaml`:

```yaml
environment: Development

chemConfig:
  maxMoleculeSize: 1000
  maxTransformationSteps: 50
  computationTimeout: 300

apiConfig:
  socketPath: "/tmp/gurnemanz.sock"

logConfig:
  logLevel: Info
  enableConsoleLog: true
```

## Features

- **Type-safe molecule operations** with phantom types
- **In-memory storage** for fast performance
- **Multi-agent collaboration** support
- **Transformation chains** with parent-child relationships
- **Socket API** for external integration
- **Structured logging** and error handling

## Development

```bash
cabal build
cabal run gurnemanz
cabal test
cabal repl
```

## Documentation

See [DOCUMENTATION.md](DOCUMENTATION.md) for detailed architecture, API reference, and implementation details.
