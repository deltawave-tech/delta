# DELTA: AI-Powered Scientific Research Platform

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

DELTA is a multi-agent platform designed for AI-driven scientific research, with a strong focus on computational drug discovery. It empowers teams of specialized AI agents to collaborate on complex scientific workflows, from literature review to molecular docking and analysis.

## Features

-   **Multi-Agent Collaboration**: Orchestrate teams of AI agents with different roles and expertise.
-   **Extensible Toolset**: A rich set of tools for scientific research, which can be easily extended.
-   **Computational Chemistry & Biology Focus**: Includes tools for:
    -   Literature search (PubMed, UniProt)
    -   Database queries (PDB, ChEMBL)
    -   Molecular docking and analysis (Vina, DiffDock, PLIP)
    -   De novo molecule generation
    -   Drug-likeness prediction (QED, SA Score)
-   **Configurable LLM Providers**: Supports various LLM providers like Anthropic (Claude) and OpenAI (GPT).
-   **Reproducible Workflows**: Saves conversation history and results for each run.

## Repository Structure

```
multi_agent_paper/
├── src/
│   ├── agent.py              # Core Agent class definition
│   ├── demo/                 # Example pipelines and scripts
│   │   └── multi_agent_pipeline.py       # Main script to run a drug discovery pipeline
│   └── tools/                # Directory for all agent tools
│       ├── tool_definitions.py # Schema and definitions for all tools
│       ├── api_based_tools/  # Tools that wrap around external APIs
│       └── ...               # Individual tool implementation files
├── autodock/                 # AutoDock Vina docking service
├── gurnemanz/                # Molecular transformation tracking service
├── ipc/                      # Inter-process communication client for Gurnemanz
├── plip/                     # Protein-ligand interaction profiler
├── server/                   # API server endpoints
├── tests/                    # Tests for the codebase
└── README.md                 # This file
```

## Installation

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd multi_agent_paper
    ```

2.  **Install dependencies using uv:**
    ```bash
    uv sync
    ```

    To activate the virtual environment:
    ```bash
    source .venv/bin/activate
    ```

    If you prefer using pip and virtual environments manually:
    ```bash
    python3 -m venv venv
    source venv/bin/activate
    pip install --upgrade pip
    pip install -e ".[dev]"
    ```

## AutoDock Setup

To use the molecular docking capabilities:

```bash
cd autodock
docker compose up --build
```

The AutoDock service will be available for molecular docking operations. For detailed configuration and usage, see the [AutoDock README](autodock/README.md).

## Configuration

The platform requires API keys for the LLM providers and potentially other services.

1.  Create a `.env` file in the root of the project.
2.  Add your API keys to the `.env` file:

    ```env
    # For OpenAI
    OPENAI_API_KEY="your-openai-key"

    # For Anthropic/Claude
    ANTHROPIC_API_KEY="your-anthropic-key"
    ```
    The application will load these environment variables automatically.

## How to Run the Pipeline

### Prerequisites

Before running the pipeline, you need to start the required services:

1. **Start the Gurnemanz molecular transformation service:**
   ```bash
   cd gurnemanz
   cabal build
   cabal run gurnemanz
   ```

2. **Start the PLIP API service:**
   ```bash
   uv run python plip/plip/plip_api.py
   ```

### Running the Main Pipeline

The main entry point for running a research workflow is `src/demo/multi_agent_pipeline.py`.

You can run it from the command line and specify the LLM provider and number of iterations:

```bash
uv run python src/demo/multi_agent_pipeline.py --llm_provider sonnet-4 --iteration_num 3
```

-   `--llm_provider`: Choose from `sonnet-4`, `sonnet-3.7`, `o3`, `gpt4.1`, `gemini`. Default is `sonnet-4`.
-   `--iteration_num`: Number of research cycles to run. Default is `3`.

Each run will create a new directory in `runs/` containing the conversation logs and any generated files.

## How to Add New Tools

The platform is designed to be extensible. You can add new tools for the agents to use by following these steps:

1.  **Implement the tool logic**: Create a new Python file in the `src/tools/` directory (or an appropriate subdirectory). This file should contain the function that performs the tool's action.

2.  **Define the tool schema**: Open `src/tools/tool_definitions.py`.
    -   Add a new entry to the `ToolNames` enum for your new tool.
    -   Add a new dictionary entry to `ANTHROPIC_TOOLS`. This entry defines the tool's `name`, `description`, and `input_schema` (using JSON Schema). This information is crucial for the LLM to understand how to use the tool.

    Example for a new tool `my_new_tool`:
    ```python
    // ... in ToolNames enum
    MY_NEW_TOOL = "my_new_tool"

    // ... in ANTHROPIC_TOOLS dict
    [ToolNames.MY_NEW_TOOL.value]: {
        "name": "my_new_tool",
        "description": "A brief but clear description of what this tool does.",
        "input_schema": {
            "type": "object",
            "properties": {
                "parameter1": {
                    "type": "string",
                    "description": "Description of the first parameter."
                },
                "parameter2": {
                    "type": "boolean",
                    "description": "Description of the second parameter."
                }
            },
            "required": ["parameter1"]
        }
    }
    ```

3.  **Register the tool**: In `src/tools/tool_definitions.py`, use the `@tool_registry.register` decorator to link the tool name to your implementation function.

    ```python
    from src.tools.my_new_tool_file import my_new_tool_function

    // ... at the end of tool_definitions.py
    @tool_registry.register("my_new_tool")
    def my_new_tool(arguments: dict, output_dir_id: str) -> dict:
        # Call your tool's logic here
        param1 = arguments.get("parameter1")
        param2 = arguments.get("parameter2")
        result = my_new_tool_function(param1, param2)
        return {"result": result}
    ```

4.  **Add the tool to an agent**: In `src/demo/multi_agent_pipeline.py` (or your custom pipeline script), add the new tool to the `tools` list of the agent that should have this capability.

    ```python
    // ... in multi_agent_pipeline.py
    database_agent = Agent(
        # ...
        tools=[
            // ... existing tools
            run_map[run][ToolNames.MY_NEW_TOOL.value],
        ]
    )
    ```

## Available Tools

The platform comes with a variety of pre-built tools for scientific research:

### Database and Literature Search
-   **search_uniprot**: Search UniProt for protein information.
-   **search_pubmed_agent**: Search PubMed for scientific articles.
-   **get_pdb_file**: Download PDB files for proteins.
-   **search_chembl_activity**: Search the ChEMBL database for bioactivity data.

### Molecular Modeling and Docking
-   **vina_docker**: Perform molecular docking using AutoDock Vina.
-   **diffdock**: Perform molecular docking using DiffDock.
-   **get_plip_report**: Analyze protein-ligand interactions using PLIP.

### Molecule Generation and Analysis
-   **molecule_generation**: Generate novel molecules.
-   **vina_mol_gen**: A comprehensive tool to generate, filter, and dock molecules against a target.
-   **QED**: Calculate the Quantitative Estimate of Drug-likeness.
-   **SA**: Calculate the Synthetic Accessibility score.
-   **molecular_similarity**: Calculate the similarity between molecules.
-   **logP** / **molecular_weight**: Calculate molecular properties.

For more details on each tool's parameters, please refer to the `input_schema` in `src/tools/tool_definitions.py`.

## References

-   **Prot2Mol**: Ünlü, A., Çevrim, E., & Doğan, T. (2024). *Prot2Mol: Target based molecule generation using protein embeddings and SELFIES molecule representation*. GitHub Repository. [https://github.com/HUBioDataLab/Prot2Mol](https://github.com/HUBioDataLab/Prot2Mol)

-   **PaperQA**: Lála, J., et al. (2023). PaperQA: Retrieval-Augmented Generative Agent for Scientific Research. *arXiv preprint arXiv:2312.07559*. [https://arxiv.org/abs/2312.07559](https://arxiv.org/abs/2312.07559)

-   **Simple and Effective Masked Diffusion Language Models**: Sahoo, S. S., et al. (2024). Simple and Effective Masked Diffusion Language Models. *arXiv:2406.07524*. [https://arxiv.org/abs/2406.07524](https://arxiv.org/abs/2406.07524)

-   **AutoDock-GPU**: The AutoDock-GPU developers. *AutoDock-GPU*. GitHub Repository. [https://github.com/ccsb-scripps/AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU)

-   **DiffDock**: Corso, G., et al. (2024). Deep Confident Steps to New Pockets: Strategies for Docking Generalization. *International Conference on Learning Representations (ICLR)*.

-   **PLIP**: Adasme, M., et al. (2021). PLIP 2021: expanding the scope of the protein-ligand interaction profiler to DNA and RNA. *Nucleic Acids Research*, gkab294. [https://doi.org/10.1093/nar/gkab294](https://doi.org/10.1093/nar/gkab294)

-   **plipinteractor**: Olgaç, A. *plipinteractor*. GitHub Repository. [https://github.com/aolgac/plipinteractor](https://github.com/aolgac/plipinteractor)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. 