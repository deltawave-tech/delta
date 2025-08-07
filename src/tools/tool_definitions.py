import os
import sys
import copy
from enum import Enum
from pathlib import Path
from rdkit.Chem import RDConfig
from typing import Dict, Any, Callable
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))

class ToolNames(Enum):
    DIFFDOCK = "diffdock"
    UNIPROT_SEARCH = "search_uniprot"
    MOLECULE_GENERATION = "molecule_generation"
    PDB_FILE = "pdb_file"
    AF_PDB_FILE = "af_pdb_file"
    PUBMED_SEARCH_AGENT = "search_pubmed_agent"
    CHEMBL_SEARCH_ACTIVITY = "search_chembl_activity",
    RDIFFUSION = "rdiffusion"
    PLIP_RUN = "get_plip_report"
    AGENT_CALL = "call_agent"
    MOCK_TOOL = "mock_tool"
    VALIDITY = "validity"
    QED = "QED"
    SA = "SA"
    READ_FILE = "read_file"
    VISUALIZE_MOLECULES = "visualize_molecules"
    TEST_MOLECULE_VIZ = "test_molecule_visualization"
    VINA_DOCKER = "vina_docker"
    MOLECULAR_SIMILARITY = "molecular_similarity"
    MOLECULAR_WEIGHT = "molecular_weight"
    LOG_P = "logP"
    MOL_GEN = "mol_gen"
    VINA_MOL_GEN = "vina_mol_gen" 
    GURNEMANZ_APPLY = "gurnemanz_apply"
    VINA_REPORT = "vina_report"
    VINA_GEN_MOL = "vina_genmol"

ANTHROPIC_TOOLS: Dict[str, dict] = {
    ToolNames.VINA_MOL_GEN.value: {
        "name": "vina_mol_gen",
        "description": """Generate drug-like molecules based on a protein target, filter them by molecular properties, and perform docking.
This comprehensive tool has TWO PRIMARY MODES:
1. DE NOVO GENERATION: Generate new molecules based on protein sequence. Set generate_de_novo to True.


This comprehensive tool:
- Generates potential binding molecules using the provided protein sequence OR evaluates existing molecules
- Filters molecules based on drug-likeness properties (QED, synthetic accessibility, logP, molecular weight)
- Performs molecular docking using Autodock Vina to evaluate binding potential
- Analyzes protein-ligand interactions using PLIP
- Returns detailed information including:
  - SMILES strings of the molecules
  - Drug-likeness scores (QED, SA, logP, MW)
  - Docking scores (binding energy in kcal/mol, lower/more negative values indicate stronger binding)
  - Protein-ligand interaction analysis
  - Paths to generated complex structures
The entire process is automated with smart filtering to prioritize promising candidates.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "protein_path": {
                    "type": "string",
                    "description": "Path to protein PDB file for docking"
                },
                "protein_sequence": {
                    "type": "string",
                    "description": "Amino acid sequence of the protein for molecule generation"
                },
                "num_molecules": {
                    "type": "integer",
                    "description": "Number of molecules to initially generate (only used when generate_de_novo is True)",
                    "default": 50
                },
                "generate_de_novo": {
                    "type": "boolean",
                    "description": "If set to True, generates new molecules based on the protein sequence. Required if molecules are not provided."
                },
                "thresholds": {
                    "type": "object",
                    "description": "Filtering criteria for molecule properties, only applied to de novo generated molecules",
                    "properties": {
                        "qed": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum QED value (0-1, higher is more drug-like)"},
                                "max": {"type": "number", "description": "Maximum QED value"}
                            }
                        },
                        "sa": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum synthetic accessibility score"},
                                "max": {"type": "number", "description": "Maximum synthetic accessibility score (lower is easier to synthesize)"}
                            }
                        },
                        "logp": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum logP value"},
                                "max": {"type": "number", "description": "Maximum logP value"}
                            }
                        },
                        "molecular_weight": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum molecular weight (Da)"},
                                "max": {"type": "number", "description": "Maximum molecular weight (Da)"}
                            }
                        }
                    }
                }
            },
            "required": ["protein_path", "protein_sequence", "generate_de_novo", "thresholds", "num_molecules"]
        }
    },
    ToolNames.VINA_REPORT.value: {
        "name": "vina_report",
        "description": """Generate drug-like molecules based on a protein target, filter them by molecular properties, and perform docking.
1. EVALUATE EXISTING MOLECULES:
This comprehensive tool:
- Performs molecular docking using Autodock Vina to evaluate binding potential
- Analyzes protein-ligand interactions using PLIP
- Returns detailed information including:
  - SMILES strings of the molecules
  - Drug-likeness scores (QED, SA, logP, MW)
  - Docking scores (binding energy in kcal/mol, lower/more negative values indicate stronger binding)
  - Protein-ligand interaction analysis
  - Paths to generated complex structures
The entire process is automated with smart filtering to prioritize promising candidates.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "protein_path": {
                    "type": "string",
                    "description": "Path to protein PDB file for docking"
                },
                "molecules": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings to evaluate"
                }
            },
            "required": ["protein_path", "molecules"]
        }
    },
    ToolNames.GURNEMANZ_APPLY.value: {
        "name": "gurnemanz_apply",
        "description": "Submit transformation and molecule data to GURNEMANZ and retrieve the response with IDs. The Haskell service returns the full data structure with embedded IDs, providing better traceability between molecules and their IDs.",
        "input_schema": {
            "type": "object",
            "properties": {
                "session_id": {
                    "type": "string",
                    "description": "Session UUID from creation (required for tracking transformations across sessions)"
                },
                "parent_transformation_id": {
                    "type": "string",
                    "description": "ID of the parent transformation (required for tracking transformation chains)"
                },
                "client_name": {
                    "type": "string",
                    "description": "Use claude for client name",
                    "default": "claude"
                },
                "transformation": {
                    "type": "object",
                    "description": "Details of the transformation being applied",
                    "properties": {
                        "type": {
                            "type": "string",
                            "description": "Required: Type of transformation (e.g., 'molecule-generation', 'property-prediction')"
                        },
                        "agent": {
                            "type": "string",
                            "description": "Required: Name of the agent performing the transformation"
                        },
                        "user_message": {
                            "type": "string",
                            "description": "Optional: Human-readable message about the transformation"
                        },
                        "method_details": {
                            "type": "object",
                            "description": "Optional: Details about the method used in the transformation",
                            "properties": {
                                "method_id": {
                                    "type": "string",
                                    "description": "Required if methodDetails is present: Identifier of the method used"
                                },
                            },
                            "required": ["method_id"]
                        },
                        "rationale": {
                            "type": "string",
                            "description": "Optional: General rationale for the transformation"
                        },
                        "protein_sequence": {
                            "type": "string",
                            "description": "Optional: Protein sequence associated with the transformation"
                        },
                        "protein_path": {
                            "type": "string",
                            "description": "Optional: Path to the protein pdb file associated with the transformation"
                        },
                        "iteration": {
                            "type": "integer",
                            "description": "Required: Iteration number for this transformation"
                        }
                    },
                    "required": ["type", "agent", "iteration"]
                },
                "molecules": {
                    "type": "array",
                    "description": "Array of molecule objects generated or modified in this transformation",
                    "items": {
                        "type": "object",
                        "properties": {
                            "structure": {
                                "type": "object",
                                "description": "Required: Molecular structure information",
                                "properties": {
                                    "smiles": {
                                        "type": "string",
                                        "description": "Required: SMILES string representation"
                                    },
                                },
                                "required": ["smiles"]
                            },
                            "properties": {
                                "type": "object",
                                "description": "Optional: Molecule properties",
                                "properties": {
                                    "parent_molecule_id": {
                                        "type": "string",
                                        "description": "Optional: ID of parent molecule"
                                    },
                                    "molecular_weight": {
                                        "type": "number",
                                        "description": "Optional: Molecular weight"
                                    },
                                    "status": {
                                        "type": "string",
                                        "description": "Optional: Molecule status (e.g., 'modified', 'de novo')"
                                    }
                                }
                            },
                            "computed_properties": {
                                "type": "object",
                                "description": "Optional: Computed or predicted properties",
                                "properties": {
                                    "logP": {"type": "number", "description": "Optional: Lipophilicity"},
                                    "qed": {"type": "number", "description": "Optional: Quantitative estimate of drug-likeness"},
                                    "sas": {"type": "number", "description": "Optional: Synthetic accessibility score"},
                                    "lipinski_count": {"type": "integer", "description": "Optional: Number of Lipinski's Rule of Five violations"},
                                    "docking_score": {"type": "number", "description": "Optional: Docking score (if applicable)"},
                                    "plip_interactions": {"type": "string", "description": "Optional: Protein-ligand interaction information"}
                                }
                            },
                            "rationale": {
                                "type": "string",
                                "description": "Required: Reasoning specific to this molecule"
                            },
                            "activity": {
                                "type": "string",
                                "description": "Optional: Activity status of the molecule (e.g., 'active', 'inactive')"
                            },
                            "friendly_id": {
                                    "type": "string",
                                    "description": "Friendly ID of the molecule"
                            },
                            "parent_friendly_id": {
                                    "type": "string",
                                    "description": "Optional: Friendly ID of the parent molecule, only molecules created by Medicinal Chemist can have parents"
                            }
                        },
                        "required": ["structure", "rationale", "friendly_id"]
                    }
                }
            },
            "required": ["session_id", "parent_transformation_id", "transformation", "molecules", "client_name"]
        }
    },
    ToolNames.DIFFDOCK.value: {
        "name": "run_diffdock",
        "description": """Run molecular docking simulation using DiffDock. For each input SMILES:
- Generates 10 ranked conformations of ligand-protein complexes
- Returns .sdf files named 'rankN_confidence-X.sdf' where:
  - N (1-10): Rank of the conformation
  - X: Confidence score (lower is better)
- Best conformations have lower rank numbers and confidence scores
- Example: rank1_confidence-1.44.sdf is better than rank2_confidence-1.53.sdf
Note: Some SMILES may fail to dock, resulting in fewer than 10 conformations.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "protein_path": {
                    "type": "string",
                    "description": "Path to protein PDB file"
                },
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings for ligands"
                }
            },
            "required": ["protein_path", "smiles"]
        }
    },
    ToolNames.VINA_DOCKER.value: {
        "name": "run_vina_docker",
        "description": """Run molecular docking simulation using AutoDock Vina. For a protein and ligand:
- Performs docking to predict binding modes and affinities
- Returns multiple files including:
  - docking_result.pdbqt: Docked ligand poses ranked by binding energy
  - protein.pdbqt: Prepared protein structure in PDBQT format
  - ligand.pdbqt: Prepared ligand structure in PDBQT format
  - vina.conf: Configuration file used for docking
  - vina.log: Detailed log with energy scores (lower scores = stronger binding)
- Binding energies are reported in kcal/mol, with more negative values indicating better binding
- The tool accepts either a ligand file path or SMILES string(s)
- If multiple SMILES are provided, all will be processed in a batch""",
        "input_schema": {
            "type": "object",
            "properties": {
                "protein_path": {
                    "type": "string",
                    "description": "Path to protein PDB file"
                },
                "ligand_path": {
                    "type": ["string", "array"],
                    "description": "Path to ligand SDF file or list of paths (use either this OR ligand_smiles OR smiles)"
                },
                "ligand_smiles": {
                    "type": ["string", "array"],
                    "description": "SMILES string or list of SMILES strings (use either this OR ligand_path OR smiles)"
                },
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings to be processed in batch for docking"
                }
            },
            "required": ["protein_path"]
        }
    },
    ToolNames.MOLECULE_GENERATION.value: {
        "name": "generate_molecules",
        "description": "Generate potential binding molecules based on a protein sequence using Prot2Mol",
        "input_schema": {
            "type": "object",
            "properties": {
                "protein_sequence": {
                    "type": "string",
                    "description": "Amino acid sequence of the protein"
                },
                "num_molecules": {
                    "type": "integer",
                    "description": "Number of molecules to generate",
                    "default": 1
                },
                "encoder_max_length": {
                    "type": "integer",
                    "description": "Maximum length for protein encoder",
                    "default": 1000
                },
                "max_length": {
                    "type": "integer",
                    "description": "Maximum length of generated molecule sequences",
                    "default": 202
                }
            },
            "required": ["protein_sequence", "num_molecules"]
        }
    },
    ToolNames.UNIPROT_SEARCH.value: {
        "name": "search_uniprot",
        "description": "Search for protein information in UniProt database and access specific fields from the results, if no PDB entries are found, it gets the AlphaFold structure",
        "input_schema": {
            "type": "object",
            "properties": {
                "query": {
                    "type": "string",
                    "description": "Gene name (e.g., 'AKT1') or UniProt accession ID (e.g., 'P31750')"
                },
                "fields": {
                    "type": "array",
                    "items": {
                        "type": "string",
                        "enum": ["sequence", "function", "domain", "activity_regulation",
                                "interpro_ids", "pfam_ids", "bindingdb_ids", "chembl_ids",
                                "drugbank_ids", "pdb_entries"]
                    },
                    "description": "List of fields to extract from the search results"
                }
            },
            "required": ["query"]
        }
    },
    ToolNames.PDB_FILE.value: {
        "name": "get_pdb_file",
        "description": "Download PDB structure file and save it to the specified directory, share the full path of the file with the user.",
        "input_schema": {
            "type": "object",
            "properties": {
                "pdb_id": {
                    "type": "string",
                    "description": "PDB ID of the structure to download"
                },
                "file_format": {
                    "type": "string",
                    "enum": ["pdb", "cif", "xml"],
                    "description": "Format of the structure file",
                    "default": "pdb"
                }
            },
            "required": ["pdb_id"]
        }
    },
    ToolNames.PUBMED_SEARCH_AGENT.value:{
        "name": "search_pubmed_agent",
        "description": """Use this tool when you need to find relevant information from PubMed.
This tool first downloads papers from PubMed using a broad 'search_term' and then asks your 'query' to a specialized pub_med_agent (PaperQA) that reads through the papers and attempts to answer your query.
This tool should outputs a summary of relevant chunks. If you receive this summaries DO NOT recall this tool but answer the question you received.

**Recommended Usage**:
1. Keep the 'search_term' broad enough (e.g., ("CDK2" OR "Cyclin-Dependent Kinase 2")) so you retrieve a sufficient set of PDFs. This ensures the PaperQA agent can refine its internal search on those downloaded PDFs.
2. If no answer is provided or is found initially, consider expanding the 'search_term' (adding synonyms, removing overly specific terms, etc.).
3. if pub_med_agent indicates it can't answer the question (insufficient or irrelevant chunks), you can:
- Refine your 'search_term' to locate different/more relevant PDFs,
- Refine your 'query' to be more explicit, or
- Do both.
4. You can use this tool sequentially, up to three attempts. If, after three attempts, no satisfactory answer is found, consider that no relevant information is available.
5. WARNING: You should only REPEATLY use this tool if it does return either "No PDF found" or "Can't answer the question". Never call this tool more than 1 time otherwise.

        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "search_term": {
                    "type": "string",
                    "description": (
                    "A broad PubMed search string to locate relevant articles. "
                    "For example: (\"CDK2\" OR \"Cyclin-Dependent Kinase 2\"). "
                    "Avoid using too many narrowing terms. "
                    "If you receive 'No PDF found', consider adding synonyms or related terms to broaden the scope. "
                    "This term should be as Broad as Possible"
                    )
                },
                "query": {
                    "type": "string",
                    "description": (
                    "The specific information you're looking for within the retrieved PDFs. "
                    "The tool will return relevant chunks of text from the downloaded papers that match your query. "
                    "Example: 'What recent advances have been reported for CDK2 inhibitors?' "
                    "If insufficient relevant information is found within the PDFs, you'll receive a notification."
                    )
                },
                "max_pdf_count": {
                    "type": "integer",
                    "description": "Maximum number of PDFs to download, in case of no pdf found, increase this number",
                    "default": 10
                }
            },
            "required": ["search_term", "query"]
        }
    },
    ToolNames.CHEMBL_SEARCH_ACTIVITY.value: {
        "name": "search_chembl_activity",
        "description": """Use this tool to retrieve compound activity data (assay type = 'B', i.e., binding assays) from ChEMBL for a given target.
It returns up to 5 'active' (pChEMBL ≥ 6) and 5 'inactive' (pChEMBL < 6) compounds, sorted by descending pChEMBL value, in a three-column Markdown table (ID/Name, SMILES, pChEMBL).
If no SMILES or no activity data are found, a corresponding message is returned.
Use this tool after you have the ChEMBL target ID (e.g., 'CHEMBL240') to fetch example active/inactive molecules for that target.
You must return this smiles and their activity right after you have the results for the rests of your teams members to see them (return the innactive and active ones)
""",
        "input_schema": {
            "type": "object",
            "properties": {
            "target_chembl_id": {
                "type": "string",
                "description": "The ChEMBL target identifier (e.g., 'CHEMBL240') for which to retrieve up to 5 active and 5 inactive compounds."
            },
            "protein_path": {
                "type": "string",
                "description": "Path to the protein protein .pdb file downloaded via the get_pdb_file tool, this is required !"
            },
            "protein_sequence": {
                "type": "string",
                "description": "The protein sequence"
            }
            },
            "required": ["target_chembl_id", "protein_path", "protein_sequence"]
        }
    },
    ToolNames.RDIFFUSION.value: {
        "name": "rdiffusion",
        "description":"""Generate molecular structures using RDiffusion you can use this without any arguments the argument are hardcoded you just need to call this function if needed without arguments.

        RFdiffusion has been developed to design proteins for a wide range of challenges, such as creating protein binders, symmetric architectures, and enzyme scaffolds. RFdiffusion applies this diffusion to protein structures. It learns to take noisy or
        incomplete protein structures and refine them into realistic and functional designs.

        You can use only this tool so if you have hesitation use it""",
        "input_schema": {
            "type": "object",
            "properties": {},
            "required": []
        }
    },
    'get_weather_city': {
        "name": "get_weather_city",
        "description": "Get weather information for a specific city",
        "input_schema": {
            "type": "object",
            "properties": {
                "city": {
                    "type": "string",
                    "description": "Name of the city"
                },
            },
            "required": ["city"]
        }
    },
    ToolNames.PLIP_RUN.value: {
        "name": "get_plip_report",
        "description": """Get PLIP (Protein-Ligand Interaction Profiler) report for a given protein structure or structures.
                            PLIP analyzes and visualizes non-covalent protein-ligand interactions including:
                            - Hydrogen bonds
                            - Hydrophobic contacts
                            - π-stacking
                            - Salt bridges
                            - Water bridges
                            - Metal complexes
                            - Halogen bonds
                            Returns a detailed report of all detected interactions.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "input_structure": {
                    "type": "string",
                    "description": "Path or list of paths to the input PDB file(s) containing protein-ligand complex"
                },
                "wait_for_result": {
                    "type": "boolean",
                    "description": "Whether to wait for analysis completion",
                    "default": True
                }
            },
            "required": ["input_structure"]
        }
    },
    ToolNames.VALIDITY.value: {
        "name": "validity",
        "description": "Check the validity of list of molecules",
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings of the molecules"
                }
            },
            "required": ["smiles"]
        }
    },
    ToolNames.QED.value: {
        "name": "QED",
        "description": "Check the QED of list of molecules",
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings of the molecules"
                }
            },
            "required": ["smiles"]
        }
    },
    ToolNames.SA.value: {
        "name": "SA",
        "description": "Check the SA of list of molecules",
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings of the molecules"
                }
            },
            "required": ["smiles"]
        }
    },
    ToolNames.READ_FILE.value: {
        "name": "read_file",
        "description": "Read a file and return the content",
        "input_schema": {
            "type": "object",
            "properties": {
                "file_path": {
                    "type": "string",
                    "description": "Path to the file to read"
                }
            },
        }
    },
    ToolNames.AGENT_CALL.value: {
        "name": "call_agent",
        "description": """Call another agent to get their expertise on a specific task. Each agent has specific expertise and tools:
        - Database Agent: Expert in data mining (UniProt, PDB, ChEMBL)
        - Modelling Agent: Expert in AI/ML approaches for drug discovery
        - Medicinal Chemist: Expert in drug-like molecule design
        - Structural Biologist: Expert in protein structure analysis
        - Ranking Agent: Expert in prioritizing candidates
        - Critic Agent: Expert in validating scientific reasoning

        The called agent will maintain their role and use their specific tools to help with the task.
        To prevent infinite loops, there is a maximum call depth of 2 (an agent can call another agent, who can call one more agent).
        """,
        "input_schema": {
            "type": "object",
            "properties": {
                "agent_title": {
                    "type": "string",
                    "description": "The title of the agent to call (e.g., 'Database Agent', 'Literature Agent')"
                },
                "query": {
                    "type": "string",
                    "description": "The question or task for the agent"
                },
                "relevant_history": {
                    "type": "array",
                    "items": {"type": "object"},
                    "description": "Relevant messages from conversation history to provide context",
                    "default": []
                },
                "call_depth": {
                    "type": "integer",
                    "description": "Current depth of agent calls (max 2)",
                    "default": 0
                }
            },
            "required": ["agent_title", "query"]
        }
    },
    ToolNames.MOCK_TOOL.value: {
        "name": "mock_tool",
        "description": "DON'T CALL THIS TOOL. A mock tool that doesn't do anything except return the input arguments. No need to use this tool",
        "input_schema": {
            "type": "object",
            "properties": {
                "any_input": {
                    "type": "string",
                    "description": "Any input string"
                },
                "optional_number": {
                    "type": "integer",
                    "description": "An optional number parameter",
                    "default": 42
                }
            },
            "required": ["any_input"]
        }
    },
    ToolNames.VISUALIZE_MOLECULES.value: {
        "name": "visualize_molecules",
        "description": "Convert SMILES strings to molecular structure images for visualizatio. Can create multiple images of multiple molecules at once.",
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings to visualize"
                },
                "labels": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Optional labels for each molecule (same length as smiles array)"
                }
            },
            "required": ["smiles"]
        }
    },
    ToolNames.TEST_MOLECULE_VIZ.value: {
        "name": "test_molecule_visualization",
        "description": "Test the molecule visualization system by returning some example molecules",
        "input_schema": {
            "type": "object",
            "properties": {
                "num_molecules": {
                    "type": "integer",
                    "description": "Number of example molecules to return",
                    "default": 3
                }
            }
        }
    },
    ToolNames.MOLECULAR_SIMILARITY.value: {
        "name": "molecular_similarity",
        "description": "Calculate the Tanimoto similarity between two molecules",
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles1": {
                    "type": "list",
                    "description": "List of SMILES strings of the first molecule"
                },
                "smiles2": {
                    "type": "list",
                    "description": "List of SMILES strings of the second molecule"
                }
            },
            "required": ["smiles1", "smiles2"]
        }
    },
    ToolNames.MOLECULAR_WEIGHT.value: {
        "name": "molecular_weight",
        "description": "Calculate the molecular weight of list of molecules",
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings of the molecules"
                }
            },
            "required": ["smiles"]
        }
    },
    ToolNames.LOG_P.value: {
        "name": "logP",
        "description": "Calculate the logP of list of molecules",
        "input_schema": {
            "type": "object",
            "properties": {
                "smiles": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "List of SMILES strings of the molecules"
                }
            },
            "required": ["smiles"]
        }
    },
    ToolNames.MOL_GEN.value: {
        "name": "mol_gen",
        "description": """Generate drug-like molecules based on a protein target, filter them by molecular properties, and perform docking.
This comprehensive tool:
- Generates potential binding molecules using the provided protein sequence
- Filters molecules based on drug-likeness properties (QED, synthetic accessibility, logP, molecular weight)
- Performs molecular docking using DiffDock to evaluate binding potential
- Analyzes protein-ligand interactions using PLIP
- Returns detailed information including:
  - SMILES strings of the generated molecules
  - Drug-likeness scores (QED, SA, logP, MW)
  - Docking confidence (confidence higher is better, we are looking for positive values)
  - Protein-ligand interaction analysis
  - Paths to generated complex structures
The entire process is automated with smart filtering to prioritize promising candidates.""",
        "input_schema": {
            "type": "object",
            "properties": {
                "protein_path": {
                    "type": "string",
                    "description": "Path to protein PDB file for docking"
                },
                "protein_sequence": {
                    "type": "string",
                    "description": "Amino acid sequence of the protein for molecule generation"
                },
                "num_molecules": {
                    "type": "integer",
                    "description": "Number of molecules to initially generate",
                    "default": 50
                },
                "molecules": {
                    "type": "array",
                    "items": {"type": "string"},
                    "description": "Optional list of SMILES strings (skips generation of de novo molecules if provided, and returns the provided molecules (not de novo molecules), you can use this to provide known active/inactive molecules if you want to evaluate them"
                },
                "thresholds": {
                    "type": "object",
                    "description": "Filtering criteria for molecule properties, you can base these criteria based on the known active/inactive molecules in the ChEMBL database",
                    "properties": {
                        "qed": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum QED value (0-1, higher is more drug-like)"},
                                "max": {"type": "number", "description": "Maximum QED value"}
                            }
                        },
                        "sa": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum synthetic accessibility score"},
                                "max": {"type": "number", "description": "Maximum synthetic accessibility score (lower is easier to synthesize)"}
                            }
                        },
                        "logp": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum logP value"},
                                "max": {"type": "number", "description": "Maximum logP value"}
                            }
                        },
                        "molecular_weight": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum molecular weight (Da)"},
                                "max": {"type": "number", "description": "Maximum molecular weight (Da)"}
                            }
                        }
                    }
                }
            },
            "required": ["protein_path", "protein_sequence", "thresholds", "num_molecules"]
        }
    },
    ToolNames.VINA_GEN_MOL.value: {
        "name": "vina_genmol",
        "description": """Generate drug-like molecules based on known ligands, filter them by molecular properties, and perform docking.
This comprehensive tool:
- Generates potential de novo molecules using NVIDIA GenMol API. 
- Filters molecules based on drug-likeness properties (QED, synthetic accessibility, logP, molecular weight)
- Performs molecular docking using Vina to evaluate binding potential
- Returns detailed information including:
  - SMILES strings of the generated molecules
  - Drug-likeness scores (QED, SA, logP, MW)
  - Docking score 
  - Protein-ligand interaction analysis
""",
        "input_schema": {
            "type": "object",
            "properties": {
                "protein_path": {
                    "type": "string",
                    "description": "Path to protein PDB file for docking"
                },
                "num_molecules": {
                    "type": "integer",
                    "description": "Number of molecules to generate",
                    "default": 10
                },
                "temperature": {
                    "type": "number",
                    "description": "Temperature parameter controlling randomness in generation",
                    "default": 1.0
                },
                "noise": {
                    "type": "number",
                    "description": "Noise parameter for molecule generation",
                    "default": 0.0
                },
                "step_size": {
                    "type": "number",
                    "description": "Step size parameter for molecule generation",
                    "default": 1.0
                },
                "scoring": {
                    "type": "string",
                    "description": "Scoring method for molecule evaluation",
                    "default": "QED"
                },
                "unique_only": {
                    "type": "boolean",
                    "description": "Whether to return only unique molecules",
                    "default": True
                },
                "thresholds": {
                    "type": "object",
                    "description": "Filtering criteria for molecule properties, you can base these criteria based on the known active/inactive molecules in the ChEMBL database",
                    "properties": {
                        "qed": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum QED value (0-1, higher is more drug-like)"},
                                "max": {"type": "number", "description": "Maximum QED value"}
                            }
                        },
                        "sa": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum synthetic accessibility score"},
                                "max": {"type": "number", "description": "Maximum synthetic accessibility score (lower is easier to synthesize)"}
                            }
                        },
                        "logp": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum logP value"},
                                "max": {"type": "number", "description": "Maximum logP value"}
                            }
                        },
                        "molecular_weight": {
                            "type": "object",
                            "properties": {
                                "min": {"type": "number", "description": "Minimum molecular weight (Da)"},
                                "max": {"type": "number", "description": "Maximum molecular weight (Da)"}
                            }
                        }
                    }
                }
            },
            "required": ["num_molecules", "protein_path", "thresholds"]
        }
    },
}

from functools import wraps
import logging
from src.tools.gurnemanz_utils import MoleculeIDGenerator

class ToolRegistry:
    """Registry for all available tools and their implementations."""

    def __init__(self):
        self._tools: Dict[str, Callable] = {}
        self._previous_sequence: str | None = None
        parent_dir = Path(__file__).resolve().parents[2]
        self.output_dir = parent_dir / "runs_metadata"
        self._id_generator: MoleculeIDGenerator = MoleculeIDGenerator()
        print(f'Output dir: {self.output_dir}')

    def register(self, name: str) -> Callable:
        """Decorator to register a new tool."""
        def decorator(func: Callable) -> Callable:
            @wraps(func)
            def wrapper(*args, **kwargs) -> Dict[str, Any]:
                logging.info(f"Executing tool: {name}")
                try:
                    return func(*args, **kwargs)
                except Exception as e:
                    logging.error(f"Error executing {name}: {str(e)}")
                    return {"error": str(e)}

            self._tools[name] = wrapper
            return wrapper
        return decorator

    def execute(self, tool_name: str, arguments: dict, run_id: str) -> dict:
        """Execute a registered tool by name."""
        if tool_name not in self._tools:
            return {"error": f"Tool '{tool_name}' not found"}

        return self._tools[tool_name](arguments, run_id)

    @property
    def available_tools(self) -> list[str]:
        """List all registered tools."""
        return list(self._tools.keys())

    @property
    def previous_sequence(self) -> str | None:
        return self._previous_sequence

    @previous_sequence.setter
    def previous_sequence(self, value: str):
        self._previous_sequence = value

tool_registry = ToolRegistry()

@tool_registry.register("call_agent")
def call_agent(arguments: dict, output_dir_id: str) -> dict:
    """Execute a call to another agent."""


    # Check call depth to prevent infinite loops
    call_depth = arguments.get('call_depth', 0)
    if call_depth >= 2:
        return {
            "error": "Maximum agent call depth exceeded",
            "message": "Cannot make further agent calls as maximum depth (2) has been reached."
        }

    # Get the agent mapping from the global scope
    try:
        from __main__ import agents_mapping
    except ImportError:
        return {"error": "Agent mapping not found in global scope"}

    agent_title = arguments['agent_title']
    query = arguments['query']
    relevant_history = arguments.get('relevant_history', [])

    # Find the requested agent
    agent = None
    for a in agents_mapping.values():
        if a.title == agent_title:
            agent = a
            break

    if not agent:
        return {"error": f"Agent '{agent_title}' not found"}

    # Create a focused conversation history
    focused_history = []
    if relevant_history:
        focused_history.extend(relevant_history)

    # Add the current query
    focused_history.append({
        "role": "user",
        "content": query
    })

    # Increment call depth for the next agent
    next_call_depth = call_depth + 1

    # Create a copy of the agent's tools with updated call_depth
    agent_tools = copy.deepcopy(agent.tools) if agent.tools else []
    for tool in agent_tools:
        if tool['name'] == 'call_agent':
            tool['input_schema']['properties']['call_depth']['default'] = next_call_depth

    # Temporarily update agent's tools
    original_tools = agent.tools
    agent.tools = agent_tools

    try:
        # Run the conversation with the agent
        response = agent.run_conversation(
            user_message=query,
            history=focused_history,
            run_id=output_dir_id
        )

        # Format the response to include both the answer and any tool results
        return {
            "agent": agent_title,
            "response": response.content,
            "results": response.results if response.results else [],
            "call_depth": next_call_depth
        }
    finally:
        # Restore original tools
        agent.tools = original_tools

@tool_registry.register("mock_tool")
def mock_tool(arguments: dict, output_dir_id: str) -> dict:
    print("Mock tool is being used")
    return {
        "message": "This is a mock tool response",
        "received_arguments": arguments,
        "output_dir_id": output_dir_id
    }

@tool_registry.register("read_file")
def read_file(arguments: dict, output_dir_id: str) -> dict:
    print("Read file is being used")
    return {"file_content": open(arguments["file_path"], "r").read()}
