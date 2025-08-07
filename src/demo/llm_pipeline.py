from dotenv import load_dotenv
load_dotenv()
import sys
import os
import json
import argparse

# Get the base directory from the current file location (to make paths relative)
# --- Add project root to sys.path ---
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
print(project_root)
sys.path.insert(0, project_root)
# --- End of sys.path modification ---
print(project_root)
run = "litellm"

from src.agent import Agent
from src.tools.tool_definitions import ToolNames, ANTHROPIC_TOOLS
run_map = {
    "litellm": ANTHROPIC_TOOLS,
    "litellm": ANTHROPIC_TOOLS
}

import traceback
import nest_asyncio
nest_asyncio.apply()
from datetime import datetime
import random
import string
import signal
import copy
from src.demo.multi_agent_prompt import PARSER_PROMPT_TEMPLATE
from src.demo.utils import (setup_logger, 
                            save_conversation_to_markdown, 
                            generate_run_id, 
                            timeout_handler,
                            format_agent_query)
# --- Gurnemanz tools ---
from src.tools.gurnemanz import gurnemanz_apply
from src.tools.gurnemanz_utils import initialize_session, submit_to_gurnemanz, format_gurnemanz_for_llm


def create_parser():
    parser = argparse.ArgumentParser(description="Simple LLM drug discovery pipeline (no tools)")
    parser.add_argument(
        "--llm_provider", 
        type=str, 
        default="sonnet-3.7",
        choices=["sonnet-4", "sonnet-3.7", "o3", "gpt4.1", "gemini", "qwen", "deepseek", "mistral"],
        help="LLM provider to use (default: sonnet-3.7)"
    )
    parser.add_argument(
        "--iteration_num",
        type=int,
        default=3,
        help="Number of iterations to run (default: 3)"
    )
    return parser

# Parse command line arguments
parser = create_parser()
args = parser.parse_args()

llm_provider = args.llm_provider
max_iterations = args.iteration_num

# Configuration
TIMEOUT_SECONDS = 36000  # 10 hour timeout

# LLM Prompts for direct molecule generation
llm_system_prompt = """You are an expert medicinal chemist and computational drug discovery scientist. Your task is to design novel small molecule inhibitors for AKT1 (Protein Kinase B alpha) using your knowledge of structure-activity relationships, medicinal chemistry principles, and drug design.

**Target Information:**
- AKT1 is a serine/threonine kinase in the PI3K/AKT/mTOR pathway
- Critical for cell survival, proliferation, and metabolism
- Key oncology target with multiple known inhibitors
- ATP-competitive and allosteric binding sites available

**Your Goal:**
Generate novel, drug-like small molecule AKT1 inhibitors with the following criteria:
- High predicted binding affinity to AKT1
- Drug-like properties (Lipinski's Rule of Five compliance)
- Synthetic accessibility
- Novel scaffolds not directly copied from known inhibitors

**Output Format:**
For each molecule, provide:
1. **Molecule ID**: A unique identifier (e.g., LLM_MOL_001)
2. **SMILES**: Valid SMILES string
3. **Rationale**: Brief explanation of your design strategy and expected properties
4. **Predicted Properties**: Your assessment of drug-likeness, potency, selectivity

Focus on innovative chemical designs while maintaining scientific rigor."""

llm_iteration_prompts = {
    1: """**ITERATION 1: Initial Discovery & Exploration**

Design 10-15 novel AKT1 inhibitor molecules using your medicinal chemistry expertise. Consider:

1. **ATP-competitive inhibitors**: Target the ATP-binding pocket with novel scaffolds
2. **Allosteric inhibitors**: Design molecules that bind to regulatory sites
3. **Hybrid approaches**: Molecules that may interact with multiple binding sites

Use your knowledge of:
- Known AKT1 pharmacophores (avoid direct copying)
- Kinase selectivity determinants
- Drug-like chemical space
- Structure-activity relationships

Generate diverse scaffolds including heterocycles, aromatic systems, and appropriate linkers. Ensure synthetic feasibility while prioritizing novelty.""",

    2: """**ITERATION 2: Lead Optimization**

Based on your previous designs, now optimize and refine your best candidates. Focus on:

1. **Scaffold hopping**: Create variants of your most promising designs
2. **Property optimization**: Improve drug-likeness, solubility, permeability
3. **Selectivity enhancement**: Design features for AKT1 vs AKT2/AKT3 selectivity
4. **ADMET optimization**: Consider metabolic stability and safety

Generate 8-12 optimized molecules that build upon your iteration 1 insights. Include:
- Modified versions of your best scaffolds
- New analogs with improved properties
- Rationale for each modification based on medicinal chemistry principles""",

    3: """**ITERATION 3: Final Lead Optimization & Selection**

Create your final set of 8-10 candidate molecules representing your best designs for AKT1 inhibition. Focus on:

1. **Best-in-class designs**: Your most promising scaffolds with optimal substitution patterns
2. **Diverse mechanisms**: Include both ATP-competitive and allosteric candidates
3. **Clinical viability**: Emphasize drug-like properties and synthetic accessibility
4. **Innovation**: Ensure your designs represent novel chemical space

For each final candidate, provide:
- Comprehensive rationale for inclusion
- Predicted binding mode and key interactions
- Assessment of development potential
- Comparison to known AKT1 inhibitors (without direct copying)
- Don't forget to use the molecule entry format"""
}

# Initialize the LLM agent (no tools)
llm_agent = Agent(
    title="LLM Drug Designer",
    expertise=(
        "Medicinal chemistry and drug design",
        "Structure-activity relationships",
        "Kinase inhibitor development",
        "ADMET optimization",
        "Computational drug discovery"
    ),
    goal=(
        "Design novel, drug-like small molecule AKT1 inhibitors using medicinal chemistry "
        "knowledge and computational drug design principles, without relying on external tools."
    ),
    role=(
        "Expert medicinal chemist responsible for generating innovative AKT1 inhibitor designs "
        "based on scientific knowledge, literature understanding, and drug design principles. "
        "Provide novel chemical structures with clear design rationales."
    ),
    llm_provider=llm_provider,
    tools=[]  # No tools - direct LLM only
)

# Parser agent to extract molecules from LLM responses
summary_parser_agent = Agent(
    title="Summary Parser",
    expertise=(
        "Interpreting unstructured research outputs",
        "Generating structured JSON data",
        "Maintaining data consistency across transformations"
    ),
    goal=(
        "Format LLM outputs into GURNEMANZ JSON schema and submit them while maintaining "
        "all important relationships between data elements."
    ),
    role=(
        "Parse LLM outputs to construct a clean formatted JSON that preserves "
        "all relevant information about molecules and their design rationales."
    ),
    llm_provider="sonnet-3.7",
    tools=[run_map[run][ToolNames.GURNEMANZ_APPLY.value]]
)

def run_llm_pipeline(max_iterations=3, run_id="999"):
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(TIMEOUT_SECONDS)
    run_id = f'{run_id}_llm_{llm_provider}_{max_iterations}'
    run_dir = os.path.join(project_root, 'runs_metadata', run_id)
    os.makedirs(run_dir, exist_ok=True)
    total_price = 0.0
    
    # Set up logger
    logger = setup_logger(run_dir, run_id)
    logger.info(f"Starting LLM pipeline with run_id: {run_id}")
    
    # Initialize molecules list for tracking across iterations
    all_molecules = []
    
    # Initialize GURNEMANZ session
    session_id, initial_transform_id = initialize_session(f"user-{run_id}")
    print(f"GURNEMANZ session initialized: {session_id}")
    gurnemanz_session = {
        "session_id": session_id,
        "last_transform_id": initial_transform_id
    }

    # Save the session info for the tool to use
    session_file = os.path.join(run_dir, "gurnemanz_session.json")
    try:
        with open(session_file, 'w') as f:
            json.dump(gurnemanz_session, f)
        logger.info(f"Saved GURNEMANZ session info to {session_file}")
    except Exception as e:
        logger.error(f"Error saving session file: {e}")
        session_file = None

    # Initialize conversation history
    main_history = [{"role": "user", "content": llm_system_prompt}]

    try:
        for iteration in range(max_iterations):
            logger.info('--------------------------------')
            logger.info(f'Starting Iteration {iteration + 1}')
            logger.info('--------------------------------')
            
            # Get the prompt for this iteration
            iteration_prompt = llm_iteration_prompts.get(iteration + 1, 
                f"Continue designing AKT1 inhibitors for iteration {iteration + 1}")
            
            logger.info(f'LLM Agent: {llm_agent.title}')
            
            # Direct LLM call - no tools, just conversation
            response = llm_agent.run_conversation(
                iteration_prompt, 
                history=main_history, 
                run_id=run_id, 
                logger=logger
            )
            total_price += response.price if response and hasattr(response, 'price') else 0

            # Parse the LLM response to extract molecules
            if response is not None and hasattr(response, 'content'):
                logger.info(f"main_history before parser: {json.dumps(main_history, indent=2)}")
                parser_prompt = PARSER_PROMPT_TEMPLATE.format(
                    agent_output=response.content,
                    iteration=iteration + 1,
                    agent_title=llm_agent.title
                )

                parser_response = summary_parser_agent.run_conversation(
                    parser_prompt,
                    history=copy.deepcopy(main_history),
                    run_id=run_id,
                    logger=logger
                )
                print("Summary parser response: ", parser_response.content)
                logger.info(f"main_history after parser: {json.dumps(main_history, indent=2)}")
                parser_price = parser_response.price if parser_response and hasattr(parser_response, 'price') else 0
                total_price += parser_price

                # Process Gurnemanz session data
                if session_file:
                    try:
                        with open(session_file, 'r') as f:
                            session_data = json.load(f)
                        raw_json = session_data.get("latest_raw_json")
                        if raw_json and 'molecules' in raw_json:
                            for molecule in raw_json['molecules']:
                                if 'moleculeId' in molecule:
                                    del molecule['moleculeId']
                                if molecule not in all_molecules:
                                    all_molecules.append(molecule)
                    except Exception as e:
                        print(f"Error loading session file: {str(e)}")

            # Save the current iteration's history
            save_conversation_to_markdown(main_history, run_dir, f"iteration_{iteration + 1}_history.md")
            
            logger.info(f"Completed Iteration {iteration + 1}")

        # Final summary
        logger.info("LLM pipeline completed successfully")
        final_summary = f"""
# LLM Pipeline Final Summary

## Project Overview
- **Target**: AKT1 inhibitor discovery
- **Approach**: Direct LLM-based molecule generation (no tools)
- **Iterations**: {max_iterations}
- **Total molecules generated**: {len(all_molecules)}

## Key Insights
This pipeline demonstrates the capability of large language models to:
1. Generate novel chemical structures based on medicinal chemistry knowledge
2. Provide scientific rationales for design decisions
3. Iterate and optimize molecular designs across multiple rounds

## Generated Molecules
Total unique molecules extracted: {len(all_molecules)}

This represents a pure LLM-based approach for comparison in ablation studies against tool-augmented pipelines.
"""
        
        with open(os.path.join(run_dir, "final_summary.md"), 'w') as f:
            f.write(final_summary)
        
        logger.info(final_summary)

    except TimeoutError:
        print("LLM pipeline stopped due to timeout")
        logger.warning("LLM pipeline stopped due to timeout")
    except Exception as e:
        print(f"Unexpected error in LLM pipeline: {str(e)}")
        print(f"Traceback: {traceback.format_exc()}")
        logger.error(f"Unexpected error in LLM pipeline: {str(e)}", exc_info=True)
    finally:
        signal.alarm(0)

    return main_history, total_price

# Run the LLM pipeline
if __name__ == "__main__":
    history, total_price = run_llm_pipeline(
        max_iterations=max_iterations,
        run_id=generate_run_id()
    )
    print(f"Total price of the LLM pipeline: {total_price:.2f} USD")
