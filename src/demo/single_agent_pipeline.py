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
from src.demo.single_agent_prompt import (
    single_agent_system_prompt,
    single_agent_iteration_1_kickoff,
    single_agent_next_iteration_system_context,
    single_agent_next_iteration_kickoff,
    single_agent_summary_prompt,
    single_agent_final_summary_prompt
)
from src.demo.multi_agent_prompt import PARSER_PROMPT_TEMPLATE
from src.demo.utils import (setup_logger, 
                            save_conversation_to_markdown, 
                            generate_run_id, 
                            timeout_handler,
                            format_agent_query)
# --- Database agent tools ---
from src.tools.api_based_tools.fetch_chembl import search_chembl_activity
from src.tools.api_based_tools.pdb_search import get_pdb_file
from src.tools.api_based_tools.literature_api import search_uniprot, search_pubmed_agent
# --- AI Expert tools ---
from src.tools.prot2mol_vina import get_vina_mol_gen_report
from src.tools.vina_plip import get_vina_report
# --- Gurnemanz tools ---
from src.tools.gurnemanz import gurnemanz_apply
from src.tools.gurnemanz_utils import initialize_session, submit_to_gurnemanz, format_gurnemanz_for_llm


def create_parser():
    parser = argparse.ArgumentParser(description="Multi-agent drug discovery pipeline")
    parser.add_argument(
        "--llm_provider", 
        type=str, 
        default="sonnet-3.7",
        choices=["sonnet-4", "sonnet-3.7", "o3", "gpt4.1", "gemini", "qwen", "deepseek", "mistral"],
        help="LLM provider to use (default: sonnet-4)"
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

# Initialize main_history with structured GURNEMANZ JSON
initial_gurnemanz = {
    "transformation": {
        "type": "initial_prompt",
        "agent": "System",
        "userMessage": "", #TODO meeting_1_prompt_new["content"],
        "iteration": 1
    },
    "molecules": []
}

# Format initial GURNEMANZ data for readability
formatted_initial = format_gurnemanz_for_llm(initial_gurnemanz)
main_history = [{"role": "user", "content": [{'type':"text", "text":single_agent_system_prompt["content"].format(max_iterations=max_iterations)}]}]


single_agent = Agent(
    title="One Person Biotech",
    expertise=(
        "Leading *in silico* drug discovery projects and strategic planning",
        "Integrating computational data from various drug discovery phases (database mining, AI generation, molecular modeling)",
        "Setting objectives for iterative research cycles and synthesizing progress",
        "Accessing molecular and protein databases (UniProt, PDB, ChEMBL)",
        "Retrieving structural files, known inhibitors, protein sequences, and activity data",
        "AI and computational chemistry",
        "AI-driven molecule generation and structure-based design", 
        "Small-molecule chemistry and drug-likeness assessment",
        "Designing feasible synthetic routes and optimizing molecules for potency and selectivity",
        "Molecular docking and binding affinity prediction",
        "Critical evaluation of scientific ideas and ensuring research rigor"
    ),
    goal=(
        "Orchestrate a comprehensive *in silico* drug discovery project towards the identification "
        "of 10 promising drug candidate molecules. This involves strategic planning, database research, "
        "AI-driven molecule generation, and medicinal chemistry optimization. Set iteration-specific goals, "
        "search molecular databases for target information, generate novel drug candidates using AI models, "
        "refine molecules based on chemical feasibility and drug-like properties, and synthesize progress "
        "across iterations while presenting optimized molecules with supporting computational data."
    ),
    role=(
        "As a one-person biotech researcher, you are the strategic leader and executor of this "
        "*in silico* drug discovery project. You must combine strategic oversight with hands-on "
        "computational work. Your comprehensive responsibilities include: "
        ""
        "**STRATEGIC LEADERSHIP:** "
        "- Define research strategy and specific objectives for each iteration "
        "- Synthesize progress across iterations and adapt strategy based on findings "
        "- Ensure project continuity and maintain focus on the goal of identifying 10 promising candidates "
        "- This project is strictly *in silico* - base all analyses on computational methods only "
        ""
        "**HANDS-ON COMPUTATIONAL WORKFLOW:** "
        "1. **Database Research:** Use UniProt, PDB (select 4EJN), and ChEMBL tools to gather target protein "
        "information, structural data, and known active/inactive compounds. "
        "2. **AI-Driven Generation (Iteration 1 only):** Generate novel molecules using the "
        "VINA_MOL_GEN tool ONCE based on the target information and structural insights. "
        "3. **Manual Medicinal Chemistry Modifications:** In subsequent iterations, propose "
        "and design structural modifications to existing molecules yourself WITHOUT using "
        "generation tools. Apply medicinal chemistry principles to improve ADMET properties, "
        "selectivity, and binding affinity. "
        "4. **Computational Evaluation:** Use the VINA_REPORT tool to evaluate all molecules "
        "(both AI-generated and manually modified) for binding affinity and drug-like properties. "
        "5. **Iterative Improvement:** Based on computational results, continue manual "
        "optimization and re-evaluation. "
        ""
        "**CRITICAL FORMATTING REQUIREMENTS:** "
        "- Always include **PDB file path**: <path> and **Protein sequence**: <sequence> "
        "- When presenting molecules, use this exact format: "
        "  * **Molecule Name/ID**: [friendly_id] "
        "  * **SMILES**: [smiles_string] "
        "  * **Docking Score**: [score] (if available) "
        "  * **Rationale**: [modification reasoning] "
        ""
        "**TOOL USAGE STRATEGY:** "
        "- Use VINA_MOL_GEN only in iteration 1 for initial molecule generation "
        "- For iterations 2+, manually design modifications based on structure-activity relationships "
        "- Always use VINA_REPORT to evaluate both generated and modified molecules "
        "- You can evaluate multiple modified molecules in a single VINA_REPORT call "
        ""
        "**SCIENTIFIC VALIDATION:** "
        "- Critically evaluate your own scientific reasoning and computational results "
        "- Identify potential gaps, limitations, or inconsistencies in your approach "
        "- Ensure scientific rigor and alignment with project objectives throughout the process "
    ),
    llm_provider=llm_provider,
    tools=[
        run_map[run][ToolNames.UNIPROT_SEARCH.value],
        run_map[run][ToolNames.PDB_FILE.value],
        run_map[run][ToolNames.CHEMBL_SEARCH_ACTIVITY.value],
        run_map[run][ToolNames.VINA_MOL_GEN.value],
        run_map[run][ToolNames.VINA_REPORT.value],
    ]
)

summary_parser_agent = Agent(
    title="Summary Parser",
    expertise=(
        "Interpreting unstructured research outputs",
        "Generating structured JSON data",
        "Maintaining data consistency across transformations"
    ),
    goal=(
        "Format agent outputs into GURNEMANZ JSON schema and submit them while maintaining "
        "all important relationships between data elements."
    ),
    role=(
        "Parse agent outputs and history to construct a clean formatted JSON that preserves "
        "all relevant information and relationships between molecules and transformations."
    ),
    llm_provider="sonnet-3.7",
    tools=[run_map[run][ToolNames.GURNEMANZ_APPLY.value]]
)

agents_mapping = {
    'single': single_agent,
    'summary_parser': summary_parser_agent
}

# This is simplified for the single-agent pipeline
cycle_agents = [{'single': single_agent}] * max_iterations

def run_research_cycle(initial_query, agent: Agent, main_history=main_history, max_iterations=3, run_id="999"):
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(TIMEOUT_SECONDS)
    run_id = f'{run_id}_single_agent_{llm_provider}_{max_iterations}'
    run_dir = os.path.join(project_root, 'runs_metadata', run_id)
    os.makedirs(run_dir, exist_ok=True)
    total_price = 0.0
    history_for_saving = []  # Maintain a separate copy for saving history across iterations
    
    # Set up logger
    logger = setup_logger(run_dir, run_id)
    logger.info(f"Starting research cycle with run_id: {run_id}")
    
    # Initialize molecules list and protein_data dictionary for tracking across iterations
    data2passIteration = {'molecules': [], 'protein_data': {
        "sequence": None,
        "pdb_path": None
    }}
    protein_data = {
        "sequence": None,
        "pdb_path": None
    }
    # Initialize GURNEMANZ session
    session_id, initial_transform_id = initialize_session(f"user-{run_id}")
    print(f"GURNEMANZ session initialized: {session_id}")
    gurnemanz_session = {
        "session_id": session_id,
        "last_transform_id": initial_transform_id  # Track the latest transform ID
    }

    # Save the session info for the tool to use
    session_file = os.path.join(run_dir, "gurnemanz_session.json")
    try:
        with open(session_file, 'w') as f:
            json.dump(gurnemanz_session, f)
        logger.info(f"Saved GURNEMANZ session info to {session_file}")
    except Exception as e:
        logger.error(f"Error saving session file: {e}")
        session_file = None  # Mark session_file as None to indicate failure

    try:
        for iteration in range(max_iterations):
            logger.info('--------------------------------')
            logger.info(f'Starting Iteration {iteration + 1}')
            logger.info('--------------------------------')
            
            # Dynamically set tools based on iteration to enforce the project phase
            if iteration == 0:  # Iteration 1: Discovery & Exploration
                # All tools are available
                agent.tools = [
                    run_map[run][ToolNames.UNIPROT_SEARCH.value],
                    run_map[run][ToolNames.PDB_FILE.value],
                    run_map[run][ToolNames.CHEMBL_SEARCH_ACTIVITY.value],
                    run_map[run][ToolNames.VINA_MOL_GEN.value],
                    run_map[run][ToolNames.VINA_REPORT.value],
                ]
                logger.info("Iteration 1: All tools are available for discovery and exploration.")
            else:  # Iterations 2 & 3: Lead Optimization & Final Selection
                # Tools are restricted to evaluation only
                agent.tools = [
                    run_map[run][ToolNames.VINA_REPORT.value],
                ]
                logger.info(f"Iteration {iteration + 1}: Tools restricted to VINA_REPORT for focused optimization.")
            
            logger.info(f'Agent title: {agent.title}')
            
            # The first iteration has a unique kickoff prompt
            current_query = initial_query if iteration == 0 else single_agent_next_iteration_kickoff.format(current_iteration_number=iteration + 1)
            
            logger.info(f"main_history beginning of iteration {iteration+1}: {json.dumps(main_history, indent=2)}")
            
            response = agent.run_conversation(current_query, history=main_history, run_id=run_id, logger=logger)
            total_price += response.price if response and hasattr(response, 'price') else 0

            # Post-process to add file paths and sequences if they were generated
            if response is not None and hasattr(response, 'content'):
                try:
                    pdb_dir = os.path.join(run_dir, "pdb_files")
                    if os.path.exists(pdb_dir):
                        pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
                        if pdb_files:
                            pdb_file_path = os.path.join(pdb_dir, pdb_files[0])
                            if pdb_file_path not in response.content:
                                response.content += f"\n\n**PDB file path**: {pdb_file_path}"
                            # Store PDB path in protein_data for passing between iterations
                            if not protein_data['pdb_path']:
                                protein_data['pdb_path'] = pdb_file_path
                    
                    uniprot_dir = os.path.join(run_dir, "uniprot_entries")
                    if os.path.exists(uniprot_dir):
                        uniprot_files = [f for f in os.listdir(uniprot_dir) if f.endswith('.json')]
                        if uniprot_files:
                            uniprot_file_path = os.path.join(uniprot_dir, uniprot_files[0])
                            with open(uniprot_file_path, 'r') as f:
                                uniprot_data = json.load(f)
                                if 'extracted_features' in uniprot_data and 'sequence' in uniprot_data['extracted_features']:
                                    protein_sequence = uniprot_data['extracted_features']['sequence']
                                    if protein_sequence and protein_sequence not in response.content:
                                         response.content += f"\n\n**Protein sequence**: {protein_sequence}"
                                    # Store protein sequence in protein_data for passing between iterations
                                    if not protein_data['sequence']:
                                        protein_data['sequence'] = protein_sequence
                except Exception as e:
                    logger.warning(f"Could not append file paths or sequence to response: {e}")

            # After agent runs, parse its output to update Gurnemanz, unless it's a summary-only turn
            if response is not None and hasattr(response, 'content'):
                logger.info(f"main_history before parser: {json.dumps(main_history, indent=2)}")
                parser_prompt = PARSER_PROMPT_TEMPLATE.format(
                    agent_output=response.content,
                    iteration=iteration + 1,
                    agent_title=agent.title
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
                                if molecule not in data2passIteration['molecules']:
                                    data2passIteration['molecules'].append(molecule)
                        
                        if raw_json and 'transformation' in raw_json:
                            trans = raw_json['transformation']
                            if 'proteinSequence' in trans and not protein_data['sequence']:
                                protein_data['sequence'] = trans['proteinSequence']
                            if 'proteinPath' in trans and not protein_data['pdb_path']:
                                protein_data['pdb_path'] = trans['proteinPath']
                        
                        # Always update protein_data in data2passIteration
                        data2passIteration['protein_data'] = protein_data
                    except Exception as e:
                        print(f"Error loading session file: {str(e)}")

            # Save the current iteration's history
            save_conversation_to_markdown(main_history, run_dir, f"iteration_{iteration + 1}_history.md")
            
            # After each iteration (except the last), generate a summary and plan for the next
            if iteration < max_iterations - 1:
                tool_restrictions = "AI Generation and broad database searches are disabled." if iteration + 1 >= 1 else "All tools are available."
                summary_query = single_agent_summary_prompt.format(
                    current_iteration_number=iteration + 1,
                    max_iterations=max_iterations,
                    next_iteration_number=iteration + 2,
                    tool_restrictions_for_next_iteration=tool_restrictions
                )
                response = agent.run_conversation(
                    summary_query,
                    history=main_history,
                    run_id=run_id,
                    logger=logger
                )
                summary = response.content
                
                # Prepare for the next iteration
                iteration_kickoff_prompt = single_agent_next_iteration_system_context.format(
                    current_iteration_number=iteration + 2,
                    max_iterations=max_iterations,
                    pr_summary_from_previous_iteration=summary,
                    molecule_highlights=json.dumps(data2passIteration, indent=3)
                )
                total_price += response.price if response and hasattr(response, 'price') else 0
                logger.info(f"Summary for Iteration {iteration+1}:\n{summary}")
                
                if response is None or not hasattr(response, 'content'):
                    print("Warning: Received invalid response for summary. Halting.")
                    break
                
                # Reset history for the next iteration, starting with the new context
                main_history = [{
                    "role": "user",
                    "content": iteration_kickoff_prompt
                }]
                logger.info(f"History reset for Iteration {iteration+2}. New context loaded.")
            else:
                # For the final iteration, generate the final report
                logger.info(f"End of research cycle. Generating final summary.")
                summary_query = single_agent_final_summary_prompt.format(max_iterations=max_iterations)
                response = agent.run_conversation(
                    summary_query,
                    history=main_history,
                    run_id=run_id,
                    logger=logger
                )
                logger.info(f"Final Report:\n{response.content}")

    except TimeoutError:
        print("Research cycle stopped due to timeout")
        logger.warning("Research cycle stopped due to timeout")
    except Exception as e:
        print(f"Unexpected error in research cycle: {str(e)}")
        print(f"Traceback: {traceback.format_exc()}") # Add this line for detailed traceback
        logger.error(f"Unexpected error in research cycle: {str(e)}", exc_info=True)
    except BaseException as e:
        print(f"A critical error occurred, forcing exit: {str(e)}")
        print(f"Traceback: {traceback.format_exc()}")
        logger.critical(f"A critical error occurred, forcing exit: {str(e)}", exc_info=True)
    finally:
        signal.alarm(0)

    return main_history, total_price

# Run the research cycle
if __name__ == "__main__":
    history_2, total_price = run_research_cycle(
        initial_query=single_agent_iteration_1_kickoff["content"].format(max_iterations=max_iterations, current_iteration=1),
        agent=single_agent,
        main_history=main_history,
        max_iterations=max_iterations,
        run_id=generate_run_id()
    )
    print(f"Total price of the research cycle: {total_price:.2f} USD")
