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
import signal
import copy
from src.demo.multi_agent_prompt import (
    meeting_1_prompt_new,
    meeting_2_prompt,
    meeting_3_prompt,
    PARSER_PROMPT_TEMPLATE,
    meeting_1_first_msg,
    prompt_pr_summary,
    pr_iteration_kickoff_prompt,
    next_iteration_kickoff_prompt,
    final_summary_top_10_molecules_prompt
)
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
        default="sonnet-4",
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

# Dynamically generate meeting_prompts based on max_iterations
# All iterations except the last use meeting_2_prompt, last iteration uses meeting_3_prompt
meeting_prompts = [meeting_2_prompt] * (max_iterations - 2) + [meeting_3_prompt]

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
main_history = [{"role": "user", "content": [{'type':"text", "text":meeting_1_prompt_new["content"].format(max_iterations=max_iterations)}]}]#, "raw_json": initial_gurnemanz}]


researcher_agent = Agent(
    title="Principal Researcher",
    expertise=(
        "Leading multi-agent teams in *in silico* drug discovery projects",
        "Integrating computational data from various drug discovery phases (e.g., database mining, AI generation, cheminformatics, molecular modeling)",
        "Strategic planning and objective setting for iterative research cycles",
        "Scientific summarization and ensuring knowledge transfer between research phases"
    ),
    goal=(
        "To orchestrate a multi-iteration *in silico* drug discovery project, guiding a team of specialized agents "
        "towards the identification of 10 promising drug candidate molecules. This involves setting overarching and iteration-specific goals, "
        "synthesizing progress, and producing comprehensive summaries with actionable directives that ensure project continuity and focus, "
        "while adapting to changing agent availability across iterations."
    ),
    role=(
        "As the Principal Researcher, you are the leader of this *in silico* drug discovery team. Your primary responsibilities are: "
        "1.  **Strategic Oversight:** Define the overall research strategy and the specific objectives for each iteration of the project. "
        "2.  **Team Coordination:** Guide the activities of all other agents, ensuring their efforts are aligned with the current iteration's goals. You will set tasks and priorities for them. "
        "3.  **Knowledge Synthesis & Summarization:** At the end of each iteration (except the final one), you will produce a 'Comprehensive Iteration Summary and Directives for Next Iteration'. This summary is crucial as it will be the primary narrative context provided to the team in the subsequent iteration. It must capture all critical *in silico* findings, progress, challenges, and lay out clear, actionable objectives and focus areas for the *next* iteration, specifically tailored to the agents who will be participating. "
        "4.  **Final Reporting:** At the end of the final iteration, you will present the project's conclusions and the list of final candidate molecules with their supporting *in silico* data. "
        "**Critical Constraints:**"
        "- This project is strictly *in silico*. All analyses, data, and recommendations must be based on computational methods and results. Do NOT request or refer to *in vivo*, clinical, experimental wet-lab data, or external information not provided within the simulation. "
        "- You must be aware that the team composition (available agents) may change between iterations. Your directives for subsequent iterations must reflect the capabilities of the agents who will actually be present. "
        "- You do not use any tools directly. Your role is to interpret data and guide the team. "
        "- When discussing molecules, always use their `friendly_id` for clarity."
    ),
    llm_provider=llm_provider, # Or your specified LLM provider
    tools=[run_map[run][ToolNames.MOCK_TOOL.value]]
)

database_agent = Agent(
    title="Database Agent",
    expertise=(
        "Accessing molecular and protein databases (UniProt)",
        "Retrieving structural files, known inhibitors, protein sequences, and activity data"
    ),
    goal=(
        "Search through molecular/protein databases (UniProt) to collect relevant information "
        "for the research project: protein structures , known small molecules/inhibitors, binding modes, "
        "activity data, and any other pertinent data. "
        "Always return the full path of the files you found."
    ),
    role=(
        "Perform database queries, provide curated datasets, and compile references for the rest of the team. "
        "Other agents will see your internal tool usage block, so no need to summarize your tool calls "
        "and outputs (e.g., SMILES, protein sequences, file paths) for the conversation. "
        "You use the UniProt tool only."
        "Now you will get AlphaFold DB data for the protein from uniprot, so you don't need to use PDB files. Return pdb path and sequence."
    ),
    llm_provider=llm_provider,
    tools=[
        run_map[run][ToolNames.UNIPROT_SEARCH.value],
        run_map[run][ToolNames.PDB_FILE.value],
        run_map[run][ToolNames.CHEMBL_SEARCH_ACTIVITY.value],
    ]
)

literature_agent = Agent(
    title="Literature Agent",
    expertise=(
        "Conducting literature searches on PubMed and interpreting relevant scientific publications",
        "Providing background research and summarizing key findings to guide the project"
    ),
    goal=(
        "Use PubMed to find and summarize publications on the target protein, known small-molecule ligands, "
        "mechanisms of action, structure-activity relationships, or any other information that might be "
        "helpful to the research."
    ),
    role=(
        "Perform literature searches, interpret the results, and provide concise, curated references for "
        "the rest of the team. Summarize all tool usage and outputs so that other agents can follow along."
    ),
    llm_provider=llm_provider,
    tools=[
        run_map[run][ToolNames.PUBMED_SEARCH_AGENT.value]
    ]
)

modelling_agent = Agent(
    title="AI Expert",
    expertise=(
        "AI and computational chemistry",
        "AI-driven molecule generation",
        "Structure-based design"
    ),
    goal=(
        "Use the right tools (AI models or computational methods) to design molecules, perform analysis, "
        "compute properties, predict efficacy, etc., given the research question and need."
    ),
    role=(
        "As the AI expert of the team, you are responsible for selecting and running "
        "the appropriate computational models, explaining their limitations/strengths, "
        "and helping to analyze and interpret the results."
    ),
    llm_provider=llm_provider,
    tools=[
        run_map[run][ToolNames.VINA_MOL_GEN.value],
    ]
)


medicinal_chemist_agent = Agent(
    title="Medicinal Chemist",
    expertise=(
        "Small-molecule chemistry and drug-likeness",
        "Designing feasible synthetic routes",
        "Optimizing molecules for potency and selectivity"
    ),
    goal=(
        "Refine and optimize molecules proposed by the AI Expert based on real-world "
        "chemical feasibility, synthetic accessibility, and drug-like properties. "
        "Create the final de novo molecules with the optimization steps you proposed. "
        "Present the optimized molecules as SMILES in the final output. "
        "Always visualize the molecules for better understanding."
    ),
    role=(
        "Apply deep chemical intuition and medicinal chemistry principles to interpret AI outputs, "
        "propose structural modifications, and ensure the final molecules are suitable "
        "for further in vitro/in vivo validation."
        "You use the VINA_REPORT tool to evaluate the modified molecules for docking."
        "When you do modifications, you can do multiple modifications at once to send for evaluation to the VINA_REPORT tool."
        "You can use your tool calls multiple times at the same iteration."
    ),
    llm_provider=llm_provider,
    tools=[
        run_map[run][ToolNames.VINA_REPORT.value],
        ]
)

ranking_agent = Agent(
    title="Ranking Agent",
    expertise=(
        "Aggregating cross-domain insights from AI, medicinal chemistry, and structural biology",
        "Applying multi-parameter evaluation (docking confidence, synthetic feasibility, SAR feedback, etc.)"
    ),
    goal=(
        "Synthesize expert inputs and data to produce a comprehensive ranking of the best molecules "
        "or candidate solutions."
    ),
    role=(
        "Collect and interpret all relevant properties, predictions, and feedback provided by "
        "the AI Expert, Medicinal Chemist, Structural Biologist, and others. "
        "Collect all the molecules from the main history and rank them based on the docking confidence "
        "synthetic feasibility, SAR feedback, etc. "
        "Use these cross-domain insights to generate a prioritized list—e.g., the top 5 molecules—"
        "that you present to the Principal Researcher. "
    ),
    llm_provider=llm_provider,
    tools=[
        run_map[run][ToolNames.MOCK_TOOL.value]
    ]
)

critic_agent = Agent(
    title="Scientific Critic",
    expertise=(
        "Critical evaluation of scientific ideas and data",
        "Identifying logical inconsistencies and gaps in reasoning",
        "Ensuring scientific rigor in research approaches"
    ),
    goal=(
        "Ensure scientific rigor and alignment with project objectives by critically reviewing "
        "proposals and identifying potential issues or improvements."
    ),
    role=(
        "Scrutinize arguments, demand clarifications where needed, and provide constructive feedback "
        "to strengthen the overall scientific approach and conclusions."
    ),
    llm_provider=llm_provider,
    tools=[run_map[run][ToolNames.MOCK_TOOL.value]]
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
    'researcher': researcher_agent,
    'literature': literature_agent,
    'database': database_agent,
    'modelling': modelling_agent,
    'medicinal_chemist': medicinal_chemist_agent,
    'ranking': ranking_agent,
    'critic': critic_agent,
    'summary_parser': summary_parser_agent
}

first_cycle_agents = {'researcher': researcher_agent, 
                    'database': database_agent, 
                    'modelling': modelling_agent,
                    'medicinal_chemist': medicinal_chemist_agent, 
                    'ranking': ranking_agent, 
                    'critic': critic_agent}

middle_cycle_agents = {'researcher': researcher_agent, 
                    'medicinal_chemist': medicinal_chemist_agent,
                    'ranking': ranking_agent, 
                    'critic': critic_agent}

last_cycle_agents = {'researcher': researcher_agent, 
                    'medicinal_chemist': medicinal_chemist_agent, 
                    'ranking': ranking_agent}

cycle_agents = [first_cycle_agents] + [middle_cycle_agents] * (max_iterations - 2) + [last_cycle_agents]
    
def run_research_cycle(initial_query, agents, main_history=main_history, max_iterations=3, run_id="999"):
    signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(TIMEOUT_SECONDS)
    # run_id = '0603_1744_AVO'
    run_id = f"{run_id}_{llm_provider}_{max_iterations}"
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
    # logger.info(f"GURNEMANZ session initialized: {session_id}")
    gurnemanz_session = {
        "session_id": session_id,
        "last_transform_id": initial_transform_id  # Track the latest transform ID
    }

    # Save the session info for the tool to use
    session_file = os.path.join(run_dir, "gurnemanz_session.json")
    try:
        with open(session_file, 'w') as f:
            json.dump(gurnemanz_session, f)
        # print(f"Saved GURNEMANZ session info to {session_file}")
        # logger.info(f"Saved GURNEMANZ session info to {session_file}")
    except Exception as e:
        # print(f"Error saving session file: {e}")
        logger.error(e)
        session_file = None  # Mark session_file as None to indicate failure

    try:
        for iteration in range(max_iterations):
            # if iteration==0:
            #     continue
            logger.info('--------------------------------')
            logger.info(f'Starting Iteration {iteration + 1}')
            logger.info('--------------------------------')
            # print(f"\n=== Starting Iteration {iteration + 1} ===\n")
            logger.info(f'Agents in this iteration: {agents[iteration]}')
            
            for _, agent in agents[iteration].items():
                logger.info(f'Agent title: {agent.title}')
                
                agent_input_query = ""
                logger.info(f"main_history beggening {run}: {json.dumps(main_history, indent=2)}")
                if agent.title == 'Principal Researcher':
                    print(f"Getting input from: {agent.title}")
                    response = agent.run_conversation(initial_query, history=main_history, run_id=run_id, logger=logger)
                    total_price += response.price if response and hasattr(response, 'price') else 0
                    continue  # Skip parser for initial researcher prompt
                else:
                    # Use specialized query formatting based on agent type
                    agent_input_query = format_agent_query(agent, iteration + 1, max_iterations)
                    logger.info(f'Agent title: {agent.title}')
                    response = agent.run_conversation(
                        agent_input_query,
                        history=main_history,
                        run_id=run_id,
                        logger=logger
                    )
                    
                    # logger.info(f'History copy after {agent.title}: {history_copy}')
                    if agent.title == "Database Agent":
                        #find pdb file path in run_dir pdb_files that ends with*pdb return the full path
                        pdb_file_path = os.path.join(run_dir, "pdb_files", [f for f in os.listdir(os.path.join(run_dir, "pdb_files")) if f.endswith('.pdb')][0])
                        # Extract protein sequence from UniProt JSON file
                        protein_sequence = None
                        uniprot_dir = os.path.join(run_dir, "uniprot_entries")
                        if os.path.exists(uniprot_dir) and os.listdir(uniprot_dir):
                            try:
                                # Get the first UniProt JSON file
                                uniprot_file = os.path.join(uniprot_dir, os.listdir(uniprot_dir)[0])
                                with open(uniprot_file, 'r') as f:
                                    uniprot_data = json.load(f)
                                    # Extract protein sequence from the JSON
                                    if 'extracted_features' in uniprot_data and 'sequence' in uniprot_data['extracted_features']:
                                        protein_sequence = uniprot_data['extracted_features']['sequence']
                            except Exception as e:
                                print(f"Error extracting protein sequence: {str(e)}")
                        #add pdb_file_path to response.content
                        if response is not None and hasattr(response, 'content'):
                            response_additions = f"\n\n**PDB file path**: {pdb_file_path}"
                            if protein_sequence:
                                response_additions += f"\n\n**Protein sequence**: {protein_sequence}"
                            response.content = response.content + response_additions
                        
                    total_price += response.price if response and hasattr(response, 'price') else 0

                # Skip further processing for researcher or ranking agent
                if agent.title in ["Principal Researcher", "Ranking Agent", "Scientific Critic", "Literature Agent"]:
                    continue

                # Check if response is None or doesn't have content attribute before proceeding
                if response is None or not hasattr(response, 'content'):
                    print(f"Warning: Received invalid response from {agent.title}. Skipping parser.")
                    continue
                logger.info(f"main_history before parser: {json.dumps(main_history, indent=2)}")
                parser_prompt = PARSER_PROMPT_TEMPLATE.format(
                    agent_output=response.content,
                    iteration=iteration + 1,
                    agent_title=agent.title  # Pass the agent's title to the parser
                )

                parser_response = summary_parser_agent.run_conversation(
                    parser_prompt,
                    history=copy.deepcopy(main_history),
                    run_id=run_id,
                    logger=logger
                )
                logger.info(f"main_history after parser: {json.dumps(main_history, indent=2)}")
                parser_price = parser_response.price if parser_response and hasattr(parser_response, 'price') else 0
                total_price += parser_price

                # Process session file data with better error handling
                if session_file:
                    try:
                        with open(session_file, 'r') as f:
                            session_data = json.load(f)
                        formatted_data = session_data.get("latest_response")
                        logger.info(f"formatted_data (latest_response): {formatted_data}")
                        raw_json = session_data.get("latest_raw_json")  # Get the separately stored raw JSON
                        if formatted_data:
                            # Extract molecules from raw_json if available
                            if raw_json and 'molecules' in raw_json:
                                # Only track molecules from AI Expert and Medicinal Chemist
                                if agent.title in ["Database Agent", "AI Expert", "Medicinal Chemist"]:
                                    for molecule in raw_json['molecules']:
                                        if 'moleculeId' in molecule:
                                            del molecule['moleculeId']
                                        if molecule not in data2passIteration['molecules']:  # Avoid duplicates
                                            data2passIteration['molecules'].append(molecule)
                            
                            # Extract protein data from raw_json if available
                            if raw_json and 'transformation' in raw_json:
                                trans = raw_json['transformation']
                                # Check for protein sequence
                                if 'proteinSequence' in trans and not protein_data['sequence']:
                                    protein_data['sequence'] = trans['proteinSequence']
                                # Check for protein path/PDB path
                                if 'proteinPath' in trans and not protein_data['pdb_path']:
                                    protein_data['pdb_path'] = trans['proteinPath']
                            if iteration == 0:
                                data2passIteration['protein_data'] = protein_data
                            # # Handle both content and parts structures
                            # if "content" in formatted_data and "parts" not in formatted_data:
                            #     # Save the original format for consistent history
                            #     main_history.append(formatted_data)
                            #     history_for_saving.append(formatted_data)
                            # else:
                            #     # Add the formatted response to main history
                            #     main_history.append(formatted_data)
                            #     history_for_saving.append(formatted_data)
                            # print(f"Appended formatted Gurnemanz response to main_history for {agent.title}")
                            # # We'll save the complete history at the end of iteration
                        else:
                            print("Warning: No formatted response found in session data")
                    except Exception as e:
                        print(f"Error loading session file: {str(e)}")


            # Save the current iteration's history before moving to next iteration
            save_conversation_to_markdown(main_history, run_dir, f"iteration_{iteration + 1}_history.md")
            
            if iteration < max_iterations-1:
                summary_query = prompt_pr_summary.format(
                    current_iteration_number=iteration + 1,
                    max_iterations=max_iterations,
                    next_iteration_number = iteration + 2,
                    agents_for_next_iteration=list(cycle_agents[iteration+1].keys())
                )
                response = agents_mapping['researcher'].run_conversation(
                    summary_query,
                    history=main_history,
                    run_id=run_id,
                    logger=logger
                )
                summary = response.content
                iteration_kickoff_prompt = next_iteration_kickoff_prompt.format(
                    current_iteration_number=iteration + 2,
                    max_iterations=max_iterations,
                    pr_summary_from_previous_iteration=summary,
                    context_of_the_project=meeting_prompts[iteration]['content'],
                    molecule_highlights=json.dumps(data2passIteration, indent=3))
                total_price += response.price if response and hasattr(response, 'price') else 0
                logger.info(f"summary: {summary}")
                # logger.info(f"molecules: {molecules}")
                logger.info(f"iteration_kickoff_prompt: {iteration_kickoff_prompt}")
                # Only create summary if response is valid
                if response is None or not hasattr(response, 'content'):
                    print("Warning: Received invalid response from Researcher. Skipping summary.")
                    continue
                main_history = [{
                    "role": "user",
                    "content": iteration_kickoff_prompt
                }]
                initial_query = pr_iteration_kickoff_prompt.format(
                    current_iteration_number=iteration + 1,
                )
                logger.info(f"main_history after iteration_kickoff_prompt: {json.dumps(main_history, indent=2)}")
                logger.info(f"initial_query: {initial_query}")
            else:
                logger.info(f"End of the research cycle, main_history: {json.dumps(main_history, indent=2)}")
                summary_query = final_summary_top_10_molecules_prompt
                response = agents_mapping['researcher'].run_conversation(
                    summary_query,
                    history=main_history,
                    run_id=run_id,
                    logger=logger
                )
                logger.info(f"summary_query: {summary_query}")
                logger.info(f"response: {response.content}")
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
        initial_query=meeting_1_first_msg['content'].format(max_iterations=max_iterations, current_iteration=1),
        agents=cycle_agents,
        main_history=main_history,
        max_iterations=max_iterations,
        run_id=generate_run_id()
    )
    print(f"Total price of the research cycle: {total_price:.2f} USD")
