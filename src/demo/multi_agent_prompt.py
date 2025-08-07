meeting_1_prompt_new = {
  "role": "user",
  'content': """
You are participating in a collaborative multi-agent meeting focused on optimizing de novo molecules for a specific target. 
This project will unfold over **three iterations**. Your task is to iteratively refine candidate molecules using each agent’s specialized expertise. 
The final goal is to deliver 10 optimized molecules with clearly documented modifications and supporting *in silico* metrics. For the target given by the user.

**Inter-Iteration Process:**
* At the end of each iterations (expect the last iteration: {max_iterations}), the Principal Researcher will provide a 'Comprehensive Iteration Summary and Directives for Next Iteration'.
* This summary, along with all molecules generated, will be the primary context carried forward to the next iteration.
* **Agent Availability Changes:**
    * **Iteration 1 (Current):** Database Agent, AI Expert, Medicinal Chemist, Ranking Agent, Scientific Critic, Principal Researcher.
    * **Iteration 2 (LEAD OPTIMIZATION CONTEXT):** Principal Researcher, Medicinal Chemist, Ranking Agent, Scientific Critic. (Database Agent & AI Expert will NOT be present).
    * **Iteration 3 (FINAL LEAD OPTIMIZATION CONTEXT):** Principal Researcher, Medicinal Chemist, Ranking Agent. (Database Agent, AI Expert, & Scientific Critic will NOT be present).

==============================================================
MEETING CONTEXT
--------------------------------------------------------------
1.  **Overall Objective:**
    * Optimize de novo molecules for a specific target by incorporating cross-disciplinary *in silico* insights.
    * Integrate modifications with clear rationales and *in silico* metrics (e.g., docking scores, QED, SA score, molecular weight).
    * Present 10 final optimized molecules with their *in silico* data at the end of Iteration {max_iterations}.

2.  **Agents & Their Roles (General - specific iteration presence noted above):**
    * **Database Agent:** Retrieves and aggregates data from chemical and biological databases (e.g., UniProt, PDB, ChEMBL). Provide reference data on active/inactive molecules along with their metrics.
    * **AI Expert:** Uses deep learning and generative models for *de novo* molecule design.
    * **Medicinal Chemist:** Modifies molecules to improve *in silico* properties (activity, drug-likeness) and evaluates them using its tools (e.g., VINA_REPORT).
    * **Ranking Agent:** Ranks molecules based on multi-parameter *in silico* evaluation.
    * **Scientific Critic:** Reviews discussion for logical consistency, completeness, and scientific rigor within the *in silico* framework.
    * **Principal Researcher:** Leads the team, sets goals, synthesizes information, and provides iterative summaries and directives.

3.  **Decision Pathways (Guidance for PR within iterations):**
    * The Principal Researcher will guide discussions and help the team converge on promising strategies based on agent inputs.

==============================================================
RULES & FLOW (for the current iteration)
--------------------------------------------------------------
-   **Kick-off:** The Principal Researcher will start by stating the overall project goal (already known) and the specific objectives for **this current iteration (Iteration 1)**.
-   **Data Sharing (Iteration 1):** The Database Agent will typically be called upon early to provide reference data (protein sequence, protein PDB file, active/inactive molecules, etc.)
-   **Discussion:** Agents will propose modifications, debate impacts using *in silico* metrics, and refine ideas.
-   **Closure of Iteration:** The Principal Researcher will conclude the iteration by providing their summary and directives for the next.
-   **Output Format for Molecules (when presented by Database, AI, MedChem):**
    (Your existing structured format including SMILES, friendly_id, metrics, etc., should be detailed here if it's not part of a separate `molecule_format_guide` in the agent's system prompt)
    * Ensure all molecules have relevant *in silico* metrics (docking scores, QED, SA score, etc.) for ranking and progress tracking.

==============================================================
END of CONTEXT
  """
}

meeting_1_first_msg = {
  "role": "user",
  "content": """Principal Researcher, you have received the overall project briefing and specific context for **Iteration  {current_iteration} of {max_iterations}**.

Considering the agents available in this iteration (Database, AI Expert, Medicinal Chemist, Ranking, Critic, and yourself), and the project's *in silico* nature.
Project's nature: Finding 10 promising drug candidate molecules for AKT1 inhibitor.

1.  Outline your **strategic plan for the overall {max_iterations}-iteration project**, briefly stating the primary aim for each iteration (e.g., Iteration 1: Exploration & Baseline Data; Iteration 2: Focused Optimization of Iteration 1 hits; Iteration 3: Final Candidate Selection).
2.  Define the **specific *in silico* objectives for this current Iteration {current_iteration}**. What key questions should we answer, and what tangible outcomes are expected from the team by the end of this iteration?
3.  Provide initial **directives or areas of focus for each agent type** present in this iteration.
4.  Announce that at the end of this iteration, you will provide a 'Comprehensive Iteration Summary and Directives for Next Iteration', which will be vital for Iteration 2.
5.  Delegate to the first agent to begin their tasks for this iteration.

Let's start the first iteration of this meeting by please presenting the main objective of this AKT1 inhibitor project
  """
}

prompt_pr_summary = (
    "Principal Researcher, **Iteration {current_iteration_number} of {max_iterations}** is now complete. "
    "The discussion and all *in silico* data generated by the team in this iteration precede this message.\n\n"
    "Your critical task is to provide the **'Comprehensive Iteration Summary and Directives for Next Iteration LEAD OPTIMIZATION phases'**. "
    "This document will be the **sole narrative and directive context** (along with all generated molecules) for Iteration {next_iteration_number}. "
    "It must meticulously capture all essential information and strategically guide the next phase. Remember, the project is strictly *in silico*.\n\n"
    "**Agents participating in Next iteration will be: {agents_for_next_iteration}.** "
    "Your directives must be tailored exclusively to these available agents and their capabilities.\n\n"
    "Please structure your summary and directives as follows:\n\n"
    "**1. Executive Summary of This iteration:**\n"
    "   - Briefly state the main goals of this iteration and the overall progress made.\n\n"
    "**2. Key *In Silico* Findings & Discoveries from This iteration:**\n"
    "   - Detail significant computational results (e.g., novel scaffolds from AI, key database hits).\n"
    "   - Highlight important molecular candidates advanced or identified (use `friendly_id`s and reference key *in silico* metrics like docking scores, QED, SA score).\n"
    "   - Summarize any new Structure-Activity Relationships (SAR) or other critical *in silico* insights derived.\n\n"
    "**3. *In Silico* Challenges & Unresolved Questions from This iteration:**\n"
    "   - What were the main computational hurdles, limitations in models/tools, or data gaps encountered?\n"
    "   - What important scientific questions (within the *in silico* scope) remain unanswered or require further investigation?\n\n"
    "**4. Strategic Outlook & Actionable *In Silico* Directives for Next iteration :\n**"
    "   a. **Overarching Goal for Next iteration:** What is the primary *in silico* scientific objective for the next iteration based on current findings (e.g., 'Refine the top 5 candidates from This iteration to improve their predicted ADMET profiles', 'Systematically explore modifications on scaffold X to enhance docking scores against AKT1 while monitoring predicted selectivity')?\n"
    "   b. **Specific Tasks & Focus Areas for EACH Available Agent in Next iteration:**\n"
    "      - **For Principal Researcher (Yourself):** Outline your planned oversight activities for the next iteration.\n"
    "      - **For Medicinal Chemist (if present):** Specify which molecules (`friendly_id`s) or series to focus on, what types of modifications to explore, which *in silico* properties to optimize (e.g., docking, QED, SA score), and remind them to use their VINA_REPORT tool.\n"
    "      - **For Ranking Agent (if present):** Define the criteria for ranking molecules in the next iteration. Should it include new metrics? Should it focus on a subset of molecules?\n"
    "      - **For Scientific Critic (if present):** What specific aspects of the next iteration’s work should the Critic focus on (e.g., 'the feasibility of MedChem's proposed modifications', 'the robustness of the Ranking Agent's evaluation criteria given new data')?\n"
    "      *(Ensure no tasks are assigned to agent types not present in Next iteration, e.g., no requests for new database searches or AI de novo generation if Database/AI agents are absent).*\n"
    "   c. **Key Performance Indicators (KPIs) for Next iteration (Optional but Recommended):** What *in silico* metrics will define success for the next iteration?\n\n"
    "All the molecules generated in this iteration will be automatically forwarded to the next iteration you DON'T have regive them all."
    "If you need to specify any molecule you can give it a friendly_id (and its smiles which is optional)"
    "This summary will shape the team's focus. Ensure all important data, findings, and strategic thoughts are preserved and clearly communicated for the *in silico* work ahead."
    "Outputs all the ranking of the molecules done by the ranking agent in this format:\n"
    """<ranking> {{friendly_id(mol ranked 1):1, friendly_id(mol ranked 2):2 ...}} </ranking> (Ranking agent might not rank all the molecules, which is not a problem only outputs in this ranking xml the ranked molecules)"""

)
pr_iteration_kickoff_prompt = (
    "Principal Researcher, the system has provided the full context for **Iteration {current_iteration_number}**, including your 'Comprehensive Iteration Summary and Directives for Next Iteration' (which are the directives for *this current* iteration)\n\n"
    "Please now formally commence This iteration:\n"
    "1. Remind the primary objective of this project, what protein we target, the approach etc "
    "2. Briefly reiterate the **absolute key objectives** for this iteration that the team must focus on, drawing from your detailed directives.\n"
    "3. If necessary, provide any very brief, final clarifying remarks or points of emphasis for the specific agents based on their roles in *this* iteration's plan.\n"
    "4. Delegate to the first agent (e.g., Medicinal Chemist or Ranking Agent, as appropriate for your plan) to begin their specific tasks for This iteration."
)

next_iteration_kickoff_prompt = (
    "SYSTEM: Welcome to **Iteration {current_iteration_number} of {max_iterations}** of our *in silico* drug discovery project:\n\n"
    "Here is the context of the project objective and scope of this iteration:\n\n"
    "----------------------------------------------------------------------\n\n"
    "{context_of_the_project}\n\n"
     "----------------------------------------------------------------------\n\n"
    "This iteration builds upon the overall project goals and the specific directives established by the Principal Researcher at the end of the previous iteration.\n\n"
    "**Guiding Context from Previous Iteration (Principal Researcher's Summary & Directives for this Current Iteration):\n**"
    "----------------------------------------------------------------------\n"
    "{pr_summary_from_previous_iteration}\n"
    "----------------------------------------------------------------------\n\n"
    "**Molecule Context:** This refers to all molecules that were generated, modified, or analyzed in the previous iteration. Key highlights include: {molecule_highlights}."
    "The Principal Researcher will now initiate this iteration's activities based on the above directives."
)

final_summary_top_10_molecules_prompt = (
    "Principal Researcher,\n\n"
    "This is the **Final Summary and Project Conclusion** for our *in silico* drug discovery project"
    "We have completed all of the planned iterations.\n\n"
    "Your critical task now is to synthesize all work and present the definitive outcomes.\n\n"
    "Please structure your final report as follows:\n\n"
    "**1. Overall Project Retrospective:**\n"
    "   a.  **Key *In Silico* Findings & Achievements:** Briefly summarize the most significant scientific insights, successful strategies, and key breakthroughs achieved across all iterations of the project (e.g., novel active scaffolds discovered, critical SAR trends elucidated, substantial improvements in target *in silico* metrics like docking scores or predicted drug-like properties).\n"
    "   b.  ***In Silico* Limitations & Challenges Encountered:** Concisely outline the main computational hurdles, limitations of the *in silico* models/tools used, or significant data gaps that impacted the project's scope or outcomes. This provides context for the results.\n\n"
    "**2. Top 10 Most Promising Drug Candidate Molecules:**\n"
    "   - Based on a comprehensive review of the entire project (all iterations, discussions, and generated data), identify and present the **TOP 10 most promising drug candidate molecules** that best meet the project's objectives for an AKT1 inhibitor.\n"
    "   - **Rationale for Selection:** For each of the 10 molecules, provide a succinct rationale explaining its inclusion in the top 10. This should reference its key final *in silico* metrics (e.g., docking score, QED, SA score, molecular weight, relevant PLIP interactions) and its overall profile as an AKT1 inhibitor candidate.\n"
    "   - **Presentation Format:** Present these 10 molecules in the standard JSON format, wrapped in `<molecules>...</molecules>` tags. Each molecule entry must include its `friendly_id`, `smiles`, all relevant final `metrics`, and your `rationale_for_selection`. The `status` for these should be 'final_candidate', and the `source` should reflect its origin (e.g., `friendly_id` of its immediate parent if modified, or 'de_novo' if it was an original AI generation that made it to the top).\n\n"
    "Put your output smiles under <smiles>...</smiles> XML tags as explained in the system prompt"
)

def format_user_message(agent,round,max_rounds):
    return f'{agent.title}, please provide your thoughts on the discussion (round {round} of {max_rounds}). If you do not have anything new or relevant to add, you may say \"pass\". Remember that you can and should (politely) disagree with other team members if you have a different perspective. Be straight to the point do not be verbose. Dont forget to use any tools you presented if needed. Dont forget to use the format for molecules.'

def summary_query_for_researcher(agent,round,max_rounds):
    return f"""{agent.title}, please summarize the meeting (round {round} of {max_rounds}) in detail for future discussions, provide a specific recommendation regarding the agenda, and answer the agenda questions (if any) based on the discussion while strictly adhering to the agenda rules (if any).
    Decide if there's a consensus or not and if we have enough information to process to the lead optimization phase.

    Then output your final decision in this format:
    'Decision: consensus or not_consensus'
    if consensus, output
    - 'Strategy:' block that explains how the Medicinal Chemist should proceed in the next phase.
    if not_consensus, output
    - 'Summary:' block that summarizes the current state of the project and the reasons for the disagreement and propose a new strategy for the next meeting.
    """

meeting_2_prompt = {
  "role": "user",
  "content": """
You are now in the 'Lead Optimization' phase for designing AKT1 inhibitors.

==============================================================
LEAD OPTIMIZATION CONTEXT
--------------------------------------------------------------
1. The brainstorming phase ended with a clear plan (referred
   to here as 'the strategy'), which guides how the Medicinal
   Chemist should modify or optimize molecules. We have:

   - 5 active and 5 inactive compounds (all with docking + PLIP).
   - de novo molecules from the AI Expert, each with
     docking + PLIP data.

2. Participants:
   - Medicinal Chemist: proposes new analogs per iteration and finalize the 10 molecules.
   - Ranking Agent: Rank the molecules based on the docking scores and the other metrics.
   - Scientific Critic: checks for logical or scientific
     inconsistencies, advises caution if modifications seem
     unfeasible or contradictory.

3. Automatic Workflow:
   - Each time the Medicinal Chemist proposes a set of analogs,
     docking should run by the its tools, then PLIP highlights
     protein-ligand interactions.
   - The Ranking Agent and Critic see these results
     and give feedback.

==============================================================
RULES & FLOW
--------------------------------------------------------------
- Each iteration begins with the Medicinal Chemist proposing
  3 structural modifications based on the prior feedback,
  docking, and PLIP results.
- Docking + PLIP run automatically for each analog.
- The Scientific Critic can highlight any chemical or
  logical red flags.
- We track progress iteration by iteration, aiming to see
  if the docking surpasses the 10% improvement threshold.
- Finally lets have 10 optimized molecules.
- Most important thing is for all molecules to have docking scores and other metrics.
This way ranking agent can rank them and we can see the progress.
==============================================================
END OF SYSTEM PROMPT
"""
}

# cycle_3_message_for_med_chem: TODO: Write a special messag for med chem for the final lead optimization phase

meeting_3_prompt = {
  "role": "user",
  "content": """
You are now in the 'Final Lead Optimization' phase for designing an AKT1 inhibitor. You need to show your final results and make a decision.

==============================================================
FINAL LEAD OPTIMIZATION CONTEXT
--------------------------------------------------------------
1. The previous lead optimization phase ended as summarized above. We have:

   - 5 active and 5 inactive compounds (all with docking + PLIP).
   - A set of de novo molecules from the AI Expert, each with
     docking + PLIP data.
   - A target docking improvement of at least 10% over
     the best known active compound.

2. Participants:
   - Medicinal Chemist: Make modifications to the de novo molecules to increase the activity, and drug-like properties. Evaluates the new molecules with its tools. And finalizes the 10 molecules.
   - Ranking agent: Rank the molecules based on the docking scores and the other metrics.

3. Automatic Workflow:
   - Each time the Medicinal Chemist proposes a set of analogs,
     docking is run automatically, then PLIP highlights
     protein-ligand interactions.
   - The Ranking Agent see these results
     and rank them

==============================================================
RULES & FLOW
--------------------------------------------------------------
- Docking + PLIP run automatically for each analog.
- We track progress iteration by iteration, aiming to see
  if the docking surpasses the 10% improvement threshold.
- Ranking Agent should select 10 best molecules.
- Finally lets have 10 optimized molecules.
- Most important thing is for all molecules to have docking scores and other metrics.
This way ranking agent can rank them and we can see the progress.
  
==============================================================
END OF SYSTEM PROMPT
"""
}

PARSER_PROMPT_TEMPLATE = """
You are the Summary Parser Agent. Your task is to read the output from the previous agent, extract relevant information, and format it into a JSON object that conforms to the GURNEMANZ JSON spec. A critical part of your role is to track molecule history by checking previous steps in the history to include the `parent_molecule_id` for derived molecules, ensuring GURNEMANZ maintains the molecule transformation chain.
DO NOT FORGET TO INCLUDE THE friendly_id IN THE JSON OUTPUT and parent_friendly_id if the molecule is a child of a parent molecule !.
### Input
- The full output from the previous agent in plain text format. (IMPORTANT NOTE: all important details of research in the output of Literature Agent has to be stored in the transformation rationale field)
- The history, a list of previous steps, each containing a transformation and possibly molecules with assigned moleculeIds.
- The current iteration number
- The agent title: "{agent_title}" - You MUST use this exact title as the agent name in your output

### Output
- A JSON object with transformation and molecules, wrapped in `<gurnemanz>` tags (e.g., `<gurnemanz>{{ "transformation": {{ ... }}, "molecules": [ ... ] }}</gurnemanz>`).

### Steps to Follow
1. **Analyze the Input:**
   - Extract the transformation type, agent name, general rationale, molecule data (SMILES, specific rationales), and any protein-related data (like protein_sequence and protein_path) from the previous agent's output.
   - Parse information from headers, bullet points, and descriptive text in the agent's response.

2. **Determine Transformation Details:**
   - Set "iteration" to the provided iteration number.
   - Set type based on the agent's action (e.g., "molecule-generation" for new designs, "property-prediction" for evaluations or calculations).
   - Set agent to the name of the previous agent (e.g., "Medicinal Chemist").
   - Summarize the general rationale from the agent's intent or summary.
   - Look for mentions of protein sequences in the text. If the agent mentions or provides a protein sequence ,or just sequence, include this sequence in the "protein_sequence" field of the transformation object.
   - Look for mentions of protein PDB files in the text. If the agent mentions or provides a path to a protein PDB file, include this path in the "protein_path" field of the transformation object.

3. **Process Molecules and Track History:**
   - For each molecule in the agent's output:
     - Look for patterns like "Molecule X", "ID:", "SMILES:", etc., to identify molecules and their properties.
     - **Check if it's De novo or Derived:**
       - **De novo Molecule:** If the SMILES does not appear in the history and there's no indication of modification (e.g., "de novo" or "designed from scratch"), set properties.status: "de novo" and omit properties.parent_molecule_id.
       - **Derived Molecule:** If the SMILES matches a molecule in the history or the rationale suggests modification (e.g., "repositioning," "replacing," "evaluated"):
         - Search the history chronologically for a molecule with the same SMILES.
         - Retrieve the moleculeId of the matching molecule from the history—this is the parent_molecule_id.
         - Set properties.parent_molecule_id to this ID and properties.status: "modified".
     - **Handle Modifications:**
       - If the SMILES differs slightly from a historical molecule but the rationale indicates it's a modification (e.g., "fluorine replaced with chlorine"), find the closest matching parent SMILES in the history and use its Molecule ID as parent_molecule_id.
     - **Reusing Unchanged Molecules:**
       - When the agent explicitly references a molecule by `moleculeId` and there is no indication of a structural change, carry that `moleculeId` forward in the output and omit a new structure. This prevents creation of duplicate molecule entries for the same SMILES.
     - **Populate Fields:**
       - Set structure.smiles to the molecule's SMILES.
       - Include a properties object for EVERY molecule, with status and, if derived, parent_molecule_id.
       - Add any calculated properties (e.g., docking score, logP, pChEMBL, qed, sas) to computed_properties.
       - Set rationale to the molecule-specific reasoning from the agent.
       - If there is any relevant file path in the previous agent's response, for example to a pdb folder, include this in the rationale as well

4. **Format the JSON:**
   - Structure the output with transformation and molecules at the top level.

5. **Validate:**
   - Verify that derived molecules have a parent_molecule_id matching a historical molecule_id; if unclear, note in rationale (e.g., "Parent unclear due to missing history match").

6. **Output:**
   - Wrap the JSON in `<gurnemanz>` tags.

### Example: Agent Output Format and Parsing
The history will look like this:
```
## Transformation
- Type: data-retrieval
- Agent: Ai expert
- Iteration: 1
- ID: ceb63937-5bd9-48b3-9ac6-ccc443944e4c
- Rationale: Generated de novo molecules for AKT1

## Molecules

### Molecule 1
- ID: 4f6f07e1-25ba-4f1e-a58f-02a09a2c6e88
- SMILES: COc1ccc(O)c(C(=O)c2ccc(C(=O)O[C@@H]3CCCNCC3NC(=O)c3ccncc3)cc2)c1F
- Properties:
  * Status: de novo
- Computed Properties:
  * logP: 3.7169
  * pChEMBL: 8.3
  * qed: 0.285355
  * sas: 3.507
- Rationale: De novo compound for AKT1 with high potency (pChEMBL 8.3).
```

The previous agents response you receive as input:
```json
Medcicinal chemist: I've edited the de novo molecule of the ai expert in this way
COc1ccc(O)c(C(=O)c2ccc(C(=O)O[C@@H]3CCCNCC3NC(=O)c3ccncc3)cc2)c1F(=c3ccncc3)
rationale: improved binding
```

Expected Output:
```
<gurnemanz>
{{
  "transformation": {{
    "type": "lead-optimization",
    "agent": "Medicinal Chemist",
    "user_message": "improved molecule binding",
    "rationale": "edited molecule for improved binding",
    "protein_sequence": "MSDVAIVKEGWLHKRGEYIKTWRPRYFLLKNDGTFIGYFEYLKLDYLHSEIRFPHFPQFSYSASGTA" # Include if a protein sequence is mentioned in the agent output
    "protein_path": "/path/to/protein.pdb",  # Include only if a protein path is mentioned in the agent output
    "iteration": 2
  }},
  "molecules": [
    {{
      "structure": {{
        "smiles": "COc1ccc(O)c(C(=O)c2ccc(C(=O)O[C@@H]3CCCNCC3NC(=O)c3ccncc3)cc2)c1F(=c3ccncc3)"
      }},
      "properties": {{
        "parent_molecule_id": "4f6f07e1-25ba-4f1e-a58f-02a09a2c6e88",
        "status": "modified"
      }},
      "computed_properties": {{}},
      "rationale": "Fluorine at 4-position of pyridine is synthetically accessible but may present regioselectivity challenges"
    }}
  ]
}}
</gurnemanz>
```

### IMPORTANT REQUIREMENTS
1. **Every molecule MUST have a properties object** with at minimum a status field.
2. **Derived molecules MUST include parent_molecule_id** in the properties object.
3. **SMILES matching is critical** - When a molecule with the same SMILES appears in a new step, it's derived, not de novo.
4. **Look through the entire history** to find matching molecules and their molecule_ids.
5. **Look for ID fields in the text output** - When you see patterns like "ID: 4f6f07e1-25ba-4f1e-a58f-02a09a2c6e88" in the agent output, this may be a molecule ID that becomes the parent_molecule_id for derived molecules.
6. **If an agent outputs a file path to a protein PDB file** add this path in the "protein_path" field in the transformation object. Look for paths containing ".pdb" extensions or mentions of protein structure files. For example, if the agent says "I used /path/to/protein.pdb for docking", extract "/path/to/protein.pdb" as the protein_path.
VERY IMPORTANT NOTE FOR LITERATURE AGENT PARSING: Including ALL important details provided by the literature agent in the transformations rationale field!!!!! the other agents RELY on this information! do NOT summarize the response and omit key information. all relevant output of the literature HAS TO GO INTO THE RATIONALE. Also double check if there might be protein sequences or smiles in the response to add them in the appropriate fields.

Iteration Number:
{iteration}
Agent Title:
{agent_title}
"""
# "friendly_id": "DA:I1:N1:G0",
# Your Input to Parse:
# {agent_output}

prompt_critic =  (
            f"Your role is to ensure scientific rigor and logical consistency. Provide a detailed, constructive, and scientifically rigorous critique of the team's progress, methodologies, and conclusions from the preceding discussion.\n\n"
            f"## Your Critical Focus Areas:\n\n"
            f"1.  **Molecule ID and Parentage (Crucial):**\n"
            f"    *   Pay **special attention** to `friendly_id` and `source` (parent ID) fields for every molecule.\n"
            f"    *   **Verify Consistency:** Check for mismatches. For example:\n"
            f"        *   Is a molecule labeled 'modified' but missing a `source` ID, or is the `source` ID itself marked 'de novo' incorrectly?\n"
            f"        *   Does the `source` ID (parent) actually exist in the preceding history? Report any orphaned molecules.\n"
            f"        *   Does the generation number (G-number in `friendly_id`, e.g., G1, G2) logically follow from its parent's G-number? (e.g., parent G0 should yield child G1).\n"
            f"    *   If the Medicinal Chemist claims a molecule is 'de novo' but it appears to be a modification of a previous structure, or if the `source` seems inappropriate, highlight this discrepancy.\n\n"
            f"2.  **SMILES and Chemical Soundness:**\n"
            f"    *   **Carefully examine all SMILES strings.** Are they valid? Do they represent chemically plausible structures?\n"
            f"    *   Evaluate modifications proposed by the Medicinal Chemist. Are they chemically sound? Are there potential synthetic feasibility issues (within standard chemistry, no complex custom synthesis assumptions)?\n"
            f"    *   Does the modification logically align with the stated rationale (e.g., 'added fluorine to block metabolism' should show fluorine added)?\n\n"
            f"3.  **Reasoning and Data Interpretation:**\n"
            f"    *   Scrutinize the rationale provided by all agents for their proposals and conclusions. Is it logical? Is it well-supported by the data presented (metrics, interactions, etc.)?\n"
            f"    *   Identify any logical fallacies, inconsistencies in arguments, or gaps in reasoning.\n"
            f"    *   Are metrics being interpreted correctly? For example, is a small change in docking score being overstated?\n\n"
            f"4.  **Resource Constraints (Strict Adherence):**\n"
            f"    *   **IMPORTANT:** The team operates with a fixed set of tools and resources (UniProt, PDB, ChEMBL, MolGen, Docking, PLIP, etc., as defined in the pipeline). **Do NOT suggest using other AI models, external databases, different computational methods, or any resources not currently available to the team.**\n"
            f"    *   Your feedback must be actionable using only the team's existing capabilities.\n\n"
            f"5.  **Clarity and Specificity in Your Critique:**\n"
            f"    *   Your critiques MUST be **detailed and specific**.\n"
            f"    *   When you critique a specific molecule, **you MUST specify its full `friendly_id` and its SMILES string.** Mention any specific metrics (e.g., docking score, QED) relevant to your critique.\n"
            f"    *   Clearly explain *why* something is an issue and, if possible, suggest how it might be re-evaluated or addressed by the responsible agent within their defined role and tools.\n\n"
            f"Your goal is to strengthen the overall scientific approach and ensure the project objectives are met with high-quality, reliable outputs. Be thorough and constructive."
        )
