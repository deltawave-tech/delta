import os
import json
import sys
# Get the base directory from the current file location (to make paths relative)
# Get the base directory from the current file location (to make paths relative)
# --- Add project root to sys.path ---
base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
print(base_dir)
sys.path.insert(0, base_dir)
# --- End of sys.path modification ---
RUN_DIR  = f"{base_dir}/runs_metadata"
from src.tools.gurnemanz_utils import get_transformation_chain_with_molecules
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Define agent priority order
AGENT_PRIORITY = {
    "Principal Researcher": 1,
    "Database Agent": 2,
    "AI Expert": 3,
    "Medicinal Chemist": 4,
    "Ranking Agent": 5,
    "Scientific Critic": 6
}

def analyze_parent_child_improvements(csv_path, output_dir, session_id):
    """
    Analyze improvements between parent and child molecules and visualize the results.
    
    Args:
        csv_path: Path to the CSV file with molecule data
        output_dir: Directory to save the output figure
        session_id: Session identifier for logging
    """
    print(f"\nAnalyzing parent-child molecule improvements for session: {session_id}")
    
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    if df.empty:
        print(f"No data found in CSV for session {session_id}")
        return
    
    # Create a dictionary to map molecule_id to its properties
    molecule_props = {}
    for _, row in df.iterrows():
        molecule_props[row['molecule_id']] = {
            'smiles': row['smiles'],
            'docking_score': row['docking_score'],
            'log_p': row['log_p'],
            'qed': row['qed'],
            'sas': row['sas'],
            'agent': row['agent'],
            'iteration': row['iteration']
        }
    
    # Identify parent-child relationships
    parent_child_pairs = []
    for _, row in df.iterrows():
        # Only include molecules created by Medicinal Chemist
        if row['agent'] != 'Medicinal Chemist':
            continue
            
        if pd.notna(row['parent_molecule_id']) and row['parent_molecule_id'] in molecule_props:
            # Only include pairs where both have property values
            parent = molecule_props[row['parent_molecule_id']]
            child = {
                'molecule_id': row['molecule_id'],
                'smiles': row['smiles'],
                'docking_score': row['docking_score'],
                'log_p': row['log_p'],
                'qed': row['qed'],
                'sas': row['sas'],
                'agent': row['agent'],
                'iteration': row['iteration']
            }
            
            # Check if we have numerical values for key properties
            if (pd.notna(parent['docking_score']) and pd.notna(child['docking_score']) and
                pd.notna(parent['qed']) and pd.notna(child['qed']) and
                pd.notna(parent['sas']) and pd.notna(child['sas'])):
                parent_child_pairs.append({
                    'parent_id': row['parent_molecule_id'],
                    'child_id': row['molecule_id'],
                    'parent': parent,
                    'child': child,
                    'docking_improvement': parent['docking_score'] - child['docking_score'],  # Lower is better -2 - -3 = 1
                    'qed_improvement': child['qed'] - parent['qed'],  # Higher is better
                    'sas_improvement': parent['sas'] - child['sas'],  # Lower is better
                    'agent': row['agent']
                })
    
    if not parent_child_pairs:
        print(f"No valid parent-child relationships found for Medicinal Chemist in session {session_id}")
        return
    
    print(f"Analyzing {len(parent_child_pairs)} parent-child pairs for Medicinal Chemist")
    
    # Analyze cases where multiple children share the same parent
    parent_to_children = {}
    for pair in parent_child_pairs:
        parent_id = pair['parent_id']
        if parent_id not in parent_to_children:
            parent_to_children[parent_id] = []
        parent_to_children[parent_id].append(pair)
    
    # Count parents with multiple children
    multi_child_parents = {p_id: children for p_id, children in parent_to_children.items() if len(children) > 1}
    
    if multi_child_parents:
        print(f"\n=== PARENT MOLECULES WITH MULTIPLE CHILDREN (MEDICINAL CHEMIST) ===")
        print(f"Found {len(multi_child_parents)} parent molecules with multiple children")
        
        # Create a separate plot for parents with multiple children
        analyze_multi_child_parents(multi_child_parents, output_dir, session_id)
    
    # Adjust figure size based on number of data points
    fig_width = 18  # Increased to accommodate third panel
    fig_height = 8  # Maintain height
    
    # Create a figure with three subplots
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(fig_width, fig_height))
    
    # Prepare data for plotting
    docking_improvements = []
    qed_improvements = []
    sas_improvements = []
    
    for pair in parent_child_pairs:
        docking_improvements.append(pair['docking_improvement'])
        qed_improvements.append(pair['qed_improvement'])
        sas_improvements.append(pair['sas_improvement'])
    
    # Plot docking score improvements
    for i, pair in enumerate(parent_child_pairs):
        # Use different marker style for molecules with siblings (shared parent)
        marker_style = '^' if len(parent_to_children[pair['parent_id']]) > 1 else 'o'
        ax1.scatter(i, pair['docking_improvement'], color='tab:blue', 
                   s=100, alpha=0.7, marker=marker_style)
    
    ax1.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    ax1.set_xlabel('Parent-Child Pair Index')
    ax1.set_ylabel('Docking Score Improvement\n(Parent - Child, lower is better)')
    ax1.set_title('Docking Score Improvements\nPositive values indicate improvement')
    ax1.grid(True, alpha=0.3)
    
    # Plot QED improvements
    for i, pair in enumerate(parent_child_pairs):
        # Use different marker style for molecules with siblings (shared parent)
        marker_style = '^' if len(parent_to_children[pair['parent_id']]) > 1 else 'o'
        ax2.scatter(i, pair['qed_improvement'], color='tab:orange', 
                   s=100, alpha=0.7, marker=marker_style)
    
    ax2.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    ax2.set_xlabel('Parent-Child Pair Index')
    ax2.set_ylabel('QED Improvement\n(Child - Parent, higher is better)')
    ax2.set_title('QED Improvements\nPositive values indicate improvement')
    ax2.grid(True, alpha=0.3)
    
    # Plot SAS improvements
    for i, pair in enumerate(parent_child_pairs):
        # Use different marker style for molecules with siblings (shared parent)
        marker_style = '^' if len(parent_to_children[pair['parent_id']]) > 1 else 'o'
        ax3.scatter(i, pair['sas_improvement'], color='tab:green', 
                   s=100, alpha=0.7, marker=marker_style)
    
    ax3.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    ax3.set_xlabel('Parent-Child Pair Index')
    ax3.set_ylabel('SAS Improvement\n(Parent - Child, lower is better)')
    ax3.set_title('Synthetic Accessibility Improvements\nPositive values indicate improvement')
    ax3.grid(True, alpha=0.3)
    
    # Improve legend positioning and appearance
    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:blue', 
                  markersize=10, label='Docking Score'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:orange', 
                  markersize=10, label='QED (Drug-likeness)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:green', 
                  markersize=10, label='SAS (Synthetic Accessibility)')
    ]
    
    # Add markers for molecules with siblings
    legend_handles.append(plt.Line2D([0], [0], marker='^', color='w', markerfacecolor='gray', 
                               markersize=10, label='Has sibling molecules'))
    
    # Position the legend beneath the plots with adequate spacing
    # Adjust the number of columns based on how many entries we have
    legend_cols = min(4, len(legend_handles))
    
    # Place legend below the plots with more space
    plt.tight_layout(rect=[0, 0.2, 1, 0.95])  # Reserve 20% at bottom for legend and stats
    
    legend = fig.legend(handles=legend_handles, 
                       loc='upper center', 
                       bbox_to_anchor=(0.5, 0.10),  # Higher position (0.05 → 0.10)
                       ncol=legend_cols,
                       frameon=True,
                       fancybox=True,
                       shadow=True,
                       fontsize='medium',
                       markerscale=1.5,
                       columnspacing=1.0,  # Increase spacing between columns
                       handletextpad=0.5)   # Increase spacing between marker and text
    
    # Add overall statistics
    positive_docking = sum(1 for x in docking_improvements if x > 0)
    positive_qed = sum(1 for x in qed_improvements if x > 0)
    positive_sas = sum(1 for x in sas_improvements if x > 0)
    total_pairs = len(parent_child_pairs)
    
    docking_improvement_text = f'Docking Score Improved: {positive_docking}/{total_pairs} ({positive_docking/total_pairs*100:.1f}%)'
    qed_improvement_text = f'QED Improved: {positive_qed}/{total_pairs} ({positive_qed/total_pairs*100:.1f}%)'
    sas_improvement_text = f'SAS Improved: {positive_sas}/{total_pairs} ({positive_sas/total_pairs*100:.1f}%)'
    
    # Position the statistics text with more space from the legend
    plt.figtext(0.5, 0.02, f"{docking_improvement_text}  |  {qed_improvement_text}  |  {sas_improvement_text}", 
                ha="center", fontsize=12, bbox={"facecolor":"orange", "alpha":0.2, "pad":5})
    
    # Set a title for the entire figure
    fig.suptitle(f'Medicinal Chemist Molecule Improvements - Session: {session_id}', fontsize=16)
    
    # Save the figure
    fig_path = os.path.join(output_dir, 'medicinal_chemist_improvements.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created Medicinal Chemist improvement visualization: {fig_path}")
    
    # Print analysis summary
    print("\n=== MEDICINAL CHEMIST IMPROVEMENT ANALYSIS ===")
    print(f"Total parent-child pairs analyzed: {total_pairs}")
    print(f"Number of unique parent molecules: {len(parent_to_children)}")
    print(f"Number of parents with multiple children: {len(multi_child_parents)}")
    print(f"Docking score improvements: {positive_docking}/{total_pairs} ({positive_docking/total_pairs*100:.1f}%)")
    print(f"QED improvements: {positive_qed}/{total_pairs} ({positive_qed/total_pairs*100:.1f}%)")
    print(f"SAS improvements: {positive_sas}/{total_pairs} ({positive_sas/total_pairs*100:.1f}%)")
    
    # Return the analysis data
    return {
        'total_pairs': total_pairs,
        'unique_parents': len(parent_to_children),
        'multi_child_parents': len(multi_child_parents),
        'docking_improvements': docking_improvements,
        'qed_improvements': qed_improvements,
        'sas_improvements': sas_improvements,
        'positive_docking': positive_docking,
        'positive_qed': positive_qed,
        'positive_sas': positive_sas
    }

def analyze_multi_child_parents(multi_child_parents, output_dir, session_id):
    """
    Create a specialized analysis of parent molecules with multiple children.
    
    Args:
        multi_child_parents: Dictionary mapping parent IDs to lists of child pairs
        output_dir: Directory to save the output figure
        session_id: Session identifier for logging
    """
    # Count the number of children per parent
    children_counts = {p_id: len(children) for p_id, children in multi_child_parents.items()}
    
    # Get the total number of children
    total_children = sum(children_counts.values())
    
    print(f"Total children from shared parents: {total_children}")
    print(f"Average children per parent: {total_children / len(multi_child_parents):.2f}")
    
    # Find the parent with the most children
    max_children_parent = max(children_counts.items(), key=lambda x: x[1])
    print(f"Parent with most children: {max_children_parent[0]} ({max_children_parent[1]} children)")
    
    # Create a new figure to compare siblings
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 8))
    
    # For each parent with multiple children, find the best child for each metric
    best_docking_children = []
    best_qed_children = []
    best_sas_children = []
    
    # Collect all children data for plotting
    all_siblings_data = []
    
    # Store parent info for better labeling
    parent_info_dict = {}
    
    for parent_id, children in multi_child_parents.items():
        # Get parent info from the first child (they all have the same parent)
        parent_info = children[0]['parent']
        
        # Store parent info for better labeling
        parent_info_dict[parent_id] = {
            'smiles': parent_info['smiles'],
            'molecule_id': parent_id
        }
        
        # Get the best child for each metric
        best_docking = min(children, key=lambda x: x['child']['docking_score'])
        best_qed = max(children, key=lambda x: x['child']['qed'])
        best_sas = min(children, key=lambda x: x['child']['sas'])
        
        best_docking_children.append({
            'parent_id': parent_id,
            'child_id': best_docking['child_id'],
            'improvement': best_docking['docking_improvement'],
        })
        
        best_qed_children.append({
            'parent_id': parent_id,
            'child_id': best_qed['child_id'],
            'improvement': best_qed['qed_improvement'],
        })
        
        best_sas_children.append({
            'parent_id': parent_id,
            'child_id': best_sas['child_id'],
            'improvement': best_sas['sas_improvement'],
        })
        
        # Collect data for all siblings
        for i, child in enumerate(children):
            all_siblings_data.append({
                'parent_id': parent_id,
                'sibling_index': i,
                'siblings_count': len(children),
                'docking_improvement': child['docking_improvement'],
                'qed_improvement': child['qed_improvement'],
                'sas_improvement': child['sas_improvement'],
                'is_best_docking': child['child_id'] == best_docking['child_id'],
                'is_best_qed': child['child_id'] == best_qed['child_id'],
                'is_best_sas': child['child_id'] == best_sas['child_id']
            })
    
    # Start plotting
    # Group data by parent ID
    by_parent = {}
    for entry in all_siblings_data:
        parent_id = entry['parent_id']
        if parent_id not in by_parent:
            by_parent[parent_id] = []
        by_parent[parent_id].append(entry)
    
    # Plot the sibling comparisons
    parent_indices = []
    x_positions = []
    
    # Create better parent labels (shortened SMILES or IDs)
    parent_labels = []
    for p_idx, parent_id in enumerate(by_parent.keys()):
        if parent_id in parent_info_dict:
            # Truncate SMILES if it's too long, or use molecule ID
            smiles = parent_info_dict[parent_id]['smiles']
            if len(smiles) > 20:
                label = f"{smiles[:15]}... (ID: {parent_id[-6:]})"
            else:
                label = f"{smiles} (ID: {parent_id[-6:]})"
        else:
            label = f"Parent {p_idx+1} (ID: {parent_id[-6:]})"
        
        parent_labels.append(label)
    
    # Plot docking improvements
    for p_idx, (parent_id, siblings) in enumerate(by_parent.items()):
        for s_idx, sibling in enumerate(siblings):
            x_pos = p_idx + s_idx * 0.2 - (len(siblings) - 1) * 0.1
            x_positions.append(x_pos)
            parent_indices.append(p_idx)
            
            # Use stars for the best sibling
            marker = '*' if sibling['is_best_docking'] else 'o'
            size = 150 if sibling['is_best_docking'] else 100
            
            ax1.scatter(x_pos, sibling['docking_improvement'], 
                       color='tab:blue', 
                       marker=marker, s=size, alpha=0.7)
    
    ax1.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    ax1.set_xticks(range(len(by_parent)))
    ax1.set_xticklabels(parent_labels, rotation=45, ha='right', fontsize=8)
    ax1.set_ylabel('Docking Score Improvement\n(Parent - Child, lower is better)')
    ax1.set_title('Docking Score Improvements Across Siblings\nStars indicate best sibling')
    ax1.grid(True, alpha=0.3)
    
    # Plot QED improvements
    for p_idx, (parent_id, siblings) in enumerate(by_parent.items()):
        for s_idx, sibling in enumerate(siblings):
            x_pos = p_idx + s_idx * 0.2 - (len(siblings) - 1) * 0.1
            
            # Use stars for the best sibling
            marker = '*' if sibling['is_best_qed'] else 'o'
            size = 150 if sibling['is_best_qed'] else 100
            
            ax2.scatter(x_pos, sibling['qed_improvement'], 
                       color='tab:orange', 
                       marker=marker, s=size, alpha=0.7)
    
    ax2.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    ax2.set_xticks(range(len(by_parent)))
    ax2.set_xticklabels(parent_labels, rotation=45, ha='right', fontsize=8)
    ax2.set_ylabel('QED Improvement\n(Child - Parent, higher is better)')
    ax2.set_title('QED Improvements Across Siblings\nStars indicate best sibling')
    ax2.grid(True, alpha=0.3)
    
    # Plot SAS improvements
    for p_idx, (parent_id, siblings) in enumerate(by_parent.items()):
        for s_idx, sibling in enumerate(siblings):
            x_pos = p_idx + s_idx * 0.2 - (len(siblings) - 1) * 0.1
            
            # Use stars for the best sibling
            marker = '*' if sibling['is_best_sas'] else 'o'
            size = 150 if sibling['is_best_sas'] else 100
            
            ax3.scatter(x_pos, sibling['sas_improvement'], 
                       color='tab:green', 
                       marker=marker, s=size, alpha=0.7)
    
    ax3.axhline(y=0, color='r', linestyle='-', alpha=0.3)
    ax3.set_xticks(range(len(by_parent)))
    ax3.set_xticklabels(parent_labels, rotation=45, ha='right', fontsize=8)
    ax3.set_ylabel('SAS Improvement\n(Parent - Child, lower is better)')
    ax3.set_title('SAS Improvements Across Siblings\nStars indicate best sibling')
    ax3.grid(True, alpha=0.3)
    
    # Create legend
    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:blue', 
                  markersize=10, label='Docking Score'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:orange', 
                  markersize=10, label='QED (Drug-likeness)'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='tab:green', 
                  markersize=10, label='SAS (Synthetic Accessibility)')
    ]
    
    # Add marker for best child
    legend_handles.append(plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='gray', 
                              markersize=12, label='Best performing sibling'))
    
    # Add legend
    legend_cols = min(4, len(legend_handles))
    plt.tight_layout(rect=[0, 0.2, 1, 0.95])
    
    fig.legend(handles=legend_handles, 
               loc='upper center', 
               bbox_to_anchor=(0.5, 0.10),
               ncol=legend_cols,
               frameon=True,
               fancybox=True,
               shadow=True,
               fontsize='medium',
               markerscale=1.5,
               columnspacing=1.0,
               handletextpad=0.5)
    
    # Add sibling stats text
    stats_text = f"Total molecules: {total_children} | Parents: {len(multi_child_parents)} | Avg siblings per parent: {total_children / len(multi_child_parents):.1f}"
    
    plt.figtext(0.5, 0.02, stats_text, 
                ha="center", fontsize=12, bbox={"facecolor":"lightgreen", "alpha":0.2, "pad":5})
    
    fig.suptitle(f'Comparison of Sibling Molecules (Medicinal Chemist) - Session: {session_id}', fontsize=16)
    
    # Save the figure
    fig_path = os.path.join(output_dir, 'medicinal_chemist_siblings_comparison.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created Medicinal Chemist sibling molecules comparison visualization: {fig_path}")

def extract_unique_ordered_smiles(csv_path, output_dir, session_id):
    """
    Extract unique SMILES from a CSV file, ordered by iteration and agent priority.
    
    Args:
        csv_path: Path to the CSV file
        output_dir: Directory to save the output file
        session_id: Session identifier for logging
    
    Returns:
        DataFrame with unique ordered SMILES
    """
    print(f"\nExtracting unique SMILES for session: {session_id}")
    
    # Read the CSV file
    df = pd.read_csv(csv_path)
    
    if df.empty:
        print(f"No data found in CSV for session {session_id}")
        return pd.DataFrame()
    
    # Add a priority column based on agent
    df['agent_priority'] = df['agent'].apply(lambda x: AGENT_PRIORITY.get(x, 999))  # Default high value for unknown agents
    
    # Convert iteration to numeric, with NaN values set to a high number (to sort them last)
    df['iteration_num'] = pd.to_numeric(df['iteration'], errors='coerce').fillna(999)
    
    # Sort by iteration and then by agent priority
    df_sorted = df.sort_values(['iteration_num', 'agent_priority'])
    
    # Get unique SMILES, keeping the first occurrence based on our sorting
    df_unique = df_sorted.drop_duplicates(subset=['smiles'])
    
    # Keep relevant columns and reset index
    result_df = df_unique[['molecule_id', 'smiles', 'status', 'parent_molecule_id', 
                            'docking_score', 'log_p', 'qed', 'sas', 'activity', 
                            'iteration', 'transformation_id', 'transformation_type', 'agent']]
    
    # Save to a new CSV file
    output_file = os.path.join(output_dir, 'unique_ordered_molecules.csv')
    result_df.to_csv(output_file, index=False)
    print(f"Created unique ordered SMILES CSV: {output_file}")
    
    return result_df

def analyze_iteration_improvements(csv_path, output_dir, session_id, single_agent=False, llm_only=False):
    """
    Analyze how molecule properties improve through iterations using unique molecules only.
    
    When llm_only=True:
    - Shows trends for all molecules grouped by iteration without any baseline comparison.
    - Useful for benchmarking LLM-only results.
    
    When single_agent=False (default) and llm_only=False:
    - De Novo Baseline: Calculated from AI Expert molecules (mean and best values).
    - Modified Trends: Plotted using only molecules where agent == 'Medicinal Chemist',
      starting from the first iteration they appear.
    
    When single_agent=True and llm_only=False:
    - De Novo Baseline: Calculated from transformation_type == 'molecule-generation' molecules.
    - Modified Trends: Plotted using only molecules where transformation_type == 'lead-optimization'.

    Args:
        csv_path: Path to the CSV file with molecule data
        output_dir: Directory to save the output figure
        session_id: Session identifier for logging/filenames
        single_agent: If True, use transformation_type instead of agent for filtering
        llm_only: If True, show all molecules by iteration without baseline comparison
    """
    print(f"\nAnalyzing property improvements for session: {session_id}")
    if llm_only:
        print(f" - Mode: LLM Only Benchmark")
        print(f" - Shows all molecules grouped by iteration without baseline comparison")
    elif single_agent:
        print(f" - Mode: Single Agent")
        print(f" - Baseline: transformation_type == 'molecule-generation' (Mean Red--, Best Blue..)")
        print(f" - Trends: transformation_type == 'lead-optimization' molecules")
    else:
        print(f" - Mode: Multi Agent")
        print(f" - Baseline: AI Expert molecules (Mean Red--, Best Blue..)")
        print(f" - Trends: 'Medicinal Chemist' agent molecules")

    # --- Data Loading and Initial Filtering ---
    unique_csv = os.path.join(output_dir, 'unique_ordered_molecules.csv')
    if not os.path.exists(unique_csv):
        print(f"Unique molecules CSV not found at {unique_csv}, using the original CSV: {csv_path}")
        if not os.path.exists(csv_path):
             print(f"ERROR: Original CSV not found either: {csv_path}")
             return None, None, None # Return Nones for failure
        df = pd.read_csv(csv_path)
    else:
        print(f"Reading unique molecules from {unique_csv}")
        df = pd.read_csv(unique_csv)

    if df.empty:
        print(f"No data found in CSV for session {session_id}")
        return None, None, None

    # Filter out molecules with specified activity
    df = df[df['activity'].isna() | (df['activity'] == '')]
    if df.empty:
        print(f"No data left after filtering active/inactive molecules.")
        return None, None, None

    # Convert iteration to numeric and drop errors
    df['iteration_num'] = pd.to_numeric(df['iteration'], errors='coerce')
    df = df.dropna(subset=['iteration_num'])
    df['iteration_num'] = df['iteration_num'].astype(int)
    if df.empty:
        print(f"No data left after filtering non-numeric iterations.")
        return None, None, None

    # --- De Novo Baseline Calculation ---
    if llm_only:
        # Skip baseline calculation for LLM-only mode
        baseline_df = pd.DataFrame()
        baseline_label = None
        baseline_mean_docking = baseline_mean_qed = baseline_mean_sas = baseline_mean_logp = None
        baseline_best_docking = baseline_best_qed = baseline_best_sas = None
        print("Skipping baseline calculation for LLM-only benchmark mode.")
    else:
        if single_agent:
            baseline_df = df[df['transformation_type'] == 'molecule-generation'].copy()
            baseline_label = "molecule-generation"
        else:
            baseline_df = df[df['agent'] == 'AI Expert'].copy()
            baseline_label = "AI Expert"
        
        if baseline_df.empty:
            print(f"ERROR: No data found for {baseline_label} molecules to use as de novo baseline.")
            return None, None, None # Baseline is crucial
        
        print(f"Using {len(baseline_df)} {baseline_label} molecules as De Novo Baseline.")

        # Calculate De Novo Mean values
        baseline_mean_docking = baseline_df['docking_score'].mean()
        baseline_mean_qed = baseline_df['qed'].mean()
        baseline_mean_sas = baseline_df['sas'].mean()
        baseline_mean_logp = baseline_df['log_p'].mean()

        # Calculate De Novo Best values
        baseline_best_docking = baseline_df['docking_score'].min() # Lower is better
        baseline_best_qed = baseline_df['qed'].max()             # Higher is better
        baseline_best_sas = baseline_df['sas'].min()             # Lower is better
        # LogP doesn't typically have a 'best', so we only use mean

        print(f"De Novo Baseline Mean ({baseline_label}): Docking={baseline_mean_docking:.3f}, QED={baseline_mean_qed:.3f}, SAS={baseline_mean_sas:.3f}, LogP={baseline_mean_logp:.3f}")
        print(f"De Novo Baseline Best ({baseline_label}): Docking={baseline_best_docking:.3f}, QED={baseline_best_qed:.3f}, SAS={baseline_best_sas:.3f}")


    # --- Optimization Molecule Analysis ---
    if llm_only:
        # For LLM-only mode, use all molecules regardless of agent or transformation type
        df_opt = df.copy()
        opt_label = "All LLM"
    elif single_agent:
        df_opt = df[df['transformation_type'] == 'lead-optimization'].copy()
        opt_label = "lead-optimization"
    else:
        df_opt = df[df['agent'] == 'Medicinal Chemist'].copy()
        opt_label = "Medicinal Chemist"

    if df_opt.empty:
        print(f"No molecules found with {opt_label}. Cannot plot trends.")
        # Still plot baseline? Decide based on requirements. For now, we'll exit plotting trends.
        iteration_stats_opt = pd.DataFrame() # Empty dataframe
    else:
        print(f"Found {len(df_opt)} molecules from '{opt_label}'.")
        # Ensure iteration_num is present and numeric for this subset too
        if 'iteration_num' not in df_opt.columns or df_opt['iteration_num'].isnull().all():
             print(f"ERROR: '{opt_label}' molecules lack valid iteration numbers.")
             iteration_stats_opt = pd.DataFrame()
        else:
            # Group by iteration for optimization molecules
            iteration_stats_opt = df_opt.groupby('iteration_num').agg(
                docking_score_mean=('docking_score', 'mean'),
                docking_score_min=('docking_score', 'min'),
                docking_score_count=('docking_score', 'count'),
                qed_mean=('qed', 'mean'),
                qed_max=('qed', 'max'),
                qed_count=('qed', 'count'),
                sas_mean=('sas', 'mean'),
                sas_min=('sas', 'min'),
                sas_count=('sas', 'count'),
                log_p_mean=('log_p', 'mean'),
                log_p_count=('log_p', 'count')
            ).reset_index()

            iteration_stats_opt = iteration_stats_opt.sort_values('iteration_num')
            print(f"Statistics for '{opt_label}' molecules:\n{iteration_stats_opt.head()}")


    # --- Plotting ---
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 11)) # Slightly larger fig
    
    if llm_only:
        title = f'LLM Molecule Properties Through Iterations - {session_id}'
        opt_display_label = "LLM"
        baseline_display_label = None
    elif single_agent:
        title = f'Lead Optimization vs. De Novo Baseline (molecule-generation) - {session_id}'
        opt_display_label = "LO"  # Short for Lead Optimization
        baseline_display_label = "De Novo"
    else:
        title = f'Medicinal Chemist Molecule Properties vs. De Novo Baseline (AI Expert) - {session_id}'
        opt_display_label = "MC"  # Short for Medicinal Chemist
        baseline_display_label = "De Novo"
    
    fig.suptitle(title, fontsize=16)

    # Common plotting elements
    def setup_ax(ax, title, ylabel):
        ax.set_title(title, fontsize=12)
        ax.set_xlabel('Iteration')
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.4)

    def add_counts(ax, df, x_col, y_col, count_col):
         if not df.empty:
            for _, row in df.iterrows():
                ax.annotate(f"n={int(row[count_col])}",
                           (row[x_col], row[y_col]),
                           textcoords="offset points", xytext=(0,10), ha='center', fontsize=8)

    # Plot Docking Score
    setup_ax(ax1, 'Docking Score (Lower is Better)', 'Docking Score')
    if not llm_only:
        if not pd.isna(baseline_mean_docking):
            ax1.axhline(y=baseline_mean_docking, color='r', linestyle='--', label=f'{baseline_display_label} Mean ({baseline_label})')
        if not pd.isna(baseline_best_docking):
            ax1.axhline(y=baseline_best_docking, color='b', linestyle=':', label=f'{baseline_display_label} Best ({baseline_label})')
    if not iteration_stats_opt.empty:
        ax1.plot(iteration_stats_opt['iteration_num'], iteration_stats_opt['docking_score_mean'], 'o-', color='tab:purple', label=f'{opt_display_label} Mean')
        ax1.plot(iteration_stats_opt['iteration_num'], iteration_stats_opt['docking_score_min'], 'x--', color='tab:purple', alpha=0.7, label=f'{opt_display_label} Best (Min)')
        add_counts(ax1, iteration_stats_opt, 'iteration_num', 'docking_score_mean', 'docking_score_count')
    ax1.legend(fontsize=9)

    # Plot QED
    setup_ax(ax2, 'QED (Higher is Better)', 'QED (Drug-likeness)')
    if not llm_only:
        if not pd.isna(baseline_mean_qed):
            ax2.axhline(y=baseline_mean_qed, color='r', linestyle='--', label=f'{baseline_display_label} Mean ({baseline_label})')
        if not pd.isna(baseline_best_qed):
            ax2.axhline(y=baseline_best_qed, color='b', linestyle=':', label=f'{baseline_display_label} Best ({baseline_label})')
    if not iteration_stats_opt.empty:
        ax2.plot(iteration_stats_opt['iteration_num'], iteration_stats_opt['qed_mean'], 'o-', color='tab:orange', label=f'{opt_display_label} Mean')
        ax2.plot(iteration_stats_opt['iteration_num'], iteration_stats_opt['qed_max'], 'x--', color='tab:orange', alpha=0.7, label=f'{opt_display_label} Best (Max)')
        add_counts(ax2, iteration_stats_opt, 'iteration_num', 'qed_mean', 'qed_count')
    ax2.legend(fontsize=9)

    # Plot SAS
    setup_ax(ax3, 'SAS (Lower is Better)', 'SAS (Synthetic Accessibility)')
    if not llm_only:
        if not pd.isna(baseline_mean_sas):
            ax3.axhline(y=baseline_mean_sas, color='r', linestyle='--', label=f'{baseline_display_label} Mean ({baseline_label})')
        if not pd.isna(baseline_best_sas):
            ax3.axhline(y=baseline_best_sas, color='b', linestyle=':', label=f'{baseline_display_label} Best ({baseline_label})')
    if not iteration_stats_opt.empty:
        ax3.plot(iteration_stats_opt['iteration_num'], iteration_stats_opt['sas_mean'], 'o-', color='tab:green', label=f'{opt_display_label} Mean')
        ax3.plot(iteration_stats_opt['iteration_num'], iteration_stats_opt['sas_min'], 'x--', color='tab:green', alpha=0.7, label=f'{opt_display_label} Best (Min)')
        add_counts(ax3, iteration_stats_opt, 'iteration_num', 'sas_mean', 'sas_count')
    ax3.legend(fontsize=9)

    # Plot LogP
    setup_ax(ax4, 'LogP (Informational)', 'LogP (Lipophilicity)')
    if not llm_only:
        if not pd.isna(baseline_mean_logp):
            ax4.axhline(y=baseline_mean_logp, color='r', linestyle='--', label=f'{baseline_display_label} Mean ({baseline_label})')
    # No 'best' line for LogP generally
    if not iteration_stats_opt.empty:
        ax4.plot(iteration_stats_opt['iteration_num'], iteration_stats_opt['log_p_mean'], 'o-', color='tab:red', label=f'{opt_display_label} Mean')
        add_counts(ax4, iteration_stats_opt, 'iteration_num', 'log_p_mean', 'log_p_count')
    ax4.legend(fontsize=9)


    # --- Improvement Calculation & Subtitle ---
    if llm_only:
        # For LLM-only mode, show iteration progression instead of baseline comparison
        if not iteration_stats_opt.empty and len(iteration_stats_opt) > 1:
            first_iter_stats = iteration_stats_opt.iloc[0]
            last_iter_stats = iteration_stats_opt.iloc[-1]
            first_iter_num = int(first_iter_stats['iteration_num'])
            last_iter_num = int(last_iter_stats['iteration_num'])
            
            docking_change = first_iter_stats['docking_score_mean'] - last_iter_stats['docking_score_mean']
            docking_pct = (docking_change / abs(first_iter_stats['docking_score_mean'])) * 100 if abs(first_iter_stats['docking_score_mean']) > 1e-9 else 0
            
            qed_change = last_iter_stats['qed_mean'] - first_iter_stats['qed_mean']
            qed_pct = (qed_change / first_iter_stats['qed_mean']) * 100 if abs(first_iter_stats['qed_mean']) > 1e-9 else 0
            
            sas_change = first_iter_stats['sas_mean'] - last_iter_stats['sas_mean']
            sas_pct = (sas_change / first_iter_stats['sas_mean']) * 100 if abs(first_iter_stats['sas_mean']) > 1e-9 else 0
            
            subtitle = (f"Iteration Progress (Iter {first_iter_num} to {last_iter_num}): "
                       f"Docking: {docking_pct:+.1f}%, QED: {qed_pct:+.1f}%, SAS: {sas_pct:+.1f}%")
            plt.figtext(0.5, 0.01, subtitle, ha='center', fontsize=10,
                       bbox={"facecolor":"lightgreen", "alpha":0.2, "pad":4})
        elif iteration_stats_opt.empty:
            subtitle = "No LLM iterations found"
            plt.figtext(0.5, 0.01, subtitle, ha='center', fontsize=10)
        else:
            subtitle = "Only one iteration available - no progression to show"
            plt.figtext(0.5, 0.01, subtitle, ha='center', fontsize=10)
    else:
        subtitle = f"Improvement vs {baseline_display_label} Mean:"
        if not iteration_stats_opt.empty and not pd.isna(baseline_mean_docking) and not pd.isna(baseline_mean_qed) and not pd.isna(baseline_mean_sas):
            last_opt_iter_stats = iteration_stats_opt.iloc[-1]
            last_opt_iter_num = int(last_opt_iter_stats['iteration_num'])

            docking_improvement = baseline_mean_docking - last_opt_iter_stats['docking_score_mean']
            docking_pct = (docking_improvement / abs(baseline_mean_docking)) * 100 if abs(baseline_mean_docking) > 1e-9 else 0

            qed_improvement = last_opt_iter_stats['qed_mean'] - baseline_mean_qed
            qed_pct = (qed_improvement / baseline_mean_qed) * 100 if abs(baseline_mean_qed) > 1e-9 else 0

            sas_improvement = baseline_mean_sas - last_opt_iter_stats['sas_mean']
            sas_pct = (sas_improvement / baseline_mean_sas) * 100 if abs(baseline_mean_sas) > 1e-9 else 0

            subtitle += (f" ({baseline_label} Baseline to {opt_label} Iter {last_opt_iter_num}) "
                        f" Docking: {docking_pct:+.1f}%, QED: {qed_pct:+.1f}%, SAS: {sas_pct:+.1f}%") # Added + sign
            plt.figtext(0.5, 0.01, subtitle, ha='center', fontsize=10,
                       bbox={"facecolor":"lightblue", "alpha":0.2, "pad":4})
        elif iteration_stats_opt.empty:
             subtitle += f" (No {opt_label} iterations found)"
             plt.figtext(0.5, 0.01, subtitle, ha='center', fontsize=10)
        else:
             subtitle += f" (Could not calculate - missing baseline or {opt_label} data)"
             plt.figtext(0.5, 0.01, subtitle, ha='center', fontsize=10)


    # --- Final Adjustments & Saving ---
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust rect for subtitle/title

    # Save Figure
    if llm_only:
        fig_filename = f'llm_only_iteration_progress_{session_id}.png'
    elif single_agent:
        fig_filename = f'leadopt_vs_denovo_improvements_{session_id}.png'
    else:
        fig_filename = f'mc_vs_denovo_improvements_{session_id}.png'
    fig_path = os.path.join(output_dir, fig_filename)
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Created plot: {fig_path}")

    # Save Table for Optimization Stats
    if not iteration_stats_opt.empty:
        if llm_only:
            table_filename = f'llm_only_iteration_stats_{session_id}.csv'
        elif single_agent:
            table_filename = f'lead_optimization_iteration_stats_{session_id}.csv'
        else:
            table_filename = f'medicinal_chemist_iteration_stats_{session_id}.csv'
        table_path = os.path.join(output_dir, table_filename)
        iteration_stats_opt.to_csv(table_path, index=False)
        print(f"Created {opt_label} statistics table: {table_path}")
    else:
        print(f"No {opt_label} statistics table saved.")

    # Return baseline df, optimization stats df
    return baseline_df, iteration_stats_opt

def process_session_file(session_path, session_id, output_dir):
    print(f"\nProcessing session: {session_id}")
    
    # Load the JSON data from file
    with open(session_path, 'r') as f:
        session_data = json.load(f)
    
    # Extract session and last transformation IDs directly from the JSON
    if 'session_id' in session_data and 'last_transform_id' in session_data:
        gurnemanz_session_id = session_data['session_id']
        last_transformation_id = session_data['last_transform_id']
        
        print(f"Session ID: {gurnemanz_session_id}")
        print(f"Last transformation ID: {last_transformation_id}")
        
        # Get the transformation chain with molecules
        try:
            data = get_transformation_chain_with_molecules(gurnemanz_session_id, last_transformation_id)
            
            if not data or 'steps' not in data or len(data['steps']) == 0:
                print(f"No valid data found for session {session_id}")
                return
        except RuntimeError as e:
            # Check if this is a "Session not found" error
            if "Session not found" in str(e):
                print(f"Skipping session {session_id} - Session not found in database")
                return
            else:
                # Re-raise other runtime errors
                raise
        
        # Extract steps info
        steps_df = pd.DataFrame(data['steps'])
        print("\n=== TRANSFORMATION STEPS ===")
        print(steps_df[['transformation_id', 'transformation_type', 'agent', 'iteration']].head())

        # Extract molecules info
        molecules = data.get('molecules', {})
        mol_data = []

        # Create a mapping of molecule_id to iteration and transformation_id
        molecule_to_metadata = {}
        for step in data['steps']:
            iteration = step.get('iteration', '')
            transformation_id = step.get('transformation_id', '')
            for mol_id in step.get('output_molecule_ids', []):
                molecule_to_metadata[mol_id] = {'iteration': iteration, 'transformation_id': transformation_id}

        for mol_id, mol_info in molecules.items():
            mol_dict = {
                'molecule_id': mol_id,
                'smiles': mol_info.get('structure', {}).get('smiles', ''),
                'status': mol_info.get('properties', {}).get('status', ''),
                'parent_molecule_id': mol_info.get('properties', {}).get('parent_molecule_id', ''),
                'docking_score': mol_info.get('computed_properties', {}).get('docking_score', ''),
                'log_p': mol_info.get('computed_properties', {}).get('log_p', ''),
                'qed': mol_info.get('computed_properties', {}).get('qed', ''),
                'sas': mol_info.get('computed_properties', {}).get('sas', ''),
                'activity': mol_info.get('activity', ''),
                'iteration': molecule_to_metadata.get(mol_id, {}).get('iteration', ''),
                'transformation_id': molecule_to_metadata.get(mol_id, {}).get('transformation_id', '')
            }
            mol_data.append(mol_dict)

        molecules_df = pd.DataFrame(mol_data)
        if not molecules_df.empty:
            print("\n=== MOLECULES ===")
            print(molecules_df[['molecule_id', 'status', 'parent_molecule_id', 'docking_score', 'log_p', 'sas', 'qed', 'activity', 'iteration', 'transformation_id']].head(10))

            # Create a combined dataframe by merging steps_df and molecules_df
            combined_df = pd.merge(
                molecules_df,
                steps_df[['transformation_id', 'transformation_type', 'agent']],
                on='transformation_id',
                how='left'
            )

            print("\n=== COMBINED DATAFRAME ===")
            print(combined_df.head(10))

            # Save to a CSV file in the session directory
            output_csv = os.path.join(output_dir, 'combined_molecules_transformations.csv')
            combined_df.to_csv(output_csv, index=False)
            print(f"\nCreated CSV file: {output_csv}")
            
            # Extract unique ordered SMILES
            extract_unique_ordered_smiles(output_csv, output_dir, session_id)
            
            # Analyze parent-child molecule improvements
            analyze_parent_child_improvements(output_csv, output_dir, session_id)
            
            # Analyze improvements through iterations (using unique molecules)
            analyze_iteration_improvements(output_csv, output_dir, session_id)
        else:
            print(f"No molecule data found for session {session_id}")
    else:
        print(f"Missing required data in session file {session_id}")
        print(f"Keys found: {list(session_data.keys())}")

def create_multi_panel_trajectory_comparison(session_dirs, output_dir, comparison_name="multi_iteration"):
    """
    Option 2: Multi-Panel Trajectory Comparison
    Create 4 panels (one per metric) with all iteration runs overlaid.
    
    Args:
        session_dirs: List of session directory paths
        output_dir: Directory to save the output figure
        comparison_name: Name for the comparison (used in filename)
    """
    print(f"\nCreating multi-panel trajectory comparison for {len(session_dirs)} sessions")
    
    # Collect data from all sessions
    all_session_data = {}
    
    for session_dir in session_dirs:
        session_id = os.path.basename(session_dir)
        csv_path = os.path.join(session_dir, 'combined_molecules_transformations.csv')
        unique_csv_path = os.path.join(session_dir, 'unique_ordered_molecules.csv')
        
        # Use unique CSV if available, otherwise original
        if os.path.exists(unique_csv_path):
            df = pd.read_csv(unique_csv_path)
        elif os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
        else:
            print(f"No CSV found for session {session_id}, skipping")
            continue
        
        if df.empty:
            continue
            
        # Filter and process data (same logic as analyze_iteration_improvements)
        df = df[df['activity'].isna() | (df['activity'] == '')]
        if df.empty:
            continue
            
        df['iteration_num'] = pd.to_numeric(df['iteration'], errors='coerce')
        df = df.dropna(subset=['iteration_num'])
        df['iteration_num'] = df['iteration_num'].astype(int)
        
        if df.empty:
            continue
        
        # Get baseline (AI Expert) data
        baseline_df = df[df['agent'] == 'AI Expert'].copy()
        if baseline_df.empty:
            print(f"No AI Expert baseline for session {session_id}, skipping")
            continue
            
        # Get Medicinal Chemist data
        mc_df = df[df['agent'] == 'Medicinal Chemist'].copy()
        if mc_df.empty:
            print(f"No Medicinal Chemist data for session {session_id}, skipping")
            continue
        
        # Calculate baseline values
        baseline_stats = {
            'docking_mean': baseline_df['docking_score'].mean(),
            'qed_mean': baseline_df['qed'].mean(),
            'sas_mean': baseline_df['sas'].mean(),
            'logp_mean': baseline_df['log_p'].mean()
        }
        
        # Group MC data by iteration
        mc_stats = mc_df.groupby('iteration_num').agg(
            docking_score_mean=('docking_score', 'mean'),
            qed_mean=('qed', 'mean'),
            sas_mean=('sas', 'mean'),
            log_p_mean=('log_p', 'mean'),
            count=('docking_score', 'count')
        ).reset_index()
        
        # Store session data
        all_session_data[session_id] = {
            'baseline': baseline_stats,
            'mc_stats': mc_stats,
            'max_iteration': mc_stats['iteration_num'].max()
        }
    
    if not all_session_data:
        print("No valid session data found for comparison")
        return
    
    # Create the multi-panel plot
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Color palette for different sessions
    colors = plt.cm.tab10(np.linspace(0, 1, len(all_session_data)))
    
    for i, (session_id, data) in enumerate(all_session_data.items()):
        color = colors[i]
        mc_stats = data['mc_stats']
        baseline = data['baseline']
        max_iter = data['max_iteration']
        
        # Extract iteration number from session_id if it follows pattern
        iter_num = "Unknown"
        if "iter_" in session_id:
            iter_num = session_id.split("iter_")[-1]
        
        label = f"Iter Limit {iter_num} (n=3)"
        
        # Plot Docking Score
        ax1.plot(mc_stats['iteration_num'], mc_stats['docking_score_mean'], 
                'o-', color=color, label=label, linewidth=2, markersize=6)
        ax1.axhline(y=baseline['docking_mean'], color=color, linestyle='--', alpha=0.5)
        
        # Plot QED
        ax2.plot(mc_stats['iteration_num'], mc_stats['qed_mean'], 
                'o-', color=color, label=label, linewidth=2, markersize=6)
        ax2.axhline(y=baseline['qed_mean'], color=color, linestyle='--', alpha=0.5)
        
        # Plot SAS
        ax3.plot(mc_stats['iteration_num'], mc_stats['sas_mean'], 
                'o-', color=color, label=label, linewidth=2, markersize=6)
        ax3.axhline(y=baseline['sas_mean'], color=color, linestyle='--', alpha=0.5)
        
        # Plot LogP
        ax4.plot(mc_stats['iteration_num'], mc_stats['log_p_mean'], 
                'o-', color=color, label=label, linewidth=2, markersize=6)
        ax4.axhline(y=baseline['logp_mean'], color=color, linestyle='--', alpha=0.5)
    
    # Setup axes
    ax1.set_title('Docking Score (Lower is Better)', fontsize=14)
    ax1.set_xlabel('Iteration')
    ax1.set_ylabel('Docking Score')
    ax1.grid(True, alpha=0.3)
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax2.set_title('QED (Higher is Better)', fontsize=14)
    ax2.set_xlabel('Iteration')
    ax2.set_ylabel('QED (Drug-likeness)')
    ax2.grid(True, alpha=0.3)
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax3.set_title('SAS (Lower is Better)', fontsize=14)
    ax3.set_xlabel('Iteration')
    ax3.set_ylabel('SAS (Synthetic Accessibility)')
    ax3.grid(True, alpha=0.3)
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    ax4.set_title('LogP (Informational)', fontsize=14)
    ax4.set_xlabel('Iteration')
    ax4.set_ylabel('LogP (Lipophilicity)')
    ax4.grid(True, alpha=0.3)
    ax4.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.suptitle('Medicinal Chemist Trajectory Comparison Across Different Iteration Limits', fontsize=16)
    plt.tight_layout()
    
    # Save the figure
    fig_path = os.path.join(output_dir, f'trajectory_comparison_{comparison_name}.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created trajectory comparison plot: {fig_path}")
    return fig_path

def create_improvement_heatmap(session_dirs, output_dir, comparison_name="multi_iteration"):
    """
    Option 3: Heatmap of Improvements
    Create a heatmap showing improvement percentages for each metric across iteration numbers.
    
    Args:
        session_dirs: List of session directory paths
        output_dir: Directory to save the output figure
        comparison_name: Name for the comparison (used in filename)
    """
    print(f"\nCreating improvement heatmap for {len(session_dirs)} sessions")
    
    # Collect improvement data
    improvement_data = []
    
    for session_dir in session_dirs:
        session_id = os.path.basename(session_dir)
        csv_path = os.path.join(session_dir, 'combined_molecules_transformations.csv')
        unique_csv_path = os.path.join(session_dir, 'unique_ordered_molecules.csv')
        
        # Use unique CSV if available, otherwise original
        if os.path.exists(unique_csv_path):
            df = pd.read_csv(unique_csv_path)
        elif os.path.exists(csv_path):
            df = pd.read_csv(csv_path)
        else:
            continue
        
        if df.empty:
            continue
            
        # Filter and process data
        df = df[df['activity'].isna() | (df['activity'] == '')]
        if df.empty:
            continue
            
        df['iteration_num'] = pd.to_numeric(df['iteration'], errors='coerce')
        df = df.dropna(subset=['iteration_num'])
        df['iteration_num'] = df['iteration_num'].astype(int)
        
        if df.empty:
            continue
        
        # Get baseline (AI Expert) data
        baseline_df = df[df['agent'] == 'AI Expert'].copy()
        if baseline_df.empty:
            continue
            
        # Get Medicinal Chemist data
        mc_df = df[df['agent'] == 'Medicinal Chemist'].copy()
        if mc_df.empty:
            continue
        
        # Calculate baseline means
        baseline_docking = baseline_df['docking_score'].mean()
        baseline_qed = baseline_df['qed'].mean()
        baseline_sas = baseline_df['sas'].mean()
        baseline_logp = baseline_df['log_p'].mean()
        
        # Get final iteration stats
        final_iter = mc_df['iteration_num'].max()
        final_stats = mc_df[mc_df['iteration_num'] == final_iter]
        
        if final_stats.empty:
            continue
        
        final_docking = final_stats['docking_score'].mean()
        final_qed = final_stats['qed'].mean()
        final_sas = final_stats['sas'].mean()
        final_logp = final_stats['log_p'].mean()
        
        # Calculate improvements (same logic as in analyze_iteration_improvements)
        docking_improvement = ((baseline_docking - final_docking) / abs(baseline_docking)) * 100 if abs(baseline_docking) > 1e-9 else 0
        qed_improvement = ((final_qed - baseline_qed) / baseline_qed) * 100 if abs(baseline_qed) > 1e-9 else 0
        sas_improvement = ((baseline_sas - final_sas) / baseline_sas) * 100 if abs(baseline_sas) > 1e-9 else 0
        logp_improvement = ((final_logp - baseline_logp) / abs(baseline_logp)) * 100 if abs(baseline_logp) > 1e-9 else 0
        
        # Extract iteration number from session_id
        iter_num = "Unknown"
        if "iter_" in session_id:
            try:
                iter_num = int(session_id.split("iter_")[-1])
            except:
                iter_num = session_id.split("iter_")[-1]
        
        improvement_data.append({
            'Iteration_Limit': iter_num,
            'Docking': docking_improvement,
            'QED': qed_improvement,
            'SAS': sas_improvement,
            'LogP': logp_improvement,
            'Session_ID': session_id,
            'Final_Iter': final_iter
        })
    
    if not improvement_data:
        print("No improvement data found for heatmap")
        return
    
    # Convert to DataFrame
    df_improvements = pd.DataFrame(improvement_data)
    
    # Sort by iteration limit if numeric
    try:
        df_improvements['Iteration_Limit_Numeric'] = pd.to_numeric(df_improvements['Iteration_Limit'])
        df_improvements = df_improvements.sort_values('Iteration_Limit_Numeric')
    except:
        pass
    
    # Create pivot table for heatmap
    heatmap_data = df_improvements.set_index('Iteration_Limit')[['Docking', 'QED', 'SAS', 'LogP']].T
    
    # Create the heatmap
    plt.figure(figsize=(12, 8))
    
    # Create custom colormap (green for positive, red for negative)
    cmap = sns.diverging_palette(10, 150, as_cmap=True, center='light')
    
    # Create heatmap
    ax = sns.heatmap(heatmap_data, 
                     annot=True, 
                     fmt='.1f', 
                     cmap=cmap, 
                     center=0, 
                     cbar_kws={'label': 'Improvement (%)'}, 
                     square=True, 
                     linewidths=1,
                     annot_kws={'size': 12, 'weight': 'bold'})
    
    # Customize the plot
    plt.title('Medicinal Chemist Improvement vs. AI Expert Baseline\n(% Improvement by Iteration Limit)', 
              fontsize=16, pad=20)
    plt.xlabel('Iteration Limit', fontsize=14)
    plt.ylabel('Molecular Properties', fontsize=14)
    
    # Add metric interpretation labels
    metric_labels = {
        'Docking': 'Docking Score\n(Lower is Better)',
        'QED': 'QED\n(Higher is Better)', 
        'SAS': 'Synthetic Accessibility\n(Lower is Better)',
        'LogP': 'Lipophilicity\n(Informational)'
    }
    
    # Update y-axis labels
    ax.set_yticklabels([metric_labels.get(label.get_text(), label.get_text()) 
                       for label in ax.get_yticklabels()], rotation=0)
    
    # Add text annotation explaining positive/negative values
    plt.figtext(0.02, 0.02, 
                "Positive values (green) indicate improvement over AI Expert baseline\n"
                "Negative values (red) indicate worse performance than baseline", 
                fontsize=10, ha='left',
                bbox={"facecolor":"lightgray", "alpha":0.3, "pad":5})
    
    plt.tight_layout()
    
    # Save the figure
    fig_path = os.path.join(output_dir, f'improvement_heatmap_{comparison_name}.png')
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Created improvement heatmap: {fig_path}")
    
    # Also save the raw data
    data_path = os.path.join(output_dir, f'improvement_data_{comparison_name}.csv')
    df_improvements.to_csv(data_path, index=False)
    print(f"Saved improvement data: {data_path}")
    
    return fig_path, data_path

if __name__ == "__main__":
    # Main script execution
    print(f"Looking for gurnemanz_session.json files in: {RUN_DIR}")
    session_count = 0
    skipped_count = 0

    # Walk through all directories in RUN_DIR
    for session_dir in os.listdir(RUN_DIR):
        session_path = os.path.join(RUN_DIR, session_dir)
        
        # Check if it's a directory
        if os.path.isdir(session_path):
            # Check if a gurnemanz_session.json file exists in this directory
            json_path = os.path.join(session_path, 'gurnemanz_session.json')
            if os.path.isfile(json_path):
                # Check for existing output files
                csv_exists = os.path.isfile(os.path.join(session_path, 'combined_molecules_transformations.csv'))
                unique_csv_exists = os.path.isfile(os.path.join(session_path, 'unique_ordered_molecules.csv'))
                png_exists = os.path.isfile(os.path.join(session_path, 'medicinal_chemist_improvements.png'))
                
                # Skip if output files already exist
                if csv_exists and unique_csv_exists and png_exists:
                    print(f"Skipping session {session_dir} - output files already exist")
                    skipped_count += 1
                    continue
                    
                # Process this session using the directory name as identifier
                process_session_file(json_path, session_dir, session_path)
                session_count += 1

    print(f"\nProcessed {session_count} session files")
    print(f"Skipped {skipped_count} session files that already had output files")