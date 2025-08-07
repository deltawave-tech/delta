import random
import string
from datetime import datetime
import json
from src.tools.gurnemanz_utils import format_gurnemanz_for_llm
import os
import logging
from logging.handlers import RotatingFileHandler
import pandas as pd
import numpy as np

# RDKit imports for molecular analysis
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, AllChem
    import torch
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    print("RDKit not available. Molecular analysis functions will not work.")

# FCD import for molecular distance calculation
try:
    from fcd_torch import FCD
    FCD_AVAILABLE = True
except ImportError:
    FCD_AVAILABLE = False
    print("FCD not available. Install fcd_torch for FCD metric calculation.")

# Utility Functions
def generate_run_id():
    timestamp = datetime.now().strftime("%m%d_%H%M")
    random_chars = ''.join(random.choices(string.ascii_uppercase + string.digits, k=3))
    return f"{timestamp}_{random_chars}"

def timeout_handler(signum, frame):
    raise TimeoutError("Research cycle stopped due to timeout")

def save_conversation_to_markdown(history, run_dir, filename):
    with open(f"{run_dir}/{filename}", 'w', encoding='utf-8') as f:
        f.write(f"# {filename.replace('.md', '').replace('_', ' ').title()}\n\n")
        for i, message in enumerate(history):
            f.write(f"## {message['role'].title()}\n\n")

            content = message.get('content', '')

            # Attempt to parse content if it's a JSON string
            try:
                # This handles strings that are valid JSON
                parsed_content = json.loads(content)
                
                # If parsing is successful, we assume it's a structured message
                if isinstance(parsed_content, list) and parsed_content:
                    # Handle formats like [{"type": "text", "text": "..."}]
                    all_text = ""
                    for item in parsed_content:
                        if isinstance(item, dict) and 'text' in item:
                            all_text += item['text']
                    if all_text:
                        f.write(f"{all_text.strip()}\n\n")
                    else: # Fallback for other list-based JSON
                         f.write(f"```json\n{json.dumps(parsed_content, indent=2)}\n```\n\n")

                elif isinstance(parsed_content, dict):
                     # Fallback for other dict-based JSON
                    f.write(f"```json\n{json.dumps(parsed_content, indent=2)}\n```\n\n")
                else:
                    # Content is a JSON-encoded string, but not a structure we handle specially.
                    # Write it out as plain text.
                    f.write(f"{content}\n\n")

            except (json.JSONDecodeError, TypeError):
                # Content is not a valid JSON string, treat as plain text.
                if isinstance(content, str):
                    f.write(f"{content.strip()}\n\n")
                else:
                    # Fallback for non-string, non-JSON content.
                    f.write(f"```\n{content}\n```\n\n")

            if i < len(history) - 1:
                f.write("---\n\n")

def validate_history(history, agent_title="Unknown"):
    """Validate history to ensure it doesn't contain 'system' role messages"""
    for i, msg in enumerate(history):
        if msg.get("role") == "system":
            print(f"WARNING: Found 'system' role in history for {agent_title} at position {i}. This will cause API errors.")
            # Replace with a user role to prevent API errors
            msg["role"] = "user"
            print(f"Automatically converted 'system' role to 'user' role to prevent API errors.")
    return history


def format_agent_query(agent, iteration, max_iterations):
    """Format a message for an agent based on its title and the current iteration."""
    if agent.title == "Database Agent":
        return (
            f"{agent.title}, for iteration {iteration}/{max_iterations}, please search the relevant databases for the current research target. "
            f"Share your findings (protein sequence, PDB path, active/inactive compounds) with the team. Follow the molecule sharing format."
            f"Keep in mind other agents don't see the outputs of your tools, so you need to share the pdb file path and protein sequence in your final outputs and the root path of the .sdf"
            f'At the end of your response this should appear: \n **PDB file path**: <pdb_file_path> \n **Protein sequence**: <protein_sequence>'
        )
    elif agent.title == "AI Expert":
        return (
            f"{agent.title}, for iteration {iteration}/{max_iterations}, generate novel molecules based on the current understanding. "
            f"Share your generated molecules and rationale. Follow the molecule sharing format."
            f"Do NOT generates molecules by your self ! You can't only generate new molecules by calling your tools !!"
        )
    elif agent.title == "Medicinal Chemist":
        return (
            f"{agent.title}, for iteration {iteration}/{max_iterations}, review the molecules presented. "
            f"Propose and execute modifications to improve them, then evaluate the modified molecules using your tools. "
            f"Share your work. Follow the molecule sharing format."
            f"When you are reference the SMILES and its Friendly ID when you are thinking a enhancement of a molecule."
            f"You can do multiple modifications at once to send for evaluation to the VINA_REPORT tool."
        )
    elif agent.title == "Ranking Agent":
        return (
            f"{agent.title}, for iteration {iteration}/{max_iterations}, rank the molecules discussed so far. "
            f"Provide your top candidates and the ranking rationale. Follow the molecule sharing format."
            f"Always keep the format for the molecules of previous agent, DO NOT MODIFY ANY METRICS OR ID."
            f"Your role is just to rank the molecules, do not modify anything."
        )
    elif agent.title == "Scientific Critic":
        return (
            f"{agent.title}, for iteration {iteration}/{max_iterations}, critically review the team's progress, "
            f"methodologies, and conclusions from the preceding discussion. Provide constructive, scientifically rigorous feedback."
            f"WHEN you critics a molecule, you MUST specify its Friendly ID and SMILES and any other specific metric important for your critics to help other agents to understand what molecule you are criticizing"
        )
    elif agent.title == "One Person Biotech":
        if iteration == 1:
            return (
                f"{agent.title}, for iteration {iteration}/{max_iterations}, as the one-person biotech researcher, "
                f"you must orchestrate the entire *in silico* drug discovery project. Start by defining your strategic "
                f"approach for all {max_iterations} iterations, then execute the following for this iteration: "
                f"1) Use database tools to gather target protein information, structures, and known compounds "
                f"2) Generate novel molecules using VINA_MOL_GEN (only in iteration 1) "
                f"3) Critically evaluate all computational results and maintain scientific rigor throughout. "
                f"Always include **PDB file path** and **Protein sequence** in your outputs. Follow the molecule sharing format."
            )
        else:
            return (
                f"{agent.title}, for iteration {iteration}/{max_iterations}, continue your strategic leadership of the project. "
                f"Focus on: 1) Manual medicinal chemistry modifications to existing molecules (do NOT use VINA_MOL_GEN) "
                f"2) Evaluate modified molecules using VINA_REPORT tool 3) Synthesize progress and adapt strategy "
                f"4) Critically review your own work for scientific rigor. When referencing molecules, use their Friendly ID and SMILES. "
                f"You can evaluate multiple modifications in a single VINA_REPORT call."
            )
    else:
        # Default message format
        return (
            f"You are {agent.title}. This is iteration {iteration} of {max_iterations}. "
            f"Please contribute your expertise to the current discussion."
        )


def setup_logger(run_dir, run_id):
    """Set up a logger for the multi-agent system to record all agent messages"""
    # Create logs directory if it doesn't exist
    logs_dir = os.path.join(run_dir, 'logs')
    os.makedirs(logs_dir, exist_ok=True)
    
    # Configure logger
    log_file = os.path.join(logs_dir, f'multi_agent_conversation_{run_id}.log')
    logger = logging.getLogger(f'multi_agent_{run_id}')  # Use unique logger name per run
    logger.setLevel(logging.INFO)
    
    # Clear any existing handlers to avoid duplicate logs
    if logger.handlers:
        logger.handlers.clear()
    
    # Create and add a file handler
    file_handler = RotatingFileHandler(log_file, maxBytes=10*1024*1024, backupCount=5)
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    # Optionally add a console handler to see logs in the terminal too
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    return logger

def calculate_fcd_scores(gen_smiles, ref_smiles):
    """
    Calculate FCD (Fréchet ChemNet Distance) scores between generated and reference molecules.
    
    Args:
        gen_smiles (list): List of generated SMILES strings
        ref_smiles (list): Reference set of SMILES strings
    
    Returns:
        dict: Dictionary containing FCD score
    """
    if not FCD_AVAILABLE:
        print("FCD calculation skipped - fcd_torch not available")
        return {}
    
    if not gen_smiles or not ref_smiles:
        print("FCD calculation skipped - insufficient molecules")
        return {}
    
    try:
        device = get_device()
        fcd = FCD(device=device)
        
        # Calculate FCD against reference set
        fcd_score = fcd(gen_smiles, ref_smiles)
        
        return {'fcd_score': float(fcd_score)}
        
    except Exception as e:
        print(f"Error calculating FCD scores: {e}")
        return {}


def calculate_molecular_metrics(molecules_df, reference_molecules=None, docking_score_col='docking_score', smiles_col='smiles'):
    """
    Calculate comprehensive molecular metrics including validity, uniqueness, novelty, 
    similarity, Lipinski properties, docking statistics, SA scores, LogP, and FCD.
    
    Args:
        molecules_df (pd.DataFrame): DataFrame containing molecules with SMILES and docking scores
        reference_molecules (list, optional): List of reference SMILES for novelty and FCD calculation
        docking_score_col (str): Column name for docking scores
        smiles_col (str): Column name for SMILES strings
    
    Returns:
        dict: Dictionary containing all calculated metrics
    """
    if not RDKIT_AVAILABLE:
        raise ImportError("RDKit is required for molecular analysis. Please install rdkit-pypi")
    
    results = {}
    smiles_list = molecules_df[smiles_col].tolist()
    
    # 1. Validity - Check if SMILES are valid
    valid_mols = []
    valid_smiles = []
    valid_indices = []
    
    for i, smiles in enumerate(smiles_list):
        if pd.isna(smiles) or smiles == '':
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            valid_mols.append(mol)
            valid_smiles.append(smiles)
            valid_indices.append(i)
    
    total_count = len(smiles_list)
    valid_count = len(valid_smiles)
    validity_rate = valid_count / total_count if total_count > 0 else 0
    
    results['validity'] = {
        'total_molecules': total_count,
        'valid_molecules': valid_count,
        'invalid_molecules': total_count - valid_count,
        'validity_rate': validity_rate
    }
    
    if valid_count == 0:
        print("No valid molecules found. Cannot calculate other metrics.")
        return results
    
    # 2. Uniqueness - Count unique valid SMILES
    unique_smiles = list(set(valid_smiles))
    uniqueness_rate = len(unique_smiles) / len(valid_smiles) if len(valid_smiles) > 0 else 0
    
    results['uniqueness'] = {
        'total_valid': len(valid_smiles),
        'unique_molecules': len(unique_smiles),
        'duplicate_molecules': len(valid_smiles) - len(unique_smiles),
        'uniqueness_rate': uniqueness_rate
    }
    
    # 3. Novelty - Compare to reference molecules
    if reference_molecules is not None:
        reference_set = set(reference_molecules)
        novel_molecules = set(unique_smiles) - reference_set
        novelty_rate = len(novel_molecules) / len(unique_smiles) if len(unique_smiles) > 0 else 0
        
        results['novelty'] = {
            'reference_molecules': len(reference_molecules),
            'unique_generated': len(unique_smiles),
            'novel_molecules': len(novel_molecules),
            'known_molecules': len(unique_smiles) - len(novel_molecules),
            'novelty_rate': novelty_rate
        }
    
    # 4. Lipinski Rule of 5 and Molecular Properties
    lipinski_results = []
    
    # Check if we have pre-computed columns
    valid_df = molecules_df.iloc[valid_indices]
    has_logp_col = 'log_p' in molecules_df.columns
    has_sas_col = 'sas' in molecules_df.columns
    has_qed_col = 'qed' in molecules_df.columns
    
    # Use existing columns if available, otherwise calculate
    if has_logp_col:
        logp_values = valid_df['log_p'].dropna().tolist()
    else:
        logp_values = []
        
    if has_sas_col:
        sa_scores = valid_df['sas'].dropna().tolist()
    else:
        sa_scores = []
        
    if has_qed_col:
        qed_values = valid_df['qed'].dropna().tolist()
    else:
        qed_values = []
    
    # Calculate missing properties and Lipinski descriptors
    for i, mol in enumerate(valid_mols):
        try:
            # Always calculate MW, HBA, HBD for Lipinski analysis
            MW = Descriptors.MolWt(mol)
            HBA = Descriptors.NOCount(mol)
            HBD = Descriptors.NHOHCount(mol)
            
            # Use existing LogP if available, otherwise calculate
            if has_logp_col and i < len(logp_values):
                LogP = logp_values[i]
            else:
                LogP = Descriptors.MolLogP(mol)
                logp_values.append(LogP)
            
            # Use existing SA score if available, otherwise calculate
  

            sa_score = sa_scores[i]

            
            # Use existing QED if available, otherwise calculate
            if has_qed_col and i < len(qed_values):
                qed_score = qed_values[i]
            else:
                from rdkit.Chem import QED
                qed_score = QED.qed(mol)
                qed_values.append(qed_score)
            
            # Check Ro5 conditions
            conditions = [MW <= 500, HBA <= 10, HBD <= 5, LogP <= 5]
            violations = 4 - conditions.count(True)
            pass_ro5 = violations <= 1
            
            lipinski_results.append({
                'molecular_weight': MW,
                'hba': HBA,
                'hbd': HBD,
                'logp': LogP,
                'sa_score': sa_score,
                'qed_score': qed_score,
                'pass_ro5': pass_ro5,
                'violations': violations
            })
        except Exception as e:
            print(f"Error calculating descriptors for molecule {i}: {e}")
            continue
    
    if lipinski_results:
        results['lipinski'] = {
            'mean_molecular_weight': np.mean([r['molecular_weight'] for r in lipinski_results]),
            'mean_logp': np.mean([r['logp'] for r in lipinski_results]),
            'mean_hbd': np.mean([r['hbd'] for r in lipinski_results]),
            'mean_hba': np.mean([r['hba'] for r in lipinski_results]),
            'mean_violations': np.mean([r['violations'] for r in lipinski_results]),
            'molecules_passing_ro5': sum(1 for r in lipinski_results if r['pass_ro5']),
            'ro5_pass_rate': sum(1 for r in lipinski_results if r['pass_ro5']) / len(lipinski_results),
        }
        
        results['synthetic_accessibility'] = {
            'mean_sa_score': np.mean(sa_scores),
            'median_sa_score': np.median(sa_scores),
            'std_sa_score': np.std(sa_scores),
            'min_sa_score': np.min(sa_scores),
            'max_sa_score': np.max(sa_scores),
            'data_source': 'existing_column' if has_sas_col else 'calculated'
        }
        
        results['logp_analysis'] = {
            'mean_logp': np.mean(logp_values),
            'median_logp': np.median(logp_values),
            'std_logp': np.std(logp_values),
            'min_logp': np.min(logp_values),
            'max_logp': np.max(logp_values),
            'data_source': 'existing_column' if has_logp_col else 'calculated'
        }
        
        results['qed_analysis'] = {
            'mean_qed': np.mean(qed_values),
            'median_qed': np.median(qed_values),
            'std_qed': np.std(qed_values),
            'min_qed': np.min(qed_values),
            'max_qed': np.max(qed_values),
            'molecules_qed_above_0_5': sum(1 for q in qed_values if q > 0.5),
            'qed_above_0_5_rate': sum(1 for q in qed_values if q > 0.5) / len(qed_values),
            'data_source': 'existing_column' if has_qed_col else 'calculated'
        }
    
    # 5. Docking Score Analysis
    if docking_score_col in molecules_df.columns:
        valid_df = molecules_df.iloc[valid_indices]
        docking_scores = valid_df[docking_score_col].dropna()
        
        if len(docking_scores) > 0:
            results['docking'] = {
                'mean_docking_score': docking_scores.mean(),
                'median_docking_score': docking_scores.median(),
                'std_docking_score': docking_scores.std(),
                'min_docking_score': docking_scores.min(),
                'max_docking_score': docking_scores.max(),
                'molecules_below_minus_8': (docking_scores < -8).sum(),
                'ratio_below_minus_8': (docking_scores < -8).mean(),
                'molecules_below_minus_7': (docking_scores < -7).sum(),
                'ratio_below_minus_7': (docking_scores < -7).mean(),
                'total_scored_molecules': len(docking_scores)
            }
    
    # 6. Similarity Analysis (if reference molecules provided)
    if reference_molecules is not None:
        try:
            # Calculate Morgan fingerprints
            ref_mols = [Chem.MolFromSmiles(smiles) for smiles in reference_molecules]
            ref_fps = []
            for mol in ref_mols:
                if mol is not None:
                    ref_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
            
            gen_fps = []
            for mol in valid_mols:
                gen_fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024))
            
            if ref_fps and gen_fps:
                # Convert to numpy arrays
                ref_fps_array = np.array([list(fp) for fp in ref_fps])
                gen_fps_array = np.array([list(fp) for fp in gen_fps])
                
                # Calculate average Tanimoto similarity
                similarities = average_agg_tanimoto(ref_fps_array, gen_fps_array, agg='max')
                
                results['similarity'] = {
                    'mean_max_similarity': np.mean(similarities),
                    'median_max_similarity': np.median(similarities),
                    'std_max_similarity': np.std(similarities),
                    'min_similarity': np.min(similarities),
                    'max_similarity': np.max(similarities)
                }
        except Exception as e:
            print(f"Error calculating similarity: {e}")
    
    # 7. FCD Analysis
    if reference_molecules:
        try:
            fcd_results = calculate_fcd_scores(valid_smiles, reference_molecules)
            results.update(fcd_results)
        except Exception as e:
            print(f"Error calculating FCD scores: {e}")
    
    return results


def average_agg_tanimoto(stock_vecs, gen_vecs, batch_size=5000, agg='max', device='cpu', p=1):
    """
    For each molecule in gen_vecs finds closest molecule in stock_vecs.
    Returns average tanimoto score between these molecules.
    
    Parameters:
        stock_vecs: numpy array <n_vectors x dim>
        gen_vecs: numpy array <n_vectors' x dim>
        agg: 'max' or 'mean'
        p: power for averaging: (mean x^p)^(1/p)
    """
    assert agg in ['max', 'mean'], "Can aggregate only max or mean"
    agg_tanimoto = np.zeros(len(gen_vecs))
    total = np.zeros(len(gen_vecs))
    
    for j in range(0, stock_vecs.shape[0], batch_size):
        x_stock = torch.tensor(stock_vecs[j:j + batch_size]).to(device).float()
        for i in range(0, gen_vecs.shape[0], batch_size):
            y_gen = torch.tensor(gen_vecs[i:i + batch_size]).to(device).float()
            y_gen = y_gen.transpose(0, 1)
            tp = torch.mm(x_stock, y_gen)
            jac = (tp / (x_stock.sum(1, keepdim=True) +
                         y_gen.sum(0, keepdim=True) - tp)).cpu().numpy()
            jac[np.isnan(jac)] = 1
            if p != 1:
                jac = jac**p
            if agg == 'max':
                agg_tanimoto[i:i + y_gen.shape[1]] = np.maximum(
                    agg_tanimoto[i:i + y_gen.shape[1]], jac.max(0))
            elif agg == 'mean':
                agg_tanimoto[i:i + y_gen.shape[1]] += jac.sum(0)
                total[i:i + y_gen.shape[1]] += jac.shape[0]
                
    if agg == 'mean':
        agg_tanimoto /= total
    if p != 1:
        agg_tanimoto = (agg_tanimoto)**(1/p)
    
    return agg_tanimoto


def print_molecular_analysis_report(results):
    """
    Print a formatted report of molecular analysis results.
    
    Args:
        results (dict): Results from calculate_molecular_metrics function
    """
    print("=" * 60)
    print("MOLECULAR ANALYSIS REPORT")
    print("=" * 60)
    
    # Validity
    if 'validity' in results:
        v = results['validity']
        print(f"\n📊 VALIDITY ANALYSIS:")
        print(f"  Total molecules: {v['total_molecules']}")
        print(f"  Valid SMILES: {v['valid_molecules']}")
        print(f"  Invalid SMILES: {v['invalid_molecules']}")
        print(f"  Validity rate: {v['validity_rate']:.2%}")
    
    # Uniqueness
    if 'uniqueness' in results:
        u = results['uniqueness']
        print(f"\n🔄 UNIQUENESS ANALYSIS:")
        print(f"  Total valid molecules: {u['total_valid']}")
        print(f"  Unique molecules: {u['unique_molecules']}")
        print(f"  Duplicate molecules: {u['duplicate_molecules']}")
        print(f"  Uniqueness rate: {u['uniqueness_rate']:.2%}")
    
    # Novelty
    if 'novelty' in results:
        n = results['novelty']
        print(f"\n🆕 NOVELTY ANALYSIS:")
        print(f"  Reference molecules: {n['reference_molecules']}")
        print(f"  Novel molecules: {n['novel_molecules']}")
        print(f"  Known molecules: {n['known_molecules']}")
        print(f"  Novelty rate: {n['novelty_rate']:.2%}")
    
    # Lipinski Rule of 5
    if 'lipinski' in results:
        l = results['lipinski']
        print(f"\n💊 LIPINSKI RULE OF 5:")
        print(f"  Mean molecular weight: {l['mean_molecular_weight']:.2f} Da")
        print(f"  Mean LogP: {l['mean_logp']:.2f}")
        print(f"  Mean HBD: {l['mean_hbd']:.2f}")
        print(f"  Mean HBA: {l['mean_hba']:.2f}")
        print(f"  Mean violations: {l['mean_violations']:.2f}")
        print(f"  Molecules passing Ro5: {l['molecules_passing_ro5']}")
        print(f"  Ro5 pass rate: {l['ro5_pass_rate']:.2%}")
    
    # Synthetic Accessibility
    if 'synthetic_accessibility' in results:
        sa = results['synthetic_accessibility']
        print(f"\n🧪 SYNTHETIC ACCESSIBILITY:")
        print(f"  Mean SA score: {sa['mean_sa_score']:.2f}")
        print(f"  Median SA score: {sa['median_sa_score']:.2f}")
        print(f"  SA score range: {sa['min_sa_score']:.2f} - {sa['max_sa_score']:.2f}")
        print(f"  Data source: {sa['data_source']}")
        print(f"  (Lower scores = easier to synthesize)")
    
    # LogP Analysis
    if 'logp_analysis' in results:
        logp = results['logp_analysis']
        print(f"\n🧊 LOGP ANALYSIS:")
        print(f"  Mean LogP: {logp['mean_logp']:.2f}")
        print(f"  Median LogP: {logp['median_logp']:.2f}")
        print(f"  LogP range: {logp['min_logp']:.2f} - {logp['max_logp']:.2f}")
        print(f"  Data source: {logp['data_source']}")
    
    # QED Analysis
    if 'qed_analysis' in results:
        qed = results['qed_analysis']
        print(f"\n💎 QED (DRUG-LIKENESS) ANALYSIS:")
        print(f"  Mean QED: {qed['mean_qed']:.3f}")
        print(f"  Median QED: {qed['median_qed']:.3f}")
        print(f"  QED range: {qed['min_qed']:.3f} - {qed['max_qed']:.3f}")
        print(f"  Molecules QED > 0.5: {qed['molecules_qed_above_0_5']} ({qed['qed_above_0_5_rate']:.2%})")
        print(f"  Data source: {qed['data_source']}")
        print(f"  (Higher QED = more drug-like)")
    
    # Docking Analysis
    if 'docking' in results:
        d = results['docking']
        print(f"\n🎯 DOCKING ANALYSIS:")
        print(f"  Mean docking score: {d['mean_docking_score']:.2f}")
        print(f"  Median docking score: {d['median_docking_score']:.2f}")
        print(f"  Score range: {d['min_docking_score']:.2f} - {d['max_docking_score']:.2f}")
        print(f"  Molecules < -8: {d['molecules_below_minus_8']} ({d['ratio_below_minus_8']:.2%})")
        print(f"  Molecules < -7: {d['molecules_below_minus_7']} ({d['ratio_below_minus_7']:.2%})")
    
    # Similarity
    if 'similarity' in results:
        s = results['similarity']
        print(f"\n🔍 SIMILARITY TO REFERENCE:")
        print(f"  Mean max similarity: {s['mean_max_similarity']:.3f}")
        print(f"  Median max similarity: {s['median_max_similarity']:.3f}")
        print(f"  Similarity range: {s['min_similarity']:.3f} - {s['max_similarity']:.3f}")
    
    # FCD Analysis
    if 'fcd_score' in results:
        print(f"\n📏 FCD (FRÉCHET CHEMNET DISTANCE):")
        print(f"  FCD score: {results['fcd_score']:.3f}")
        print(f"  (Lower FCD = more similar distribution)")
    
    print("\n" + "=" * 60)


def analyze_molecules_from_file(file_path, smiles_col='smiles', docking_col='docking_score', 
                               reference_smiles=None, save_report=True):
    """
    Convenience function to analyze molecules from a CSV file.
    
    Args:
        file_path (str): Path to CSV file containing molecules
        smiles_col (str): Column name for SMILES strings
        docking_col (str): Column name for docking scores
        reference_smiles (list): List of reference SMILES for novelty calculation
        save_report (bool): Whether to save a detailed report to file
    
    Returns:
        dict: Analysis results
    """
    try:
        df = pd.read_csv(file_path)
        results = calculate_molecular_metrics(df, reference_smiles, docking_col, smiles_col)
        
        print_molecular_analysis_report(results)
        
        if save_report:
            report_path = file_path.replace('.csv', '_analysis_report.txt')
            with open(report_path, 'w') as f:
                # Redirect print output to file
                import sys
                original_stdout = sys.stdout
                sys.stdout = f
                print_molecular_analysis_report(results)
                sys.stdout = original_stdout
            print(f"\nDetailed report saved to: {report_path}")
        
        return results
        
    except Exception as e:
        print(f"Error analyzing molecules from file: {e}")
        return None

def get_device():
    """Get the best available device for computation"""
    if torch.cuda.is_available():
        return 'cuda'
    elif hasattr(torch.backends, 'mps') and torch.backends.mps.is_available():
        return 'mps'
    else:
        return 'cpu'

