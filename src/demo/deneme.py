import os
import sys
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
print(project_root)
sys.path.insert(0, project_root)

from src.tools.genmol_vina import get_vina_genmol_report
from src.tools.tool_definitions import ANTHROPIC_TOOLS, ToolNames, tool_registry
from src.agent import Agent
from src.tools.prot2mol_vina import get_vina_mol_gen_report
from src.tools.vina_plip import get_vina_report
from src.tools.api_based_tools.literature_api import search_uniprot
from src.tools.api_based_tools.pdb_search import get_pdb_file, get_alphafold_file
from src.tools.api_based_tools.fetch_chembl import search_chembl_activity
from src.tools.api_based_tools.api_tools import run_diffdock
from src.tools.prot2mol import generate_molecules
from src.tools.genmol_vina import get_vina_genmol_report
import pandas as pd
from src.tools.genmol import genmol
import time
start_time = time.time()

vina_akt1_args = {
    "protein_path": "/Users/atabeyunlu/multi_agent_paper/runs_metadata/0618_2215_GV4_sonnet-3.7_3/pdb_files/4EJN.pdb"
}

print("Running Vina analysis...")
vina_akt1_res = tool_registry.execute(tool_name="vina_genmol", arguments=vina_akt1_args, run_id="genmol_rep_test")
