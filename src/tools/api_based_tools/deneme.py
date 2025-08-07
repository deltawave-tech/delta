import sys
import os

# Get the base directory from the current file location (to make paths relative)
# --- Add project root to sys.path ---
project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
print(project_root)
sys.path.insert(0, project_root)
# --- End of sys.path modification ---
print(project_root)

from src.tools.api_based_tools.literature_api import search_uniprot
from src.tools.tool_definitions import tool_registry


tool_registry.execute(tool_name="search_uniprot", arguments={"query": "SLC25A44"}, run_id="SLC25A44")