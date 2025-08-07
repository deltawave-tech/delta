#!/bin/bash
set -e

# Use Python 3.9 explicitly
PYTHON=/opt/conda/bin/python3.9

echo "Verifying Python environment using $($PYTHON --version)"

# Verify we can import both RDKit and Meeko from conda
$PYTHON -c "
import sys
import os
print('Python version:', sys.version)
print('Python path:', sys.path)

try:
    from rdkit import Chem
    print('✅ RDKit successfully imported from conda')
except ImportError as e:
    print('❌ RDKit import failed:', e)
    exit(1)

try:
    from meeko.preparation import PDBQTWriterLegacy, MoleculePreparation
    print('✅ Meeko successfully imported from conda')
except ImportError as e:
    print('❌ Failed to import Meeko:', e)
    exit(1)

try:
    import fastapi
    print(f'✅ FastAPI version: {fastapi.__version__}')
except ImportError as e:
    print('❌ FastAPI import failed:', e)
    exit(1)

try:
    import pydantic
    print(f'✅ Pydantic version: {pydantic.__version__}')
except ImportError as e:
    print('❌ Pydantic import failed:', e)
    exit(1)
"

# Print package info for debugging
echo "Installed packages:"
echo "Conda packages:"
conda list | grep -E 'meeko|rdkit'
echo "Pip packages:"
/opt/conda/bin/pip freeze | grep -E 'fastapi|pydantic|uvicorn|starlette'

echo "Python environment verification completed"