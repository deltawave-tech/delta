#!/bin/bash
set -e

echo "Container starting up..."

# Use Python 3.9 explicitly
export PYTHONPATH=/opt/conda/lib/python3.9/site-packages:/app
export PYTHON=/opt/conda/bin/python3.9

echo "Using Python: $($PYTHON --version)"

# Always verify that Meeko is accessible
echo "Verifying Meeko installation..."
chmod +x /app/install_meeko.sh
/app/install_meeko.sh

# Run the import test to verify all dependencies
echo "Testing all package imports..."
$PYTHON /app/test_imports.py

# Start the API server
echo "Starting API server..."
cd /app
exec /opt/conda/bin/uvicorn autodock_api:app --host 0.0.0.0 --port 8686 --reload