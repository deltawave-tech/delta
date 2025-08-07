"""
Test script to verify all required imports work correctly.
This helps diagnose any issues with the Python environment.
"""
import sys
import os
import importlib

def check_import(module_name, package=None):
    """Attempt to import a module and print its details"""
    try:
        if package:
            module = importlib.import_module(module_name, package=package)
        else:
            module = importlib.import_module(module_name)
        
        version = getattr(module, "__version__", "unknown")
        path = getattr(module, "__file__", "unknown")
        
        print(f"✅ Successfully imported {module_name} (version: {version})")
        print(f"   Path: {path}")
        return True
    except ImportError as e:
        print(f"❌ Failed to import {module_name}: {str(e)}")
        return False

def main():
    """Test all required imports"""
    print("=== Python Environment Information ===")
    print(f"Python version: {sys.version}")
    print(f"Python executable: {sys.executable}")
    print("\n=== Testing Meeko Imports ===")
    check_import("meeko")
    check_import("meeko.preparation")
    
    print("\n=== Testing RDKit Imports ===")
    check_import("rdkit")
    check_import("rdkit.Chem")
    
    print("\n=== Testing FastAPI Imports ===")
    check_import("fastapi")
    check_import("pydantic")
    check_import("uvicorn")
    check_import("starlette")
    
    print("\n=== Testing Other Required Imports ===")
    check_import("asyncio")
    check_import("pathlib")
    check_import("tempfile")
    check_import("logging")
    
    print("\n=== Python Path ===")
    for path in sys.path:
        print(f" - {path}")
    
    print("\n=== Environment Variables ===")
    for key, value in os.environ.items():
        if key.startswith("PYTHON") or key.startswith("PATH"):
            print(f"{key}={value}")

if __name__ == "__main__":
    main()