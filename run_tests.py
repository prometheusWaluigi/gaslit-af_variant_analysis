#!/usr/bin/env python3
"""
Test runner for GASLIT-AF Variant Analysis.

This script provides a convenient way to run the test suite with different options.
"""

import os
import sys
import argparse
import subprocess
from pathlib import Path

def parse_args():
    """Parse command-line arguments for the test runner."""
    parser = argparse.ArgumentParser(description="GASLIT-AF Variant Analysis Test Runner")
    
    parser.add_argument("--unit", action="store_true", help="Run unit tests only")
    parser.add_argument("--integration", action="store_true", help="Run integration tests only")
    parser.add_argument("--all", action="store_true", help="Run all tests (default)")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    parser.add_argument("--specific", type=str, help="Run a specific test module or function")
    parser.add_argument("--skip-slow", action="store_true", help="Skip slow tests")
    parser.add_argument("--clinical", action="store_true", help="Run clinical variant tests only")
    parser.add_argument("--api", action="store_true", help="Run API integration tests only")
    parser.add_argument("--mock-api", action="store_true", help="Use mocked API responses for tests")
    
    args = parser.parse_args()
    
    # If no test type is specified, run all tests
    if not (args.unit or args.integration or args.specific):
        args.all = True
    
    return args

def run_tests(args):
    """Run the tests based on the specified options."""
    # Base pytest command
    cmd = ["pytest"]
    
    # Add verbose flag if requested
    if args.verbose:
        cmd.append("-v")
    
    # Add test selection
    if args.unit:
        cmd.append("tests/unit/")
    elif args.integration:
        cmd.append("tests/integration/")
    elif args.specific:
        cmd.append(args.specific)
    elif args.clinical:
        cmd.append("tests/unit/test_clinical_variants.py")
        cmd.append("tests/unit/test_clinical_integration.py")
        cmd.append("tests/integration/test_clinical_api_integration.py")
    elif args.api:
        cmd.append("tests/unit/test_api_integration.py")
        cmd.append("tests/integration/test_clinical_api_integration.py")
    elif args.all:
        cmd.append("tests/")
    
    # Skip slow tests if requested
    if args.skip_slow:
        cmd.append("-m")
        cmd.append("not slow")
    
    # Use mocked API responses for tests
    if args.mock_api:
        os.environ["USE_MOCK_API"] = "1"
    
    # Print the command being run
    print(f"Running: {' '.join(cmd)}")
    
    # Run the tests
    result = subprocess.run(cmd)
    
    return result.returncode

def main():
    """Main entry point for the test runner."""
    args = parse_args()
    
    # Ensure we're in the project root directory
    project_root = Path(__file__).parent
    os.chdir(project_root)
    
    # Print available test modules
    if args.verbose:
        print("Available test modules:")
        unit_tests = list(Path("tests/unit").glob("test_*.py"))
        integration_tests = list(Path("tests/integration").glob("test_*.py"))
        
        print("Unit tests:")
        for test in unit_tests:
            print(f"  - {test.stem}")
        
        print("Integration tests:")
        for test in integration_tests:
            print(f"  - {test.stem}")
        
        print("")
    
    # Run the tests
    return run_tests(args)

if __name__ == "__main__":
    sys.exit(main())
