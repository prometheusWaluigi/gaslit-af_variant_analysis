#!/usr/bin/env python3
"""
GASLIT-AF Variant Analysis - Main Script

This script performs genomic variant analysis targeting GASLIT-AF gene clusters
using Intel oneAPI with SYCL for GPU acceleration.

The analysis includes:
- Memory-bounded chunking for large VCF files
- Parallel processing with GPU acceleration
- Biological system-level analysis
- Visualization and reporting
- Caching for improved performance
"""

import logging
from rich.logging import RichHandler
from rich.console import Console

# Import modular components
from src.gaslit_af.cli import parse_args
from src.gaslit_af.workflow import run_analysis_workflow

# Configure logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af")

def main():
    """Main entry point for GASLIT-AF Variant Analysis."""
    # Parse command-line arguments
    args = parse_args()
    
    # Run the analysis workflow
    run_analysis_workflow(args)

if __name__ == "__main__":
    main()
