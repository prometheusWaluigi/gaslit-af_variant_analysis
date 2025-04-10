#!/usr/bin/env python3
"""
Script to run the full GASLIT-AF analysis pipeline:
1. Run direct analysis on all VCF files in the data folder
2. Generate a combined HTML report of all results
"""

import os
import sys
import subprocess
import argparse
from pathlib import Path
import logging
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af-full-pipeline")

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="GASLIT-AF Full Analysis Pipeline")
    
    # Input/output arguments
    parser.add_argument("--data-dir", type=str, default="data", 
                        help="Directory containing VCF files")
    parser.add_argument("--output-dir", type=str, default="analysis_results", 
                        help="Base output directory for results")
    parser.add_argument("--file-pattern", type=str, default="*.vcf.gz", 
                        help="File pattern to match VCF files")
    parser.add_argument("--report-file", type=str, default="gaslit_af_combined_report.html", 
                        help="Output HTML report file")
    
    # Processing arguments
    parser.add_argument("--chunk-size", type=int, default=10000, 
                        help="Number of variants to process at once")
    parser.add_argument("--open-browser", action="store_true", 
                        help="Automatically open report in browser")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert string paths to Path objects
    args.data_dir = Path(args.data_dir)
    args.output_dir = Path(args.output_dir)
    
    return args

def run_direct_analysis(args):
    """Run direct analysis on all VCF files."""
    log.info("Running direct analysis on all VCF files...")
    
    cmd = [
        "python", "run_direct_analysis.py",
        "--data-dir", str(args.data_dir),
        "--output-base-dir", str(args.output_dir),
        "--file-pattern", args.file_pattern,
        "--chunk-size", str(args.chunk_size)
    ]
    
    try:
        subprocess.run(cmd, check=True)
        log.info("Direct analysis completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        log.error(f"Error running direct analysis: {e}")
        return False

def generate_combined_report(args):
    """Generate combined HTML report."""
    log.info("Generating combined HTML report...")
    
    cmd = [
        "python", "generate_combined_report.py",
        "--results-dir", str(args.output_dir),
        "--output-file", args.report_file,
        "--direct-only"
    ]
    
    try:
        subprocess.run(cmd, check=True)
        log.info(f"Combined report generated successfully: {args.report_file}")
        return True
    except subprocess.CalledProcessError as e:
        log.error(f"Error generating combined report: {e}")
        return False

def open_browser(report_file):
    """Open report in browser."""
    log.info(f"Opening report in browser: {report_file}")
    
    if sys.platform == 'darwin':  # macOS
        cmd = ['open', report_file]
    elif sys.platform == 'win32':  # Windows
        cmd = ['start', report_file]
    else:  # Linux
        cmd = ['xdg-open', report_file]
    
    try:
        subprocess.run(cmd, check=True)
        log.info("Browser opened successfully")
        return True
    except subprocess.CalledProcessError as e:
        log.error(f"Error opening browser: {e}")
        return False

def main():
    """Main function."""
    args = parse_args()
    
    log.info("Starting GASLIT-AF full analysis pipeline...")
    
    # Run direct analysis
    if not run_direct_analysis(args):
        log.error("Direct analysis failed, aborting pipeline")
        sys.exit(1)
    
    # Generate combined report
    if not generate_combined_report(args):
        log.error("Combined report generation failed, aborting pipeline")
        sys.exit(1)
    
    # Open browser if requested
    if args.open_browser:
        open_browser(args.report_file)
    
    log.info("GASLIT-AF full analysis pipeline completed successfully")
    log.info(f"Results saved to: {args.output_dir}")
    log.info(f"Combined report: {args.report_file}")

if __name__ == "__main__":
    main()
