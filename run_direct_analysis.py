#!/usr/bin/env python3
"""
Script to run the direct GASLIT-AF analysis on all VCF files in the data folder.
This script uses the DirectSystemAnalyzer to extract gene information directly from
the ANN field in the VCF files, which is more effective at finding GASLIT-AF genes.
"""

import os
import sys
import argparse
from pathlib import Path
import logging
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn

# Import the DirectSystemAnalyzer
try:
    from direct_system_analysis import DirectSystemAnalyzer
except ImportError:
    # If running from the same directory, try relative import
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from direct_system_analysis import DirectSystemAnalyzer

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af-batch-direct")

def parse_args():
    """Parse command-line arguments for batch processing."""
    parser = argparse.ArgumentParser(description="GASLIT-AF Batch Direct Variant Analysis")
    
    # Input/output arguments
    parser.add_argument("--data-dir", type=str, default="data", 
                        help="Directory containing VCF files")
    parser.add_argument("--output-base-dir", type=str, default="analysis_results", 
                        help="Base output directory for results")
    parser.add_argument("--file-pattern", type=str, default="*.vcf.gz", 
                        help="File pattern to match VCF files")
    
    # Processing arguments
    parser.add_argument("--chunk-size", type=int, default=10000, 
                        help="Number of variants to process at once")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert string paths to Path objects
    args.data_dir = Path(args.data_dir)
    args.output_base_dir = Path(args.output_base_dir)
    
    return args

def main():
    """Main function to run the direct analysis on all VCF files."""
    # Parse command-line arguments
    args = parse_args()
    
    # Find all VCF files in the data directory
    vcf_files = list(args.data_dir.glob(args.file_pattern))
    
    if not vcf_files:
        log.error(f"No VCF files found in {args.data_dir} matching pattern {args.file_pattern}")
        sys.exit(1)
    
    log.info(f"Found {len(vcf_files)} VCF files to analyze:")
    for vcf_file in vcf_files:
        log.info(f"  - {vcf_file}")
    
    # Create the base output directory
    args.output_base_dir.mkdir(parents=True, exist_ok=True)
    
    # Process each VCF file
    with Progress(
        TextColumn("[bold blue]{task.description}"),
        BarColumn(),
        TextColumn("[bold]{task.completed}/{task.total}"),
        TimeElapsedColumn(),
        console=console
    ) as progress:
        task = progress.add_task("[green]Processing VCF files", total=len(vcf_files))
        
        for vcf_file in vcf_files:
            # Create a unique output directory for this VCF file
            output_dir = args.output_base_dir / f"direct_{vcf_file.stem.replace('.vcf', '')}"
            output_dir.mkdir(parents=True, exist_ok=True)
            
            log.info(f"\n{'='*80}")
            log.info(f"Processing {vcf_file}")
            log.info(f"Output directory: {output_dir}")
            log.info(f"{'='*80}\n")
            
            try:
                # Create analyzer
                analyzer = DirectSystemAnalyzer()
                
                # Analyze VCF file
                success = analyzer.analyze_vcf(
                    vcf_path=str(vcf_file),
                    output_dir=str(output_dir),
                    chunk_size=args.chunk_size
                )
                
                if success:
                    log.info(f"Successfully analyzed VCF file: {vcf_file}")
                    log.info(f"Results saved to: {output_dir}")
                else:
                    log.error(f"Failed to analyze VCF file: {vcf_file}")
            except Exception as e:
                log.error(f"Error processing {vcf_file}: {e}")
                import traceback
                log.error(traceback.format_exc())
            
            # Update progress
            progress.update(task, advance=1)
    
    log.info("\n\n")
    log.info(f"{'='*80}")
    log.info(f"All analyses completed!")
    log.info(f"Results saved to: {args.output_base_dir}")
    log.info(f"{'='*80}")

if __name__ == "__main__":
    main()
