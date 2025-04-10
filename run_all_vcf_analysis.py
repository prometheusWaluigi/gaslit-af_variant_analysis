#!/usr/bin/env python3
"""
Script to run the GASLIT-AF analysis on all VCF files in the data folder.
"""

import os
import sys
import argparse
from pathlib import Path
import logging
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn

# Import the workflow module
try:
    from src.gaslit_af.workflow import run_analysis_workflow
except ImportError:
    # If running from the same directory, try relative import
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from src.gaslit_af.workflow import run_analysis_workflow

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af-batch")

def parse_args():
    """Parse command-line arguments for batch processing."""
    parser = argparse.ArgumentParser(description="GASLIT-AF Batch Variant Analysis")
    
    # Input/output arguments
    parser.add_argument("--data-dir", type=str, default="data", 
                        help="Directory containing VCF files")
    parser.add_argument("--output-base-dir", type=str, default="analysis_results", 
                        help="Base output directory for results")
    parser.add_argument("--file-pattern", type=str, default="*.vcf.gz", 
                        help="File pattern to match VCF files")
    
    # Performance tuning arguments
    parser.add_argument("--batch-size", type=int, default=2000000, 
                        help="Batch size for processing")
    parser.add_argument("--max-ram", type=int, default=64, 
                        help="Maximum RAM usage in GB")
    parser.add_argument("--ram-buffer", type=int, default=16, 
                        help="RAM buffer in GB")
    parser.add_argument("--threads", type=int, default=16, 
                        help="Number of worker threads")
    parser.add_argument("--sample-limit", type=int, 
                        help="Limit number of variants to process (for testing)")
    
    # Caching arguments
    parser.add_argument("--cache-dir", default="./cache", 
                        help="Directory for caching intermediate results")
    parser.add_argument("--cache-max-age", type=int, default=24, 
                        help="Maximum age of cache in hours")
    parser.add_argument("--no-cache", action="store_true", 
                        help="Disable caching")
    
    # Analysis options
    parser.add_argument("--system-analysis", action="store_true", default=True,
                        help="Perform biological system analysis")
    parser.add_argument("--use-pysam", action="store_true", 
                        help="Use pysam for advanced variant processing")
    parser.add_argument("--use-streaming", action="store_true",
                        help="Use streaming processor for very large VCF files (optimized memory usage)")
    parser.add_argument("--dbsnp-path", type=str, 
                        help="Path to dbSNP VCF file for rsID mapping")
    parser.add_argument("--known-variants-only", action="store_true", 
                        help="Only process known variants from the gene variant map")
    parser.add_argument("--clinical-data", type=str,
                        help="Path to clinical variant data JSON file for clinical annotations")
    
    # API annotation options
    parser.add_argument("--api-annotation", action="store_true",
                        help="Enable annotation with external APIs (Ensembl, MyVariant.info)")
    parser.add_argument("--api-sources", nargs="+", default=["ensembl", "myvariant"],
                        help="API sources to use for annotation (default: ensembl myvariant)")
    parser.add_argument("--api-cache-dir", type=str, default="./cache/api",
                        help="Directory to cache API responses")
    parser.add_argument("--api-cache-ttl", type=int, default=24,
                        help="Cache time-to-live in hours for API responses")
                        
    # ANNOVAR annotation options
    parser.add_argument("--use-annovar", action="store_true",
                        help="Enable annotation with ANNOVAR for advanced functional annotations")
    parser.add_argument("--annovar-path", type=str,
                        help="Path to ANNOVAR installation directory")
    parser.add_argument("--humandb-path", type=str,
                        help="Path to ANNOVAR humandb directory")
    
    # Variant enrichment options
    parser.add_argument("--variant-enrichment", action="store_true",
                        help="Enable advanced variant enrichment with AF-specific annotations")
    parser.add_argument("--enrichment-cache-dir", type=str, default="./cache/enrichment",
                        help="Directory to cache enrichment data")
    parser.add_argument("--enrichment-cache-ttl", type=int, default=24,
                        help="Cache time-to-live in hours for enrichment data")
    
    # Output options
    parser.add_argument("--no-visualization", action="store_true", 
                        help="Skip visualization generation")
    parser.add_argument("--no-report", action="store_true", 
                        help="Skip HTML report generation")
    parser.add_argument("--enhanced-report", action="store_true", 
                        help="Generate enhanced report with interactive visualizations and symptom correlations")
    parser.add_argument("--include-symptoms", action="store_true",
                        help="Include symptom correlations in enhanced report")
    parser.add_argument("--open-browser", action="store_true", 
                        help="Automatically open report in browser")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert string paths to Path objects
    args.data_dir = Path(args.data_dir)
    args.output_base_dir = Path(args.output_base_dir)
    args.cache_dir = Path(args.cache_dir)
    
    return args

def main():
    """Main function to run the analysis on all VCF files."""
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
            output_dir = args.output_base_dir / vcf_file.stem.replace('.vcf', '')
            output_dir.mkdir(parents=True, exist_ok=True)
            
            log.info(f"\n{'='*80}")
            log.info(f"Processing {vcf_file}")
            log.info(f"Output directory: {output_dir}")
            log.info(f"{'='*80}\n")
            
            # Create a copy of the args with the current VCF file and output directory
            class Args:
                pass
            
            current_args = Args()
            for key, value in vars(args).items():
                setattr(current_args, key, value)
            
            current_args.vcf_path = vcf_file
            current_args.output_dir = output_dir
            
            try:
                # Run the analysis workflow
                run_analysis_workflow(current_args)
                log.info(f"Analysis completed for {vcf_file}")
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
