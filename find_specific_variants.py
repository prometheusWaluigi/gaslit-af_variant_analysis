#!/usr/bin/env python3
"""
Find Specific Variants Script

This script specifically targets the variants mentioned in the gene variant map
using pysam for advanced variant processing with Intel oneAPI acceleration.
"""

import os
import sys
import logging
import pandas as pd
from pathlib import Path
import argparse
from rich.console import Console
from rich.logging import RichHandler

# Import modular components
from src.gaslit_af.device import initialize_device
from src.gaslit_af.advanced_variant_processing import process_vcf_with_pysam, VariantProcessor, KNOWN_SNPS

# Configure logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af")

# Define specific variants of interest from the gene variant map
SPECIFIC_VARIANTS = {
    # Cognition & Brain Function
    "CHRM2": ["rs8191992", "rs2350780"],
    "DRD2": ["rs6277"],
    "TFAM": ["rs1937"],
    "BCL2": ["rs956572"],
    "MAPRE1": [],  # Add specific SNPs when available
    
    # Sleep Traits
    "ADA": ["rs73598374"],
    "CHRM2": ["rs8191992"],  # Also listed under cognition
    
    # Rare & Neurological Conditions
    "ADGRV1": ["rs575602255", "rs555466095"],
    "C19orf12": ["rs146170087"],
    "PRSS1": ["rs202003805", "rs1232891794"],
    "ATM": ["rs531617441"]
}

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Find Specific Variants in VCF File")
    parser.add_argument("vcf_path", help="Path to VCF file")
    parser.add_argument("--output-dir", type=str, default="output", help="Output directory")
    parser.add_argument("--dbsnp-path", type=str, help="Path to dbSNP VCF file for rsID mapping")
    parser.add_argument("--batch-size", type=int, default=4000000, help="Batch size for processing")
    parser.add_argument("--threads", type=int, default=16, help="Number of worker threads")
    
    args = parser.parse_args()
    args.output_dir = Path(args.output_dir)
    return args

def main():
    """Main entry point for finding specific variants."""
    args = parse_args()
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize device queue
    queue = initialize_device()
    
    # Print configuration
    console.print("GASLIT-AF Specific Variant Finder")
    console.print("Configuration:")
    console.print(f"  VCF File: {args.vcf_path}")
    console.print(f"  Output Directory: {args.output_dir}")
    console.print(f"  Device: {queue.sycl_device}")
    console.print(f"  Target Variants: {sum(len(snps) for snps in SPECIFIC_VARIANTS.values())}")
    console.print("")
    
    # Process VCF file with pysam
    log.info(f"🔍 Searching for specific variants in: {args.vcf_path}")
    
    # Create variant processor
    processor = VariantProcessor(queue=queue, threads=args.threads)
    
    # Load dbSNP data if available
    if args.dbsnp_path:
        processor.load_rsid_map(args.dbsnp_path)
    
    # Get the set of target genes
    target_genes = set(SPECIFIC_VARIANTS.keys())
    
    # Process VCF file
    gene_counts, variant_df = process_vcf_with_pysam(
        vcf_path=args.vcf_path,
        target_genes=target_genes,
        dbsnp_path=args.dbsnp_path,
        queue=queue,
        threads=args.threads,
        batch_size=args.batch_size
    )
    
    if variant_df.empty:
        log.warning("No variants found matching the specific targets")
        return
    
    # Filter for specific variants
    specific_rsids = []
    for gene, snps in SPECIFIC_VARIANTS.items():
        specific_rsids.extend(snps)
    
    specific_variants = variant_df[variant_df['rsid'].isin(specific_rsids)]
    
    if specific_variants.empty:
        log.warning("None of the specific target variants were found")
    else:
        log.info(f"Found {len(specific_variants)} specific target variants")
        
        # Save specific variants to CSV
        specific_variants_path = args.output_dir / "specific_variants.csv"
        specific_variants.to_csv(specific_variants_path, index=False)
        log.info(f"Saved specific variants to: {specific_variants_path}")
        
        # Generate variant report
        variant_report_path = processor.generate_variant_report(specific_variants, args.output_dir)
        log.info(f"Generated variant report: {variant_report_path}")
        
        # Print summary
        console.print("\n[bold green]Specific Variants Found:[/]")
        for _, row in specific_variants.iterrows():
            gene = row.get('gene', '')
            rsid = row.get('rsid', '')
            genotype = row.get('genotype', '')
            chrom = row.get('chrom', '')
            pos = row.get('pos', '')
            
            console.print(f"  [bold]{gene}[/] - {rsid} ({genotype}) at {chrom}:{pos}")

if __name__ == "__main__":
    main()
