#!/usr/bin/env python3
"""
Benchmark script for GASLIT-AF Variant Analysis.

This script benchmarks the performance of the variant analysis pipeline
on different VCF files and checks for specific variants of interest.
"""

import os
import sys
import time
import json
import pandas as pd
import logging
from pathlib import Path
from datetime import datetime
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn, TimeRemainingColumn

# Add project root to Python path
sys.path.insert(0, str(Path(__file__).parent))

# Import GASLIT-AF modules
from src.gaslit_af.cli import parse_args
from src.gaslit_af.device import initialize_device
from src.gaslit_af.advanced_variant_processing import VariantProcessor, process_vcf_with_pysam
from src.gaslit_af.gene_lists import GASLIT_AF_GENES, KNOWN_SNPS

# Configure logging
logging.basicConfig(level=logging.INFO, 
                   format='[%(asctime)s] %(levelname)-8s %(message)s',
                   datefmt='%Y-%m-%d %H:%M:%S')
log = logging.getLogger("gaslit-af-benchmark")

# Rich console for pretty output
console = Console()

# Define the specific variants we're looking for
TARGET_VARIANTS = [
    # Cognition & Brain Function
    {"gene": "CHRM2", "rsid": "rs8191992", "genotype": "TT", "trait": "Executive function, memory, attention"},
    {"gene": "CHRM2", "rsid": "rs2350780", "genotype": "AA", "trait": "Executive function, memory, attention"},
    {"gene": "DRD2", "rsid": "rs6277", "genotype": "AG", "trait": "Dopamine modulation, cognitive flexibility"},
    {"gene": "TFAM", "rsid": "rs1937", "genotype": "GG", "trait": "Mitochondrial efficiency, energy for brain cells"},
    {"gene": "BCL2", "rsid": "rs956572", "genotype": "GG", "trait": "Neuroprotection, stress resilience"},
    {"gene": "ST8SIA6", "rsid": "multiple", "genotype": "TT", "trait": "Enhanced neuronal growth and connectivity"},
    {"gene": "CHRNA5", "rsid": "multiple", "genotype": "CT", "trait": "Enhanced neuronal growth and connectivity"},
    {"gene": "NRG1", "rsid": "multiple", "genotype": "varies", "trait": "Enhanced neuronal growth and connectivity"},
    {"gene": "MAPRE1", "rsid": "multiple", "genotype": "varies", "trait": "Enhanced neuronal growth and connectivity"},
    {"gene": "GYPC", "rsid": "multiple", "genotype": "varies", "trait": "Enhanced neuronal growth and connectivity"},
    {"gene": "CABP5", "rsid": "multiple", "genotype": "varies", "trait": "Enhanced neuronal growth and connectivity"},
    
    # Sleep Traits
    {"gene": "ADA", "rsid": "rs73598374", "genotype": "TC", "trait": "Deep sleep, longer delta wave cycles"},
    {"gene": "VRK1", "rsid": "multiple", "genotype": "GT", "trait": "Higher sleep quality"},
    {"gene": "CHRM2", "rsid": "multiple", "genotype": "GG", "trait": "Higher sleep quality"},
    {"gene": "RALYL", "rsid": "multiple", "genotype": "varies", "trait": "Higher sleep quality"},
    {"gene": "FOXO6", "rsid": "multiple", "genotype": "varies", "trait": "Higher sleep quality"},
    
    # Rare Conditions
    {"gene": "ADGRV1", "rsid": "rs575602255", "genotype": "AG", "trait": "Usher Syndrome II (vision + hearing)"},
    {"gene": "ADGRV1", "rsid": "rs555466095", "genotype": "CG", "trait": "Usher Syndrome II (vision + hearing)"},
    {"gene": "C19orf12", "rsid": "rs146170087", "genotype": "TC", "trait": "NBIA-4 (Neurodegeneration with motor + cognitive decline)"},
    {"gene": "PRSS1", "rsid": "rs202003805", "genotype": "CT", "trait": "Hereditary Pancreatitis"},
    {"gene": "PRSS1", "rsid": "rs1232891794", "genotype": "GC", "trait": "Hereditary Pancreatitis"},
    {"gene": "ATM", "rsid": "rs531617441", "genotype": "AG", "trait": "Increased DNA repair-related cancer susceptibility"}
]

def run_benchmark(vcf_files, thread_counts, batch_sizes, output_dir):
    """
    Run the benchmark on multiple VCF files with different thread counts and batch sizes.
    
    Args:
        vcf_files: List of VCF files to benchmark
        thread_counts: List of thread counts to test
        batch_sizes: List of batch sizes to test
        output_dir: Directory to save benchmark results
    """
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Initialize results dictionary
    results = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "system_info": {
            "cpu": os.cpu_count(),
            "device": str(initialize_device().sycl_device)
        },
        "benchmarks": []
    }
    
    # Create target genes set from the variants we're looking for
    target_genes = {variant["gene"] for variant in TARGET_VARIANTS if variant["gene"] != "multiple"}
    
    # Initialize variant processor
    processor = VariantProcessor(threads=max(thread_counts))
    
    # Run benchmarks
    total_benchmarks = len(vcf_files) * len(thread_counts) * len(batch_sizes)
    benchmark_index = 0
    
    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console
    ) as progress:
        task = progress.add_task("[cyan]Running benchmarks...", total=total_benchmarks)
        
        for vcf_file in vcf_files:
            for thread_count in thread_counts:
                for batch_size in batch_sizes:
                    benchmark_index += 1
                    vcf_path = Path(vcf_file)
                    
                    if not vcf_path.exists():
                        log.warning(f"VCF file not found: {vcf_path}")
                        progress.update(task, advance=1)
                        continue
                    
                    # Update progress description
                    progress.update(task, description=f"[cyan]Benchmark {benchmark_index}/{total_benchmarks}: {vcf_path.name} with {thread_count} threads, batch size {batch_size}")
                    
                    try:
                        # Initialize device
                        queue = initialize_device()
                        
                        # Start timing
                        start_time = time.time()
                        
                        # Process VCF file
                        gene_counts, variant_df = process_vcf_with_pysam(
                            vcf_path=str(vcf_path),
                            target_genes=target_genes,
                            dbsnp_path=None,
                            queue=queue,
                            threads=thread_count,
                            batch_size=batch_size,
                            max_ram_usage=64,
                            ram_buffer=16
                        )
                        
                        # End timing
                        end_time = time.time()
                        processing_time = end_time - start_time
                        
                        # Check for target variants
                        found_variants = []
                        if variant_df is not None and not variant_df.empty:
                            for variant in TARGET_VARIANTS:
                                if variant["rsid"] == "multiple":
                                    # For variants with multiple SNPs, just check if the gene is present
                                    gene_matches = variant_df[variant_df["gene"] == variant["gene"]]
                                    if not gene_matches.empty:
                                        found_variants.append({
                                            "gene": variant["gene"],
                                            "rsid": "multiple",
                                            "found": True,
                                            "count": len(gene_matches)
                                        })
                                else:
                                    # For specific SNPs, check if the exact variant is present
                                    matches = variant_df[(variant_df["gene"] == variant["gene"]) & 
                                                        (variant_df["rsid"] == variant["rsid"])]
                                    if not matches.empty:
                                        for _, row in matches.iterrows():
                                            found_variants.append({
                                                "gene": variant["gene"],
                                                "rsid": variant["rsid"],
                                                "found": True,
                                                "genotype": row.get("genotype", "unknown"),
                                                "expected_genotype": variant["genotype"]
                                            })
                        
                        # Store benchmark results
                        benchmark_result = {
                            "vcf_file": str(vcf_path),
                            "thread_count": thread_count,
                            "batch_size": batch_size,
                            "processing_time_seconds": processing_time,
                            "variant_count": len(variant_df) if variant_df is not None else 0,
                            "gene_count": len(gene_counts) if gene_counts else 0,
                            "found_variants": found_variants,
                            "target_variants_found": len(found_variants),
                            "target_variants_total": len(TARGET_VARIANTS)
                        }
                        
                        results["benchmarks"].append(benchmark_result)
                        
                        # Log results
                        log.info(f"Benchmark completed: {vcf_path.name} with {thread_count} threads, batch size {batch_size}")
                        log.info(f"Processing time: {processing_time:.2f} seconds")
                        log.info(f"Found {len(found_variants)} target variants out of {len(TARGET_VARIANTS)}")
                        
                    except Exception as e:
                        log.error(f"Error processing {vcf_path.name}: {e}")
                        # Store error in results
                        results["benchmarks"].append({
                            "vcf_file": str(vcf_path),
                            "thread_count": thread_count,
                            "batch_size": batch_size,
                            "error": str(e)
                        })
                    
                    # Update progress
                    progress.update(task, advance=1)
    
    # Save results to JSON
    results_file = output_dir / f"benchmark_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
    with open(results_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    log.info(f"Benchmark results saved to {results_file}")
    
    # Display summary table
    display_benchmark_summary(results)
    
    return results

def display_benchmark_summary(results):
    """
    Display a summary table of benchmark results.
    
    Args:
        results: Benchmark results dictionary
    """
    table = Table(title="GASLIT-AF Variant Analysis Benchmark Summary")
    
    table.add_column("VCF File", style="cyan")
    table.add_column("Threads", style="magenta")
    table.add_column("Batch Size", style="blue")
    table.add_column("Time (s)", style="green")
    table.add_column("Variants Found", style="yellow")
    table.add_column("Target Variants", style="red")
    
    for benchmark in results["benchmarks"]:
        if "error" in benchmark:
            table.add_row(
                Path(benchmark["vcf_file"]).name,
                str(benchmark["thread_count"]),
                str(benchmark["batch_size"]),
                "ERROR",
                "ERROR",
                "ERROR"
            )
        else:
            table.add_row(
                Path(benchmark["vcf_file"]).name,
                str(benchmark["thread_count"]),
                str(benchmark["batch_size"]),
                f"{benchmark['processing_time_seconds']:.2f}",
                str(benchmark["variant_count"]),
                f"{benchmark['target_variants_found']}/{benchmark['target_variants_total']}"
            )
    
    console.print(table)

def main():
    """Main entry point for the benchmark script."""
    # Define VCF files to benchmark
    data_dir = Path("/home/k10/dev/windsage/gaslit-af_variant_analysis/data")
    
    # Focus on the main SNP-indel VCF file for primary benchmark
    vcf_files = [data_dir / "KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.snp-indel.genome.vcf.gz"]
    
    # Add other VCF files if they exist and are small enough for quick testing
    cnv_file = data_dir / "KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.cnv.vcf.gz"
    sv_file = data_dir / "KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.sv.vcf.gz"
    
    if cnv_file.exists() and cnv_file.stat().st_size < 50_000_000:  # Only include if < 50MB
        vcf_files.append(cnv_file)
    
    if sv_file.exists() and sv_file.stat().st_size < 50_000_000:  # Only include if < 50MB
        vcf_files.append(sv_file)
    
    # Define thread counts to test - focus on a smaller range for quicker results
    thread_counts = [16, 32]
    
    # Define batch sizes to test - focus on a smaller range for quicker results
    batch_sizes = [1000000, 2000000]
    
    # Define output directory
    output_dir = Path("./benchmark_results")
    
    # Run benchmark
    run_benchmark(vcf_files, thread_counts, batch_sizes, output_dir)

if __name__ == "__main__":
    main()
