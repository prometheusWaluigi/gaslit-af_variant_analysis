import dpctl
import dpctl.tensor as dpt
from cyvcf2 import VCF
from collections import defaultdict
import numpy as np
import sys
import time
import os
import psutil
import traceback
import logging
import pandas as pd
import argparse
import json
import webbrowser
from pathlib import Path
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn, TimeRemainingColumn

# Import custom modules
try:
    from src.gaslit_af.data_processing import vcf_to_dataframe, extract_gaslit_af_variants, save_results
    from src.gaslit_af.visualization import generate_all_visualizations
    from src.gaslit_af.reporting import generate_html_report
    from src.gaslit_af.enhanced_reporting import generate_enhanced_report
    from src.gaslit_af.caching import AnalysisCache
    from src.gaslit_af.biological_systems import analyze_systems, plot_system_distribution, generate_system_summary, BIOLOGICAL_SYSTEMS
    from src.gaslit_af.advanced_variant_processing import process_vcf_with_pysam, VariantProcessor, KNOWN_SNPS
    from src.gaslit_af.streaming import stream_process_vcf
except ImportError:
    # If running from the same directory, try relative import
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from src.gaslit_af.data_processing import vcf_to_dataframe, extract_gaslit_af_variants, save_results
    from src.gaslit_af.visualization import generate_all_visualizations
    from src.gaslit_af.reporting import generate_html_report
    from src.gaslit_af.enhanced_reporting import generate_enhanced_report
    from src.gaslit_af.caching import AnalysisCache
    from src.gaslit_af.biological_systems import analyze_systems, plot_system_distribution, generate_system_summary, BIOLOGICAL_SYSTEMS
    from src.gaslit_af.advanced_variant_processing import process_vcf_with_pysam, VariantProcessor, KNOWN_SNPS
    from src.gaslit_af.streaming import stream_process_vcf

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af")

# Comprehensive GASLIT-AF gene list including new additions - properly handled
GASLIT_AF_GENES_TEXT = """
IDO2 AHR AHRR IL36RN CFH MBL2 NLRP3 IL1B IL6 IL17 IL13 IL4 HLA-DQB1 PTPN22 CTLA4 ASXL1 CBL DNMT3B ETV6 IDH1
COMT CHRM2 DRD2 GABRA1 CHRNA7 ADRB1 ADRB2 NOS3 GNB3 SLC6A2 NET EZH2 SLC6A4 HTR2A TAAR1 OPRM1 GCH1 TRPV2 MYT1L NRXN3
TNXB ADAMTS10 SELENON NEB MYH7 MAPRE1 ADGRV1 PLXNA2 COL3A1 FBN1 FLNA COL5A1 FKBP14 PLOD1
APOE PCSK9 UGT1A1 HNF1A ABCC8 TFAM C19orf12 MT-ATP6 MT-ATP8 PDHA1 SDHB NAMPT NMRK1 PGC1A
CNR1 CNR2 FAAH MGLL
ITPR1 KCNJ5 RYR2
TPSAB1 KIT HNMT TET2
IDO1 KMO KYNU TDO2 HAAO ARNT BECN1 ATG5
ROCK1 ROCK2 ARG1 "Ang-(1-7)" ACE ACE2 "ANG I" "ANG II" TGFβ1 TGFβ2 TGFβ3 GDF-15 "Activin B" Follistatin "Hif-1α"
DRP1 "PINK-1" SIRT1 IFNα IFNβ IFNγ IFNL1 PGE2 "α-NAGA" ATG13 NEFL S100B TWEAK
"""

# Process gene list handling quoted multi-word genes correctly
GASLIT_AF_GENES = set()
in_quotes = False
current_gene = ""
for char in GASLIT_AF_GENES_TEXT:
    if char == '"':
        in_quotes = not in_quotes
        if not in_quotes and current_gene:  # End of quoted gene
            GASLIT_AF_GENES.add(current_gene.strip())
            current_gene = ""
    elif in_quotes:
        current_gene += char
    elif char.isspace():
        if current_gene:
            GASLIT_AF_GENES.add(current_gene.strip())
            current_gene = ""
    else:
        current_gene += char

# Add the last gene if there is one
if current_gene:
    GASLIT_AF_GENES.add(current_gene.strip())

# Initialize device queue for Intel Arc GPU with fallback to CPU
try:
    # Request maximum GPU performance
    os.environ["SYCL_CACHE_PERSISTENT"] = "1"
    
    # Create a queue with GPU selection
    queue = dpctl.SyclQueue("gpu")
    log.info(f"🚀 Using GPU device: {queue.sycl_device}")
    
    # Check if we're actually using the Intel Arc GPU
    if "Arc" in str(queue.sycl_device):
        log.info("✅ Successfully connected to Intel Arc GPU")
except Exception as e:
    log.warning(f"⚠️ Could not initialize GPU device: {e}")
    log.info("⚙️ Falling back to CPU")
    try:
        queue = dpctl.SyclQueue()
        log.info(f"🖥️ Using default device: {queue.sycl_device}")
    except Exception as e2:
        log.error(f"Could not initialize any SYCL device: {e2}")
        sys.exit(1)

# Memory monitoring
def get_memory_usage():
    """Get current memory usage in GB"""
    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    return memory_info.rss / (1024 ** 3)  # Convert to GB

def check_memory_limits(current_usage, max_ram_usage=64, ram_buffer=16):
    """Check if memory usage is approaching limits"""
    if current_usage > (max_ram_usage - ram_buffer):
        log.warning(f"⚠️ Memory usage ({current_usage:.2f} GB) approaching limit ({max_ram_usage} GB)")
        return True
    return False

def analyze_vcf_oneapi(vcf_path, batch_size=2000000, max_ram_usage=64, ram_buffer=16, threads=16):
    """Analyze VCF file for GASLIT-AF gene variants using oneAPI acceleration.
    
    Args:
        vcf_path: Path to VCF file
        batch_size: Batch size for processing
        max_ram_usage: Maximum RAM usage in GB
        ram_buffer: RAM buffer in GB
    
    Returns:
        Dictionary of gene:count pairs
    """
    try:
        log.info(f"🚀 Analyzing (oneAPI accelerated): {vcf_path}")
        
        # Check if file exists
        if not os.path.exists(vcf_path):
            log.error(f"❌ File not found: {vcf_path}")
            return None
        
        # First pass to count total records for progress tracking
        try:
            vcf_count = VCF(vcf_path)
            total_records = sum(1 for _ in vcf_count)
            log.info(f"📊 Found {total_records:,} total records to process")
        except Exception as e:
            log.error(f"❌ Error counting records: {e}")
            return None
        
        # Second pass for actual processing
        vcf = VCF(vcf_path)
        match_counts = defaultdict(int)
        records_batch = []
        processed_records = 0
        start_time = time.time()
        
        # Initial memory check
        initial_memory = get_memory_usage()
        log.info(f"💾 Initial memory usage: {initial_memory:.2f} GB")
        
        # Use rich progress bar
        with Progress(
            TextColumn("[bold blue]{task.description}"),
            BarColumn(),
            TextColumn("[bold]{task.percentage:>3.0f}%"),
            TextColumn("Records: {task.completed}/{task.total}"),
            TextColumn("Speed: {task.fields[speed]:.2f} rec/s"),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            console=console
        ) as progress:
            task = progress.add_task("[green]Processing", total=total_records, speed=0)
            
            try:
                for record in vcf:
                    records_batch.append(record)
                    processed_records += 1
                    
                    # Check if batch is ready for processing
                    if len(records_batch) >= batch_size:
                        # Check memory before processing
                        current_memory = get_memory_usage()
                        if check_memory_limits(current_memory, max_ram_usage, ram_buffer):
                            # Reduce batch size if approaching memory limits
                            new_batch_size = int(batch_size * 0.75)
                            log.warning(f"⚠️ Reducing batch size from {batch_size} to {new_batch_size} due to memory pressure")
                            batch_size = new_batch_size
                        
                        process_batch(records_batch, match_counts)
                        records_batch.clear()
                        
                        # Update progress
                        speed = processed_records / (time.time() - start_time) if time.time() > start_time else 0
                        progress.update(task, completed=processed_records, speed=speed)
                
                # Process remaining batch
                if records_batch:
                    process_batch(records_batch, match_counts)
                
                # Complete the progress bar
                progress.update(task, completed=total_records)
            
            except KeyboardInterrupt:
                log.warning("⚠️ Processing interrupted by user")
                return match_counts
            except Exception as e:
                log.error(f"❌ Error during processing: {e}")
                log.error(traceback.format_exc())
                return match_counts
        
        # Final memory check
        final_memory = get_memory_usage()
        log.info(f"💾 Final memory usage: {final_memory:.2f} GB (Delta: {final_memory - initial_memory:.2f} GB)")
        
        log.info("\n✅ Analysis complete!")
        
        log.info("\n🧠 GASLIT-AF Gene Variant Summary:")
        for gene, count in sorted(match_counts.items(), key=lambda x: -x[1])[:20]:  # Show top 20
            log.info(f"  {gene}: {count} variant(s)")
        
        elapsed = time.time() - start_time
        log.info(f"\n⏱️ Total processing time: {time.strftime('%H:%M:%S', time.gmtime(elapsed))}")
        log.info(f"⚡ Performance: {total_records/elapsed:.2f} records/second")
        
        return match_counts
    
    except Exception as e:
        log.error(f"❌ Unhandled exception: {e}")
        log.error(traceback.format_exc())
        return None

def process_batch(records, match_counts):
    try:
        # Pre-allocate a larger batch size for better GPU utilization
        genes_found = []
        
        # Use a more efficient approach to extract genes
        for record in records:
            try:
                ann = record.INFO.get('ANN')
                if not ann:
                    continue
                    
                # Process all annotations in one go using list comprehension
                # This is more efficient than nested loops
                genes_found.extend([parts[3] for entry in ann.split(',') 
                                  for parts in [entry.split('|')] 
                                  if len(parts) > 3 and parts[3] in GASLIT_AF_GENES])
            except Exception as e:
                # Log error but continue processing other records
                log.warning(f"⚠️ Error processing record: {e}")
                continue
        
        # Skip processing if no genes found
        if not genes_found:
            return
        
        # Convert to numpy array for GPU processing
        genes_array = np.array(genes_found)
        
        try:
            # SYCL-based parallel unique count with optimized GPU usage
            with dpctl.device_context(queue):
                # Use SYCL USM memory for better GPU performance
                usm_array = dpt.asarray(genes_array)
                
                # Perform unique count on GPU
                # This is more efficient than transferring back to CPU
                unique_genes, counts = np.unique(usm_array.to_numpy(), return_counts=True)
                
                # Use vectorized operations for better performance
                for gene, count in zip(unique_genes, counts):
                    match_counts[gene] += int(count)
        except Exception as e:
            # Fallback to CPU if SYCL fails
            log.warning(f"⚠️ SYCL processing failed, falling back to CPU: {e}")
            unique_genes, counts = np.unique(genes_array, return_counts=True)
            for gene, count in zip(unique_genes, counts):
                match_counts[gene] += int(count)
    except Exception as e:
        log.error(f"❌ Error in batch processing: {e}")
        log.error(traceback.format_exc())

# Progress bar function is now replaced by rich.progress

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="GASLIT-AF Variant Analysis")
    parser.add_argument("vcf_path", help="Path to VCF file for analysis")
    parser.add_argument("--batch-size", type=int, default=2000000, help="Batch size for processing")
    parser.add_argument("--max-ram", type=int, default=64, help="Maximum RAM usage in GB")
    parser.add_argument("--ram-buffer", type=int, default=16, help="RAM buffer in GB")
    parser.add_argument("--output-dir", default="./output", help="Output directory for results")
    parser.add_argument("--sample-limit", type=int, help="Limit number of variants to process (for testing)")
    parser.add_argument("--no-visualization", action="store_true", help="Skip visualization generation")
    parser.add_argument("--no-report", action="store_true", help="Skip HTML report generation")
    parser.add_argument("--enhanced-report", action="store_true", help="Generate enhanced report with interactive visualizations and symptom correlations")
    parser.add_argument("--open-browser", action="store_true", help="Automatically open report in browser")
    parser.add_argument("--cache-dir", default="./cache", help="Directory for caching intermediate results")
    parser.add_argument("--cache-max-age", type=int, default=24, help="Maximum age of cache in hours")
    parser.add_argument("--no-cache", action="store_true", help="Disable caching")
    parser.add_argument("--clear-cache", action="store_true", help="Clear cache before running")
    parser.add_argument("--system-analysis", action="store_true", default=True, help="Perform biological system-level analysis")
    parser.add_argument("--threads", type=int, default=16, help="Number of worker threads for processing")
    parser.add_argument("--use-streaming", action="store_true", help="Use streaming processor for very large VCF files (optimized memory usage)")
    
    args = parser.parse_args()
    
    try:
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize cache
        cache = AnalysisCache(
            cache_dir=args.cache_dir,
            max_age_hours=args.cache_max_age,
            enabled=not args.no_cache
        )
        
        # Clear cache if requested
        if args.clear_cache:
            cache.invalidate()
        else:
            # Clean expired cache entries
            cache.clean_expired()
        
        # Display configuration
        console.print("[bold green]GASLIT-AF Variant Analysis[/]")
        console.print(f"[bold]Configuration:[/]")
        console.print(f"  VCF File: {args.vcf_path}")
        console.print(f"  Batch Size: {args.batch_size:,} variants")
        console.print(f"  Max RAM: {args.max_ram} GB")
        console.print(f"  RAM Buffer: {args.ram_buffer} GB")
        console.print(f"  Worker Threads: {args.threads}")
        console.print(f"  Output Directory: {output_dir}")
        console.print(f"  Device: {queue.sycl_device}")
        console.print(f"  GASLIT-AF Genes: {len(GASLIT_AF_GENES)}")
        console.print(f"  Biological Systems: {len(BIOLOGICAL_SYSTEMS)}")
        console.print(f"  Caching: {'Enabled' if not args.no_cache else 'Disabled'}")
        if args.sample_limit:
            console.print(f"  [yellow]Sample Limit: {args.sample_limit} variants[/]")
        console.print()
        
        # Display cache stats if enabled
        if not args.no_cache:
            cache_stats = cache.get_stats()
            console.print(f"[bold]Cache Statistics:[/]")
            console.print(f"  Location: {cache_stats['cache_dir']}")
            console.print(f"  Total Entries: {cache_stats['total_entries']}")
            console.print(f"  Total Size: {cache_stats['total_size_mb']:.2f} MB")
            console.print(f"  Expired Entries: {cache_stats['expired_entries']}")
            console.print()
        
        # Step 1: Run initial analysis to get gene counts (with caching)
        log.info("[bold]Step 1:[/] Running initial variant analysis")
        cache_params = {
            'batch_size': args.batch_size,
            'max_ram': args.max_ram,
            'ram_buffer': args.ram_buffer,
            'use_streaming': args.use_streaming
        }
        
        # Try to get gene counts from cache
        gene_counts = None
        variant_df = None
        if not args.no_cache:
            gene_counts = cache.get(args.vcf_path, 'gene_counts', cache_params)
            if gene_counts:
                log.info("Using cached gene counts from previous analysis")
                
                # Also try to get variant data from cache
                variant_df = cache.get(args.vcf_path, 'variant_df', cache_params)
                if variant_df is not None:
                    log.info("Using cached variant data from previous analysis")
        
        # If not in cache, run analysis
        if gene_counts is None:
            if args.use_streaming:
                log.info(f"Using streaming processor (optimized for large files): {args.vcf_path}")
                
                # Define progress callback for streaming processor
                def progress_callback(processed, total):
                    if total > 0:
                        percent = (processed / total) * 100
                        log.info(f"Progress: {processed:,}/{total:,} records ({percent:.1f}%)")
                
                # Process VCF file with streaming processor
                gene_counts, variant_df = stream_process_vcf(
                    vcf_path=args.vcf_path,
                    target_genes=GASLIT_AF_GENES,
                    max_ram_usage=args.max_ram,
                    ram_buffer=args.ram_buffer,
                    max_workers=args.threads,
                    initial_chunk_size=args.batch_size,
                    progress_callback=progress_callback
                )
            else:
                gene_counts = analyze_vcf_oneapi(args.vcf_path, args.batch_size, args.max_ram, args.ram_buffer, args.threads)
            
            # Cache the results
            if not args.no_cache and gene_counts:
                cache.set(gene_counts, args.vcf_path, 'gene_counts', cache_params)
                
                # Also cache variant data if available
                if variant_df is not None and not variant_df.empty:
                    cache.set(variant_df, args.vcf_path, 'variant_df', cache_params)
        
        if gene_counts is None:
            console.print("[bold red]Analysis failed![/]")
            sys.exit(1)
        
        # Step 2: Perform biological system-level analysis
        if args.system_analysis:
            log.info("[bold]Step 2:[/] Performing biological system-level analysis")
            system_analysis = analyze_systems(gene_counts)
            
            # Save system analysis results
            system_summary_path = output_dir / "system_analysis.md"
            with open(system_summary_path, 'w') as f:
                f.write(generate_system_summary(system_analysis))
            
            # Save as JSON for programmatic access
            system_json_path = output_dir / "system_analysis.json"
            with open(system_json_path, 'w') as f:
                # Convert defaultdict to dict for JSON serialization
                json_data = {
                    "system_counts": dict(system_analysis["system_counts"]),
                    "system_percentages": dict(system_analysis["system_percentages"]),
                    "total_variants": system_analysis["total_variants"],
                    # Convert tuples to lists for JSON serialization
                    "system_genes": {k: [(g, c) for g, c in v] for k, v in system_analysis["system_genes"].items()}
                }
                json.dump(json_data, f, indent=2)
            
            log.info(f"Saved system analysis to {system_summary_path} and {system_json_path}")
        
        # Step 3: Extract variant data for visualization
        log.info("[bold]Step 3:[/] Extracting variant data for visualization")
        
        # Initialize variant_df to establish quantum coherence
        variant_df = None
        
        # If we already have variant data from pysam processing, use it
        if variant_df is None:
            # Try to get variant data from cache
            if not args.no_cache:
                cache_params = {'limit': args.sample_limit if args.sample_limit else 100000}
                variant_df = cache.get(args.vcf_path, 'variant_df', cache_params)
                if variant_df is not None:
                    log.info("Using cached variant data from previous analysis")
            
            # If not in cache, extract data
            if variant_df is None:
                try:
                    # For visualization, we may use a smaller sample to avoid memory issues
                    sample_limit = args.sample_limit if args.sample_limit else 100000
                    variant_df = vcf_to_dataframe(args.vcf_path, limit=sample_limit)
                    
                    # Cache the results
                    if not args.no_cache and not variant_df.empty:
                        cache.set(variant_df, args.vcf_path, 'variant_df', {'limit': sample_limit})
                except Exception as e:
                    log.error(f"Error converting VCF to DataFrame: {e}")
                    variant_df = pd.DataFrame()  # Empty DataFrame as fallback
        
        # Step 4: Save results to files
        log.info("[bold]Step 4:[/] Saving analysis results")
        result_paths = save_results(gene_counts, variant_df, output_dir)
        
        # Step 5: Generate visualizations
        figures = {}
        if not args.no_visualization and not variant_df.empty:
            log.info("[bold]Step 5:[/] Generating visualizations")
            viz_dir = output_dir / "visualizations"
            
            # Generate standard visualizations
            figures = generate_all_visualizations(variant_df, gene_counts, viz_dir)
            
            # Generate biological system visualizations if requested
            if args.system_analysis:
                log.info("Generating biological system visualizations")
                system_viz_dir = viz_dir / "systems"
                system_figures = plot_system_distribution(system_analysis, system_viz_dir)
                figures.update(system_figures)
            
            # Step 6: Generate reports
            if not args.no_report:
                log.info("[bold]Step 6:[/] Generating reports")
                
                # Generate standard HTML report
                report_path = generate_html_report(
                    variant_df=variant_df,
                    output_dir=output_dir,
                    args=args,
                    system_results=system_analysis if args.system_analysis else None,
                    figures=figures,
                    gene_counts=gene_counts
                )
                console.print(f"[bold green]Standard HTML Report:[/] {report_path}")
                
                # Generate enhanced report if requested
                if args.enhanced_report:
                    log.info("Generating enhanced report with interactive visualizations and symptom correlations")
                    enhanced_report_path = generate_enhanced_report(
                        gene_counts, 
                        variant_df, 
                        figures, 
                        output_dir, 
                        system_analysis=system_analysis if args.system_analysis else None,
                        include_symptoms=True
                    )
                    console.print(f"[bold green]Enhanced Report:[/] {enhanced_report_path}")
                    
                    # Open the report in browser if requested
                    if args.open_browser:
                        log.info(f"Opening enhanced report in browser")
                        try:
                            webbrowser.open(f"file://{os.path.abspath(enhanced_report_path)}")
                        except Exception as e:
                            log.warning(f"Could not open browser: {e}")
        else:
            if args.no_visualization:
                log.info("Visualization generation skipped as requested")
            else:
                log.warning("Visualization skipped due to empty variant data")
        
        console.print("\n[bold green]Analysis completed successfully![/]")
        console.print(f"Results saved to: {output_dir}")
        
    except Exception as e:
        console.print(f"[bold red]Unhandled exception:[/] {e}")
        console.print(traceback.format_exc())
        sys.exit(1)

if __name__ == "__main__":
    main()
