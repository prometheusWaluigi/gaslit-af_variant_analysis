"""
Workflow module for GASLIT-AF Variant Analysis.
Handles the main analysis workflow and orchestration.

This module has been refactored for improved maintainability with smaller,
more focused functions and better error handling.
"""

import os
import sys
import time
import logging
import traceback
import webbrowser
import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any, Union
from rich.console import Console
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn, TimeRemainingColumn

# Import internal modules
from src.gaslit_af.device import (
    initialize_device, 
    get_memory_usage, 
    get_memory_pressure,
    get_optimal_batch_size,
    throttle_processing,
    manage_memory_proactively
)
from src.gaslit_af.gene_lists import GASLIT_AF_GENES, KNOWN_SNPS
from src.gaslit_af.data_processing import vcf_to_dataframe, extract_gaslit_af_variants, save_results, process_batch_with_data, update_gene_counts
from src.gaslit_af.streaming import stream_process_vcf
from src.gaslit_af.visualization import generate_all_visualizations
from src.gaslit_af.reporting import generate_html_report
from src.gaslit_af.enhanced_reporting import generate_enhanced_report
from src.gaslit_af.caching import AnalysisCache
from src.gaslit_af.biological_systems import analyze_systems, plot_system_distribution, generate_system_summary
from src.gaslit_af.advanced_variant_processing import process_vcf_with_pysam, VariantProcessor
from src.gaslit_af.clinical_integration import ClinicalIntegration
from src.gaslit_af.api_integration import VariantAPIIntegration
from src.gaslit_af.variant_enrichment import VariantEnricher, enrich_variants_with_af_data
from src.gaslit_af.annovar_integration import AnnovarIntegration, enrich_variants_with_annovar, is_annovar_available
from src.gaslit_af.rccx_analysis import analyze_rccx_cnv_from_vcf, RCCX_REGION_GRCH38
from src.gaslit_af.exceptions import (
    GaslitAFError, 
    DataProcessingError, 
    AnnotationError,
    DeviceError,
    MemoryError,
    CacheError,
    APIError,
    FileError,
    ConfigurationError,
    VisualizationError,
    ReportingError,
    retry_operation
)

# Configure logging
log = logging.getLogger("gaslit-af")
console = Console()

def setup_environment(args) -> Tuple[AnalysisCache, Optional[ClinicalIntegration], Optional[VariantAPIIntegration], Optional[VariantEnricher], Any]:
    """
    Set up the analysis environment.
    
    Args:
        args: Command-line arguments namespace
        
    Returns:
        Tuple containing cache, clinical_integration, api_integration, variant_enricher, queue
    """
    try:
        # Create output directory
        args.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize cache
        cache = AnalysisCache(
            cache_dir=args.cache_dir,
            max_age_hours=args.cache_max_age,
            enabled=not args.no_cache
        )
        
        # Initialize clinical integration if clinical data is provided
        clinical_integration = None
        if hasattr(args, 'clinical_data') and args.clinical_data:
            log.info(f"Initializing clinical integration with data from {args.clinical_data}")
            clinical_integration = ClinicalIntegration(args.clinical_data)
        
        # Initialize API integration if enabled
        api_integration = None
        if hasattr(args, 'api_annotation') and args.api_annotation:
            log.info(f"Initializing API integration with sources: {', '.join(args.api_sources)}")
            api_integration = VariantAPIIntegration(
                cache_dir=Path(args.api_cache_dir),
                cache_ttl=args.api_cache_ttl
            )
            
        # Initialize variant enricher if enabled
        variant_enricher = None
        if hasattr(args, 'variant_enrichment') and args.variant_enrichment:
            log.info("Initializing variant enrichment for AF-specific analysis")
            variant_enricher = VariantEnricher(
                cache_dir=Path(args.enrichment_cache_dir) if hasattr(args, 'enrichment_cache_dir') else None,
                cache_ttl=args.enrichment_cache_ttl if hasattr(args, 'enrichment_cache_ttl') else 24
            )
        
        # Initialize device queue
        queue = initialize_device()
        
        return cache, clinical_integration, api_integration, variant_enricher, queue
    
    except Exception as e:
        error_msg = f"Error setting up environment: {e}"
        log.error(error_msg)
        raise ConfigurationError(error_msg, details=str(e))

def print_configuration(args, queue, cache):
    """
    Print configuration information.
    
    Args:
        args: Command-line arguments namespace
        queue: SYCL queue
        cache: AnalysisCache instance
    """
    console.print("GASLIT-AF Variant Analysis")
    console.print("Configuration:")
    console.print(f"  VCF File: {args.vcf_path}")
    console.print(f"  Batch Size: {args.batch_size:,} variants")
    console.print(f"  Max RAM: {args.max_ram} GB")
    console.print(f"  RAM Buffer: {args.ram_buffer} GB")
    console.print(f"  Worker Threads: {args.threads}")
    console.print(f"  Output Directory: {args.output_dir}")
    console.print(f"  Device: {queue.sycl_device}")
    console.print(f"  GASLIT-AF Genes: {len(GASLIT_AF_GENES)}")
    console.print(f"  Biological Systems: 8")
    console.print(f"  Processing Mode: {'Streaming' if hasattr(args, 'use_streaming') and args.use_streaming else 'Standard'}")
    console.print(f"  Caching: {'Enabled' if not args.no_cache else 'Disabled'}")
    console.print("")
    
    # Print cache statistics
    if not args.no_cache:
        console.print("Cache Statistics:")
        console.print(f"  Location: {args.cache_dir}")
        console.print(f"  Total Entries: {cache.count_entries()}")
        console.print(f"  Total Size: {cache.get_total_size():.2f} MB")
        console.print(f"  Expired Entries: {cache.count_expired_entries()}")
        console.print("")

def perform_variant_analysis(args, cache, queue) -> Tuple[Dict, pd.DataFrame]:
    """
    Perform initial variant analysis.
    
    Args:
        args: Command-line arguments namespace
        cache: AnalysisCache instance
        queue: SYCL queue
        
    Returns:
        Tuple containing gene_counts and variant_df
    """
    log.info("[bold]Step 1:[/] Running initial variant analysis")
    
    # Check if results are in cache
    gene_counts = None
    variant_df = None
    if not args.no_cache:
        cache_params = {'use_pysam': args.use_pysam, 'use_streaming': getattr(args, 'use_streaming', False)}
        gene_counts = cache.get(args.vcf_path, 'gene_counts', cache_params)
        if gene_counts is not None:
            log.info("Using cached gene counts from previous analysis")
            
            # Also try to get variant data from cache
            variant_df = cache.get(args.vcf_path, 'variant_df', cache_params)
            if variant_df is not None:
                log.info("Using cached variant data from previous analysis")
    
    # If not in cache, analyze VCF file
    if gene_counts is None:
        try:
            # Get the set of target genes - either all GASLIT-AF genes or just the ones in the known variants map
            target_genes = set(KNOWN_SNPS.keys()) if args.known_variants_only else GASLIT_AF_GENES
            
            # Choose processing method based on args
            if hasattr(args, 'use_streaming') and args.use_streaming:
                log.info(f" Analyzing with streaming processor (optimized for large files): {args.vcf_path}")
                
                # Define progress callback for streaming processor
                def progress_callback(processed, total):
                    if total > 0:
                        percent = (processed / total) * 100
                        log.info(f" Progress: {processed:,}/{total:,} records ({percent:.1f}%)")
                
                # Process VCF file with streaming processor
                gene_counts, variant_df = stream_process_vcf(
                    vcf_path=args.vcf_path,
                    target_genes=target_genes,
                    max_ram_usage=args.max_ram,
                    ram_buffer=args.ram_buffer,
                    max_workers=args.threads,
                    initial_chunk_size=args.batch_size,
                    progress_callback=progress_callback
                )
            elif args.use_pysam:
                log.info(f" Analyzing with pysam (oneAPI accelerated): {args.vcf_path}")
                
                # Process VCF file with pysam
                gene_counts, variant_df = process_vcf_with_pysam(
                    vcf_path=args.vcf_path,
                    target_genes=target_genes,
                    dbsnp_path=args.dbsnp_path,
                    queue=queue,
                    threads=args.threads,
                    batch_size=args.batch_size,
                    max_ram_usage=args.max_ram,
                    ram_buffer=args.ram_buffer
                )
                
                # Generate variant report
                if not variant_df.empty:
                    processor = VariantProcessor(queue=queue, threads=args.threads)
                    variant_report_path = processor.generate_variant_report(variant_df, args.output_dir)
                    log.info(f"Generated variant report: {variant_report_path}")
            else:
                log.info(f" Analyzing (oneAPI accelerated): {args.vcf_path}")
                gene_counts = analyze_vcf_oneapi(
                    args.vcf_path,
                    batch_size=args.batch_size,
                    max_ram_usage=args.max_ram,
                    ram_buffer=args.ram_buffer,
                    threads=args.threads
                )
        except Exception as e:
            error_msg = f"Error analyzing VCF file: {e}"
            log.error(error_msg)
            raise DataProcessingError(error_msg, details=str(e))
        
        # Cache the results
        if not args.no_cache and gene_counts:
            cache_params = {'use_pysam': args.use_pysam, 'use_streaming': getattr(args, 'use_streaming', False)}
            cache.set(gene_counts, args.vcf_path, 'gene_counts', cache_params)
            
            # Also cache variant data if available
            if variant_df is not None and not variant_df.empty:
                cache.set(variant_df, args.vcf_path, 'variant_df', cache_params)
    
    if gene_counts is None:
        error_msg = "Analysis failed to produce gene counts"
        log.error(error_msg)
        raise DataProcessingError(error_msg)
    
    return gene_counts, variant_df

def perform_system_analysis(args, variant_df) -> Tuple[Optional[Dict], List]:
    """
    Perform biological system-level analysis.
    
    Args:
        args: Command-line arguments namespace
        variant_df: DataFrame with variant information
        
    Returns:
        Tuple containing system_results and rccx_results
    """
    system_results = None
    rccx_results = []  # Initialize RCCX results
    
    if variant_df is None or variant_df.empty:
        log.warning("Biological system analysis skipped due to empty variant data")
        return None, []
    
    try:
        log.info("Performing biological system-level analysis")
        system_results = analyze_systems(variant_df, GASLIT_AF_GENES)
        
        # Save system analysis results
        system_summary = generate_system_summary(system_results)
        system_output_path = args.output_dir / "system_analysis.md"
        with open(system_output_path, 'w') as f:
            f.write(system_summary)
        
        # Save as JSON for programmatic access
        system_json_path = args.output_dir / "system_analysis.json"
        with open(system_json_path, 'w') as f:
            # Convert defaultdict to dict for JSON serialization
            json_data = {
                "system_counts": dict(system_results["system_counts"]),
                "system_percentages": dict(system_results["system_percentages"]),
                "total_variants": system_results["total_variants"],
                # Convert tuples to lists for JSON serialization
                "system_genes": {k: [(g, c) for g, c in v] for k, v in system_results["system_genes"].items()}
            }
            json.dump(json_data, f, indent=2)
        console.print(f"[bold green]System Analysis Results:[/] {system_output_path}")
        
        # --- RCCX Structural Variant Analysis (from VCF) ---
        # Ensure VCF path is valid before proceeding
        if Path(args.vcf_path).is_file():
            log.info(f"Scanning VCF for potential RCCX structural variants ({RCCX_REGION_GRCH38['name']})...")
            # Assuming GRCh38 for now, could pass reference from args if needed
            rccx_results = analyze_rccx_cnv_from_vcf(args.vcf_path, reference="GRCh38")
            if rccx_results:
                console.print(f"[bold cyan]RCCX Scan:[/] Found {len(rccx_results)} potential SV/CNV records in VCF.")
            else:
                console.print("[bold cyan]RCCX Scan:[/] No potential SV/CNV records found in VCF.")
        else:
            log.warning(f"RCCX Scan skipped: VCF file not found at {args.vcf_path}")
    
    except Exception as e:
        error_msg = f"Error performing system analysis: {e}"
        log.error(error_msg)
        # Don't raise an exception here, just return None for system_results
        system_results = None
    
    return system_results, rccx_results

def extract_variant_data_for_visualization(args, cache, variant_df):
    """
    Extract variant data for visualization if not already available.
    
    Args:
        args: Command-line arguments namespace
        cache: AnalysisCache instance
        variant_df: Existing variant DataFrame or None
        
    Returns:
        DataFrame with variant information
    """
    log.info("Extracting variant data for visualization")
    
    # If we already have variant data from pysam processing, use it
    if variant_df is not None and not variant_df.empty:
        return variant_df
    
    # Try to get variant data from cache
    if not args.no_cache:
        cache_params = {'limit': args.sample_limit if args.sample_limit else 100000, 'use_pysam': args.use_pysam}
        variant_df = cache.get(args.vcf_path, 'variant_df', cache_params)
        if variant_df is not None:
            log.info("Using cached variant data from previous analysis")
            return variant_df
    
    # If not in cache, extract data
    try:
        # For visualization, we may use a smaller sample to avoid memory issues
        sample_limit = args.sample_limit if args.sample_limit else 100000
        variant_df = vcf_to_dataframe(args.vcf_path, limit=sample_limit)
        
        # Cache the results
        if not args.no_cache and not variant_df.empty:
            cache.set(variant_df, args.vcf_path, 'variant_df', {'limit': sample_limit, 'use_pysam': False})
        
        return variant_df
    except Exception as e:
        error_msg = f"Error converting VCF to DataFrame: {e}"
        log.error(error_msg)
        # Return empty DataFrame as fallback
        return pd.DataFrame()

def annotate_variants_with_api(args, api_integration, variant_df):
    """
    Annotate variants with external APIs.
    
    Args:
        args: Command-line arguments namespace
        api_integration: VariantAPIIntegration instance
        variant_df: DataFrame with variant information
        
    Returns:
        Annotated DataFrame
    """
    if api_integration is None or variant_df.empty or not args.api_annotation:
        log.info("Skipping API annotation (not enabled or no variants to annotate)")
        return variant_df
    
    try:
        log.info("Annotating variants with external APIs")
        console.print(f"Annotating variants with sources: {', '.join(args.api_sources)}")
        
        # Annotate variants
        annotated_df = api_integration.annotate_variants(variant_df, sources=args.api_sources)
        
        # Log annotation summary
        api_cols = [col for col in annotated_df.columns if col.startswith('ensembl_') 
                  or col in ['clinvar_significance', 'cadd_phred', 'gnomad_af', 'sift_pred', 'polyphen_pred']]
        
        for col in api_cols:
            non_null = annotated_df[col].notna().sum()
            if non_null > 0:
                log.info(f"  - {col}: {non_null} variants annotated")
        
        return annotated_df
    except Exception as e:
        error_msg = f"Error annotating variants with API: {e}"
        log.error(error_msg)
        # Return original DataFrame if annotation fails
        return variant_df

def annotate_variants_with_annovar(args, variant_df):
    """
    Annotate variants with ANNOVAR.
    
    Args:
        args: Command-line arguments namespace
        variant_df: DataFrame with variant information
        
    Returns:
        Annotated DataFrame
    """
    if variant_df.empty or not hasattr(args, 'use_annovar') or not args.use_annovar:
        if not hasattr(args, 'use_annovar') or not args.use_annovar:
            log.info("Skipping ANNOVAR annotation (not enabled)")
        elif variant_df.empty:
            log.info("Skipping ANNOVAR annotation (no variants to annotate)")
        return variant_df
    
    try:
        log.info("Annotating variants with ANNOVAR")
        
        # Check if ANNOVAR is available
        if not is_annovar_available():
            log.warning("ANNOVAR not available, skipping annotation")
            console.print("[bold yellow]ANNOVAR not found.[/] To enable ANNOVAR annotation, install ANNOVAR and provide paths using --annovar-path and --humandb-path")
            return variant_df
        
        # Create ANNOVAR directory
        annovar_dir = args.output_dir / "annovar"
        annovar_dir.mkdir(exist_ok=True)
        
        # Get ANNOVAR paths from args if provided
        annovar_path = getattr(args, 'annovar_path', None)
        humandb_path = getattr(args, 'humandb_path', None)
        
        # Perform ANNOVAR annotation
        log.info(f"Enriching variants with ANNOVAR annotations")
        enriched_df = enrich_variants_with_annovar(
            variant_df,
            output_dir=annovar_dir,
            annovar_path=annovar_path,
            humandb_path=humandb_path
        )
        
        if not enriched_df.equals(variant_df):
            # Save ANNOVAR-annotated variants
            annovar_output = args.output_dir / "annovar_variants.csv"
            enriched_df.to_csv(annovar_output, index=False)
            console.print(f"[bold green]ANNOVAR-annotated variants saved to:[/] {annovar_output}")
            
            # Log annotation summary
            annovar_cols = [col for col in enriched_df.columns if col.startswith('Func.') 
                          or col.startswith('ExonicFunc.') or col.startswith('Gene.') 
                          or col in ['CLNSIG', 'CLNDN', 'ExAC_ALL', 'gnomAD_exome_ALL']]
            
            for col in annovar_cols:
                non_null = enriched_df[col].notna().sum()
                if non_null > 0:
                    log.info(f"  - {col}: {non_null} variants annotated")
            
            return enriched_df
        else:
            log.warning("ANNOVAR annotation did not add any new information")
            return variant_df
    except Exception as e:
        error_msg = f"Error annotating variants with ANNOVAR: {e}"
        log.error(error_msg)
        # Return original DataFrame if annotation fails
        return variant_df

def enrich_variants(args, variant_df):
    """
    Enrich variants with AF-specific annotations.
    
    Args:
        args: Command-line arguments namespace
        variant_df: DataFrame with variant information
        
    Returns:
        Enriched DataFrame
    """
    if not hasattr(args, 'variant_enrichment') or not args.variant_enrichment or variant_df.empty:
        return variant_df
    
    try:
        log.info("Enriching variants with AF-specific annotations")
        
        # Create enrichment directory
        enrichment_dir = args.output_dir / "enrichment"
        enrichment_dir.mkdir(exist_ok=True)
        
        # Perform enrichment
        enriched_df, enrichment_report = enrich_variants_with_af_data(
            variant_df,
            cache_dir=Path(args.enrichment_cache_dir) if hasattr(args, 'enrichment_cache_dir') else None,
            output_dir=enrichment_dir
        )
        
        if not enriched_df.empty:
            # Save enriched variants
            enriched_output = args.output_dir / "enriched_variants.csv"
            enriched_df.to_csv(enriched_output, index=False)
            console.print(f"[bold green]Enriched variants saved to:[/] {enriched_output}")
            
            if enrichment_report:
                console.print(f"[bold green]AF Enrichment Report:[/] {enrichment_report}")
            
            return enriched_df
        else:
            return variant_df
    except Exception as e:
        error_msg = f"Error enriching variants: {e}"
        log.error(error_msg)
        # Return original DataFrame if enrichment fails
        return variant_df

def generate_visualizations_and_reports(args, variant_df, gene_counts, system_results, rccx_results):
    """
    Generate visualizations and reports.
    
    Args:
        args: Command-line arguments namespace
        variant_df: DataFrame with variant information
        gene_counts: Dictionary of gene:count pairs
        system_results: System analysis results
        rccx_results: RCCX analysis results
        
    Returns:
        Dictionary of generated figures
    """
    figures = {}
    
    if args.no_visualization or variant_df.empty:
        if args.no_visualization:
            log.info("Visualization generation skipped as requested")
        else:
            log.warning("Visualization skipped due to empty variant data")
        return figures
    
    try:
        log.info("Generating visualizations")
        viz_dir = args.output_dir / "visualizations"
        
        # Generate standard visualizations
        figures = generate_all_visualizations(variant_df, gene_counts, viz_dir)
        
        # Generate biological system visualizations if requested
        if args.system_analysis and system_results:
            log.info("Generating biological system visualizations")
            system_viz_dir = viz_dir / "systems"
            system_figures = plot_system_distribution(system_results, system_viz_dir)
            figures.update(system_figures)
        
        # Generate reports
        if not args.no_report:
            log.info("Generating standard HTML report")
            report_path = generate_html_report(
                variant_df=variant_df,
                output_dir=args.output_dir,
                args=args,
                system_results=system_results,
                rccx_results=rccx_results
            )
            console.print(f"[bold green]Standard Report:[/] {report_path}")
            
            # Generate enhanced report if requested
            if args.enhanced_report:
                log.info(f"Generating enhanced report in {args.output_dir} (including RCCX results if available)")
                enhanced_report_path = generate_enhanced_report(
                    variant_df=variant_df,
                    output_dir=args.output_dir,
                    system_results=system_results,
                    rccx_results=rccx_results,
                    figures=figures,
                    include_symptoms=args.include_symptoms
                )
                log.info(f"Enhanced HTML report saved to {enhanced_report_path}")
        
        return figures
    except Exception as e:
        error_msg = f"Error generating visualizations and reports: {e}"
        log.error(error_msg)
        # Return empty dictionary if visualization fails
        return {}

def analyze_vcf_oneapi(vcf_path, batch_size=2000000, max_ram_usage=64, ram_buffer=16, threads=16):
    """
    Analyze VCF file for GASLIT-AF gene variants using oneAPI acceleration.
    
    Args:
        vcf_path: Path to VCF file
        batch_size: Batch size for processing
        max_ram_usage: Maximum RAM usage in GB
        ram_buffer: RAM buffer in GB
        threads: Number of worker threads
        
    Returns:
        Dictionary of gene:count pairs
    """
    from src.gaslit_af.data_processing import process_batch_with_data
    
    # Initialize device queue
    queue = initialize_device()
    
    try:
        # Open VCF file
        from cyvcf2 import VCF
        vcf = VCF(vcf_path)
        
        # Get total number of records for progress tracking
        log.info("Counting total records in VCF file...")
        total_records = 0
        for _ in vcf:
            total_records += 1
            # Break if we've counted enough records to estimate
            if total_records >= 100000:
                # Estimate total based on file size
                import os
                file_size = os.path.getsize(vcf_path)
                records_per_byte = total_records / file_size
                total_records = int(file_size * records_per_byte)
                log.info(f"Estimated total records: {total_records:,} (based on first 100,000 records)")
                break
        
        # Reset VCF iterator
        vcf = VCF(vcf_path)
        
        log.info(f" Found {total_records:,} total records to process")
        
        # Initialize progress bar
        progress_bar = Progress(
            TextColumn("Processing"),
            BarColumn(),
            TextColumn("{task.percentage:>3.0f}%"),
            TextColumn("Records: {task.completed}/{task.total}"),
            # Handle potential None value for speed
            TextColumn("Speed: {task.speed if task.speed is not None else 0.0:.2f} rec/s"),
            TimeElapsedColumn(),
            TimeRemainingColumn()
        )
        
        # Initialize counters and data structures
        match_counts = {}
        variant_data = []
        processed = 0
        
        # Monitor memory usage
        initial_memory = get_memory_usage()
        log.info(f" Initial memory usage: {initial_memory:.2f} GB")
        
        # Process in batches
        with progress_bar as progress:
            task = progress.add_task("", total=total_records)
            
            batch = []
            start_time = time.time()
            
            for record in vcf:
                batch.append(record)
                
                if len(batch) >= batch_size:
                    # Process batch
                    genes_found, batch_variant_data = process_batch_with_data(batch, GASLIT_AF_GENES)
                    update_gene_counts(genes_found, match_counts, queue)
                    variant_data.extend(batch_variant_data)
                    
                    # Update progress
                    processed += len(batch)
                    progress.update(task, completed=processed)
                    
                    # Proactive memory management
                    current_memory = get_memory_usage()
                    continue_processing, new_batch_size = manage_memory_proactively(
                        current_memory, 
                        max_ram_usage, 
                        batch_size,
                        min_batch_size=10000
                    )
                    
                    # Update batch size if needed
                    if new_batch_size != batch_size:
                        log.info(f" Adjusting batch size: {batch_size:,} → {new_batch_size:,} variants")
                        batch_size = new_batch_size
                    
                    # Check if we should stop processing due to memory pressure
                    if not continue_processing:
                        log.error(f" Stopping processing due to critical memory pressure")
                        break
                    
                    # Clear batch
                    batch = []
            
            # Process final batch if any
            if batch:
                genes_found, batch_variant_data = process_batch_with_data(batch, GASLIT_AF_GENES)
                update_gene_counts(genes_found, match_counts, queue)
                variant_data.extend(batch_variant_data)
                
                processed += len(batch)
                progress.update(task, completed=processed)
        
        # Final memory usage
        final_memory = get_memory_usage()
        log.info(f" Final memory usage: {final_memory:.2f} GB (Delta: {final_memory - initial_memory:.2f} GB)")
        
        log.info("\n Analysis complete!")
        
        log.info("\n GASLIT-AF Gene Variant Summary:")
        
        # Calculate processing time and speed
        end_time = time.time()
        processing_time = end_time - start_time
        processing_speed = processed / processing_time if processing_time > 0 else 0
        
        # Format processing time as HH:MM:SS
        hours, remainder = divmod(processing_time, 3600)
        minutes, seconds = divmod(remainder, 60)
        formatted_time = f"{int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}"
        
        log.info(f" Total processing time: {formatted_time}")
        log.info(f" Performance: {processing_speed:.2f} records/second")
        
        return match_counts
    
    except Exception as e:
        error_msg = f"Error analyzing VCF file: {e}"
        log.error(error_msg)
        log.error(traceback.format_exc())
        raise DataProcessingError(error_msg, details=str(e))

def run_analysis_workflow(args):
    """
    Run the main analysis workflow.
    
    Args:
        args: Command-line arguments namespace
    """
    try:
        # Step 1: Set up environment
        cache, clinical_integration, api_integration, variant_enricher, queue = setup_environment(args)
        
        # Step 2: Print configuration
        print_configuration(args, queue, cache)
        
        # Step 3: Run initial variant analysis
        gene_counts, variant_df = perform_variant_analysis(args, cache, queue)
        
        # Step 4: Perform biological system-level analysis
        system_results, rccx_results = perform_system_analysis(args, variant_df)
        
        # Step 5: Extract variant data for visualization if needed
        if variant_df is None or variant_df.empty:
            variant_df = extract_variant_data_for_visualization(args, cache, variant_df)
        
        # Step 6: Annotate variants with external APIs if enabled
        variant_df = annotate_variants_with_api(args, api_integration, variant_df)
        
        # Step 7: Annotate variants with ANNOVAR if available
        variant_df = annotate_variants_with_annovar(args, variant_df)
        
        # Step 8: Save results to files
        log.info("Saving analysis results")
        result_paths = save_results(gene_counts, variant_df, args.output_dir)
        
        # Step 9: Perform variant enrichment if enabled
        variant_df = enrich_variants(args, variant_df)
        
        # Step 10: Generate visualizations and reports
        figures = generate_visualizations_and_reports(args, variant_df, gene_counts, system_results, rccx_results)
        
        console.print("\n[bold green]Analysis completed successfully![/]")
        console.print(f"Results saved to: {args.output_dir}")
        
    except GaslitAFError as e:
        console.print(f"[bold red]GASLIT-AF Error:[/] {e}")
        if hasattr(e, 'details') and e.details:
            console.print(f"[bold red]Details:[/] {e.details}")
        sys.exit(1)
    except Exception as e:
        console.print(f"[bold red]Unhandled exception:[/] {e}")
        console.print(traceback.format_exc())
        sys.exit(1)
