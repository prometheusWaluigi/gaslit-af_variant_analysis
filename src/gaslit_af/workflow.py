"""
Workflow module for GASLIT-AF Variant Analysis.
Handles the main analysis workflow and orchestration.
"""

import os
import sys
import time
import logging
import traceback
import webbrowser
import pandas as pd
from pathlib import Path
from rich.console import Console
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn, TimeRemainingColumn

# Import internal modules
from src.gaslit_af.device import initialize_device, get_memory_usage
from src.gaslit_af.gene_lists import GASLIT_AF_GENES, KNOWN_SNPS
from src.gaslit_af.data_processing import vcf_to_dataframe, extract_gaslit_af_variants, save_results, process_batch_with_data
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

# Configure logging
log = logging.getLogger("gaslit-af")
console = Console()

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
        
        # Initialize counters
        match_counts = {}
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
                    
                    # Check memory usage
                    current_memory = get_memory_usage()
                    if current_memory > (max_ram_usage - ram_buffer):
                        log.warning(f" Memory usage ({current_memory:.2f} GB) is approaching limit ({max_ram_usage} GB)")
                    
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
        log.error(f"Error analyzing VCF file: {e}")
        log.error(traceback.format_exc())
        return {}

def run_analysis_workflow(args):
    """
    Run the main analysis workflow.
    
    Args:
        args: Command-line arguments namespace
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
        
        # Print configuration
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
        
        # Step 1: Run initial variant analysis
        log.info("[bold]Step 1:[/] Running initial variant analysis")
        
        # Check if results are in cache
        gene_counts = None
        variant_df = None
        if not args.no_cache:
            cache_params = {'use_pysam': args.use_pysam}
            gene_counts = cache.get(args.vcf_path, 'gene_counts', cache_params)
            if gene_counts is not None:
                log.info("Using cached gene counts from previous analysis")
                
                # Also try to get variant data from cache
                variant_df = cache.get(args.vcf_path, 'variant_df', cache_params)
                if variant_df is not None:
                    log.info("Using cached variant data from previous analysis")
        
        # If not in cache, analyze VCF file
        if gene_counts is None:
            if args.use_pysam:
                log.info(f" Analyzing with pysam (oneAPI accelerated): {args.vcf_path}")
                
                # Get the set of target genes - either all GASLIT-AF genes or just the ones in the known variants map
                target_genes = set(KNOWN_SNPS.keys()) if args.known_variants_only else GASLIT_AF_GENES
                
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
            
            # Cache the results
            if not args.no_cache and gene_counts:
                cache_params = {'use_pysam': args.use_pysam}
                cache.set(gene_counts, args.vcf_path, 'gene_counts', cache_params)
                
                # Also cache variant data if available
                if variant_df is not None and not variant_df.empty:
                    cache.set(variant_df, args.vcf_path, 'variant_df', cache_params)
        
        if gene_counts is None:
            console.print("[bold red]Analysis failed![/]")
            sys.exit(1)
        
        # Step 2: Perform biological system-level analysis
        system_analysis = None
        if args.system_analysis:
            log.info("[bold]Step 2:[/] Performing biological system-level analysis")
            system_analysis = analyze_systems(gene_counts)
            
            # Save system analysis results
            system_summary = generate_system_summary(system_analysis)
            system_summary_path = args.output_dir / "system_analysis.md"
            with open(system_summary_path, 'w') as f:
                f.write(system_summary)
            
            # Save as JSON for programmatic access
            system_json_path = args.output_dir / "system_analysis.json"
            with open(system_json_path, 'w') as f:
                # Convert defaultdict to dict for JSON serialization
                import json
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
        
        # If we already have variant data from pysam processing, use it
        if variant_df is None:
            # Try to get variant data from cache
            if not args.no_cache:
                cache_params = {'limit': args.sample_limit if args.sample_limit else 100000, 'use_pysam': args.use_pysam}
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
                        cache.set(variant_df, args.vcf_path, 'variant_df', {'limit': sample_limit, 'use_pysam': False})
                except Exception as e:
                    log.error(f"Error converting VCF to DataFrame: {e}")
                    variant_df = pd.DataFrame()  # Empty DataFrame as fallback
        
        # Step 4: Annotate variants with external APIs if enabled
        if api_integration and not variant_df.empty and args.api_annotation:
            log.info("[bold]Step 4:[/] Annotating variants with external APIs")
            console.print(f"Annotating variants with sources: {', '.join(args.api_sources)}")
            
            # Annotate variants
            variant_df = api_integration.annotate_variants(variant_df, sources=args.api_sources)
            
            # Log annotation summary
            api_cols = [col for col in variant_df.columns if col.startswith('ensembl_') 
                      or col in ['clinvar_significance', 'cadd_phred', 'gnomad_af', 'sift_pred', 'polyphen_pred']]
            
            for col in api_cols:
                non_null = variant_df[col].notna().sum()
                if non_null > 0:
                    log.info(f"  - {col}: {non_null} variants annotated")
        else:
            log.info("Skipping API annotation (not enabled or no variants to annotate)")
            
        # Step 4.5: Annotate variants with ANNOVAR if available
        if not variant_df.empty and hasattr(args, 'use_annovar') and args.use_annovar:
            log.info("[bold]Step 4.5:[/] Annotating variants with ANNOVAR")
            
            # Check if ANNOVAR is available
            if is_annovar_available():
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
                    # Update variant_df with ANNOVAR annotations
                    variant_df = enriched_df
                    
                    # Save ANNOVAR-annotated variants
                    annovar_output = args.output_dir / "annovar_variants.csv"
                    enriched_df.to_csv(annovar_output, index=False)
                    console.print(f"[bold green]ANNOVAR-annotated variants saved to:[/] {annovar_output}")
                    
                    # Log annotation summary
                    annovar_cols = [col for col in variant_df.columns if col.startswith('Func.') 
                                  or col.startswith('ExonicFunc.') or col.startswith('Gene.') 
                                  or col in ['CLNSIG', 'CLNDN', 'ExAC_ALL', 'gnomAD_exome_ALL']]
                    
                    for col in annovar_cols:
                        non_null = variant_df[col].notna().sum()
                        if non_null > 0:
                            log.info(f"  - {col}: {non_null} variants annotated")
                else:
                    log.warning("ANNOVAR annotation did not add any new information")
            else:
                log.warning("ANNOVAR not available, skipping annotation")
                console.print("[bold yellow]ANNOVAR not found.[/] To enable ANNOVAR annotation, install ANNOVAR and provide paths using --annovar-path and --humandb-path")
        else:
            if not hasattr(args, 'use_annovar') or not args.use_annovar:
                log.info("Skipping ANNOVAR annotation (not enabled)")
            elif variant_df.empty:
                log.info("Skipping ANNOVAR annotation (no variants to annotate)")
        
        # Step 5: Save results to files
        log.info("[bold]Step 5:[/] Saving analysis results")
        result_paths = save_results(gene_counts, variant_df, args.output_dir)
        
        # Step 6: Perform variant enrichment if enabled
        if hasattr(args, 'variant_enrichment') and args.variant_enrichment and not variant_df.empty:
            log.info("[bold]Step 6:[/] Enriching variants with AF-specific annotations")
            
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
                # Update variant_df with enriched data
                variant_df = enriched_df
                
                # Save enriched variants
                enriched_output = args.output_dir / "enriched_variants.csv"
                enriched_df.to_csv(enriched_output, index=False)
                console.print(f"[bold green]Enriched variants saved to:[/] {enriched_output}")
                
                if enrichment_report:
                    console.print(f"[bold green]AF Enrichment Report:[/] {enrichment_report}")
        
        # Step 7: Generate visualizations
        figures = {}
        if not args.no_visualization and not variant_df.empty:
            log.info("[bold]Step 6:[/] Generating visualizations")
            viz_dir = args.output_dir / "visualizations"
            
            # Generate standard visualizations
            figures = generate_all_visualizations(variant_df, gene_counts, viz_dir)
            
            # Generate biological system visualizations if requested
            if args.system_analysis and system_analysis:
                log.info("Generating biological system visualizations")
                system_viz_dir = viz_dir / "systems"
                system_figures = plot_system_distribution(system_analysis, system_viz_dir)
                figures.update(system_figures)
            
            # Step 8: Generate reports
            if not args.no_report:
                log.info("[bold]Step 7:[/] Generating reports")
                
                # Generate standard HTML report
                report_path = generate_html_report(gene_counts, variant_df, figures, args.output_dir, 
                                                  system_analysis if args.system_analysis else None)
                console.print(f"[bold green]Standard HTML Report:[/] {report_path}")
                
                # Generate enhanced report if requested
                if args.enhanced_report:
                    log.info("Generating enhanced report with interactive visualizations and symptom correlations")
                    enhanced_report_path = generate_enhanced_report(
                        gene_counts, 
                        variant_df, 
                        figures, 
                        args.output_dir, 
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
                
                # Generate clinical report if clinical data is available
                if clinical_integration and clinical_integration.clinical_data_loaded:
                    log.info("Generating clinical variant report")
                    clinical_report_path = clinical_integration.generate_clinical_report(variant_df, args.output_dir)
                    if clinical_report_path:
                        console.print(f"[bold green]Clinical Report:[/] {clinical_report_path}")
                        
                        # Add clinical annotations to the variant data
                        annotated_df = clinical_integration.annotate_variants(variant_df)
                        annotated_output = args.output_dir / "clinical_variants.csv"
                        annotated_df.to_csv(annotated_output, index=False)
                        console.print(f"[bold green]Clinical Variant Annotations:[/] {annotated_output}")
                        
                        # Generate clinical summary
                        clinical_summary = clinical_integration.get_clinical_summary(variant_df)
                        summary_path = args.output_dir / "clinical_summary.json"
                        with open(summary_path, 'w') as f:
                            import json
                            json.dump(clinical_summary, f, indent=2)
                        console.print(f"[bold green]Clinical Summary:[/] {summary_path}")
                
                # Generate integrated AF report if enrichment was performed
                if hasattr(args, 'variant_enrichment') and args.variant_enrichment and 'af_pathogenic' in variant_df.columns:
                    log.info("Generating integrated AF variant report")
                    # Count pathogenic AF variants
                    pathogenic_count = variant_df['af_pathogenic'].sum() if 'af_pathogenic' in variant_df.columns else 0
                    af_related_count = variant_df['is_af_related'].sum() if 'is_af_related' in variant_df.columns else 0
                    
                    if pathogenic_count > 0:
                        console.print(f"[bold cyan]AF-Related Findings:[/] {af_related_count} variants related to atrial fibrillation")
                        console.print(f"[bold cyan]Potentially Pathogenic:[/] {pathogenic_count} variants with potential AF pathogenicity")
                        
                        # Add summary to the report
                        af_summary = {
                            "af_related_variants": int(af_related_count),
                            "af_pathogenic_variants": int(pathogenic_count),
                            "af_categories": variant_df['af_category'].value_counts().to_dict() if 'af_category' in variant_df.columns else {}
                        }
                        
                        af_summary_path = args.output_dir / "af_enrichment_summary.json"
                        with open(af_summary_path, 'w') as f:
                            import json
                            json.dump(af_summary, f, indent=2)
                        console.print(f"[bold cyan]AF Enrichment Summary:[/] {af_summary_path}")
        else:
            if args.no_visualization:
                log.info("Visualization generation skipped as requested")
            else:
                log.warning("Visualization skipped due to empty variant data")
        
        console.print("\n[bold green]Analysis completed successfully![/]")
        console.print(f"Results saved to: {args.output_dir}")
        
    except Exception as e:
        console.print(f"[bold red]Unhandled exception:[/] {e}")
        console.print(traceback.format_exc())
        sys.exit(1)
