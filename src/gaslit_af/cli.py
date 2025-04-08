"""
Command-line interface module for GASLIT-AF Variant Analysis.
Handles argument parsing and configuration.
"""

import argparse
from pathlib import Path
from src.gaslit_af.workflow import run_analysis_workflow


def parse_args():
    """
    Parse command-line arguments for GASLIT-AF Variant Analysis.
    
    Returns:
        Parsed arguments namespace
    """
    parser = argparse.ArgumentParser(description="GASLIT-AF Variant Analysis")
    
    # Input/output arguments
    parser.add_argument("vcf_path", help="Path to VCF file for analysis")
    parser.add_argument("--output-dir", type=str, default="output", 
                        help="Output directory for results")
    
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
    parser.add_argument("--system-analysis", action="store_true", 
                        help="Perform biological system analysis")
    parser.add_argument("--use-pysam", action="store_true", 
                        help="Use pysam for advanced variant processing")
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
    parser.add_argument("--annovar-protocols", nargs="+", 
                        default=["refGene", "exac03", "gnomad211_exome", "clinvar_20220320", "dbnsfp42a"],
                        help="ANNOVAR annotation protocols to use")
    parser.add_argument("--annovar-operations", nargs="+",
                        help="ANNOVAR operations for each protocol (must match number of protocols)")
    parser.add_argument("--download-annovar-dbs", action="store_true",
                        help="Download required ANNOVAR databases if not present")
    
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
    parser.add_argument("--open-browser", action="store_true", 
                        help="Automatically open report in browser")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert string paths to Path objects
    args.output_dir = Path(args.output_dir)
    args.cache_dir = Path(args.cache_dir)
    
    return args


if __name__ == "__main__":
    # Parse command-line arguments
    parsed_args = parse_args()
    
    # Run the main analysis workflow
    run_analysis_workflow(parsed_args)
