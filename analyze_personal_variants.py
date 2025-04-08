#!/usr/bin/env python3
"""
Personal Variant Analysis for GASLIT-AF Framework

This script analyzes personal variant data in JSON format and integrates it with
the GASLIT-AF framework, creating a quantum coherence bridge between personal
genomic architecture and the theoretical model parameters.
"""

import os
import json
import logging
import argparse
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional
from rich.console import Console
from rich.logging import RichHandler
from rich.table import Table

# Import GASLIT-AF modules
from src.gaslit_af.gene_lists import GASLIT_AF_GENES, KNOWN_SNPS
from src.gaslit_af.biological_systems import get_system_for_gene
from src.gaslit_af.visualization import generate_all_visualizations
from src.gaslit_af.enhanced_reporting import generate_enhanced_report
from src.gaslit_af.gene_mapping import GeneMapper, enhance_gene_mapping, create_gene_mapping_database

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af")

def parse_args():
    """Parse command-line arguments for personal variant analysis."""
    parser = argparse.ArgumentParser(description="Personal Variant Analysis for GASLIT-AF Framework")
    
    parser.add_argument("--variant-file", type=str, default="data/personal_variants.json",
                      help="Input JSON file with personal variants (default: data/personal_variants.json)")
    parser.add_argument("--output-dir", type=str, default="results/personal_analysis",
                      help="Output directory for analysis results (default: results/personal_analysis)")
    parser.add_argument("--format", type=str, choices=["standard", "recursive"], default="recursive",
                      help="Input format (default: recursive)")
    parser.add_argument("--enhanced-report", action="store_true",
                      help="Generate enhanced report with interactive visualizations")
    parser.add_argument("--open-browser", action="store_true",
                      help="Automatically open report in browser")
    parser.add_argument("--gene-db", type=str,
                      help="Path to gene mapping database (JSON)")
    parser.add_argument("--create-gene-db", action="store_true",
                      help="Create gene mapping database")
    parser.add_argument("--external-mappings", type=str,
                      help="Path to external gene mappings (JSON)")
    
    return parser.parse_args()

def load_personal_variants(variant_file: Path, format_type: str) -> List[Dict[str, Any]]:
    """
    Load personal variants from JSON file.
    
    Args:
        variant_file: Path to JSON file
        format_type: Format type (standard or recursive)
        
    Returns:
        List of variant dictionaries
    """
    with open(variant_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    if format_type == "standard":
        return data.get("variants", [])
    else:  # recursive
        # Extract variants from recursive structure
        variants = []
        
        # Extract from parameters
        for param, param_data in data.get("gaslit_af_model", {}).get("parameters", {}).items():
            for variant in param_data.get("variants", []):
                if variant not in variants:
                    variants.append(variant)
        
        # Extract from biological systems
        for system, system_data in data.get("gaslit_af_model", {}).get("biological_systems", {}).items():
            for variant in system_data.get("variants", []):
                if variant not in variants:
                    variants.append(variant)
        
        return variants

def convert_to_dataframe(variants: List[Dict[str, Any]]) -> pd.DataFrame:
    """
    Convert variant list to DataFrame.
    
    Args:
        variants: List of variant dictionaries
        
    Returns:
        DataFrame with variant data
    """
    # Extract key fields for DataFrame
    variant_data = []
    for variant in variants:
        row = {
            "chrom": "",  # Would need to extract from variant_id or additional data
            "pos": 0,     # Would need to extract from variant_id or additional data
            "ref": variant.get("genotype", "")[0] if variant.get("genotype", "") else "",
            "alt": variant.get("risk_allele", ""),
            "gene": variant.get("gene", ""),
            "rsid": variant.get("variant_id", ""),
            "genotype": variant.get("genotype", ""),
            "quality": 100,  # Default quality
            "impact": variant.get("classification", ""),
            "trait": variant.get("condition", "")
        }
        
        # Add GASLIT-AF parameters if available
        gaslit_params = variant.get("gaslit_parameters", {})
        for param, value in gaslit_params.items():
            row[f"gaslit_{param}"] = value
        
        variant_data.append(row)
    
    return pd.DataFrame(variant_data)

def analyze_variants(variant_df: pd.DataFrame, gene_db_path: Optional[str] = None, external_mappings_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Analyze variants and generate summary statistics with enhanced gene mapping.
    
    Args:
        variant_df: DataFrame with variant data
        gene_db_path: Path to gene database file (optional)
        external_mappings_path: Path to external mappings file (optional)
        
    Returns:
        Dictionary with analysis results
    """
    # Convert DataFrame to list of dictionaries for gene mapper
    variants_list = variant_df.to_dict('records')
    
    # Create gene mapper
    mapper = GeneMapper(gene_db_path, external_mappings_path)
    
    # Enhance variants with improved gene mapping
    log.info("Enhancing gene mapping for variants")
    enhanced_variants = mapper.enhance_variants(variants_list)
    
    # Convert enhanced variants back to DataFrame for analysis
    enhanced_df = pd.DataFrame(enhanced_variants)
    
    # Count variants by gene with enhanced mapping
    gene_counts = variant_df["gene"].value_counts().to_dict()
    
    # Count variants by classification
    classification_counts = {}
    if "impact" in variant_df.columns:
        classification_counts = variant_df["impact"].value_counts().to_dict()
    elif "classification" in variant_df.columns:
        classification_counts = variant_df["classification"].value_counts().to_dict()
    
    # Count variants by condition
    condition_counts = {}
    if "trait" in variant_df.columns:
        condition_counts = variant_df["trait"].value_counts().to_dict()
    
    # Extract pathway information from enhanced variants
    pathway_counts = {}
    for variant in enhanced_variants:
        pathways = variant.get("pathways", [])
        if isinstance(pathways, list):
            for pathway in pathways:
                if pathway:
                    pathway_counts[pathway] = pathway_counts.get(pathway, 0) + 1
    
    # Calculate GASLIT-AF parameters from enhanced variants
    # First check for existing parameters in the DataFrame
    gaslit_params = {}
    for param in ["gaslit_γ", "gaslit_Λ", "gaslit_Ω", "gaslit_Χ", "gaslit_σ"]:
        if param in variant_df.columns:
            gaslit_params[param] = variant_df[param].mean()
    
    # If no parameters in DataFrame, calculate from enhanced variants
    if not gaslit_params:
        param_sums = {
            "gaslit_γ": 0.0,  # Genetic fragility
            "gaslit_Λ": 0.0,  # Allostatic load
            "gaslit_Ω": 0.0,  # Endocannabinoid buffering
            "gaslit_Χ": 0.0,  # Physiological coherence
            "gaslit_σ": 0.0   # Entropy production
        }
        param_counts = {k: 0 for k in param_sums.keys()}
        
        # Extract parameter scores from enhanced variants
        for variant in enhanced_variants:
            param_scores = variant.get('parameter_scores', {})
            if param_scores:
                if 'γ' in param_scores and param_scores['γ'] > 0:
                    param_sums["gaslit_γ"] += param_scores['γ']
                    param_counts["gaslit_γ"] += 1
                if 'Λ' in param_scores and param_scores['Λ'] > 0:
                    param_sums["gaslit_Λ"] += param_scores['Λ']
                    param_counts["gaslit_Λ"] += 1
                if 'Ω' in param_scores and param_scores['Ω'] > 0:
                    param_sums["gaslit_Ω"] += param_scores['Ω']
                    param_counts["gaslit_Ω"] += 1
                if 'Χ' in param_scores and param_scores['Χ'] > 0:
                    param_sums["gaslit_Χ"] += param_scores['Χ']
                    param_counts["gaslit_Χ"] += 1
                if 'σ' in param_scores and param_scores['σ'] > 0:
                    param_sums["gaslit_σ"] += param_scores['σ']
                    param_counts["gaslit_σ"] += 1
        
        # Calculate average parameter values
        for param in param_sums:
            if param_counts[param] > 0:
                gaslit_params[param] = param_sums[param] / param_counts[param]
            else:
                # Default values if no enhanced data
                if param == "gaslit_γ":
                    gaslit_params[param] = 0.01
                elif param == "gaslit_Λ":
                    gaslit_params[param] = 0.00
                elif param == "gaslit_Ω":
                    gaslit_params[param] = 0.50
                elif param == "gaslit_Χ":
                    gaslit_params[param] = 0.50
                elif param == "gaslit_σ":
                    gaslit_params[param] = 0.00
    
    # Group variants by biological system using enhanced mapping
    system_variants = {}
    
    for variant in enhanced_variants:
        system = variant.get('biological_system', None)
        if not system:
            gene = variant.get('gene', '')
            if gene:
                system = get_system_for_gene(gene)
            else:
                system = "Other"
        
        if system not in system_variants:
            system_variants[system] = []
        
        system_variants[system].append(variant)
    
    # Count variants by system
    system_counts = {system: len(variants) for system, variants in system_variants.items()}
    
    # Return comprehensive analysis results
    return {
        "total_variants": len(variant_df),
        "gene_counts": gene_counts,
        "classification_counts": classification_counts,
        "condition_counts": condition_counts,
        "system_counts": system_counts,
        "pathway_counts": pathway_counts,
        "gaslit_params": gaslit_params,
        "system_variants": system_variants,
        "enhanced_variants": enhanced_variants
    }

def display_analysis_results(analysis_results: Dict[str, Any]):
    """
    Display analysis results in the console.
    
    Args:
        analysis_results: Dictionary with analysis results
    """
    console.print("\n[bold cyan]GASLIT-AF Personal Variant Analysis[/]")
    console.print(f"Total Variants: {analysis_results['total_variants']}")
    
    # Display gene counts
    gene_table = Table(title="Variants by Gene")
    gene_table.add_column("Gene", style="cyan")
    gene_table.add_column("Count", style="green")
    
    for gene, count in sorted(analysis_results["gene_counts"].items(), key=lambda x: x[1], reverse=True):
        gene_table.add_row(gene, str(count))
    
    console.print(gene_table)
    
    # Display system counts
    system_table = Table(title="Variants by Biological System")
    system_table.add_column("System", style="cyan")
    system_table.add_column("Count", style="green")
    
    for system, count in sorted(analysis_results["system_counts"].items(), key=lambda x: x[1], reverse=True):
        system_table.add_row(system, str(count))
    
    console.print(system_table)
    
    # Display GASLIT-AF parameters
    if analysis_results["gaslit_params"]:
        param_table = Table(title="GASLIT-AF Parameters (Average)")
        param_table.add_column("Parameter", style="cyan")
        param_table.add_column("Value", style="green")
        
        for param, value in analysis_results["gaslit_params"].items():
            param_name = param.replace("gaslit_", "")
            param_table.add_row(param_name, f"{value:.2f}")
        
        console.print(param_table)
    
    # Display classification counts
    class_table = Table(title="Variants by Classification")
    class_table.add_column("Classification", style="cyan")
    class_table.add_column("Count", style="green")
    
    for classification, count in sorted(analysis_results["classification_counts"].items(), key=lambda x: x[1], reverse=True):
        class_table.add_row(classification, str(count))
    
    console.print(class_table)
    
    # Display pathway counts if available
    if "pathway_counts" in analysis_results and analysis_results["pathway_counts"]:
        pathway_table = Table(title="Top Biological Pathways")
        pathway_table.add_column("Pathway", style="cyan")
        pathway_table.add_column("Count", style="green")
        
        for pathway, count in sorted(analysis_results["pathway_counts"].items(), key=lambda x: x[1], reverse=True)[:10]:  # Show top 10
            pathway_table.add_row(pathway, str(count))
        
        console.print(pathway_table)

def main():
    """Main entry point for personal variant analysis."""
    args = parse_args()
    
    # Ensure variant file exists
    variant_file = Path(args.variant_file)
    if not variant_file.exists():
        log.error(f"Variant file not found: {variant_file}")
        return
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load personal variants
    log.info(f"Loading personal variants from {variant_file}")
    variants = load_personal_variants(variant_file, args.format)
    log.info(f"Loaded {len(variants)} variants")
    
    # Convert to DataFrame
    variant_df = convert_to_dataframe(variants)
    
    # Create gene mapping database if requested
    gene_db_path = None
    if args.create_gene_db:
        gene_db_path = output_dir / "gene_mapping_database.json"
        log.info(f"Creating gene mapping database at {gene_db_path}")
        create_gene_mapping_database(str(gene_db_path), variant_df)
    elif args.gene_db:
        gene_db_path = args.gene_db
    
    # Save variant DataFrame
    variant_csv = output_dir / "personal_variants.csv"
    variant_df.to_csv(variant_csv, index=False)
    log.info(f"Saved variant data to {variant_csv}")
    
    # Analyze variants with enhanced gene mapping
    log.info("Analyzing variants with enhanced gene mapping")
    analysis_results = analyze_variants(
        variant_df=variant_df,
        gene_db_path=gene_db_path,
        external_mappings_path=args.external_mappings
    )
    
    # Save analysis results
    results_json = output_dir / "analysis_results.json"
    with open(results_json, 'w', encoding='utf-8') as f:
        json.dump(analysis_results, f, indent=2)
    log.info(f"Saved analysis results to {results_json}")
    
    # Save enhanced variants if available
    if "enhanced_variants" in analysis_results:
        enhanced_variants_path = output_dir / "enhanced_variants.json"
        with open(enhanced_variants_path, 'w', encoding='utf-8') as f:
            json.dump(analysis_results["enhanced_variants"], f, indent=2)
        log.info(f"Saved enhanced variants to {enhanced_variants_path}")
        
        # Convert enhanced variants to DataFrame
        enhanced_df = pd.DataFrame(analysis_results["enhanced_variants"])
        enhanced_csv_path = output_dir / "enhanced_variants.csv"
        enhanced_df.to_csv(enhanced_csv_path, index=False)
        log.info(f"Saved enhanced variants to {enhanced_csv_path}")
    
    # Display analysis results
    display_analysis_results(analysis_results)
    
    # Generate visualizations
    log.info("Generating visualizations")
    vis_dir = output_dir / "visualizations"
    vis_dir.mkdir(exist_ok=True)
    
    # Use enhanced variant data for visualizations if available
    visualization_df = variant_df
    if "enhanced_variants" in analysis_results:
        enhanced_df = pd.DataFrame(analysis_results["enhanced_variants"])
        if not enhanced_df.empty:
            visualization_df = enhanced_df
    
    figures = generate_all_visualizations(
        variant_data=visualization_df,
        gene_counts=analysis_results["gene_counts"],
        output_dir=vis_dir
    )
    
    # Generate report
    if args.enhanced_report:
        log.info("Generating enhanced report")
        report_path = generate_enhanced_report(
            gene_counts=analysis_results["gene_counts"],
            variant_data=visualization_df,
            figures=figures,
            output_dir=output_dir,
            system_analysis=analysis_results["system_variants"],
            include_symptoms=True
        )
        
        log.info(f"Enhanced report generated: {report_path}")
        
        # Open report in browser if requested
        if args.open_browser and report_path:
            import webbrowser
            webbrowser.open(f"file://{report_path}")
    
    log.info("[bold green]Personal variant analysis completed successfully with enhanced gene mapping![/]")

if __name__ == "__main__":
    main()
