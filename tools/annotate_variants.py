#!/usr/bin/env python3
"""
Variant Annotation Tool for GASLIT-AF

This tool demonstrates the API integration capabilities by annotating variants
from a CSV file or directly by rsID.
"""

import os
import sys
import json
import argparse
import pandas as pd
from pathlib import Path
from rich.console import Console
from rich.table import Table
from rich.progress import Progress

# Add project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import API integration module
from src.gaslit_af.api_integration import VariantAPIIntegration

# Initialize console
console = Console()

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Variant Annotation Tool for GASLIT-AF")
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest="command", help="Command to run")
    
    # Annotate file command
    file_parser = subparsers.add_parser("file", help="Annotate variants from a CSV file")
    file_parser.add_argument("input_file", help="Input CSV file with variants")
    file_parser.add_argument("--output-file", help="Output CSV file for annotated variants")
    file_parser.add_argument("--id-column", default="rsid", help="Column name for variant IDs")
    file_parser.add_argument("--sources", nargs="+", default=["ensembl", "myvariant"],
                            help="Sources to use for annotation")
    
    # Annotate variant command
    variant_parser = subparsers.add_parser("variant", help="Annotate a single variant by rsID")
    variant_parser.add_argument("variant_id", help="Variant ID (e.g., rs429358)")
    variant_parser.add_argument("--output-file", help="Output JSON file for variant details")
    variant_parser.add_argument("--format", choices=["json", "table"], default="table",
                               help="Output format (json or table)")
    
    # Batch annotate variants command
    batch_parser = subparsers.add_parser("batch", help="Annotate multiple variants by rsID")
    batch_parser.add_argument("variant_ids", nargs="+", help="Variant IDs (e.g., rs429358 rs7412)")
    batch_parser.add_argument("--output-file", help="Output JSON file for variant details")
    batch_parser.add_argument("--format", choices=["json", "table"], default="table",
                             help="Output format (json or table)")
    
    # Common options
    for p in [file_parser, variant_parser, batch_parser]:
        p.add_argument("--cache-dir", help="Directory to cache API responses")
        p.add_argument("--cache-ttl", type=int, default=24, 
                      help="Cache time-to-live in hours")
    
    return parser.parse_args()

def annotate_file(args):
    """Annotate variants from a CSV file."""
    console.print(f"[bold]Annotating variants from file:[/] {args.input_file}")
    
    try:
        # Load variants from CSV
        df = pd.read_csv(args.input_file)
        console.print(f"Loaded {len(df)} variants from file")
        
        # Check if ID column exists
        if args.id_column not in df.columns:
            console.print(f"[bold red]Error:[/] Column '{args.id_column}' not found in file")
            console.print(f"Available columns: {', '.join(df.columns)}")
            return
        
        # Initialize API integration
        api = VariantAPIIntegration(cache_dir=args.cache_dir, cache_ttl=args.cache_ttl)
        
        # Annotate variants
        with Progress() as progress:
            task = progress.add_task("[cyan]Annotating variants...", total=1)
            
            annotated_df = api.annotate_variants(df, sources=args.sources)
            
            progress.update(task, completed=1)
        
        # Save annotated variants
        output_file = args.output_file or Path(args.input_file).with_suffix('.annotated.csv')
        annotated_df.to_csv(output_file, index=False)
        
        console.print(f"[bold green]Annotated variants saved to:[/] {output_file}")
        
        # Print summary of annotations
        table = Table(title="Annotation Summary")
        table.add_column("Source", style="cyan")
        table.add_column("Column", style="green")
        table.add_column("Non-null Values", style="magenta")
        table.add_column("Example", style="yellow")
        
        # Ensembl columns
        ensembl_cols = [col for col in annotated_df.columns if col.startswith('ensembl_')]
        for col in ensembl_cols:
            non_null = annotated_df[col].notna().sum()
            example = annotated_df[col].dropna().iloc[0] if non_null > 0 else ""
            table.add_row("Ensembl", col, str(non_null), str(example)[:50])
        
        # MyVariant columns
        myvariant_cols = ['clinvar_significance', 'cadd_phred', 'gnomad_af', 'sift_pred', 'polyphen_pred']
        for col in myvariant_cols:
            if col in annotated_df.columns:
                non_null = annotated_df[col].notna().sum()
                example = annotated_df[col].dropna().iloc[0] if non_null > 0 else ""
                table.add_row("MyVariant", col, str(non_null), str(example)[:50])
        
        console.print(table)
        
    except Exception as e:
        console.print(f"[bold red]Error:[/] {e}")

def annotate_variant(args):
    """Annotate a single variant by rsID."""
    console.print(f"[bold]Annotating variant:[/] {args.variant_id}")
    
    try:
        # Initialize API integration
        api = VariantAPIIntegration(cache_dir=args.cache_dir, cache_ttl=args.cache_ttl)
        
        # Get variant details
        details = api.get_variant_details(args.variant_id)
        
        # Output format
        if args.format == "json":
            # Save to file if specified
            if args.output_file:
                with open(args.output_file, 'w') as f:
                    json.dump(details, f, indent=2)
                console.print(f"[bold green]Variant details saved to:[/] {args.output_file}")
            else:
                console.print(json.dumps(details, indent=2))
        else:
            # Print as table
            table = Table(title=f"Variant: {args.variant_id}")
            table.add_column("Attribute", style="cyan")
            table.add_column("Value", style="green")
            
            # Add basic details
            table.add_row("Gene", str(details.get("gene")))
            table.add_row("Consequence", str(details.get("consequence")))
            table.add_row("Impact", str(details.get("impact")))
            table.add_row("Clinical Significance", str(details.get("clinical_significance")))
            table.add_row("Allele Frequency", str(details.get("allele_frequency")))
            
            # Add pathogenicity scores
            for score, value in details.get("pathogenicity_scores", {}).items():
                table.add_row(f"Score: {score}", str(value))
            
            console.print(table)
            
            # Save to file if specified
            if args.output_file:
                with open(args.output_file, 'w') as f:
                    json.dump(details, f, indent=2)
                console.print(f"[bold green]Variant details saved to:[/] {args.output_file}")
        
    except Exception as e:
        console.print(f"[bold red]Error:[/] {e}")

def annotate_batch(args):
    """Annotate multiple variants by rsID."""
    console.print(f"[bold]Annotating {len(args.variant_ids)} variants:[/] {', '.join(args.variant_ids)}")
    
    try:
        # Initialize API integration
        api = VariantAPIIntegration(cache_dir=args.cache_dir, cache_ttl=args.cache_ttl)
        
        # Get variant details for each variant
        all_details = []
        
        with Progress() as progress:
            task = progress.add_task("[cyan]Annotating variants...", total=len(args.variant_ids))
            
            for variant_id in args.variant_ids:
                details = api.get_variant_details(variant_id)
                all_details.append(details)
                progress.update(task, advance=1)
        
        # Output format
        if args.format == "json":
            # Save to file if specified
            if args.output_file:
                with open(args.output_file, 'w') as f:
                    json.dump(all_details, f, indent=2)
                console.print(f"[bold green]Variant details saved to:[/] {args.output_file}")
            else:
                console.print(json.dumps(all_details, indent=2))
        else:
            # Print as table
            for details in all_details:
                variant_id = details.get("variant_id")
                table = Table(title=f"Variant: {variant_id}")
                table.add_column("Attribute", style="cyan")
                table.add_column("Value", style="green")
                
                # Add basic details
                table.add_row("Gene", str(details.get("gene")))
                table.add_row("Consequence", str(details.get("consequence")))
                table.add_row("Impact", str(details.get("impact")))
                table.add_row("Clinical Significance", str(details.get("clinical_significance")))
                table.add_row("Allele Frequency", str(details.get("allele_frequency")))
                
                # Add pathogenicity scores
                for score, value in details.get("pathogenicity_scores", {}).items():
                    table.add_row(f"Score: {score}", str(value))
                
                console.print(table)
                console.print("")
            
            # Save to file if specified
            if args.output_file:
                with open(args.output_file, 'w') as f:
                    json.dump(all_details, f, indent=2)
                console.print(f"[bold green]Variant details saved to:[/] {args.output_file}")
        
    except Exception as e:
        console.print(f"[bold red]Error:[/] {e}")

def main():
    """Main function."""
    args = parse_args()
    
    if args.command == "file":
        annotate_file(args)
    elif args.command == "variant":
        annotate_variant(args)
    elif args.command == "batch":
        annotate_batch(args)
    else:
        console.print("[bold red]Error:[/] No command specified")
        console.print("Use --help for usage information")

if __name__ == "__main__":
    main()
