#!/usr/bin/env python3
"""
Script to generate a combined HTML report for all GASLIT-AF direct analysis results.
"""

import os
import sys
import json
import argparse
import pandas as pd
from pathlib import Path
from datetime import datetime
import shutil
from rich.console import Console
from rich.logging import RichHandler
import logging

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af-report")

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="GASLIT-AF Combined Report Generator")
    
    parser.add_argument("--results-dir", type=str, default="analysis_results", 
                        help="Directory containing analysis results")
    parser.add_argument("--output-file", type=str, default="gaslit_af_combined_report.html", 
                        help="Output HTML file")
    parser.add_argument("--direct-only", action="store_true", default=True,
                        help="Only include direct analysis results")
    
    args = parser.parse_args()
    
    # Convert string paths to Path objects
    args.results_dir = Path(args.results_dir)
    
    return args

def find_direct_analysis_dirs(results_dir, direct_only=True):
    """Find all direct analysis directories."""
    analysis_dirs = []
    
    for item in results_dir.iterdir():
        if item.is_dir():
            if direct_only and not item.name.startswith("direct_"):
                continue
            
            # Check if this is a direct analysis directory
            system_analysis_file = item / "system_analysis.json"
            gene_counts_file = item / "gene_counts.csv"
            
            if system_analysis_file.exists() and gene_counts_file.exists():
                analysis_dirs.append(item)
    
    return analysis_dirs

def load_analysis_results(analysis_dir):
    """Load analysis results from a directory."""
    system_analysis_file = analysis_dir / "system_analysis.json"
    gene_counts_file = analysis_dir / "gene_counts.csv"
    
    # Load system analysis
    with open(system_analysis_file, 'r') as f:
        system_analysis = json.load(f)
    
    # Load gene counts
    gene_counts = pd.read_csv(gene_counts_file)
    
    # Get visualization paths
    system_distribution_png = analysis_dir / "system_distribution.png"
    system_distribution_pie_png = analysis_dir / "system_distribution_pie.png"
    top_genes_png = analysis_dir / "top_genes.png"
    
    # Check if visualization files exist
    visualizations = {
        "system_distribution": system_distribution_png if system_distribution_png.exists() else None,
        "system_distribution_pie": system_distribution_pie_png if system_distribution_pie_png.exists() else None,
        "top_genes": top_genes_png if top_genes_png.exists() else None
    }
    
    return {
        "name": analysis_dir.name,
        "system_analysis": system_analysis,
        "gene_counts": gene_counts,
        "visualizations": visualizations
    }

def generate_html_report(analysis_results, output_file):
    """Generate HTML report from analysis results."""
    # Create output directory if it doesn't exist
    output_dir = Path(output_file).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create images directory
    images_dir = output_dir / "images"
    images_dir.mkdir(exist_ok=True)
    
    # Copy visualization images
    for result in analysis_results:
        for viz_name, viz_path in result["visualizations"].items():
            if viz_path:
                dest_path = images_dir / f"{result['name']}_{viz_name}.png"
                shutil.copy(viz_path, dest_path)
    
    # Generate HTML
    html = f"""<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>GASLIT-AF Combined Analysis Report</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                margin: 0;
                padding: 20px;
                background-color: #f8f9fa;
            }}
            .container {{
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 30px;
                border-radius: 10px;
                box-shadow: 0 0 10px rgba(0, 0, 0, 0.1);
            }}
            h1, h2, h3, h4 {{
                color: #343a40;
            }}
            .section {{
                margin-bottom: 30px;
                padding-bottom: 20px;
                border-bottom: 1px solid #e9ecef;
            }}
            .system-card {{
                margin-bottom: 20px;
                border-radius: 5px;
                overflow: hidden;
            }}
            .system-card-header {{
                padding: 10px 15px;
                background-color: #007bff;
                color: white;
                font-weight: bold;
            }}
            .system-card-body {{
                padding: 15px;
                background-color: #f8f9fa;
            }}
            .gene-table {{
                width: 100%;
                margin-top: 20px;
            }}
            .gene-table th {{
                background-color: #007bff;
                color: white;
            }}
            .visualization {{
                margin-top: 20px;
                text-align: center;
            }}
            .visualization img {{
                max-width: 100%;
                height: auto;
                border-radius: 5px;
                box-shadow: 0 0 5px rgba(0, 0, 0, 0.2);
            }}
            .summary-table {{
                width: 100%;
                margin-top: 20px;
            }}
            .summary-table th {{
                background-color: #007bff;
                color: white;
            }}
            .nav-tabs {{
                margin-bottom: 20px;
            }}
            .tab-content {{
                padding: 20px;
                background-color: white;
                border: 1px solid #dee2e6;
                border-top: none;
                border-radius: 0 0 5px 5px;
            }}
        </style>
    </head>
    <body>
        <div class="container">
            <h1 class="text-center mb-4">GASLIT-AF Combined Analysis Report</h1>
            <p class="text-center text-muted">Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            
            <div class="section">
                <h2>Summary</h2>
                <p>This report combines the results of GASLIT-AF direct analysis on multiple VCF files.</p>
                
                <table class="table table-striped summary-table">
                    <thead>
                        <tr>
                            <th>VCF File</th>
                            <th>Total Variants</th>
                            <th>Top System</th>
                            <th>Top Gene</th>
                        </tr>
                    </thead>
                    <tbody>
    """
    
    # Add summary rows
    for result in analysis_results:
        # Get top system
        system_counts = result["system_analysis"]["system_counts"]
        top_system = max(system_counts.items(), key=lambda x: x[1])[0]
        
        # Get top gene
        top_gene_row = result["gene_counts"].iloc[0]
        top_gene = f"{top_gene_row['gene']} ({top_gene_row['count']})"
        
        html += f"""
                        <tr>
                            <td>{result["name"]}</td>
                            <td>{result["system_analysis"]["total_variants"]}</td>
                            <td>{top_system}</td>
                            <td>{top_gene}</td>
                        </tr>
        """
    
    html += """
                    </tbody>
                </table>
            </div>
            
            <div class="section">
                <h2>Detailed Results</h2>
                
                <ul class="nav nav-tabs" id="resultTabs" role="tablist">
    """
    
    # Add tab headers
    for i, result in enumerate(analysis_results):
        active = "active" if i == 0 else ""
        html += f"""
                    <li class="nav-item" role="presentation">
                        <button class="nav-link {active}" id="tab-{i}" data-bs-toggle="tab" data-bs-target="#content-{i}" type="button" role="tab" aria-controls="content-{i}" aria-selected="{str(i == 0).lower()}">{result["name"]}</button>
                    </li>
        """
    
    html += """
                </ul>
                
                <div class="tab-content" id="resultTabsContent">
    """
    
    # Add tab content
    for i, result in enumerate(analysis_results):
        active = "show active" if i == 0 else ""
        
        html += f"""
                    <div class="tab-pane fade {active}" id="content-{i}" role="tabpanel" aria-labelledby="tab-{i}">
                        <h3>{result["name"]}</h3>
                        <p>Total variants: {result["system_analysis"]["total_variants"]}</p>
                        
                        <div class="row">
                            <div class="col-md-6">
                                <h4>System Distribution</h4>
                                <div class="visualization">
                                    <img src="images/{result['name']}_system_distribution.png" alt="System Distribution">
                                </div>
                            </div>
                            <div class="col-md-6">
                                <h4>System Distribution (Pie Chart)</h4>
                                <div class="visualization">
                                    <img src="images/{result['name']}_system_distribution_pie.png" alt="System Distribution Pie Chart">
                                </div>
                            </div>
                        </div>
                        
                        <div class="row mt-4">
                            <div class="col-md-12">
                                <h4>Top Genes</h4>
                                <div class="visualization">
                                    <img src="images/{result['name']}_top_genes.png" alt="Top Genes">
                                </div>
                            </div>
                        </div>
                        
                        <div class="row mt-4">
                            <div class="col-md-12">
                                <h4>Top 20 Genes by Variant Count</h4>
                                <table class="table table-striped gene-table">
                                    <thead>
                                        <tr>
                                            <th>Gene</th>
                                            <th>Count</th>
                                            <th>System</th>
                                        </tr>
                                    </thead>
                                    <tbody>
        """
        
        # Add top 20 genes
        for _, row in result["gene_counts"].head(20).iterrows():
            html += f"""
                                        <tr>
                                            <td>{row['gene']}</td>
                                            <td>{row['count']}</td>
                                            <td>{row['system']}</td>
                                        </tr>
            """
        
        html += """
                                    </tbody>
                                </table>
                            </div>
                        </div>
                        
                        <div class="row mt-4">
                            <div class="col-md-12">
                                <h4>Biological Systems</h4>
        """
        
        # Add system cards
        for system, genes in result["system_analysis"]["system_genes"].items():
            count = result["system_analysis"]["system_counts"][system]
            percentage = result["system_analysis"]["system_percentages"][system]
            
            html += f"""
                                <div class="system-card">
                                    <div class="system-card-header">
                                        {system} ({count} variants, {percentage:.2f}%)
                                    </div>
                                    <div class="system-card-body">
                                        <p>Genes: {", ".join(genes)}</p>
                                    </div>
                                </div>
            """
        
        html += """
                            </div>
                        </div>
                    </div>
        """
    
    html += """
                </div>
            </div>
        </div>
        
        <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/js/bootstrap.bundle.min.js"></script>
    </body>
    </html>
    """
    
    # Write HTML to file
    with open(output_file, 'w') as f:
        f.write(html)
    
    return output_file

def main():
    """Main function."""
    args = parse_args()
    
    # Find direct analysis directories
    analysis_dirs = find_direct_analysis_dirs(args.results_dir, args.direct_only)
    
    if not analysis_dirs:
        log.error(f"No direct analysis results found in {args.results_dir}")
        sys.exit(1)
    
    log.info(f"Found {len(analysis_dirs)} direct analysis results:")
    for analysis_dir in analysis_dirs:
        log.info(f"  - {analysis_dir}")
    
    # Load analysis results
    analysis_results = []
    for analysis_dir in analysis_dirs:
        log.info(f"Loading results from {analysis_dir}")
        results = load_analysis_results(analysis_dir)
        analysis_results.append(results)
    
    # Generate HTML report
    output_file = generate_html_report(analysis_results, args.output_file)
    
    log.info(f"Generated combined report: {output_file}")
    log.info(f"Open the report in a browser to view the results")

if __name__ == "__main__":
    main()
