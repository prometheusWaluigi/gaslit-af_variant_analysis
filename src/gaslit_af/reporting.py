"""
Reporting module for GASLIT-AF Variant Analysis.
Provides functions to generate HTML reports with visualizations and analysis results.
"""

import os
import time
import json
import pandas as pd
import plotly.io as pio
from jinja2 import Template
import logging
from pathlib import Path
from typing import Optional, Dict, List

# Configure logging
log = logging.getLogger("gaslit-af")

def generate_html_report(
    variant_df: Optional[pd.DataFrame] = None,
    output_dir: str = ".",
    args: Optional[object] = None, 
    system_results: Optional[Dict] = None,
    rccx_results: Optional[List[Dict]] = None, 
    figures: Optional[Dict] = None, 
    gene_counts: Optional[Dict] = None 
):
    """
    Generate an interactive HTML report with analysis results and visualizations.
    
    Args:
        variant_df: DataFrame with processed variant information.
        output_dir: Directory to save the report.
        args: Namespace object containing command line arguments.
        system_results: Dictionary containing results from biological system analysis.
        rccx_results: List of dictionaries containing RCCX SV/CNV findings.
        figures: (Optional/Legacy) Dictionary of plotly figures.
        gene_counts: (Optional/Legacy) Dictionary of gene:count pairs.
    
    Returns:
        Path to the generated HTML report
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    report_path = output_dir / f"gaslit_af_report_{timestamp}.html"

    log.info(f"Generating HTML report: {report_path}")

    # Calculate gene counts if variant_df is provided and gene_counts isn't
    if gene_counts is None and variant_df is not None and not variant_df.empty and 'GENE' in variant_df.columns:
        gene_counts = variant_df['GENE'].value_counts().to_dict()
    elif gene_counts is None:
        gene_counts = {}

    if variant_df is None:
        variant_df = pd.DataFrame() 

    # Prepare top genes from counts
    if gene_counts:
        top_genes = pd.DataFrame(list(gene_counts.items()),
                               columns=['Gene', 'VariantCount'])
        top_genes = top_genes.sort_values('VariantCount', ascending=False).head(20)
    else:
        top_genes = pd.DataFrame(columns=['Gene', 'VariantCount'])

    # Convert figures to HTML (handle potential None)
    plot_html = {}
    if figures:
        for name, fig in figures.items():
            if fig is not None:
                try:
                    plot_html[name] = pio.to_html(fig, full_html=False, include_plotlyjs='cdn')
                except Exception as e:
                    log.warning(f"Failed to convert plotly figure '{name}' to HTML: {e}")
                    plot_html[name] = "<p>Error generating plot.</p>"

    # Basic statistics
    stats = {
        'total_genes': len(gene_counts),
        'total_variants': len(variant_df) if not variant_df.empty else 0,
        'analysis_time': time.strftime("%Y-%m-%d %H:%M:%S"),
        'top_gene': top_genes.iloc[0]['Gene'] if not top_genes.empty else 'N/A',
        'top_gene_count': int(top_genes.iloc[0]['VariantCount']) if not top_genes.empty else 0,
        'vcf_path': args.vcf_path if args else 'N/A'
    }

    # Add system analysis stats if available
    if system_results and 'system_counts' in system_results:
        stats['total_systems'] = len(system_results['system_counts'])
        top_system = max(system_results['system_counts'].items(), key=lambda item: item[1]) if system_results['system_counts'] else ('N/A', 0)
        stats['top_system'] = top_system[0]
        stats['top_system_count'] = top_system[1]
    else:
        stats['total_systems'] = 0
        stats['top_system'] = 'N/A'
        stats['top_system_count'] = 0

    # HTML template
    template_str = """
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>GASLIT-AF Variant Analysis Report</title>
        <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0-alpha1/dist/css/bootstrap.min.css" rel="stylesheet">
        <style>
            body { 
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                margin: 0;
                padding: 20px;
                background-color: #f8f9fa;
            }
            .container {
                max-width: 1200px;
                margin: 0 auto;
                background-color: white;
                padding: 30px;
                border-radius: 10px;
                box-shadow: 0 0 20px rgba(0,0,0,0.1);
            }
            h1 { 
                color: #2c3e50;
                border-bottom: 2px solid #3498db;
                padding-bottom: 10px;
                margin-bottom: 30px;
            }
            h2 {
                color: #2c3e50;
                margin-top: 40px;
                margin-bottom: 20px;
                border-left: 5px solid #3498db;
                padding-left: 15px;
            }
            .stats-card {
                background-color: #f1f8fe;
                border-radius: 10px;
                padding: 20px;
                margin-bottom: 30px;
                box-shadow: 0 0 10px rgba(0,0,0,0.05);
            }
            .stat-item {
                margin-bottom: 15px;
            }
            .stat-label {
                font-weight: bold;
                color: #3498db;
            }
            .stat-value {
                font-size: 1.2em;
                color: #2c3e50;
            }
            .plot-container {
                margin: 30px 0;
                border: 1px solid #e0e0e0;
                border-radius: 10px;
                padding: 20px;
                background-color: white;
            }
            .table-container {
                overflow-x: auto;
                margin-top: 20px;
            }
            table {
                width: 100%;
                border-collapse: collapse;
                margin-bottom: 20px;
            }
            th, td {
                padding: 12px 15px;
                border: 1px solid #dee2e6;
                text-align: left;
            }
            th {
                background-color: #e9ecef;
                font-weight: bold;
            }
            tbody tr:nth-child(even) {
                background-color: #f8f9fa;
            }
            tbody tr:hover {
                background-color: #e2e6ea;
            }
            pre {
                background-color: #eee;
                border: 1px solid #ccc;
                padding: 10px;
                border-radius: 5px;
                overflow-x: auto;
                white-space: pre-wrap; /* Wrap long lines */
                word-wrap: break-word; /* Break words if necessary */
            }
        </style>
    </head>
    <body>
        <div class="container">
            <h1>GASLIT-AF Variant Analysis Report</h1>
            
            <div class="stats-card">
                <h2>Analysis Summary</h2>
                <div class="row">
                    <div class="col-md-4">
                        <div class="stat-item">
                            <div class="stat-label">Total GASLIT-AF Genes</div>
                            <div class="stat-value">{{ stats.total_genes }}</div>
                        </div>
                    </div>
                    <div class="col-md-4">
                        <div class="stat-item">
                            <div class="stat-label">Total Variants</div>
                            <div class="stat-value">{{ stats.total_variants }}</div>
                        </div>
                    </div>
                    <div class="col-md-4">
                        <div class="stat-item">
                            <div class="stat-label">Analysis Time</div>
                            <div class="stat-value">{{ stats.analysis_time }}</div>
                        </div>
                    </div>
                </div>
                <div class="row mt-3">
                    <div class="col-md-6">
                        <div class="stat-item">
                            <div class="stat-label">Top Gene</div>
                            <div class="stat-value">{{ stats.top_gene }} ({{ stats.top_gene_count }} variants)</div>
                        </div>
                    </div>
                    <div class="col-md-6">
                        <div class="stat-item">
                            <div class="stat-label">VCF File</div>
                            <div class="stat-value">{{ stats.vcf_path }}</div>
                        </div>
                    </div>
                    {% if stats.total_systems > 0 %}
                    <div class="col-md-6">
                        <div class="stat-item">
                            <div class="stat-label">Top Biological System</div>
                            <div class="stat-value">{{ stats.top_system }} ({{ stats.top_system_count }} variants)</div>
                        </div>
                    </div>
                    {% endif %}
                </div>
            </div>

            <h2>Chromosome Distribution</h2>
            <div class="plot-container">
                {{ plots.chromosome_distribution | safe if plots.chromosome_distribution else "No data available" }}
            </div>
            
            <h2>Variant Type Distribution</h2>
            <div class="plot-container">
                {{ plots.variant_type_distribution | safe if plots.variant_type_distribution else "No data available" }}
            </div>
            
            <h2>Transition/Transversion Analysis</h2>
            <div class="plot-container">
                {{ plots.transition_transversion | safe if plots.transition_transversion else "No data available" }}
            </div>
            
            <h2>Top Variant-Enriched Genes</h2>
            <div class="plot-container">
                {{ plots.gene_counts | safe if plots.gene_counts else "No data available" }}
            </div>
            
            <h2>Gene Variant Network</h2>
            <div class="plot-container">
                {{ plots.gene_network | safe if plots.gene_network else "No data available" }}
            </div>
            
            {% if plots.system_pie or plots.system_bar or plots.system_sunburst %}
            <h2>Biological Systems Analysis</h2>
            
            {% if plots.system_pie %}
            <div class="plot-container">
                <h3>Variant Distribution by Biological System</h3>
                {{ plots.system_pie | safe }}
            </div>
            {% endif %}
            
            {% if plots.system_bar %}
            <div class="plot-container">
                <h3>Variant Counts by Biological System</h3>
                {{ plots.system_bar | safe }}
            </div>
            {% endif %}
            
            {% if plots.system_sunburst %}
            <div class="plot-container">
                <h3>Hierarchical View of Variants by System and Gene</h3>
                {{ plots.system_sunburst | safe }}
            </div>
            {% endif %}
            
            {% if plots.system_heatmap %}
            <div class="plot-container">
                <h3>Top Genes by Biological System</h3>
                {{ plots.system_heatmap | safe }}
            </div>
            {% endif %}
            {% endif %}
            
            {% if rccx_results %}
            <h2>RCCX Locus Structural Variant Scan (from VCF)</h2>
            <p>Potential structural variants identified within the RCCX region ({{ RCCX_REGION_GRCH38.chrom }}:{{ "{:,}".format(RCCX_REGION_GRCH38.start) }}-{{ "{:,}".format(RCCX_REGION_GRCH38.end) }}) based on VCF annotations.</p>
            <div class="table-container">
                <table class="table table-striped table-bordered">
                    <thead>
                        <tr>
                            <th>Type</th>
                            <th>Chromosome</th>
                            <th>Start</th>
                            <th>End</th>
                            <th>Approx. Length</th>
                            <th>VCF Record Details</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for sv in rccx_results %}
                        <tr>
                            <td>{{ sv.type }}</td>
                            <td>{{ sv.chrom }}</td>
                            <td>{{ "{:,}".format(sv.start) }}</td>
                            <td>{{ "{:,}".format(sv.end) }}</td>
                            <td>{{ "{:,}".format(sv.length) if sv.length is not none else 'N/A' }}</td>
                            <td><pre>{{ sv.record_details }}</pre></td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
            {% endif %}
            
            <h2>Top 20 Genes by Variant Count</h2>
            <div class="table-container">
                <table class="table">
                    <thead>
                        <tr>
                            <th>Gene</th>
                            <th>Variant Count</th>
                        </tr>
                    </thead>
                    <tbody>
                        {% for _, row in top_genes.iterrows() %}
                        <tr>
                            <td>{{ row.Gene }}</td>
                            <td>{{ row.VariantCount }}</td>
                        </tr>
                        {% endfor %}
                    </tbody>
                </table>
            </div>
            
            <div class="footer">
                <p>GASLIT-AF Variant Analysis | Generated on {{ stats.analysis_time }}</p>
            </div>
        </div>
    </body>
    </html>
    """

    # Render template
    template = Template(template_str)
    html_content = template.render(
        stats=stats,
        plots=plot_html,
        top_genes=top_genes,
        rccx_results=rccx_results,
        RCCX_REGION_GRCH38={"chrom": "chr6", "start": 31950000, "end": 32050000}  # Default RCCX region coordinates
    )

    # Write to file
    with open(report_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    log.info(f"HTML report generated successfully: {report_path}")
    return report_path
