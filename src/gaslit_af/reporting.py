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

# Configure logging
log = logging.getLogger("gaslit-af")

def generate_html_report(gene_counts, variant_data, figures, output_dir, system_analysis=None):
    """
    Generate an interactive HTML report with analysis results and visualizations.
    
    Args:
        gene_counts: Dictionary of gene:count pairs
        variant_data: DataFrame with variant information
        figures: Dictionary of plotly figures
        output_dir: Directory to save the report
    
    Returns:
        Path to the generated HTML report
    """
    os.makedirs(output_dir, exist_ok=True)
    timestamp = time.strftime("%Y%m%d-%H%M%S")
    report_path = os.path.join(output_dir, f"gaslit_af_report_{timestamp}.html")
    
    log.info(f"Generating HTML report: {report_path}")
    
    # Convert figures to HTML
    plot_html = {}
    for name, fig in figures.items():
        if fig is not None:
            plot_html[name] = pio.to_html(fig, full_html=False, include_plotlyjs='cdn')
    
    # Prepare data for the report
    top_genes = pd.DataFrame(list(gene_counts.items()), 
                           columns=['Gene', 'VariantCount'])
    top_genes = top_genes.sort_values('VariantCount', ascending=False).head(20)
    
    # Basic statistics
    stats = {
        'total_genes': len(gene_counts),
        'total_variants': len(variant_data) if variant_data is not None else 0,
        'analysis_time': time.strftime("%Y-%m-%d %H:%M:%S"),
        'top_gene': top_genes.iloc[0]['Gene'] if not top_genes.empty else 'None',
        'top_gene_count': top_genes.iloc[0]['VariantCount'] if not top_genes.empty else 0
    }
    
    # Add system analysis stats if available
    if system_analysis:
        stats['total_systems'] = len(system_analysis['system_counts'])
        # Get top system
        top_system = max(system_analysis['system_counts'].items(), key=lambda x: x[1]) if system_analysis['system_counts'] else ('None', 0)
        stats['top_system'] = top_system[0]
        stats['top_system_count'] = top_system[1]
    
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
                margin: 30px 0;
                overflow-x: auto;
            }
            table {
                width: 100%;
                border-collapse: collapse;
            }
            th {
                background-color: #3498db;
                color: white;
                padding: 12px;
                text-align: left;
            }
            td {
                padding: 10px;
                border-bottom: 1px solid #e0e0e0;
            }
            tr:nth-child(even) {
                background-color: #f9f9f9;
            }
            .footer {
                margin-top: 50px;
                text-align: center;
                color: #7f8c8d;
                font-size: 0.9em;
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
                    {% if stats.top_system %}
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
        top_genes=top_genes
    )
    
    # Write to file
    with open(report_path, 'w') as f:
        f.write(html_content)
    
    log.info(f"HTML report generated successfully: {report_path}")
    return report_path
