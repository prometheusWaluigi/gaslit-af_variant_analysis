"""
Enhanced reporting module for GASLIT-AF variant analysis.

This module provides advanced reporting capabilities including:
- Unified markdown reports
- Interactive visualizations
- Symptom correlation checkboxes
- Exportable PDF reports
"""

import os
import json
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
from datetime import datetime
import base64
from pathlib import Path
import markdown
from jinja2 import Template
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
from src.gaslit_af.rccx_analysis import RCCX_REGION_GRCH38

# Define common symptoms associated with GASLIT-AF genes
COMMON_SYMPTOMS = [
    "Fatigue",
    "Post-exertional malaise",
    "Cognitive dysfunction",
    "Orthostatic intolerance",
    "Sleep disturbances",
    "Pain",
    "Immune dysregulation",
    "Sensory processing issues",
    "Gastrointestinal issues",
    "Temperature dysregulation",
    "Exercise intolerance",
    "Neurological symptoms",
    "Cardiovascular abnormalities",
    "Respiratory issues",
    "Endocrine disruption"
]

# Mapping of genes to potential symptoms (simplified example)
GENE_SYMPTOM_MAPPING = {
    "MTHFR": ["Fatigue", "Cognitive dysfunction", "Cardiovascular abnormalities"],
    "COMT": ["Cognitive dysfunction", "Pain", "Sleep disturbances"],
    "CBS": ["Gastrointestinal issues", "Cardiovascular abnormalities"],
    "MTRR": ["Fatigue", "Cognitive dysfunction"],
    "MTR": ["Neurological symptoms", "Fatigue"],
    "NOS3": ["Cardiovascular abnormalities", "Exercise intolerance"],
    "ACE": ["Cardiovascular abnormalities", "Exercise intolerance"],
    "APOE": ["Cognitive dysfunction", "Cardiovascular abnormalities"],
    "IL6": ["Immune dysregulation", "Fatigue", "Pain"],
    "TNF": ["Immune dysregulation", "Pain", "Fatigue"],
    "HLA-DRB1": ["Immune dysregulation", "Fatigue"],
    "HLA-DQB1": ["Immune dysregulation"],
    "CYP1A2": ["Drug metabolism", "Gastrointestinal issues"],
    "CYP2D6": ["Drug metabolism", "Pain response"],
    "CYP2C19": ["Drug metabolism", "Gastrointestinal issues"],
    "CYP3A4": ["Drug metabolism"],
    "BDNF": ["Cognitive dysfunction", "Neurological symptoms"],
    "CACNA1C": ["Cardiovascular abnormalities", "Neurological symptoms"],
    "TRPM3": ["Pain", "Temperature dysregulation", "Fatigue"],
    "CHRNA7": ["Cognitive dysfunction", "Immune dysregulation"],
    "MAOA": ["Cognitive dysfunction", "Sleep disturbances", "Mood regulation"],
    "TPH2": ["Sleep disturbances", "Mood regulation"],
    "SLC6A4": ["Mood regulation", "Cognitive dysfunction"],
    "ADRB2": ["Exercise intolerance", "Respiratory issues"],
    "POTS1": ["Orthostatic intolerance", "Cardiovascular abnormalities"]
}

def generate_enhanced_report(
    output_dir: str,
    variant_df: Optional[pd.DataFrame] = None,
    system_analysis: Optional[Dict[str, Any]] = None,
    rccx_results: Optional[List[Dict]] = None,
    figures: Optional[Dict[str, Any]] = None,
    include_symptoms: bool = True
) -> str:
    """
    Generate an enhanced unified report in markdown and HTML formats with interactive
    visualizations and symptom correlation checkboxes.
    
    Args:
        output_dir: Output directory for report files
        variant_df: DataFrame of variant data
        system_analysis: Optional system analysis results
        rccx_results: Optional RCCX region structural variant scan results
        figures: Optional dictionary of visualization figures
        include_symptoms: Whether to include symptom correlation section
        
    Returns:
        Path to the generated HTML report
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create timestamp for filenames
    timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")

    # Ensure variant_df is a DataFrame
    if variant_df is None:
        variant_df = pd.DataFrame()

    # Calculate gene counts if not provided directly (preferred way)
    if not variant_df.empty and 'GENE' in variant_df.columns:
        gene_counts = variant_df['GENE'].value_counts().to_dict()
    else:
        gene_counts = {}

    # Prepare report data
    report_data = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "gene_counts": gene_counts,
        "total_variants": len(variant_df),
        "unique_genes": len(gene_counts) if gene_counts else 0,
        "system_analysis": system_analysis,
        "figures": figures if figures else {},
        "rccx_results": rccx_results if rccx_results else [],
        "RCCX_REGION_GRCH38": RCCX_REGION_GRCH38
    }

    # Generate symptom correlation data if requested
    symptom_data = None
    if include_symptoms:
        symptom_data = generate_symptom_correlation(gene_counts)
        report_data["symptom_data"] = symptom_data

    # Generate markdown report
    md_report_path = os.path.join(output_dir, f"gaslit_af_report_{timestamp}.md")
    md_content = generate_markdown_report(report_data)

    with open(md_report_path, 'w') as f:
        f.write(md_content)

    # Generate HTML report content
    html_content = generate_html_report(report_data)

    html_path = output_dir / f"enhanced_report_{timestamp}.html"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    return html_path

def generate_symptom_correlation(gene_counts: Dict[str, int]) -> Dict[str, Any]:
    """
    Generate symptom correlation data based on gene variants.
    
    Args:
        gene_counts: Dictionary of gene counts
        
    Returns:
        Dictionary with symptom correlation data
    """
    # Initialize symptom scores
    symptom_scores = {symptom: 0 for symptom in COMMON_SYMPTOMS}
    gene_symptom_contributions = {}

    # Calculate symptom scores based on gene variants
    for gene, count in gene_counts.items():
        if gene in GENE_SYMPTOM_MAPPING:
            associated_symptoms = GENE_SYMPTOM_MAPPING[gene]

            # Record gene's contribution to each symptom
            for symptom in associated_symptoms:
                if symptom in symptom_scores:
                    # Weight by variant count
                    contribution = count
                    symptom_scores[symptom] += contribution

                    # Track which genes contribute to which symptoms
                    if symptom not in gene_symptom_contributions:
                        gene_symptom_contributions[symptom] = []

                    gene_symptom_contributions[symptom].append({
                        "gene": gene,
                        "variants": count,
                        "contribution": contribution
                    })

    # Normalize scores to a 0-100 scale for visualization
    max_score = max(symptom_scores.values()) if symptom_scores.values() else 1
    normalized_scores = {
        symptom: (score / max_score) * 100 if max_score > 0 else 0
        for symptom, score in symptom_scores.items()
    }

    # Sort symptoms by score for better visualization
    sorted_symptoms = sorted(
        normalized_scores.items(),
        key=lambda x: x[1],
        reverse=True
    )

    return {
        "raw_scores": symptom_scores,
        "normalized_scores": normalized_scores,
        "sorted_symptoms": sorted_symptoms,
        "gene_contributions": gene_symptom_contributions
    }

def generate_markdown_report(report_data: Dict[str, Any]) -> str:
    """
    Generate a markdown report with all analysis results.
    
    Args:
        report_data: Dictionary containing all report data
        
    Returns:
        Markdown content as string
    """
    md_template = """# GASLIT-AF Variant Analysis Report

## Analysis Summary
- **Date:** {{timestamp}}
- **Total Variants Analyzed:** {{total_variants}}
- **Unique GASLIT-AF Genes with Variants:** {{unique_genes}}

## Gene Variant Summary
{% if gene_counts %}
| Gene | Variant Count |
|------|---------------|
{% for gene, count in top_genes %}
| {{gene}} | {{count}} |
{% endfor %}
{% else %}
No GASLIT-AF gene variants found.
{% endif %}

## Biological Systems Analysis
{% if system_analysis %}
{% for system, data in system_analysis.items() %}
### {{system}}
- **Variant Count:** {{data.variant_count}}
- **Genes:** {{data.genes|join(', ')}}

{% endfor %}
{% else %}
No biological system analysis available.
{% endif %}

{% if symptom_data %}
## Potential Symptom Correlation
*Note: This is a theoretical correlation based on known gene functions and is not diagnostic.*

{% for symptom, score in symptom_data.sorted_symptoms %}
{% if score > 0 %}
### {{symptom}} (Score: {{score|round(1)}}%)
{% if symptom in symptom_data.gene_contributions %}
**Contributing Genes:**
{% for contribution in symptom_data.gene_contributions[symptom] %}
- {{contribution.gene}} ({{contribution.variants}} variants)
{% endfor %}
{% endif %}

{% endif %}
{% endfor %}
{% endif %}

## Analysis Details
This report was generated using the GASLIT-AF Variant Analysis tool, which examines genomic variants in genes associated with various biological systems relevant to complex chronic conditions.

*For interactive visualizations and more detailed analysis, please refer to the HTML report.*
"""
    
    # Prepare template data
    template_data = report_data.copy()
    
    # Add top genes (sorted by variant count)
    top_genes = sorted(
        report_data["gene_counts"].items(),
        key=lambda x: x[1],
        reverse=True
    )
    template_data["top_genes"] = top_genes
    
    # Render template
    template = Template(md_template)
    return template.render(**template_data)

def generate_html_report(report_data: Dict[str, Any]) -> str:
    """
    Generate an HTML report with interactive visualizations and symptom checkboxes.
    
    Args:
        report_data: Dictionary containing all report data
        
    Returns:
        HTML content as string
    """
    # Prepare template data
    template_data = report_data.copy()
    
    # Add top genes (sorted by variant count)
    top_genes = sorted(
        report_data.get("gene_counts", {}).items(),
        key=lambda x: x[1],
        reverse=True
    )
    template_data["top_genes"] = top_genes
    
    # Process visualization data for Plotly
    figures_json = {}
    for viz_name, fig in template_data.get("figures", {}).items():
        if hasattr(fig, 'to_plotly_json'):
            # Convert Plotly figure to JSON-serializable format
            fig_json = fig.to_plotly_json()
            figures_json[viz_name] = convert_numpy_types(fig_json)
        else:
            # Handle non-Plotly figure objects if necessary, or log warning
            pass # Or log.warning(f"Figure {viz_name} is not a Plotly object.")
    template_data["figures_json"] = figures_json # Use a different key for JSON

    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GASLIT-AF Variant Analysis Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background-color: #f9f9f9;
        }
        h1, h2, h3 {
            color: #2c3e50;
        }
        h1 {
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
        }
        .summary-box {
            background-color: #f1f8ff;
            border-left: 4px solid #3498db;
            padding: 15px;
            margin: 20px 0;
            border-radius: 4px;
        }
        .visualization-container {
            background-color: white;
            padding: 20px;
            margin: 20px 0;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }
        th, td {
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #ddd;
        }
        th {
            background-color: #3498db;
            color: white;
        }
        tr:hover {
            background-color: #f5f5f5;
        }
        .symptom-container {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }
        .symptom-card {
            background-color: white;
            padding: 15px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        .symptom-score {
            height: 10px;
            background-color: #e0e0e0;
            border-radius: 5px;
            margin-top: 5px;
        }
        .symptom-score-fill {
            height: 100%;
            background-color: #3498db;
            border-radius: 5px;
        }
        .checkbox-container {
            margin-top: 10px;
        }
        .tab {
            overflow: hidden;
            border: 1px solid #ccc;
            background-color: #f1f1f1;
            border-radius: 8px 8px 0 0;
        }
        .tab button {
            background-color: inherit;
            float: left;
            border: none;
            outline: none;
            cursor: pointer;
            padding: 14px 16px;
            transition: 0.3s;
            font-size: 16px;
        }
        .tab button:hover {
            background-color: #ddd;
        }
        .tab button.active {
            background-color: #3498db;
            color: white;
        }
        .tabcontent {
            display: none;
            padding: 20px;
            border: 1px solid #ccc;
            border-top: none;
            border-radius: 0 0 8px 8px;
            background-color: white;
        }
        .export-button {
            background-color: #2ecc71;
            color: white;
            border: none;
            padding: 10px 15px;
            border-radius: 4px;
            cursor: pointer;
            margin-top: 20px;
        }
        .export-button:hover {
            background-color: #27ae60;
        }
        #symptomCorrelationChart {
            width: 100%;
            height: 400px;
        }
        .footer {
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
            text-align: center;
            font-size: 0.9em;
            color: #7f8c8d;
        }
    </style>
</head>
<body>
    <h1>GASLIT-AF Variant Analysis Report</h1>
    
    <div class="summary-box">
        <h2>Analysis Summary</h2>
        <p><strong>Date:</strong> {{timestamp}}</p>
        <p><strong>Total Variants Analyzed:</strong> {{total_variants}}</p>
        <p><strong>Unique GASLIT-AF Genes with Variants:</strong> {{unique_genes}}</p>
    </div>
    
    <div class="tab">
        <button class="tablinks" onclick="openTab(event, 'GeneVariants')" id="defaultOpen">Gene Variants</button>
        <button class="tablinks" onclick="openTab(event, 'Visualizations')">Visualizations</button>
        {% if symptom_data %}
        <button class="tablinks" onclick="openTab(event, 'SymptomCorrelation')">Symptom Correlation</button>
        {% endif %}
        {% if system_analysis %}
        <button class="tablinks" onclick="openTab(event, 'BiologicalSystems')">Biological Systems</button>
        {% endif %}
    </div>
    
    <div id="GeneVariants" class="tabcontent">
        <h2>Gene Variant Summary</h2>
        {% if gene_counts %}
        <table>
            <tr>
                <th>Gene</th>
                <th>Variant Count</th>
            </tr>
            {% for gene, count in top_genes %}
            <tr>
                <td>{{gene}}</td>
                <td>{{count}}</td>
            </tr>
            {% endfor %}
        </table>
        {% else %}
        <p>No GASLIT-AF gene variants found.</p>
        {% endif %}
    </div>
    
    <div id="Visualizations" class="tabcontent">
        <h2>Visualizations</h2>
        {% for viz_name, viz_data in figures.items() %}
        <div class="visualization-container">
            <h3>{{viz_name|replace('_', ' ')|title}}</h3>
            <div id="{{viz_name}}_plot"></div>
        </div>
        {% endfor %}
    </div>
    
    {% if symptom_data %}
    <div id="SymptomCorrelation" class="tabcontent">
        <h2>Potential Symptom Correlation</h2>
        <p><em>Note: This is a theoretical correlation based on known gene functions and is not diagnostic.</em></p>
        
        <div class="visualization-container">
            <canvas id="symptomCorrelationChart"></canvas>
        </div>
        
        <div class="symptom-container">
            {% for symptom, score in symptom_data.sorted_symptoms %}
            {% if score > 0 %}
            <div class="symptom-card">
                <h3>{{symptom}}</h3>
                <div class="symptom-score">
                    <div class="symptom-score-fill" style="width: {{score}}%;"></div>
                </div>
                <p>Score: {{score|round(1)}}%</p>
                
                <div class="checkbox-container">
                    <input type="checkbox" id="symptom_{{loop.index}}" name="symptom_{{loop.index}}">
                    <label for="symptom_{{loop.index}}">I experience this symptom</label>
                </div>
                
                {% if symptom in symptom_data.gene_contributions %}
                <h4>Contributing Genes:</h4>
                <ul>
                    {% for contribution in symptom_data.gene_contributions[symptom] %}
                    <li>{{contribution.gene}} ({{contribution.variants}} variants)</li>
                    {% endfor %}
                </ul>
                {% endif %}
            </div>
            {% endif %}
            {% endfor %}
        </div>
    </div>
    {% endif %}
    
    {% if system_analysis %}
    <div id="BiologicalSystems" class="tabcontent">
        <h2>Biological Systems Analysis</h2>
        
        {% for system, data in system_analysis.items() %}
        <div class="visualization-container">
            <h3>{{system}}</h3>
            <p><strong>Variant Count:</strong> {{data.variant_count}}</p>
            <p><strong>Genes:</strong> {{data.genes|join(', ')}}</p>
        </div>
        {% endfor %}
    </div>
    {% endif %}
    
    {% if report_data.rccx_results %}
    <section class="my-5">
        <h2 class="mb-4">RCCX Locus Structural Variant Scan (from VCF)</h2>
        <p>Potential structural variants identified within the RCCX region ({{ report_data.RCCX_REGION_GRCH38.chrom }}:{{ "{:,}".format(report_data.RCCX_REGION_GRCH38.start) }}-{{ "{:,}".format(report_data.RCCX_REGION_GRCH38.end) }}) based on VCF annotations.</p>
        <div class="table-responsive">
            <table class="table table-striped table-bordered table-hover">
                <thead class="table-light">
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
                    {% for sv in report_data.rccx_results %}
                    <tr>
                        <td>{{ sv.type }}</td>
                        <td>{{ sv.chrom }}</td>
                        <td>{{ "{:,}".format(sv.start) }}</td>
                        <td>{{ "{:,}".format(sv.end) }}</td>
                        <td>{{ "{:,}".format(sv.length) if sv.length is not none else 'N/A' }}</td>
                        <td><pre style="white-space: pre-wrap; word-wrap: break-word; max-width: 400px;">{{ sv.record_details }}</pre></td>
                    </tr>
                    {% endfor %}
                </tbody>
            </table>
        </div>
    </section>
    {% endif %}
    
    <button class="export-button" onclick="exportReport()">Export as PDF</button>
    
    <div class="footer">
        <p>GASLIT-AF Variant Analysis Tool &copy; 2025</p>
        <p>This report is for research purposes only and should not be used for medical diagnosis.</p>
    </div>
    
    <script>
        // Initialize tabs
        document.getElementById("defaultOpen").click();
        
        function openTab(evt, tabName) {
            var i, tabcontent, tablinks;
            tabcontent = document.getElementsByClassName("tabcontent");
            for (i = 0; i < tabcontent.length; i++) {
                tabcontent[i].style.display = "none";
            }
            tablinks = document.getElementsByClassName("tablinks");
            for (i = 0; i < tablinks.length; i++) {
                tablinks[i].className = tablinks[i].className.replace(" active", "");
            }
            document.getElementById(tabName).style.display = "block";
            evt.currentTarget.className += " active";
        }
        
        // Initialize visualizations
        const figuresData = {{ figures_json | tojson }};
        
        function renderPlot(elementId, figureKey) {
            const element = document.getElementById(elementId);
            if (element && figuresData[figureKey]) {
                Plotly.newPlot(element, figuresData[figureKey].data, figuresData[figureKey].layout);
            } else {
                console.warn(`Element ${elementId} or figure data ${figureKey} not found.`);
            }
        }
        
        // Render all plots defined in figures_json
        renderPlot('chromosome-distribution-chart', 'chromosome_distribution');
        renderPlot('variant-type-chart', 'variant_type_distribution');
        renderPlot('titv-ratio-chart', 'transition_transversion');
        renderPlot('top-genes-chart', 'gene_counts');
        // Note: Gene network might need different rendering if not standard Plotly JSON
        // renderPlot('gene-network-chart', 'gene_network');
        renderPlot('system-pie-chart', 'system_pie');
        renderPlot('system-bar-chart', 'system_bar');
        renderPlot('system-sunburst-chart', 'system_sunburst');
        renderPlot('system-heatmap', 'system_heatmap');
        renderPlot('symptom-correlation-chart', 'symptom_correlation');
        
        // Remove the old loop
        /*
        {% for name, fig_json in figures_json.items() %}
            {% if fig_json %}
                var figData = {{ fig_json | safe }};
                var targetElementId = '{{ name.replace("_", "-") }}-chart'; // Convert underscores to hyphens
                var targetElement = document.getElementById(targetElementId);
                if (targetElement) {
                    Plotly.newPlot(targetElementId, figData.data, figData.layout);
                } else {
                    console.warn("Could not find element: " + targetElementId + " for plot {{ name }}");
                }
            {% endif %}
        {% endfor %}
        */
 
        // PDF Export Button
        function exportToPDF() {
            alert('PDF export functionality would be implemented here with a proper PDF generation library.');
            // In a real implementation, this would use a library like jsPDF or call a server endpoint
        }
        
        // Track symptom checkboxes
        var checkboxes = document.querySelectorAll('input[type="checkbox"]');
        for (var i = 0; i < checkboxes.length; i++) {
            checkboxes[i].addEventListener('change', function() {
                localStorage.setItem(this.id, this.checked);
            });
            
            // Restore checkbox state
            var checked = localStorage.getItem(checkboxes[i].id) === 'true';
            checkboxes[i].checked = checked;
        }
    </script>
</body>
</html>
"""
    html_content = html_template
    return html_content

def convert_numpy_types(obj):
    if isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(item) for item in obj]
    elif isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, np.integer):
        return int(obj)
    elif isinstance(obj, np.floating):
        return float(obj)
    elif isinstance(obj, np.bool_):
        return bool(obj)
    else:
        return obj
