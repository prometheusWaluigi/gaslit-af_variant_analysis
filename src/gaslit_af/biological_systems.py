"""
Biological systems module for GASLIT-AF Variant Analysis.
Provides functions to categorize genes by biological pathways and analyze variants at the system level.
"""

import pandas as pd
import numpy as np
import logging
from collections import defaultdict
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Configure logging
log = logging.getLogger("gaslit-af")

# Define biological systems and their associated genes as specified in the PRD
BIOLOGICAL_SYSTEMS = {
    "Immune & Inflammatory Pathways": [
        "IDO2", "AHR", "AHRR", "IL36RN", "CFH", "MBL2", "NLRP3", "IL1B", "IL6", 
        "IL17", "IL13", "IL4", "HLA-DQB1", "PTPN22", "CTLA4", "ASXL1", "CBL", 
        "DNMT3B", "ETV6", "IDH1"
    ],
    "Autonomic & Neurotransmitter Pathways": [
        "COMT", "CHRM2", "DRD2", "GABRA1", "CHRNA7", "ADRB1", "ADRB2", "NOS3", 
        "GNB3", "SLC6A2", "NET", "EZH2", "SLC6A4", "HTR2A", "TAAR1", "OPRM1", 
        "GCH1", "TRPV2", "MYT1L", "NRXN3"
    ],
    "Structural & Connective Tissue Integrity": [
        "TNXB", "ADAMTS10", "SELENON", "NEB", "MYH7", "MAPRE1", "ADGRV1", 
        "PLXNA2", "COL3A1", "FBN1", "FLNA", "COL5A1", "FKBP14", "PLOD1"
    ],
    "Metabolic, Mitochondrial & Oxidative Stress": [
        "APOE", "PCSK9", "UGT1A1", "HNF1A", "ABCC8", "TFAM", "C19orf12", 
        "MT-ATP6", "MT-ATP8", "PDHA1", "SDHB", "NAMPT", "NMRK1", "PGC1A"
    ],
    "Endocannabinoid System (ECS)": [
        "CNR1", "CNR2", "FAAH", "MGLL"
    ],
    "Calcium & Ion Channels": [
        "ITPR1", "KCNJ5", "RYR2"
    ],
    "Mast Cell Activation & Histamine Metabolism": [
        "TPSAB1", "KIT", "HNMT", "TET2"
    ],
    "Kynurenine Pathway": [
        "IDO1", "KMO", "KYNU", "TDO2", "HAAO", "ARNT", "BECN1", "ATG5"
    ]
}

# Create a reverse mapping from gene to system
GENE_TO_SYSTEM = {}
for system, genes in BIOLOGICAL_SYSTEMS.items():
    for gene in genes:
        GENE_TO_SYSTEM[gene] = system

def get_system_for_gene(gene):
    """
    Get the biological system for a gene.
    
    Args:
        gene: Gene symbol
    
    Returns:
        Biological system name or "Other" if not found
    """
    return GENE_TO_SYSTEM.get(gene, "Other")

def analyze_systems(gene_counts):
    """
    Analyze variant distribution across biological systems.
    
    Args:
        gene_counts: Dictionary of gene:count pairs
    
    Returns:
        Dictionary with system-level analysis results
    """
    log.info("Analyzing variant distribution across biological systems")
    
    # Initialize system counts
    system_counts = defaultdict(int)
    system_genes = defaultdict(list)
    
    # Count variants by system
    for gene, count in gene_counts.items():
        system = get_system_for_gene(gene)
        system_counts[system] += count
        system_genes[system].append((gene, count))
    
    # Sort genes within each system by variant count
    for system in system_genes:
        system_genes[system] = sorted(system_genes[system], key=lambda x: -x[1])
    
    # Calculate percentages
    total_variants = sum(system_counts.values())
    system_percentages = {system: (count / total_variants * 100) if total_variants > 0 else 0 
                         for system, count in system_counts.items()}
    
    # Prepare results
    results = {
        "system_counts": dict(system_counts),
        "system_percentages": system_percentages,
        "system_genes": system_genes,
        "total_variants": total_variants
    }
    
    log.info(f"Found variants across {len(system_counts)} biological systems")
    
    return results

def plot_system_distribution(system_analysis, output_dir):
    """
    Create visualizations for biological system variant distribution.
    
    Args:
        system_analysis: Results from analyze_systems function
        output_dir: Directory to save visualizations
    
    Returns:
        Dictionary of plotly figures
    """
    figures = {}
    
    # Extract data
    system_counts = system_analysis["system_counts"]
    system_percentages = system_analysis["system_percentages"]
    
    # Prepare DataFrame
    df = pd.DataFrame({
        "Biological System": list(system_counts.keys()),
        "Variant Count": list(system_counts.values()),
        "Percentage": [system_percentages[sys] for sys in system_counts.keys()]
    })
    
    # Sort by variant count
    df = df.sort_values("Variant Count", ascending=False)
    
    # Create bar chart
    fig_bar = px.bar(
        df, 
        x="Biological System", 
        y="Variant Count",
        color="Variant Count",
        color_continuous_scale="Viridis",
        title="Variant Distribution by Biological System",
        labels={"Variant Count": "Number of Variants"}
    )
    
    fig_bar.update_layout(
        xaxis_title="Biological System",
        yaxis_title="Number of Variants",
        xaxis={'categoryorder': 'total descending'},
        template="plotly_white"
    )
    
    figures["system_bar"] = fig_bar
    
    # Create pie chart
    fig_pie = px.pie(
        df,
        values="Variant Count",
        names="Biological System",
        title="Proportion of Variants by Biological System",
        color_discrete_sequence=px.colors.qualitative.Bold
    )
    
    fig_pie.update_traces(textposition='inside', textinfo='percent+label')
    fig_pie.update_layout(template="plotly_white")
    
    figures["system_pie"] = fig_pie
    
    # Create heatmap of top genes per system
    system_genes = system_analysis["system_genes"]
    
    # Prepare data for heatmap
    heatmap_data = []
    for system, genes in system_genes.items():
        # Take top 5 genes per system
        for i, (gene, count) in enumerate(genes[:5]):
            heatmap_data.append({
                "Biological System": system,
                "Gene": gene,
                "Variant Count": count,
                "Rank": i + 1
            })
    
    if heatmap_data:
        heatmap_df = pd.DataFrame(heatmap_data)
        
        # Create heatmap
        fig_heatmap = px.density_heatmap(
            heatmap_df,
            x="Biological System",
            y="Gene",
            z="Variant Count",
            title="Top Genes by Biological System",
            color_continuous_scale="Viridis"
        )
        
        fig_heatmap.update_layout(
            xaxis_title="Biological System",
            yaxis_title="Gene",
            template="plotly_white"
        )
        
        figures["system_heatmap"] = fig_heatmap
    
    # Create sunburst chart
    sunburst_data = []
    for system, genes in system_genes.items():
        for gene, count in genes:
            sunburst_data.append({
                "System": system,
                "Gene": gene,
                "Count": count
            })
    
    if sunburst_data:
        sunburst_df = pd.DataFrame(sunburst_data)
        
        fig_sunburst = px.sunburst(
            sunburst_df,
            path=['System', 'Gene'],
            values='Count',
            title="Hierarchical View of Variants by System and Gene"
        )
        
        fig_sunburst.update_layout(template="plotly_white")
        
        figures["system_sunburst"] = fig_sunburst
    
    return figures

def generate_system_summary(system_analysis):
    """
    Generate a text summary of the biological system analysis.
    
    Args:
        system_analysis: Results from analyze_systems function
    
    Returns:
        String with formatted summary
    """
    system_counts = system_analysis["system_counts"]
    system_genes = system_analysis["system_genes"]
    total_variants = system_analysis["total_variants"]
    
    # Sort systems by variant count
    sorted_systems = sorted(system_counts.items(), key=lambda x: -x[1])
    
    summary = ["# Biological System Analysis Summary\n"]
    summary.append(f"Total variants analyzed across all systems: {total_variants:,}\n")
    
    summary.append("## Variant Distribution by Biological System\n")
    for system, count in sorted_systems:
        percentage = (count / total_variants * 100) if total_variants > 0 else 0
        summary.append(f"### {system}: {count:,} variants ({percentage:.1f}%)\n")
        
        # List top genes for this system
        top_genes = system_genes[system][:5]  # Top 5 genes
        if top_genes:
            summary.append("Top genes in this system:\n")
            for gene, gene_count in top_genes:
                gene_percentage = (gene_count / count * 100) if count > 0 else 0
                summary.append(f"- {gene}: {gene_count:,} variants ({gene_percentage:.1f}% of system)\n")
        
        summary.append("\n")
    
    return "".join(summary)
