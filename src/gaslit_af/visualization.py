"""
Visualization module for GASLIT-AF Variant Analysis.
Provides functions to generate various plots and visualizations from VCF data.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import networkx as nx
from collections import defaultdict
import logging

# Configure logging
log = logging.getLogger("gaslit-af")

def save_visualization(fig, output_dir, filename, formats=None):
    """
    Save visualization in multiple formats.
    
    Args:
        fig: Matplotlib or Plotly figure
        output_dir: Directory to save the visualization
        filename: Base filename without extension
        formats: List of formats to save (default: ['png', 'html', 'svg'])
    """
    if formats is None:
        formats = ['png', 'html', 'svg']
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Handle different figure types
    if 'plotly' in str(type(fig)):
        if 'html' in formats:
            fig.write_html(os.path.join(output_dir, f"{filename}.html"))
        if 'png' in formats:
            fig.write_image(os.path.join(output_dir, f"{filename}.png"))
        if 'svg' in formats:
            fig.write_image(os.path.join(output_dir, f"{filename}.svg"))
        if 'json' in formats:
            fig.write_json(os.path.join(output_dir, f"{filename}.json"))
    else:  # Matplotlib figure
        for fmt in formats:
            if fmt != 'html':  # HTML not supported for matplotlib
                plt.savefig(os.path.join(output_dir, f"{filename}.{fmt}"), 
                           bbox_inches='tight', dpi=300)
    
    log.info(f"Saved visualization: {filename} in formats: {formats}")

def plot_chromosome_distribution(variant_data, output_dir):
    """
    Plot chromosome-level variant distribution.
    
    Args:
        variant_data: DataFrame with chromosome column ('CHROM' or 'chrom')
        output_dir: Directory to save visualizations
    """
    # Handle both uppercase and lowercase column names
    chrom_col = 'CHROM' if 'CHROM' in variant_data.columns else 'chrom'
    
    # Count variants per chromosome
    chrom_counts = variant_data[chrom_col].value_counts().reset_index()
    chrom_counts.columns = ['Chromosome', 'Count']
    
    # Sort chromosomes naturally (1, 2, ..., 10, 11, ... X, Y, MT)
    def chrom_key(chrom):
        # Handle empty or None values
        if not chrom or pd.isna(chrom):
            return float('inf') + 100  # Place empty values at the end
            
        if isinstance(chrom, str) and chrom.startswith('chr'):
            chrom = chrom[3:]
        try:
            return int(chrom)
        except (ValueError, TypeError):
            # X, Y, MT, etc. come after numbers
            if chrom == 'MT' or chrom == 'M':
                return float('inf')
            elif len(str(chrom)) == 1:
                return float('inf') - ord(str(chrom)[0])
            else:
                # For other non-standard chromosomes
                return float('inf') - 1
    
    chrom_counts['SortKey'] = chrom_counts['Chromosome'].apply(chrom_key)
    chrom_counts = chrom_counts.sort_values('SortKey').drop('SortKey', axis=1)
    
    # Create plotly figure
    fig = px.bar(chrom_counts, x='Chromosome', y='Count',
                title='Variant Distribution by Chromosome',
                labels={'Count': 'Number of Variants'},
                color='Count',
                color_continuous_scale='viridis')
    
    fig.update_layout(
        xaxis_title='Chromosome',
        yaxis_title='Number of Variants',
        template='plotly_white'
    )
    
    save_visualization(fig, output_dir, 'chromosome_distribution')
    return fig

def plot_variant_type_distribution(variant_data, output_dir):
    """
    Plot distribution of variant types (SNP, insertion, deletion, etc.).
    
    Args:
        variant_data: DataFrame with variant type information
        output_dir: Directory to save visualizations
    """
    # Handle both uppercase and lowercase column names
    ref_col = 'REF' if 'REF' in variant_data.columns else 'ref'
    alt_col = 'ALT' if 'ALT' in variant_data.columns else 'alt'
    
    # Determine variant types
    def get_variant_type(row):
        ref_len = len(row[ref_col])
        alt_len = len(row[alt_col])
        
        if ref_len == 1 and alt_len == 1:
            return 'SNP'
        elif ref_len > alt_len:
            return 'Deletion'
        elif ref_len < alt_len:
            return 'Insertion'
        else:
            return 'Other'
    
    variant_data['VariantType'] = variant_data.apply(get_variant_type, axis=1)
    type_counts = variant_data['VariantType'].value_counts().reset_index()
    type_counts.columns = ['Variant Type', 'Count']
    
    # Create plotly figure
    fig = px.pie(type_counts, values='Count', names='Variant Type',
                title='Distribution of Variant Types',
                color_discrete_sequence=px.colors.qualitative.Bold)
    
    fig.update_traces(textposition='inside', textinfo='percent+label')
    fig.update_layout(template='plotly_white')
    
    save_visualization(fig, output_dir, 'variant_type_distribution')
    return fig

def plot_transition_transversion_ratio(variant_data, output_dir):
    """
    Plot transition/transversion ratios.
    
    Args:
        variant_data: DataFrame with reference and alternate allele columns
        output_dir: Directory to save visualizations
    """
    # Handle both uppercase and lowercase column names
    ref_col = 'REF' if 'REF' in variant_data.columns else 'ref'
    alt_col = 'ALT' if 'ALT' in variant_data.columns else 'alt'
    
    # Filter for SNPs only
    snps = variant_data[(variant_data[ref_col].str.len() == 1) & 
                        (variant_data[alt_col].str.len() == 1)]
    
    # Define transitions and transversions
    transitions = [('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')]
    
    def get_mutation_type(row):
        if (row[ref_col], row[alt_col]) in transitions:
            return 'Transition'
        else:
            return 'Transversion'
    
    snps['MutationType'] = snps.apply(get_mutation_type, axis=1)
    mutation_counts = snps['MutationType'].value_counts()
    
    # Calculate ratio
    ti_tv_ratio = mutation_counts.get('Transition', 0) / mutation_counts.get('Transversion', 1)
    
    # Create figure
    fig = make_subplots(rows=1, cols=2, 
                        specs=[[{"type": "pie"}, {"type": "indicator"}]],
                        subplot_titles=('Transition vs Transversion', 'Ti/Tv Ratio'))
    
    fig.add_trace(
        go.Pie(
            labels=mutation_counts.index,
            values=mutation_counts.values,
            hole=0.3,
            marker_colors=['#1f77b4', '#ff7f0e']
        ),
        row=1, col=1
    )
    
    fig.add_trace(
        go.Indicator(
            mode="number+gauge",
            value=ti_tv_ratio,
            title={'text': "Ti/Tv Ratio"},
            gauge={
                'axis': {'range': [0, 3]},
                'bar': {'color': "darkblue"},
                'steps': [
                    {'range': [0, 1], 'color': "lightgray"},
                    {'range': [1, 2], 'color': "gray"},
                    {'range': [2, 3], 'color': "darkgray"}
                ],
                'threshold': {
                    'line': {'color': "red", 'width': 4},
                    'thickness': 0.75,
                    'value': 2.0  # Expected Ti/Tv ratio for whole genome
                }
            }
        ),
        row=1, col=2
    )
    
    fig.update_layout(
        title_text="Transition/Transversion Analysis",
        template="plotly_white"
    )
    
    save_visualization(fig, output_dir, 'transition_transversion_ratio')
    return fig

def plot_gene_variant_counts(gene_counts, output_dir, top_n=20):
    """
    Plot top variant-enriched genes.
    
    Args:
        gene_counts: Dictionary of gene:count pairs
        output_dir: Directory to save visualizations
        top_n: Number of top genes to show
    """
    # Convert to DataFrame and get top N genes
    df = pd.DataFrame(list(gene_counts.items()), columns=['Gene', 'Variants'])
    df = df.sort_values('Variants', ascending=False).head(top_n)
    
    # Create horizontal bar chart
    fig = px.bar(df, y='Gene', x='Variants', 
                orientation='h',
                title=f'Top {top_n} Variant-Enriched Genes',
                color='Variants',
                color_continuous_scale='Viridis')
    
    fig.update_layout(
        yaxis={'categoryorder': 'total ascending'},
        xaxis_title='Number of Variants',
        yaxis_title='Gene',
        template='plotly_white'
    )
    
    save_visualization(fig, output_dir, 'top_variant_enriched_genes')
    return fig

def create_gene_network(gene_counts, output_dir, min_variants=5):
    """
    Create a network visualization of gene variants.
    
    Args:
        gene_counts: Dictionary of gene:count pairs
        output_dir: Directory to save visualizations
        min_variants: Minimum number of variants for a gene to be included
    """
    # Filter genes with sufficient variants
    significant_genes = {gene: count for gene, count in gene_counts.items() 
                        if count >= min_variants}
    
    if len(significant_genes) < 2:
        log.warning(f"Not enough genes with {min_variants}+ variants for network visualization")
        return None
    
    # Create network graph
    G = nx.Graph()
    
    # Add nodes with variant counts as size
    for gene, count in significant_genes.items():
        G.add_node(gene, size=count, count=count)
    
    # Add edges based on biological relationships (simplified example)
    # In a real implementation, this would use actual biological pathway data
    # For this example, we'll connect genes that share similar variant counts
    genes = list(significant_genes.keys())
    counts = np.array(list(significant_genes.values()))
    
    # Normalize counts
    normalized_counts = (counts - counts.min()) / (counts.max() - counts.min() + 1e-10)
    
    # Connect genes with similar normalized counts (simple heuristic)
    for i in range(len(genes)):
        for j in range(i+1, len(genes)):
            similarity = 1 - abs(normalized_counts[i] - normalized_counts[j])
            if similarity > 0.8:  # Arbitrary threshold
                G.add_edge(genes[i], genes[j], weight=similarity)
    
    # Create positions using a layout algorithm
    pos = nx.spring_layout(G, seed=42)
    
    # Create plotly figure
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines')
    
    node_x = []
    node_y = []
    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)
    
    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers',
        hoverinfo='text',
        marker=dict(
            showscale=True,
            colorscale='YlGnBu',
            reversescale=True,
            color=[],
            size=[],
            colorbar=dict(
                thickness=15,
                title='Variant Count',
                xanchor='left',
                titleside='right'
            ),
            line_width=2))
    
    # Set node size and color based on variant count
    node_sizes = []
    node_colors = []
    node_text = []
    for node in G.nodes():
        count = G.nodes[node]['count']
        node_sizes.append(np.sqrt(count) * 5)  # Scale for better visualization
        node_colors.append(count)
        node_text.append(f'{node}: {count} variants')
    
    node_trace.marker.size = node_sizes
    node_trace.marker.color = node_colors
    node_trace.text = node_text
    
    # Create figure
    fig = go.Figure(data=[edge_trace, node_trace],
                layout=go.Layout(
                    title='Gene Variant Network',
                    titlefont_size=16,
                    showlegend=False,
                    hovermode='closest',
                    margin=dict(b=20,l=5,r=5,t=40),
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False))
                )
    
    save_visualization(fig, output_dir, 'gene_variant_network')
    return fig

def generate_all_visualizations(variant_data, gene_counts, output_dir):
    """
    Generate all visualizations for the analysis.
    
    Args:
        variant_data: DataFrame with variant information
        gene_counts: Dictionary of gene:count pairs
        output_dir: Directory to save visualizations
    
    Returns:
        Dictionary of figure objects
    """
    os.makedirs(output_dir, exist_ok=True)
    log.info(f"Generating visualizations in {output_dir}")
    
    figures = {}
    
    # Generate each visualization
    figures['chromosome_distribution'] = plot_chromosome_distribution(variant_data, output_dir)
    figures['variant_type_distribution'] = plot_variant_type_distribution(variant_data, output_dir)
    figures['transition_transversion'] = plot_transition_transversion_ratio(variant_data, output_dir)
    figures['gene_counts'] = plot_gene_variant_counts(gene_counts, output_dir)
    figures['gene_network'] = create_gene_network(gene_counts, output_dir)
    
    log.info(f"Generated {len(figures)} visualizations")
    return figures
