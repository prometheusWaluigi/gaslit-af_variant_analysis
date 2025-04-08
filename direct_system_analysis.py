#!/usr/bin/env python3
"""
Direct System Analysis for GASLIT-AF Variant Analysis

This script establishes a direct quantum coherence bridge between annotated
VCF files and the GASLIT-AF biological systems, bypassing the standard
analysis pipeline to create a more robust recursive mapping.

The script processes annotated VCF files, extracts gene information from
the ANN field, and maps variants to biological systems defined in the
GASLIT-AF theoretical framework.
"""

import os
import sys
import gzip
import json
import logging
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from collections import defaultdict
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn

# Import visualization libraries
try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af-direct-analysis")

# Import GASLIT-AF biological systems
try:
    from src.gaslit_af.biological_systems import BIOLOGICAL_SYSTEMS
except ImportError:
    # If running from the same directory, try relative import
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    try:
        from src.gaslit_af.biological_systems import BIOLOGICAL_SYSTEMS
    except ImportError:
        # Define biological systems manually if import fails
        log.warning("Could not import BIOLOGICAL_SYSTEMS, using hardcoded definition")
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

class DirectSystemAnalyzer:
    """
    Direct System Analyzer for GASLIT-AF Variant Analysis.
    
    This class establishes a direct quantum coherence bridge between
    annotated VCF files and the GASLIT-AF biological systems.
    """
    
    def __init__(self):
        """
        Initialize the Direct System Analyzer.
        """
        # Initialize system analysis results
        self.system_counts = defaultdict(int)
        self.system_genes = defaultdict(set)
        self.gene_counts = defaultdict(int)
        self.total_variants = 0
        
        # Log initialization
        log.info(f"Initialized Direct System Analyzer with {len(BIOLOGICAL_SYSTEMS)} biological systems")
        for system, genes in BIOLOGICAL_SYSTEMS.items():
            log.info(f"  - {system}: {len(genes)} genes")
    
    def analyze_vcf(self, vcf_path: str, output_dir: str, chunk_size: int = 10000):
        """
        Analyze an annotated VCF file and map variants to biological systems.
        
        Args:
            vcf_path: Path to annotated VCF file
            output_dir: Directory to save analysis results
            chunk_size: Number of variants to process at once
        """
        log.info(f"Analyzing VCF file: {vcf_path}")
        
        # Check if input file exists
        if not os.path.exists(vcf_path):
            log.error(f"Input VCF file not found: {vcf_path}")
            return False
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Process VCF file
        try:
            # Determine if input is gzipped
            is_gzipped = vcf_path.endswith('.gz')
            
            # Open input file
            if is_gzipped:
                vcf_file = gzip.open(vcf_path, 'rt')
            else:
                vcf_file = open(vcf_path, 'r')
            
            # Skip header lines
            for line in vcf_file:
                if not line.startswith('#'):
                    break
            
            # Reopen file to count total variants
            vcf_file.close()
            if is_gzipped:
                vcf_file = gzip.open(vcf_path, 'rt')
            else:
                vcf_file = open(vcf_path, 'r')
            
            # Count total variants
            total_variants = 0
            for line in vcf_file:
                if not line.startswith('#'):
                    total_variants += 1
            
            log.info(f"Processing {total_variants} variants")
            
            # Reopen file to process variants
            vcf_file.close()
            if is_gzipped:
                vcf_file = gzip.open(vcf_path, 'rt')
            else:
                vcf_file = open(vcf_path, 'r')
            
            # Skip header lines
            for line in vcf_file:
                if not line.startswith('#'):
                    break
            
            # Process variants in chunks
            processed_variants = 0
            annotated_variants = 0
            system_variants = 0
            
            with Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("{task.completed}/{task.total}"),
                TextColumn("({task.percentage:.0f}%)"),
                TimeElapsedColumn(),
                console=console
            ) as progress:
                task = progress.add_task("Analyzing variants", total=total_variants)
                
                variant_buffer = []
                variant_buffer.append(line)  # Add the first non-header line
                
                for line in vcf_file:
                    if line.startswith('#'):
                        continue
                    
                    variant_buffer.append(line)
                    
                    if len(variant_buffer) >= chunk_size:
                        # Process buffer
                        processed, annotated, system = self._process_variant_chunk(variant_buffer)
                        processed_variants += processed
                        annotated_variants += annotated
                        system_variants += system
                        
                        # Update progress
                        progress.update(task, advance=processed)
                        
                        # Clear buffer
                        variant_buffer = []
                
                # Process remaining variants
                if variant_buffer:
                    processed, annotated, system = self._process_variant_chunk(variant_buffer)
                    processed_variants += processed
                    annotated_variants += annotated
                    system_variants += system
                    
                    # Update progress
                    progress.update(task, advance=processed)
            
            # Close file
            vcf_file.close()
            
            # Save analysis results
            self._save_analysis_results(output_dir)
            
            log.info(f"Analysis complete:")
            log.info(f"  - Processed variants: {processed_variants}")
            log.info(f"  - Annotated variants: {annotated_variants}")
            log.info(f"  - System variants: {system_variants}")
            
            return True
        
        except Exception as e:
            log.error(f"Error analyzing VCF file: {e}")
            return False
    
    def _process_variant_chunk(self, variant_lines: List[str]) -> Tuple[int, int, int]:
        """
        Process a chunk of variant lines.
        
        Args:
            variant_lines: List of variant lines
            
        Returns:
            Tuple of (processed_count, annotated_count, system_count)
        """
        processed_count = 0
        annotated_count = 0
        system_count = 0
        
        for line in variant_lines:
            processed_count += 1
            
            # Parse variant line
            fields = line.strip().split('\t')
            if len(fields) < 8:
                # Invalid line, skip
                continue
            
            # Extract variant information
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # Check for ANN field in INFO column
            if 'ANN=' not in info:
                continue
            
            annotated_count += 1
            
            # Extract gene information from ANN field
            ann_start = info.find('ANN=') + 4
            ann_end = info.find(';', ann_start)
            if ann_end == -1:
                ann_end = len(info)
            
            ann_field = info[ann_start:ann_end]
            
            # Parse ANN field
            genes = set()
            for ann_entry in ann_field.split(','):
                # Simple annotation format: Gene|GeneID
                gene_parts = ann_entry.split('|')
                if len(gene_parts) > 0:
                    gene = gene_parts[0]
                    genes.add(gene)
            
            # Map genes to biological systems
            for gene in genes:
                # Update gene counts
                self.gene_counts[gene] += 1
                
                # Check if gene is in a biological system
                if gene in GENE_TO_SYSTEM:
                    system = GENE_TO_SYSTEM[gene]
                    
                    # Update system counts
                    self.system_counts[system] += 1
                    self.system_genes[system].add(gene)
                    
                    # Increment system variant count
                    system_count += 1
                    
                    # Increment total variant count
                    self.total_variants += 1
        
        return processed_count, annotated_count, system_count
    
    def _save_analysis_results(self, output_dir: str):
        """
        Save analysis results to output directory.
        
        Args:
            output_dir: Directory to save analysis results
        """
        # Calculate system percentages
        system_percentages = {}
        if self.total_variants > 0:
            for system, count in self.system_counts.items():
                system_percentages[system] = (count / self.total_variants) * 100
        
        # Convert system genes to list for JSON serialization
        system_genes_list = {}
        for system, genes in self.system_genes.items():
            system_genes_list[system] = list(genes)
        
        # Create analysis results
        analysis_results = {
            'total_variants': self.total_variants,
            'system_counts': dict(self.system_counts),
            'system_percentages': system_percentages,
            'system_genes': system_genes_list
        }
        
        # Save analysis results to JSON file
        system_analysis_path = os.path.join(output_dir, 'system_analysis.json')
        with open(system_analysis_path, 'w') as f:
            json.dump(analysis_results, f, indent=2)
        
        log.info(f"Saved system analysis results to {system_analysis_path}")
        
        # Save gene counts to CSV file
        gene_counts_path = os.path.join(output_dir, 'gene_counts.csv')
        gene_counts_df = pd.DataFrame([
            {'gene': gene, 'count': count, 'system': GENE_TO_SYSTEM.get(gene, 'Unknown')}
            for gene, count in self.gene_counts.items()
        ])
        gene_counts_df.sort_values('count', ascending=False, inplace=True)
        gene_counts_df.to_csv(gene_counts_path, index=False)
        
        log.info(f"Saved gene counts to {gene_counts_path}")
        
        # Generate system analysis visualization
        self._generate_system_visualization(output_dir)
    
    def _generate_system_visualization(self, output_dir: str):
        """
        Generate system analysis visualization.
        
        Args:
            output_dir: Directory to save visualization
        """
        if self.total_variants == 0:
            log.warning("No variants found in biological systems, skipping visualization")
            return
        
        if not PLOTLY_AVAILABLE:
            log.warning("Plotly not available, skipping visualization")
            return
        
        try:
            # Create system counts dataframe
            system_df = pd.DataFrame([
                {'system': system, 'count': count, 'percentage': (count / self.total_variants) * 100}
                for system, count in self.system_counts.items()
            ])
            system_df.sort_values('count', ascending=False, inplace=True)
            
            # Create bar chart
            fig = px.bar(
                system_df,
                x='system',
                y='count',
                color='system',
                title=f'GASLIT-AF Biological System Variant Distribution (Total: {self.total_variants:,})',
                labels={'system': 'Biological System', 'count': 'Variant Count'},
                height=600
            )
            
            # Update layout
            fig.update_layout(
                xaxis_title='Biological System',
                yaxis_title='Variant Count',
                xaxis={'categoryorder': 'total descending'},
                showlegend=False
            )
            
            # Save visualization
            system_viz_path = os.path.join(output_dir, 'system_distribution')
            fig.write_html(f"{system_viz_path}.html")
            fig.write_image(f"{system_viz_path}.png", width=1200, height=800)
            
            log.info(f"Generated system distribution visualization: {system_viz_path}.html")
            
            # Create pie chart
            fig = px.pie(
                system_df,
                values='count',
                names='system',
                title=f'GASLIT-AF Biological System Variant Distribution (Total: {self.total_variants:,})',
                height=600
            )
            
            # Update layout
            fig.update_layout(
                showlegend=True
            )
            
            # Save visualization
            system_pie_path = os.path.join(output_dir, 'system_distribution_pie')
            fig.write_html(f"{system_pie_path}.html")
            fig.write_image(f"{system_pie_path}.png", width=1200, height=800)
            
            log.info(f"Generated system distribution pie chart: {system_pie_path}.html")
            
            # Create top genes visualization
            top_genes_df = pd.DataFrame([
                {'gene': gene, 'count': count, 'system': GENE_TO_SYSTEM.get(gene, 'Unknown')}
                for gene, count in self.gene_counts.items()
                if gene in GENE_TO_SYSTEM  # Only include genes in biological systems
            ])
            top_genes_df.sort_values('count', ascending=False, inplace=True)
            top_genes_df = top_genes_df.head(20)  # Top 20 genes
            
            # Create bar chart
            fig = px.bar(
                top_genes_df,
                x='gene',
                y='count',
                color='system',
                title=f'Top 20 GASLIT-AF Genes by Variant Count',
                labels={'gene': 'Gene', 'count': 'Variant Count', 'system': 'Biological System'},
                height=600
            )
            
            # Update layout
            fig.update_layout(
                xaxis_title='Gene',
                yaxis_title='Variant Count',
                xaxis={'categoryorder': 'total descending'}
            )
            
            # Save visualization
            top_genes_path = os.path.join(output_dir, 'top_genes')
            fig.write_html(f"{top_genes_path}.html")
            fig.write_image(f"{top_genes_path}.png", width=1200, height=800)
            
            log.info(f"Generated top genes visualization: {top_genes_path}.html")
        
        except Exception as e:
            log.error(f"Error generating system visualization: {e}")

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Direct System Analysis for GASLIT-AF Variant Analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "vcf_path",
        help="Path to annotated VCF file"
    )
    
    parser.add_argument(
        "--output-dir",
        default="analysis_results/direct_analysis",
        help="Directory to save analysis results"
    )
    
    parser.add_argument(
        "--chunk-size",
        type=int,
        default=10000,
        help="Number of variants to process at once"
    )
    
    return parser.parse_args()

def main():
    """
    Main entry point for Direct System Analysis.
    """
    args = parse_args()
    
    # Create analyzer
    analyzer = DirectSystemAnalyzer()
    
    # Analyze VCF file
    success = analyzer.analyze_vcf(
        vcf_path=args.vcf_path,
        output_dir=args.output_dir,
        chunk_size=args.chunk_size
    )
    
    if success:
        log.info(f"Successfully analyzed VCF file: {args.vcf_path}")
        log.info(f"Results saved to: {args.output_dir}")
    else:
        log.error(f"Failed to analyze VCF file")
        sys.exit(1)

if __name__ == "__main__":
    main()
