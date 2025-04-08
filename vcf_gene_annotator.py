#!/usr/bin/env python3
"""
VCF Gene Annotator for GASLIT-AF Variant Analysis

This script adds gene annotations to VCF files that don't have them,
creating a quantum coherence bridge between raw genomic data and the
GASLIT-AF theoretical framework.

The script uses the pyensembl library to map genomic coordinates to gene
symbols, establishing a recursive connection between variant positions and
the biological systems in the GASLIT-AF model.
"""

import os
import sys
import gzip
import logging
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn

# Try to import pyensembl for gene annotation
try:
    import pyensembl
    PYENSEMBL_AVAILABLE = True
except ImportError:
    PYENSEMBL_AVAILABLE = False

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af-annotator")

# Import GASLIT-AF gene list
try:
    from src.gaslit_af.gene_lists import GASLIT_AF_GENES
except ImportError:
    # If running from the same directory, try relative import
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    try:
        from src.gaslit_af.gene_lists import GASLIT_AF_GENES
    except ImportError:
        # Define a minimal set of GASLIT-AF genes if import fails
        log.warning("Could not import GASLIT_AF_GENES, using minimal set")
        GASLIT_AF_GENES_TEXT = """
        IDO2 AHR AHRR IL36RN CFH MBL2 NLRP3 IL1B IL6 IL17 IL13 IL4 HLA-DQB1 PTPN22 CTLA4 ASXL1 CBL DNMT3B ETV6 IDH1
        COMT CHRM2 DRD2 GABRA1 CHRNA7 ADRB1 ADRB2 NOS3 GNB3 SLC6A2 NET EZH2 SLC6A4 HTR2A TAAR1 OPRM1 GCH1 TRPV2 MYT1L NRXN3
        TNXB ADAMTS10 SELENON NEB MYH7 MAPRE1 ADGRV1 PLXNA2 COL3A1 FBN1 FLNA COL5A1 FKBP14 PLOD1
        APOE PCSK9 UGT1A1 HNF1A ABCC8 TFAM C19orf12 MT-ATP6 MT-ATP8 PDHA1 SDHB NAMPT NMRK1 PGC1A
        CNR1 CNR2 FAAH MGLL
        ITPR1 KCNJ5 RYR2
        TPSAB1 KIT HNMT TET2
        IDO1 KMO KYNU TDO2 HAAO ARNT BECN1 ATG5
        """
        GASLIT_AF_GENES = set()
        for gene in GASLIT_AF_GENES_TEXT.split():
            if gene.strip():
                GASLIT_AF_GENES.add(gene.strip())

class VcfGeneAnnotator:
    """
    VCF Gene Annotator for GASLIT-AF Variant Analysis.
    
    This class adds gene annotations to VCF files, creating a quantum coherence
    bridge between raw genomic data and the GASLIT-AF theoretical framework.
    """
    
    def __init__(self, reference_genome: str = 'GRCh38'):
        """
        Initialize the VCF Gene Annotator.
        
        Args:
            reference_genome: Reference genome version (GRCh38 or GRCh37)
        """
        self.reference_genome = reference_genome
        
        # Initialize gene position database
        self.gene_db_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            'data',
            f'gene_positions_{reference_genome.lower()}.csv'
        )
        
        # Check if gene position database exists
        if not os.path.exists(self.gene_db_path):
            log.warning(f"Gene position database not found: {self.gene_db_path}")
            log.warning("Please run gene_position_db.py to create the database")
            
        # Load gene positions
        self.load_gene_positions()
        
        # Initialize gene position cache
        self.gene_position_cache = {}
        
        # Load known gene positions from file if available
        self.load_gene_positions()
    
    def load_gene_positions(self, file_path: Optional[str] = None):
        """
        Load known gene positions from file.
        
        Args:
            file_path: Path to gene positions file (CSV)
        """
        if file_path is None:
            file_path = self.gene_db_path
        
        if os.path.exists(file_path):
            try:
                gene_df = pd.read_csv(file_path)
                for _, row in gene_df.iterrows():
                    gene = row['gene']
                    chrom = row['chromosome']
                    start = row['start']
                    end = row['end']
                    self.gene_position_cache[gene] = (chrom, start, end)
                
                log.info(f"Loaded {len(self.gene_position_cache)} gene positions from {file_path}")
            except Exception as e:
                log.error(f"Error loading gene positions: {e}")
    
    def save_gene_positions(self, file_path: Optional[str] = None):
        """
        Save known gene positions to file.
        
        Args:
            file_path: Path to gene positions file (CSV)
        """
        if file_path is None:
            # Create data directory if it doesn't exist
            data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
            os.makedirs(data_dir, exist_ok=True)
            
            file_path = self.gene_db_path
        
        try:
            gene_data = []
            for gene, (chrom, start, end) in self.gene_position_cache.items():
                gene_data.append({
                    'gene': gene,
                    'chromosome': chrom,
                    'start': start,
                    'end': end
                })
            
            gene_df = pd.DataFrame(gene_data)
            gene_df.to_csv(file_path, index=False)
            
            log.info(f"Saved {len(self.gene_position_cache)} gene positions to {file_path}")
        except Exception as e:
            log.error(f"Error saving gene positions: {e}")
    
    def get_gene_at_position(self, chrom: str, pos: int) -> List[str]:
        """
        Get genes at a specific genomic position.
        
        Args:
            chrom: Chromosome (e.g., '1', 'X')
            pos: Position on chromosome
            
        Returns:
            List of gene symbols at the position
        """
        genes = []
        
        # Check gene position cache for overlaps
        for gene, (gene_chrom, start, end) in self.gene_position_cache.items():
            if gene_chrom == chrom and start <= pos <= end:
                if gene not in genes:
                    genes.append(gene)
        
        return genes
    
    def is_gaslit_af_gene(self, gene: str) -> bool:
        """
        Check if a gene is in the GASLIT-AF gene list.
        
        Args:
            gene: Gene symbol
            
        Returns:
            True if the gene is in the GASLIT-AF gene list
        """
        return gene in GASLIT_AF_GENES
    
    def annotate_vcf(self, input_vcf: str, output_vcf: str, chunk_size: int = 10000):
        """
        Annotate a VCF file with gene information.
        
        Args:
            input_vcf: Path to input VCF file
            output_vcf: Path to output VCF file
            chunk_size: Number of variants to process at once
        """
        log.info(f"Annotating VCF file: {input_vcf}")
        
        # Check if input file exists
        if not os.path.exists(input_vcf):
            log.error(f"Input VCF file not found: {input_vcf}")
            return False
        
        # Create output directory if it doesn't exist
        output_dir = os.path.dirname(output_vcf)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
        
        # Open input and output files
        try:
            # Determine if input is gzipped
            is_gzipped = input_vcf.endswith('.gz')
            
            # Open input file
            if is_gzipped:
                in_file = gzip.open(input_vcf, 'rt')
            else:
                in_file = open(input_vcf, 'r')
            
            # Open output file
            if output_vcf.endswith('.gz'):
                out_file = gzip.open(output_vcf, 'wt')
            else:
                out_file = open(output_vcf, 'w')
            
            # Process header lines
            header_lines = []
            for line in in_file:
                if line.startswith('#'):
                    # Add ANN format to header if not present
                    if line.startswith('##INFO') and 'ID=ANN' in line:
                        # ANN header already exists
                        header_lines.append(line)
                    elif line.startswith('#CHROM'):
                        # Add ANN header before CHROM line
                        header_lines.append('##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: Gene|Gene_ID">\n')
                        header_lines.append(line)
                    else:
                        header_lines.append(line)
                else:
                    # First non-header line, write all headers and break
                    for header in header_lines:
                        out_file.write(header)
                    
                    # Reset file pointer to beginning of file
                    in_file.seek(0)
                    
                    # Skip header lines
                    for _ in header_lines:
                        next(in_file)
                    
                    break
            
            # Count total variants for progress tracking
            total_variants = 0
            with gzip.open(input_vcf, 'rt') if is_gzipped else open(input_vcf, 'r') as count_file:
                for line in count_file:
                    if not line.startswith('#'):
                        total_variants += 1
            
            log.info(f"Processing {total_variants} variants")
            
            # Process variants in chunks
            processed_variants = 0
            annotated_variants = 0
            gaslit_af_variants = 0
            
            with Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TextColumn("{task.completed}/{task.total}"),
                TextColumn("({task.percentage:.0f}%)"),
                TimeElapsedColumn(),
                console=console
            ) as progress:
                task = progress.add_task("Annotating variants", total=total_variants)
                
                variant_buffer = []
                for line in in_file:
                    if line.startswith('#'):
                        continue
                    
                    variant_buffer.append(line)
                    
                    if len(variant_buffer) >= chunk_size:
                        # Process buffer
                        processed, annotated, gaslit = self._process_variant_chunk(variant_buffer, out_file)
                        processed_variants += processed
                        annotated_variants += annotated
                        gaslit_af_variants += gaslit
                        
                        # Update progress
                        progress.update(task, advance=processed)
                        
                        # Clear buffer
                        variant_buffer = []
                
                # Process remaining variants
                if variant_buffer:
                    processed, annotated, gaslit = self._process_variant_chunk(variant_buffer, out_file)
                    processed_variants += processed
                    annotated_variants += annotated
                    gaslit_af_variants += gaslit
                    
                    # Update progress
                    progress.update(task, advance=processed)
            
            # Close files
            in_file.close()
            out_file.close()
            
            # Save gene positions
            self.save_gene_positions()
            
            log.info(f"Annotation complete:")
            log.info(f"  - Processed variants: {processed_variants}")
            log.info(f"  - Annotated with genes: {annotated_variants}")
            log.info(f"  - GASLIT-AF gene variants: {gaslit_af_variants}")
            
            return True
        
        except Exception as e:
            log.error(f"Error annotating VCF file: {e}")
            return False
    
    def _process_variant_chunk(self, variant_lines: List[str], out_file) -> Tuple[int, int, int]:
        """
        Process a chunk of variant lines.
        
        Args:
            variant_lines: List of variant lines
            out_file: Output file handle
            
        Returns:
            Tuple of (processed_count, annotated_count, gaslit_af_count)
        """
        processed_count = 0
        annotated_count = 0
        gaslit_af_count = 0
        
        for line in variant_lines:
            processed_count += 1
            
            # Parse variant line
            fields = line.strip().split('\t')
            if len(fields) < 8:
                # Invalid line, write as is
                out_file.write(line)
                continue
            
            # Extract variant information
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            info = fields[7]
            
            # Get genes at this position
            genes = self.get_gene_at_position(chrom, pos)
            
            # Check if any genes are in GASLIT-AF gene list
            gaslit_af_genes = [gene for gene in genes if self.is_gaslit_af_gene(gene)]
            
            if genes:
                annotated_count += 1
                
                # Create ANN field
                ann_entries = []
                for gene in genes:
                    # Simple annotation format: Gene|GeneID
                    ann_entries.append(f"{gene}|{gene}")
                
                ann_field = ','.join(ann_entries)
                
                # Add ANN field to INFO column
                if info == '.':
                    info = f"ANN={ann_field}"
                else:
                    info = f"{info};ANN={ann_field}"
                
                # Update INFO field
                fields[7] = info
                
                # Check if any GASLIT-AF genes
                if gaslit_af_genes:
                    gaslit_af_count += 1
            
            # Write updated line
            out_file.write('\t'.join(fields) + '\n')
        
        return processed_count, annotated_count, gaslit_af_count

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="VCF Gene Annotator for GASLIT-AF Variant Analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "input_vcf",
        help="Input VCF file"
    )
    
    parser.add_argument(
        "output_vcf",
        help="Output VCF file"
    )
    
    parser.add_argument(
        "--reference",
        choices=["GRCh38", "GRCh37"],
        default="GRCh38",
        help="Reference genome version"
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
    Main entry point for VCF Gene Annotator.
    """
    args = parse_args()
    
    # Create annotator
    annotator = VcfGeneAnnotator(reference_genome=args.reference)
    
    # Annotate VCF file
    success = annotator.annotate_vcf(
        input_vcf=args.input_vcf,
        output_vcf=args.output_vcf,
        chunk_size=args.chunk_size
    )
    
    if success:
        log.info(f"Successfully annotated VCF file: {args.output_vcf}")
    else:
        log.error(f"Failed to annotate VCF file")
        sys.exit(1)

if __name__ == "__main__":
    main()
