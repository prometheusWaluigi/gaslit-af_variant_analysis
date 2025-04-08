#!/usr/bin/env python3
"""
Gene Position Database Generator for GASLIT-AF Variant Analysis

This script creates a quantum coherence bridge between genomic coordinates
and gene symbols by downloading and processing gene position data from
NCBI Gene and GENCODE databases.

The script establishes a recursive mapping between chromosomal positions
and the biological systems in the GASLIT-AF model, allowing for variant
enrichment analysis without requiring external annotation tools.
"""

import os
import sys
import gzip
import json
import logging
import argparse
import requests
import pandas as pd
from pathlib import Path
from typing import Dict, List, Set, Tuple, Optional
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af-gene-db")

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

class GenePositionDatabase:
    """
    Gene Position Database for GASLIT-AF Variant Analysis.
    
    This class creates and manages a database of gene positions,
    establishing a quantum coherence bridge between genomic coordinates
    and gene symbols.
    """
    
    def __init__(self, data_dir: Optional[str] = None, reference: str = 'GRCh38'):
        """
        Initialize the Gene Position Database.
        
        Args:
            data_dir: Directory to store gene position data
            reference: Reference genome version (GRCh38 or GRCh37)
        """
        self.reference = reference
        
        # Set data directory
        if data_dir is None:
            self.data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
        else:
            self.data_dir = data_dir
        
        # Create data directory if it doesn't exist
        os.makedirs(self.data_dir, exist_ok=True)
        
        # Set file paths
        self.gene_info_path = os.path.join(self.data_dir, 'gene_info.gz')
        self.gene_history_path = os.path.join(self.data_dir, 'gene_history.gz')
        self.gencode_path = os.path.join(self.data_dir, f'gencode_{reference.lower()}.gz')
        self.gene_positions_path = os.path.join(self.data_dir, f'gene_positions_{reference.lower()}.csv')
        
        # Initialize gene position database
        self.gene_positions = {}
        
        # Load gene positions if file exists
        if os.path.exists(self.gene_positions_path):
            self.load_gene_positions()
        
        # Get GASLIT-AF genes
        self.gaslit_af_genes = set(GASLIT_AF_GENES)
        log.info(f"Loaded {len(self.gaslit_af_genes)} GASLIT-AF genes")
    
    def load_gene_positions(self):
        """
        Load gene positions from CSV file.
        """
        try:
            gene_df = pd.read_csv(self.gene_positions_path)
            for _, row in gene_df.iterrows():
                gene = row['gene']
                chrom = row['chromosome']
                start = row['start']
                end = row['end']
                self.gene_positions[gene] = (chrom, start, end)
            
            log.info(f"Loaded {len(self.gene_positions)} gene positions from {self.gene_positions_path}")
        except Exception as e:
            log.error(f"Error loading gene positions: {e}")
    
    def save_gene_positions(self):
        """
        Save gene positions to CSV file.
        """
        try:
            gene_data = []
            for gene, (chrom, start, end) in self.gene_positions.items():
                gene_data.append({
                    'gene': gene,
                    'chromosome': chrom,
                    'start': start,
                    'end': end
                })
            
            gene_df = pd.DataFrame(gene_data)
            gene_df.to_csv(self.gene_positions_path, index=False)
            
            log.info(f"Saved {len(self.gene_positions)} gene positions to {self.gene_positions_path}")
        except Exception as e:
            log.error(f"Error saving gene positions: {e}")
    
    def download_ncbi_gene_info(self):
        """
        Download gene information from NCBI Gene database.
        """
        gene_info_url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz"
        gene_history_url = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_history.gz"
        
        # Download gene_info.gz
        if not os.path.exists(self.gene_info_path):
            log.info(f"Downloading gene information from NCBI: {gene_info_url}")
            try:
                response = requests.get(gene_info_url, stream=True)
                response.raise_for_status()
                
                with open(self.gene_info_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                
                log.info(f"Downloaded gene information to {self.gene_info_path}")
            except Exception as e:
                log.error(f"Error downloading gene information: {e}")
                return False
        
        # Download gene_history.gz
        if not os.path.exists(self.gene_history_path):
            log.info(f"Downloading gene history from NCBI: {gene_history_url}")
            try:
                response = requests.get(gene_history_url, stream=True)
                response.raise_for_status()
                
                with open(self.gene_history_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                
                log.info(f"Downloaded gene history to {self.gene_history_path}")
            except Exception as e:
                log.error(f"Error downloading gene history: {e}")
                return False
        
        return True
    
    def download_gencode_data(self):
        """
        Download gene annotation data from GENCODE.
        """
        if self.reference == 'GRCh38':
            gencode_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gff3.gz"
        else:
            gencode_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/GRCh37_mapping/gencode.v44lift37.annotation.gff3.gz"
        
        if not os.path.exists(self.gencode_path):
            log.info(f"Downloading GENCODE data: {gencode_url}")
            try:
                response = requests.get(gencode_url, stream=True)
                response.raise_for_status()
                
                with open(self.gencode_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        f.write(chunk)
                
                log.info(f"Downloaded GENCODE data to {self.gencode_path}")
            except Exception as e:
                log.error(f"Error downloading GENCODE data: {e}")
                return False
        
        return True
    
    def process_ncbi_gene_info(self):
        """
        Process NCBI gene information to extract human gene positions.
        """
        log.info("Processing NCBI gene information")
        
        # Check if gene_info.gz exists
        if not os.path.exists(self.gene_info_path):
            log.error(f"Gene information file not found: {self.gene_info_path}")
            return False
        
        try:
            # Process gene_info.gz
            with gzip.open(self.gene_info_path, 'rt') as f:
                # Skip header
                next(f)
                
                # Process each line
                for line in f:
                    fields = line.strip().split('\t')
                    
                    # Check if this is a human gene
                    if fields[0] == '9606':  # Taxonomy ID for Homo sapiens
                        gene_id = fields[1]
                        symbol = fields[2]
                        
                        # Check if this is a GASLIT-AF gene
                        if symbol in self.gaslit_af_genes:
                            # Extract chromosome, start, and end positions
                            map_location = fields[7]
                            
                            # Parse map location
                            if map_location != '-':
                                # Extract chromosome
                                chrom = map_location.split('|')[0].split('q')[0].split('p')[0]
                                
                                # Clean up chromosome name
                                if chrom.startswith('chr'):
                                    chrom = chrom[3:]
                                
                                # Add to gene positions
                                # We'll update start and end positions later
                                self.gene_positions[symbol] = (chrom, 0, 0)
            
            log.info(f"Processed NCBI gene information, found {len(self.gene_positions)} GASLIT-AF genes")
            return True
        except Exception as e:
            log.error(f"Error processing NCBI gene information: {e}")
            return False
    
    def process_gencode_data(self):
        """
        Process GENCODE data to extract gene positions.
        """
        log.info("Processing GENCODE data")
        
        # Check if gencode file exists
        if not os.path.exists(self.gencode_path):
            log.error(f"GENCODE file not found: {self.gencode_path}")
            return False
        
        try:
            # Process gencode file
            with gzip.open(self.gencode_path, 'rt') as f:
                # Skip header lines
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    
                    # Check if this is a gene feature
                    if fields[2] == 'gene':
                        # Extract gene information
                        chrom = fields[0]
                        start = int(fields[3])
                        end = int(fields[4])
                        
                        # Clean up chromosome name
                        if chrom.startswith('chr'):
                            chrom = chrom[3:]
                        
                        # Extract gene symbol from attributes
                        attributes = fields[8].split(';')
                        gene_name = None
                        
                        for attr in attributes:
                            if attr.startswith('gene_name='):
                                gene_name = attr.split('=')[1]
                                break
                        
                        # Check if this is a GASLIT-AF gene
                        if gene_name in self.gaslit_af_genes:
                            # Update gene position
                            self.gene_positions[gene_name] = (chrom, start, end)
            
            log.info(f"Processed GENCODE data, updated {len(self.gene_positions)} gene positions")
            return True
        except Exception as e:
            log.error(f"Error processing GENCODE data: {e}")
            return False
    
    def create_gene_position_database(self):
        """
        Create gene position database by downloading and processing
        gene information from NCBI and GENCODE.
        """
        # Download NCBI gene information
        if not self.download_ncbi_gene_info():
            return False
        
        # Download GENCODE data
        if not self.download_gencode_data():
            return False
        
        # Process NCBI gene information
        if not self.process_ncbi_gene_info():
            return False
        
        # Process GENCODE data
        if not self.process_gencode_data():
            return False
        
        # Save gene positions
        self.save_gene_positions()
        
        return True
    
    def create_manual_gene_positions(self):
        """
        Create manual gene positions for GASLIT-AF genes.
        
        This is a fallback method when online resources are not available.
        """
        log.info("Creating manual gene positions for GASLIT-AF genes")
        
        # Define manual gene positions for key GASLIT-AF genes
        manual_positions = {
            # Immune & Inflammatory System
            "IDO2": ("8", 39754899, 39785977),
            "AHR": ("7", 17298622, 17346152),
            "AHRR": ("5", 304291, 438405),
            "IL36RN": ("2", 112750310, 112756266),
            "CFH": ("1", 196621007, 196716634),
            "MBL2": ("10", 52760791, 52770796),
            "NLRP3": ("1", 247416156, 247449108),
            "IL1B": ("2", 112829751, 112837503),
            "IL6": ("7", 22725889, 22732002),
            "IL17": ("6", 52051185, 52055436),
            "IL13": ("5", 132656263, 132658199),
            "IL4": ("5", 132673986, 132682977),
            "HLA-DQB1": ("6", 32627244, 32634466),
            "PTPN22": ("1", 113813811, 113871761),
            "CTLA4": ("2", 203867786, 203873960),
            
            # Autonomic & Neurotransmitter System
            "COMT": ("22", 19929262, 19957498),
            "CHRM2": ("7", 136553416, 136705002),
            "DRD2": ("11", 113280317, 113346413),
            "GABRA1": ("5", 161274197, 161326975),
            "CHRNA7": ("15", 32322685, 32464722),
            "ADRB1": ("10", 114044167, 114049628),
            "ADRB2": ("5", 148826593, 148828634),
            "NOS3": ("7", 150688083, 150711687),
            "GNB3": ("12", 6819635, 6826818),
            "SLC6A2": ("16", 55665421, 55717858),
            "SLC6A4": ("17", 30194318, 30235968),
            "HTR2A": ("13", 46831542, 46897509),
            "OPRM1": ("6", 154331631, 154568001),
            
            # Structural & Connective Tissue System
            "TNXB": ("6", 32041155, 32115334),
            "COL3A1": ("2", 188974758, 189010070),
            "FBN1": ("15", 48408313, 48645925),
            "FLNA": ("X", 154348529, 154374699),
            "COL5A1": ("9", 134641933, 134844885),
            
            # Metabolic System
            "APOE": ("19", 44905781, 44909393),
            "PCSK9": ("1", 55039475, 55064852),
            "UGT1A1": ("2", 233760313, 233763087),
            "HNF1A": ("12", 120978516, 121002512),
            "ABCC8": ("11", 17414432, 17498392),
            "TFAM": ("10", 58385022, 58399221),
            "PDHA1": ("X", 19343703, 19366354),
            "SDHB": ("1", 17018722, 17054170),
            
            # Endocannabinoid System
            "CNR1": ("6", 88139864, 88166359),
            "CNR2": ("1", 23870455, 23878322),
            "FAAH": ("1", 46028686, 46048410),
            "MGLL": ("3", 127407836, 127543320),
            
            # Calcium & Ion Channels
            "ITPR1": ("3", 4691693, 4889520),
            "KCNJ5": ("11", 128761823, 128786607),
            "RYR2": ("1", 236749921, 237115382),
            
            # Mast Cell Activation
            "TPSAB1": ("16", 1268786, 1270544),
            "KIT": ("4", 54657918, 54740715),
            "HNMT": ("2", 138438276, 138479452),
            "TET2": ("4", 105233932, 105312855),
            
            # Kynurenine Pathway
            "IDO1": ("8", 39779230, 39791142),
            "KMO": ("1", 241747259, 241768864),
            "KYNU": ("2", 143040812, 143078487),
            "TDO2": ("4", 156317203, 156335892),
            "HAAO": ("2", 42751192, 42766927)
        }
        
        # Add manual positions to gene positions
        for gene, (chrom, start, end) in manual_positions.items():
            if gene in self.gaslit_af_genes:
                self.gene_positions[gene] = (chrom, start, end)
        
        # Save gene positions
        self.save_gene_positions()
        
        log.info(f"Created manual gene positions for {len(self.gene_positions)} GASLIT-AF genes")
        return True
    
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
        
        # Check gene position database for overlaps
        for gene, (gene_chrom, start, end) in self.gene_positions.items():
            if gene_chrom == chrom and start <= pos <= end:
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
        return gene in self.gaslit_af_genes

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Gene Position Database Generator for GASLIT-AF Variant Analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--data-dir",
        help="Directory to store gene position data"
    )
    
    parser.add_argument(
        "--reference",
        choices=["GRCh38", "GRCh37"],
        default="GRCh38",
        help="Reference genome version"
    )
    
    parser.add_argument(
        "--manual",
        action="store_true",
        help="Create manual gene positions (offline mode)"
    )
    
    return parser.parse_args()

def main():
    """
    Main entry point for Gene Position Database Generator.
    """
    args = parse_args()
    
    # Create gene position database
    db = GenePositionDatabase(data_dir=args.data_dir, reference=args.reference)
    
    if args.manual:
        # Create manual gene positions
        success = db.create_manual_gene_positions()
    else:
        # Create gene position database
        success = db.create_gene_position_database()
    
    if success:
        log.info(f"Successfully created gene position database: {db.gene_positions_path}")
    else:
        log.error(f"Failed to create gene position database")
        sys.exit(1)

if __name__ == "__main__":
    main()
