"""
ClinVar Parser Module for GASLIT-AF Variant Analysis.

This module is responsible for parsing ClinVar data files into usable formats.
"""

import logging
import pandas as pd
import gzip
from pathlib import Path
from typing import Optional, Dict, List

from src.gaslit_af.clinvar.downloader import ClinVarDownloader

# Configure logging
log = logging.getLogger("gaslit-af")

class ClinVarParser:
    """Parses ClinVar data files."""
    
    def __init__(self, downloader: ClinVarDownloader):
        """
        Initialize the ClinVar parser.
        
        Args:
            downloader: ClinVarDownloader instance
        """
        self.downloader = downloader
        self.cache_dir = downloader.cache_dir
    
    def parse_variant_summary(self, force_download: bool = False) -> pd.DataFrame:
        """
        Parse the variant_summary.txt file.
        
        Args:
            force_download: Whether to force download even if file exists
            
        Returns:
            DataFrame containing variant summary data
        """
        # Download the file if needed
        file_path = self.downloader.download_file("variant_summary", force_download)
        
        # Parse the file
        log.info(f"Parsing variant summary file: {file_path}")
        df = pd.read_csv(file_path, sep='\t', compression='gzip', low_memory=False)
        
        # Save a processed version for faster loading next time
        processed_path = self.cache_dir / "variant_summary_processed.parquet"
        df.to_parquet(processed_path, index=False)
        
        log.info(f"Parsed {len(df)} variants from variant summary")
        return df
    
    def parse_vcf(self, assembly: str = "GRCh38", force_download: bool = False) -> pd.DataFrame:
        """
        Parse a ClinVar VCF file.
        
        Args:
            assembly: Genome assembly (GRCh37 or GRCh38)
            force_download: Whether to force download even if file exists
            
        Returns:
            DataFrame containing VCF data
        """
        # Determine file type based on assembly
        file_type = f"vcf_{assembly.lower()}"
        
        # Download the file if needed
        file_path = self.downloader.download_file(file_type, force_download)
        
        # Parse the file
        log.info(f"Parsing VCF file for {assembly}: {file_path}")
        
        # Initialize data collectors
        data = []
        
        # VCF format fields
        vcf_fields = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        
        # Open and parse VCF file
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                # Skip header lines
                if line.startswith('#'):
                    continue
                
                # Split line into fields
                fields = line.strip().split('\t')
                
                # Basic fields
                variant = {
                    "chrom": fields[0],
                    "pos": int(fields[1]),
                    "id": fields[2],
                    "ref": fields[3],
                    "alt": fields[4]
                }
                
                # Parse INFO field
                info_dict = {}
                info_fields = fields[7].split(';')
                for info in info_fields:
                    if '=' in info:
                        key, value = info.split('=', 1)
                        info_dict[key] = value
                
                # Extract clinvar ID and clinical significance
                variant["clinvar_id"] = info_dict.get("CLNID", "").split(',')[0] if "CLNID" in info_dict else ""
                variant["significance"] = info_dict.get("CLNSIG", "").replace('_', ' ') if "CLNSIG" in info_dict else ""
                
                # Extract gene
                variant["gene"] = info_dict.get("GENEINFO", "").split(':')[0] if "GENEINFO" in info_dict else ""
                
                # Add to data
                data.append(variant)
                
                # For memory efficiency, process in batches
                if len(data) >= 100000:
                    log.info(f"Processed {len(data)} variants...")
                    # Process more if needed
        
        # Create DataFrame
        df = pd.DataFrame(data)
        
        # Save a processed version for faster loading next time
        processed_path = self.cache_dir / f"clinvar_vcf_{assembly.lower()}_processed.parquet"
        df.to_parquet(processed_path, index=False)
        
        log.info(f"Parsed {len(df)} variants from VCF")
        return df
    
    def load_processed_vcf(self, assembly: str = "GRCh38", force_reprocess: bool = False) -> pd.DataFrame:
        """
        Load a processed VCF file, reprocessing if necessary.
        
        Args:
            assembly: Genome assembly (GRCh37 or GRCh38)
            force_reprocess: Whether to force reprocessing even if processed file exists
            
        Returns:
            DataFrame containing VCF data
        """
        processed_path = self.cache_dir / f"clinvar_vcf_{assembly.lower()}_processed.parquet"
        
        if not force_reprocess and processed_path.exists():
            log.info(f"Loading processed VCF for {assembly} from {processed_path}")
            return pd.read_parquet(processed_path)
        
        return self.parse_vcf(assembly=assembly, force_download=force_reprocess)
