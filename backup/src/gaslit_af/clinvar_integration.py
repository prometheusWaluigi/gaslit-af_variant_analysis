"""
ClinVar Integration Module for GASLIT-AF Variant Analysis.

This module handles downloading, parsing, and caching ClinVar data locally.
It provides utilities to integrate ClinVar annotations with variant analysis results.
"""

import os
import gzip
import json
import logging
import requests
import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from datetime import datetime

# Configure logging
log = logging.getLogger("gaslit-af")

class ClinVarDownloader:
    """Handles downloading ClinVar data files."""
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the ClinVar downloader.
        
        Args:
            cache_dir: Directory to cache downloaded files
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path("./cache/clinvar")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Base URLs
        self.ftp_base_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
        
        # Available file types
        self.file_types = {
            "vcv_xml": f"{self.ftp_base_url}/xml/ClinVarVCVRelease_00-latest.xml.gz",
            "rcv_xml": f"{self.ftp_base_url}/xml/RCV_release/ClinVarRCVRelease_00-latest.xml.gz",
            "variant_summary": f"{self.ftp_base_url}/tab_delimited/variant_summary.txt.gz",
            "vcf_grch37": f"{self.ftp_base_url}/vcf_GRCh37/clinvar.vcf.gz",
            "vcf_grch38": f"{self.ftp_base_url}/vcf_GRCh38/clinvar.vcf.gz",
            "var_citations": f"{self.ftp_base_url}/tab_delimited/var_citations.txt",
            "cross_references": f"{self.ftp_base_url}/tab_delimited/cross_references.txt"
        }
    
    def download_file(self, file_type: str, force_download: bool = False) -> Path:
        """
        Download a ClinVar file.
        
        Args:
            file_type: Type of file to download (from self.file_types)
            force_download: Whether to force download even if file exists
            
        Returns:
            Path to the downloaded file
        """
        if file_type not in self.file_types:
            raise ValueError(f"Unknown file type: {file_type}. Available types: {list(self.file_types.keys())}")
        
        url = self.file_types[file_type]
        filename = url.split("/")[-1]
        output_path = self.cache_dir / filename
        
        # Check if file already exists and is recent (less than 30 days old)
        if not force_download and output_path.exists():
            file_age_days = (datetime.now() - datetime.fromtimestamp(output_path.stat().st_mtime)).days
            if file_age_days < 30:
                log.info(f"Using cached {file_type} file (age: {file_age_days} days)")
                return output_path
        
        # Download the file
        log.info(f"Downloading {file_type} from {url}")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        log.info(f"Downloaded {file_type} to {output_path}")
        return output_path
    
    def get_latest_release_date(self) -> str:
        """
        Get the latest ClinVar release date.
        
        Returns:
            Latest release date as a string (YYYY-MM-DD)
        """
        try:
            # Check the variant_summary file for release date
            url = self.file_types["variant_summary"]
            response = requests.head(url)
            response.raise_for_status()
            
            # Extract date from Last-Modified header
            last_modified = response.headers.get('Last-Modified')
            if last_modified:
                release_date = datetime.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
                return release_date.strftime('%Y-%m-%d')
            
            return "Unknown"
        except Exception as e:
            log.error(f"Error getting latest release date: {e}")
            return "Unknown"


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
    
    def load_variant_summary(self, force_reprocess: bool = False) -> pd.DataFrame:
        """
        Load the processed variant summary data.
        
        Args:
            force_reprocess: Whether to force reprocessing even if processed file exists
            
        Returns:
            DataFrame containing variant summary data
        """
        processed_path = self.cache_dir / "variant_summary_processed.parquet"
        
        if not force_reprocess and processed_path.exists():
            log.info(f"Loading processed variant summary from {processed_path}")
            return pd.read_parquet(processed_path)
        
        return self.parse_variant_summary(force_download=force_reprocess)
    
    def parse_vcf(self, assembly: str = "GRCh38", force_download: bool = False) -> pd.DataFrame:
        """
        Parse the ClinVar VCF file.
        
        Args:
            assembly: Genome assembly version (GRCh37 or GRCh38)
            force_download: Whether to force download even if file exists
            
        Returns:
            DataFrame containing VCF data
        """
        # Determine file type based on assembly
        file_type = f"vcf_{assembly.lower()}"
        if file_type not in self.downloader.file_types:
            raise ValueError(f"Unknown assembly: {assembly}. Available: GRCh37, GRCh38")
        
        # Download the file if needed
        file_path = self.downloader.download_file(file_type, force_download)
        
        # Parse VCF file
        log.info(f"Parsing VCF file for {assembly}: {file_path}")
        
        # Extract header and column names
        header_lines = []
        with gzip.open(file_path, 'rt') as f:
            for line in f:
                if line.startswith('#'):
                    header_lines.append(line.strip())
                else:
                    break
        
        # Get column names from the last header line
        column_names = header_lines[-1].replace('#', '').split('\t')
        
        # Read the VCF data
        vcf_df = pd.read_csv(
            file_path, 
            sep='\t',
            comment='#',
            names=column_names,
            compression='gzip',
            low_memory=False
        )
        
        # Process INFO field to extract clinical significance
        def extract_clnsig(info_str):
            if 'CLNSIG=' in info_str:
                clnsig_part = info_str.split('CLNSIG=')[1].split(';')[0]
                return clnsig_part
            return None
        
        vcf_df['CLNSIG'] = vcf_df['INFO'].apply(extract_clnsig)
        
        # Save a processed version for faster loading next time
        processed_path = self.cache_dir / f"clinvar_vcf_{assembly.lower()}_processed.parquet"
        vcf_df.to_parquet(processed_path, index=False)
        
        log.info(f"Parsed {len(vcf_df)} variants from VCF for {assembly}")
        return vcf_df
    
    def load_vcf(self, assembly: str = "GRCh38", force_reprocess: bool = False) -> pd.DataFrame:
        """
        Load the processed VCF data.
        
        Args:
            assembly: Genome assembly version (GRCh37 or GRCh38)
            force_reprocess: Whether to force reprocessing even if processed file exists
            
        Returns:
            DataFrame containing VCF data
        """
        processed_path = self.cache_dir / f"clinvar_vcf_{assembly.lower()}_processed.parquet"
        
        if not force_reprocess and processed_path.exists():
            log.info(f"Loading processed VCF for {assembly} from {processed_path}")
            return pd.read_parquet(processed_path)
        
        return self.parse_vcf(assembly=assembly, force_download=force_reprocess)


class ClinVarIntegration:
    """Integrates ClinVar data with variant analysis."""
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the ClinVar integration.
        
        Args:
            cache_dir: Directory to cache ClinVar data
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path("./cache/clinvar")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        self.downloader = ClinVarDownloader(self.cache_dir)
        self.parser = ClinVarParser(self.downloader)
        
        # Cached data
        self._variant_summary = None
        self._vcf_grch37 = None
        self._vcf_grch38 = None
    
    def get_release_info(self) -> Dict[str, str]:
        """
        Get information about the latest ClinVar release.
        
        Returns:
            Dictionary with release information
        """
        release_date = self.downloader.get_latest_release_date()
        
        return {
            "release_date": release_date,
            "cache_dir": str(self.cache_dir)
        }
    
    def get_variant_summary(self, force_reload: bool = False) -> pd.DataFrame:
        """
        Get the variant summary data.
        
        Args:
            force_reload: Whether to force reload the data
            
        Returns:
            DataFrame containing variant summary data
        """
        if self._variant_summary is None or force_reload:
            self._variant_summary = self.parser.load_variant_summary(force_reprocess=force_reload)
        
        return self._variant_summary
    
    def get_vcf(self, assembly: str = "GRCh38", force_reload: bool = False) -> pd.DataFrame:
        """
        Get the VCF data.
        
        Args:
            assembly: Genome assembly version (GRCh37 or GRCh38)
            force_reload: Whether to force reload the data
            
        Returns:
            DataFrame containing VCF data
        """
        if assembly.lower() == "grch37":
            if self._vcf_grch37 is None or force_reload:
                self._vcf_grch37 = self.parser.load_vcf(assembly="GRCh37", force_reprocess=force_reload)
            return self._vcf_grch37
        else:
            if self._vcf_grch38 is None or force_reload:
                self._vcf_grch38 = self.parser.load_vcf(assembly="GRCh38", force_reprocess=force_reload)
            return self._vcf_grch38
    
    def annotate_variants(self, variants_df: pd.DataFrame, assembly: str = "GRCh38") -> pd.DataFrame:
        """
        Annotate variants with ClinVar data.
        
        Args:
            variants_df: DataFrame of variants
            assembly: Genome assembly version (GRCh37 or GRCh38)
            
        Returns:
            DataFrame with ClinVar annotations
        """
        if variants_df is None or variants_df.empty:
            return variants_df
        
        # Create a copy to avoid modifying the original
        annotated_df = variants_df.copy()
        
        # Load ClinVar data
        variant_summary = self.get_variant_summary()
        
        # Add ClinVar columns
        annotated_df['clinvar_id'] = None
        annotated_df['clinvar_significance'] = None
        annotated_df['clinvar_review_status'] = None
        annotated_df['clinvar_conditions'] = None
        
        # Determine which column to use for matching
        rsid_col = next((col for col in ["rsid", "variant_id", "id"] if col in annotated_df.columns), None)
        
        if rsid_col is not None:
            # Match by rsID
            log.info(f"Annotating variants by rsID using column: {rsid_col}")
            
            # Filter out missing or invalid rsIDs
            valid_variants = annotated_df[annotated_df[rsid_col].notna() & 
                                        annotated_df[rsid_col].str.startswith('rs', na=False)]
            
            if not valid_variants.empty:
                # Create a mapping of rsID to ClinVar data
                rsid_to_clinvar = {}
                for _, row in variant_summary[variant_summary['RS# (dbSNP)'].notna()].iterrows():
                    rsid = f"rs{int(row['RS# (dbSNP)'])}"
                    rsid_to_clinvar[rsid] = {
                        'clinvar_id': row['#AlleleID'],
                        'clinvar_significance': row['ClinicalSignificance'],
                        'clinvar_review_status': row['ReviewStatus'],
                        'clinvar_conditions': row['PhenotypeList']
                    }
                
                # Apply annotations
                for idx, row in valid_variants.iterrows():
                    rsid = row[rsid_col]
                    if rsid in rsid_to_clinvar:
                        for col, value in rsid_to_clinvar[rsid].items():
                            annotated_df.at[idx, col] = value
        
        else:
            # Try to match by position
            chrom_col = next((col for col in ["chrom", "chromosome", "chr"] if col in annotated_df.columns), None)
            pos_col = next((col for col in ["pos", "position", "start"] if col in annotated_df.columns), None)
            ref_col = next((col for col in ["ref", "reference"] if col in annotated_df.columns), None)
            alt_col = next((col for col in ["alt", "alternate"] if col in annotated_df.columns), None)
            
            if all([chrom_col, pos_col, ref_col, alt_col]):
                log.info(f"Annotating variants by position using columns: {chrom_col}, {pos_col}, {ref_col}, {alt_col}")
                
                # Load VCF data for the specified assembly
                vcf_df = self.get_vcf(assembly=assembly)
                
                # Standardize chromosome format
                def standardize_chrom(chrom):
                    if isinstance(chrom, str) and chrom.startswith('chr'):
                        return chrom[3:]
                    return str(chrom)
                
                annotated_df['_chrom'] = annotated_df[chrom_col].apply(standardize_chrom)
                vcf_df['_CHROM'] = vcf_df['#CHROM'].apply(standardize_chrom)
                
                # Match variants by position and alleles
                for idx, row in annotated_df.iterrows():
                    chrom = row['_chrom']
                    pos = row[pos_col]
                    ref = row[ref_col]
                    alt = row[alt_col]
                    
                    # Find matching variants in VCF
                    matches = vcf_df[
                        (vcf_df['_CHROM'] == chrom) & 
                        (vcf_df['POS'] == pos) & 
                        (vcf_df['REF'] == ref) & 
                        (vcf_df['ALT'] == alt)
                    ]
                    
                    if not matches.empty:
                        match = matches.iloc[0]
                        
                        # Extract ClinVar ID from INFO field
                        if 'ALLELEID=' in match['INFO']:
                            clinvar_id = match['INFO'].split('ALLELEID=')[1].split(';')[0]
                            annotated_df.at[idx, 'clinvar_id'] = clinvar_id
                        
                        # Extract clinical significance
                        if match['CLNSIG'] is not None:
                            annotated_df.at[idx, 'clinvar_significance'] = match['CLNSIG']
                        
                        # Extract review status
                        if 'CLNREVSTAT=' in match['INFO']:
                            review_status = match['INFO'].split('CLNREVSTAT=')[1].split(';')[0]
                            annotated_df.at[idx, 'clinvar_review_status'] = review_status
                        
                        # Extract conditions
                        if 'CLNDISDB=' in match['INFO']:
                            conditions = match['INFO'].split('CLNDISDB=')[1].split(';')[0]
                            annotated_df.at[idx, 'clinvar_conditions'] = conditions
                
                # Drop temporary column
                annotated_df = annotated_df.drop('_chrom', axis=1)
            
            else:
                log.warning("Could not annotate variants: no suitable columns for matching")
        
        # Count annotations
        annotation_count = annotated_df['clinvar_id'].notna().sum()
        log.info(f"Annotated {annotation_count} out of {len(annotated_df)} variants with ClinVar data")
        
        return annotated_df
    
    def get_variant_details(self, variant_id: str) -> Dict:
        """
        Get detailed information for a variant from ClinVar.
        
        Args:
            variant_id: Variant identifier (rsID or ClinVar ID)
            
        Returns:
            Dictionary with variant details
        """
        # Load ClinVar data
        variant_summary = self.get_variant_summary()
        
        # Check if variant_id is an rsID
        if variant_id.startswith('rs'):
            rs_number = int(variant_id[2:])
            matches = variant_summary[variant_summary['RS# (dbSNP)'] == rs_number]
        
        # Check if variant_id is a ClinVar ID
        elif variant_id.isdigit():
            matches = variant_summary[variant_summary['#AlleleID'] == int(variant_id)]
        
        else:
            return {"error": f"Invalid variant ID format: {variant_id}"}
        
        if matches.empty:
            return {"error": f"Variant not found in ClinVar: {variant_id}"}
        
        # Get the first match
        variant = matches.iloc[0]
        
        # Extract details
        details = {
            "variant_id": variant_id,
            "clinvar_id": variant['#AlleleID'],
            "gene": variant['GeneSymbol'],
            "chromosome": variant['Chromosome'],
            "start": variant['Start'],
            "stop": variant['Stop'],
            "reference_allele": variant['ReferenceAllele'],
            "alternate_allele": variant['AlternateAllele'],
            "clinical_significance": variant['ClinicalSignificance'],
            "review_status": variant['ReviewStatus'],
            "conditions": variant['PhenotypeList'].split(';') if isinstance(variant['PhenotypeList'], str) else [],
            "last_evaluated": variant['LastEvaluated'],
            "rs_db_snp": f"rs{int(variant['RS# (dbSNP)'])}" if not pd.isna(variant['RS# (dbSNP)']) else None,
            "nsv_db_var": variant['nsv/esv (dbVar)'],
            "origin": variant['Origin'],
            "assembly": variant['Assembly'],
            "cytogenetic": variant['Cytogenetic'],
            "review_status_stars": variant['ReviewStatus'].count('star') if isinstance(variant['ReviewStatus'], str) else 0,
            "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variant['#AlleleID']}/"
        }
        
        return details


# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Initialize ClinVar integration
    clinvar = ClinVarIntegration()
    
    # Get release info
    release_info = clinvar.get_release_info()
    print(f"ClinVar Release Date: {release_info['release_date']}")
    print(f"Cache Directory: {release_info['cache_dir']}")
    
    # Get variant details for APOE e4 (rs429358)
    variant_details = clinvar.get_variant_details("rs429358")
    print("\nVariant Details for APOE e4 (rs429358):")
    for key, value in variant_details.items():
        print(f"  {key}: {value}")
    
    # Create a test DataFrame
    test_df = pd.DataFrame({
        "rsid": ["rs429358", "rs7412", "rs6311"],
        "gene": ["APOE", "APOE", "HTR2A"],
        "chrom": ["19", "19", "13"],
        "pos": [44908684, 44908822, 47409034],
        "ref": ["T", "C", "G"],
        "alt": ["C", "T", "A"]
    })
    
    # Annotate variants
    annotated_df = clinvar.annotate_variants(test_df)
    print("\nAnnotated Variants:")
    print(annotated_df[["rsid", "gene", "clinvar_significance", "clinvar_review_status"]].to_string(index=False))
