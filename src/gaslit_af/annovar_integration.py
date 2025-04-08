"""
ANNOVAR integration module for GASLIT-AF Variant Analysis.

This module provides integration with ANNOVAR for advanced variant annotation,
creating a recursive enrichment layer that enhances the quantum coherence between
genomic architecture and the GASLIT-AF model parameters.
"""

import os
import logging
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Set, Optional, Any, Union

# Configure logging
log = logging.getLogger("gaslit-af")

class AnnovarIntegration:
    """
    Integration with ANNOVAR for advanced variant annotation.
    
    This class provides methods to convert VCF files to ANNOVAR format,
    run ANNOVAR annotation, and integrate the results back into the
    GASLIT-AF variant analysis workflow.
    """
    
    def __init__(self, annovar_path: Optional[str] = None, 
                 humandb_path: Optional[str] = None,
                 build_version: str = "hg38"):
        """
        Initialize ANNOVAR integration.
        
        Args:
            annovar_path: Path to ANNOVAR installation directory
            humandb_path: Path to ANNOVAR humandb directory
            build_version: Genome build version (default: hg38)
        """
        # Find ANNOVAR path if not provided
        self.annovar_path = annovar_path or self._find_annovar_path()
        
        # Set humandb path
        if humandb_path:
            self.humandb_path = humandb_path
        elif self.annovar_path:
            self.humandb_path = os.path.join(self.annovar_path, "humandb")
        else:
            self.humandb_path = None
            
        # Set build version
        self.build_version = build_version
        
        # Check if ANNOVAR is available
        self.is_available = self._check_annovar_availability()
        
        if self.is_available:
            log.info(f"ANNOVAR integration initialized: {self.annovar_path}")
            log.info(f"Using humandb: {self.humandb_path}")
        else:
            log.warning("ANNOVAR not found or not properly configured")
    
    def _find_annovar_path(self) -> Optional[str]:
        """
        Find ANNOVAR installation path.
        
        Returns:
            Path to ANNOVAR installation directory or None if not found
        """
        # Check common installation locations
        possible_paths = [
            "/home/k10/dev/windsage/gaslit-af_variant_analysis/tools/annovar",
            "/usr/local/annovar",
            "/opt/annovar",
            os.path.expanduser("~/annovar")
        ]
        
        for path in possible_paths:
            if os.path.exists(path) and os.path.exists(os.path.join(path, "table_annovar.pl")):
                return path
        
        # Try to find in PATH
        try:
            result = subprocess.run(["which", "table_annovar.pl"], 
                                   capture_output=True, text=True, check=False)
            if result.returncode == 0:
                # Extract directory from path
                bin_path = result.stdout.strip()
                return os.path.dirname(bin_path)
        except Exception as e:
            log.debug(f"Error finding ANNOVAR in PATH: {e}")
        
        return None
    
    def _check_annovar_availability(self) -> bool:
        """
        Check if ANNOVAR is available and properly configured.
        
        Returns:
            True if ANNOVAR is available, False otherwise
        """
        if not self.annovar_path:
            return False
        
        # Check for essential ANNOVAR scripts
        required_scripts = ["table_annovar.pl", "convert2annovar.pl"]
        for script in required_scripts:
            script_path = os.path.join(self.annovar_path, script)
            if not os.path.exists(script_path):
                log.warning(f"Required ANNOVAR script not found: {script_path}")
                return False
        
        # Check for humandb directory
        if not self.humandb_path or not os.path.exists(self.humandb_path):
            log.warning(f"ANNOVAR humandb directory not found: {self.humandb_path}")
            return False
        
        return True
    
    def download_databases(self, databases: List[str]) -> bool:
        """
        Download required ANNOVAR databases.
        
        Args:
            databases: List of database names to download
            
        Returns:
            True if successful, False otherwise
        """
        if not self.is_available:
            log.error("ANNOVAR not available, cannot download databases")
            return False
        
        success = True
        for db in databases:
            log.info(f"Downloading ANNOVAR database: {db}")
            cmd = [
                "perl", os.path.join(self.annovar_path, "annotate_variation.pl"),
                "-buildver", self.build_version,
                "-downdb", db,
                self.humandb_path
            ]
            
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=False)
                if result.returncode != 0:
                    log.error(f"Error downloading database {db}: {result.stderr}")
                    success = False
                else:
                    log.info(f"Successfully downloaded database: {db}")
            except Exception as e:
                log.error(f"Exception downloading database {db}: {e}")
                success = False
        
        return success
    
    def convert_vcf_to_annovar(self, vcf_path: str, output_dir: str) -> Optional[str]:
        """
        Convert VCF file to ANNOVAR format.
        
        Args:
            vcf_path: Path to VCF file
            output_dir: Output directory
            
        Returns:
            Path to converted ANNOVAR file or None if conversion failed
        """
        if not self.is_available:
            log.error("ANNOVAR not available, cannot convert VCF")
            return None
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Get base filename without extension
        base_name = os.path.basename(vcf_path)
        if base_name.endswith(".vcf.gz"):
            base_name = base_name[:-7]
        elif base_name.endswith(".vcf"):
            base_name = base_name[:-4]
        
        # Output file path
        output_file = os.path.join(output_dir, f"{base_name}.avinput")
        
        # Convert VCF to ANNOVAR format
        log.info(f"Converting VCF to ANNOVAR format: {vcf_path}")
        cmd = [
            "perl", os.path.join(self.annovar_path, "convert2annovar.pl"),
            "-format", "vcf4",
            vcf_path,
            "-outfile", output_file,
            "-includeinfo"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                log.error(f"Error converting VCF to ANNOVAR format: {result.stderr}")
                return None
            
            log.info(f"Successfully converted VCF to ANNOVAR format: {output_file}")
            return output_file
        except Exception as e:
            log.error(f"Exception converting VCF to ANNOVAR format: {e}")
            return None
    
    def annotate_variants(self, input_file: str, output_dir: str, 
                         protocols: List[str] = None,
                         operations: List[str] = None) -> Optional[str]:
        """
        Annotate variants using ANNOVAR.
        
        Args:
            input_file: Path to input file in ANNOVAR format
            output_dir: Output directory
            protocols: List of annotation protocols (databases)
            operations: List of operations for each protocol
            
        Returns:
            Path to annotated output file or None if annotation failed
        """
        if not self.is_available:
            log.error("ANNOVAR not available, cannot annotate variants")
            return None
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Get base filename without extension
        base_name = os.path.basename(input_file)
        if base_name.endswith(".avinput"):
            base_name = base_name[:-8]
        
        # Default protocols and operations if not provided
        if not protocols:
            protocols = ["refGene", "exac03", "gnomad211_exome", "clinvar_20220320", "dbnsfp42a"]
        
        if not operations:
            operations = ["g"] * len(protocols)
        
        if len(protocols) != len(operations):
            log.error("Number of protocols must match number of operations")
            return None
        
        # Output prefix
        output_prefix = os.path.join(output_dir, base_name)
        
        # Run ANNOVAR annotation
        log.info(f"Annotating variants with ANNOVAR: {input_file}")
        cmd = [
            "perl", os.path.join(self.annovar_path, "table_annovar.pl"),
            input_file,
            self.humandb_path,
            "-buildver", self.build_version,
            "-out", output_prefix,
            "-remove",
            "-protocol", ",".join(protocols),
            "-operation", ",".join(operations),
            "-nastring", ".",
            "-csvout"
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                log.error(f"Error annotating variants with ANNOVAR: {result.stderr}")
                return None
            
            # Expected output file
            output_file = f"{output_prefix}.{self.build_version}_multianno.csv"
            
            if os.path.exists(output_file):
                log.info(f"Successfully annotated variants with ANNOVAR: {output_file}")
                return output_file
            else:
                log.error(f"ANNOVAR output file not found: {output_file}")
                return None
        except Exception as e:
            log.error(f"Exception annotating variants with ANNOVAR: {e}")
            return None
    
    def integrate_annotations(self, variant_df: pd.DataFrame, 
                             annovar_file: str) -> pd.DataFrame:
        """
        Integrate ANNOVAR annotations into variant DataFrame.
        
        Args:
            variant_df: DataFrame of variants
            annovar_file: Path to ANNOVAR annotation file
            
        Returns:
            DataFrame with integrated annotations
        """
        if not os.path.exists(annovar_file):
            log.error(f"ANNOVAR annotation file not found: {annovar_file}")
            return variant_df
        
        try:
            # Load ANNOVAR annotations
            log.info(f"Loading ANNOVAR annotations: {annovar_file}")
            annovar_df = pd.read_csv(annovar_file)
            
            # Create position keys for matching
            variant_df['pos_key'] = variant_df['chrom'].astype(str) + ':' + variant_df['pos'].astype(str)
            annovar_df['pos_key'] = annovar_df['Chr'].astype(str) + ':' + annovar_df['Start'].astype(str)
            
            # Select key columns from ANNOVAR results
            annotation_columns = [
                'pos_key', 'Func.refGene', 'Gene.refGene', 'ExonicFunc.refGene',
                'AAChange.refGene', 'ExAC_ALL', 'gnomAD_exome_ALL',
                'CLNSIG', 'CLNDN', 'CLNREVSTAT'
            ]
            
            # Filter to only include columns that exist
            annotation_columns = [col for col in annotation_columns if col in annovar_df.columns]
            
            if len(annotation_columns) < 2:  # Need at least pos_key and one annotation
                log.warning("Not enough annotation columns found in ANNOVAR results")
                return variant_df
            
            # Create a mapping dictionary for faster lookups
            annotation_map = {}
            for _, row in annovar_df[annotation_columns].iterrows():
                pos_key = row['pos_key']
                annotation_map[pos_key] = {col: row[col] for col in annotation_columns if col != 'pos_key'}
            
            # Add annotations to variant DataFrame
            for col in annotation_columns:
                if col != 'pos_key' and col not in variant_df.columns:
                    variant_df[col] = None
            
            # Update annotations
            for i, row in variant_df.iterrows():
                pos_key = row['pos_key']
                if pos_key in annotation_map:
                    for col, value in annotation_map[pos_key].items():
                        variant_df.at[i, col] = value
            
            # Remove temporary key column
            variant_df = variant_df.drop('pos_key', axis=1)
            
            log.info(f"Successfully integrated ANNOVAR annotations into variant DataFrame")
            return variant_df
        
        except Exception as e:
            log.error(f"Error integrating ANNOVAR annotations: {e}")
            return variant_df
    
    def process_variants(self, variant_df: pd.DataFrame, 
                        output_dir: str) -> pd.DataFrame:
        """
        Process variants through ANNOVAR annotation pipeline.
        
        Args:
            variant_df: DataFrame of variants
            output_dir: Output directory
            
        Returns:
            DataFrame with ANNOVAR annotations
        """
        if not self.is_available:
            log.warning("ANNOVAR not available, skipping annotation")
            return variant_df
        
        if variant_df.empty:
            log.warning("Empty variant DataFrame, nothing to annotate")
            return variant_df
        
        try:
            # Create temporary directory for ANNOVAR files
            annovar_dir = os.path.join(output_dir, "annovar_temp")
            os.makedirs(annovar_dir, exist_ok=True)
            
            # Create temporary ANNOVAR input file
            avinput_path = os.path.join(annovar_dir, "variants.avinput")
            
            # Convert DataFrame to ANNOVAR input format
            with open(avinput_path, 'w') as f:
                for _, row in variant_df.iterrows():
                    # Format: chr start end ref alt
                    chrom = row['chrom']
                    pos = row['pos']
                    ref = row['ref']
                    alt = row['alt']
                    
                    # ANNOVAR requires end position
                    end = pos + len(ref) - 1
                    
                    f.write(f"{chrom}\t{pos}\t{end}\t{ref}\t{alt}\n")
            
            # Run ANNOVAR annotation
            annotated_file = self.annotate_variants(
                input_file=avinput_path,
                output_dir=annovar_dir
            )
            
            if not annotated_file:
                log.error("ANNOVAR annotation failed")
                return variant_df
            
            # Integrate annotations
            enriched_df = self.integrate_annotations(variant_df, annotated_file)
            
            return enriched_df
        
        except Exception as e:
            log.error(f"Error in ANNOVAR processing pipeline: {e}")
            return variant_df

# Function to check if ANNOVAR is installed
def is_annovar_available() -> bool:
    """
    Check if ANNOVAR is available on the system.
    
    Returns:
        True if ANNOVAR is available, False otherwise
    """
    integration = AnnovarIntegration()
    return integration.is_available

# Function to enrich variants with ANNOVAR annotations
def enrich_variants_with_annovar(variant_df: pd.DataFrame, 
                               output_dir: str,
                               annovar_path: Optional[str] = None,
                               humandb_path: Optional[str] = None) -> pd.DataFrame:
    """
    Enrich variant DataFrame with ANNOVAR annotations.
    
    Args:
        variant_df: DataFrame of variants
        output_dir: Output directory
        annovar_path: Path to ANNOVAR installation directory
        humandb_path: Path to ANNOVAR humandb directory
        
    Returns:
        DataFrame with ANNOVAR annotations
    """
    integration = AnnovarIntegration(
        annovar_path=annovar_path,
        humandb_path=humandb_path
    )
    
    if not integration.is_available:
        log.warning("ANNOVAR not available, skipping annotation")
        return variant_df
    
    return integration.process_variants(variant_df, output_dir)
