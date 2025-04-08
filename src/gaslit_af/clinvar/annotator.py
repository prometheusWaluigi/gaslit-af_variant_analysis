"""
ClinVar Annotator Module for GASLIT-AF Variant Analysis.

This module provides the main interface for annotating variants using ClinVar data.
"""

import logging
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List, Any, Union

from src.gaslit_af.clinvar.downloader import ClinVarDownloader
from src.gaslit_af.clinvar.parser import ClinVarParser
from src.gaslit_af.clinvar.indexer import ClinVarIndexer
from src.gaslit_af.clinvar.cache_manager import ClinVarCache

# Configure logging
log = logging.getLogger("gaslit-af")

class ClinVarAnnotator:
    """Annotates variants with ClinVar data."""
    
    def __init__(self, cache: ClinVarCache):
        """
        Initialize the ClinVar annotator.
        
        Args:
            cache: ClinVarCache instance
        """
        self.cache = cache
    
    def annotate_variant(self, variant: Dict) -> Dict:
        """
        Annotate a single variant with ClinVar information.
        
        Args:
            variant: Dictionary containing variant information
            
        Returns:
            Dictionary with variant annotated with ClinVar data
        """
        # Extract variant information
        chrom = variant.get('chrom') or variant.get('chromosome')
        pos = variant.get('pos') or variant.get('position')
        ref = variant.get('ref') or variant.get('reference')
        alt = variant.get('alt') or variant.get('alternate')
        rs_id = variant.get('rs_id') or variant.get('rsid')
        gene = variant.get('gene')
        
        # Look up variant in ClinVar
        clinvar_data = []
        
        # Try by rs_id first if available
        if rs_id:
            clinvar_data = self.cache.lookup_variant(rs_id=rs_id)
        
        # If not found by rs_id, try by genomic coordinates
        if not clinvar_data and chrom and pos:
            clinvar_data = self.cache.lookup_variant(chrom=chrom, pos=pos, ref=ref, alt=alt)
            
        # If still not found, try by gene
        if not clinvar_data and gene:
            clinvar_data = self.cache.lookup_variant(gene=gene)
        
        # Add ClinVar annotation if found
        if clinvar_data:
            clinvar_variant = clinvar_data[0]  # Use first match if multiple
            variant['clinvar_id'] = clinvar_variant.get('clinvar_id')
            variant['clinical_significance'] = clinvar_variant.get('significance')
            variant['is_pathogenic'] = 'pathogenic' in clinvar_variant.get('significance', '').lower()
            variant['is_benign'] = 'benign' in clinvar_variant.get('significance', '').lower()
            variant['is_vus'] = 'uncertain' in clinvar_variant.get('significance', '').lower()
        else:
            variant['clinvar_id'] = None
            variant['clinical_significance'] = 'Not found in ClinVar'
            variant['is_pathogenic'] = False
            variant['is_benign'] = False
            variant['is_vus'] = False
            
        return variant
    
    def annotate_variants_df(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate a DataFrame of variants with ClinVar information.
        
        Args:
            variants_df: DataFrame of variants
            
        Returns:
            Annotated DataFrame
        """
        if variants_df is None or variants_df.empty:
            return pd.DataFrame()
            
        # Convert DataFrame to list of dictionaries
        variants = variants_df.to_dict('records')
        
        # Annotate each variant
        annotated_variants = []
        for variant in variants:
            annotated_variants.append(self.annotate_variant(variant))
            
        # Convert back to DataFrame
        return pd.DataFrame(annotated_variants)


class ClinVarIntegration:
    """Main interface for ClinVar integration with variant analysis pipeline."""
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the ClinVar integration.
        
        Args:
            cache_dir: Directory to cache ClinVar data
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path("./cache/clinvar")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize components
        self.cache = ClinVarCache(self.cache_dir)
        self.downloader = self.cache.downloader
        self.parser = self.cache.parser
        self.indexer = self.cache.indexer
        self.annotator = ClinVarAnnotator(self.cache)
        
        # Cached data
        self._variant_summary = None
        self._vcf_grch37 = None
        self._vcf_grch38 = None
    
    def refresh_cache(self, force_download: bool = False) -> Dict:
        """
        Refresh all ClinVar data caches.
        
        Args:
            force_download: Whether to force download even if cache is valid
            
        Returns:
            Dictionary with cache statistics
        """
        # Refresh variant summary
        variant_summary_stats = self.cache.refresh_variant_summary(force_download)
        
        # Refresh VCF data for GRCh38
        vcf_38_stats = self.cache.refresh_vcf_data("GRCh38", force_download)
        
        # Refresh VCF data for GRCh37
        vcf_37_stats = self.cache.refresh_vcf_data("GRCh37", force_download)
        
        return self.cache.get_cache_stats()
    
    def get_cache_stats(self) -> Dict:
        """
        Get statistics about the cache.
        
        Returns:
            Dictionary with cache statistics
        """
        return self.cache.get_cache_stats()
    
    def clear_cache(self, cache_type: Optional[str] = None):
        """
        Clear the cache.
        
        Args:
            cache_type: Type of cache to clear (None for all)
        """
        self.cache.clear_cache(cache_type)
    
    def annotate_variants(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """
        Annotate variants with ClinVar information.
        
        Args:
            variants_df: DataFrame of variants
            
        Returns:
            Annotated DataFrame
        """
        return self.annotator.annotate_variants_df(variants_df)
    
    def get_pathogenic_variants(self, variants_df: pd.DataFrame) -> pd.DataFrame:
        """
        Filter variants to only include those that are pathogenic or likely pathogenic.
        
        Args:
            variants_df: DataFrame of variants
            
        Returns:
            DataFrame with only pathogenic variants
        """
        if variants_df is None or variants_df.empty:
            return pd.DataFrame()
        
        # Annotate variants first
        annotated_df = self.annotate_variants(variants_df)
        
        # Filter for pathogenic variants
        return annotated_df[annotated_df['is_pathogenic'] == True]
    
    def lookup_variant(self, **kwargs) -> List[Dict]:
        """
        Look up variants in the ClinVar index.
        
        Args:
            **kwargs: Query parameters (rs_id, clinvar_id, chrom, pos, ref, alt, gene)
            
        Returns:
            List of matching variants
        """
        return self.cache.lookup_variant(**kwargs)
    
    def get_variant_by_rsid(self, rs_id: str) -> Optional[Dict]:
        """
        Get a variant by its rsID.
        
        Args:
            rs_id: The rsID of the variant
            
        Returns:
            Dictionary with variant information, or None if not found
        """
        variants = self.cache.lookup_variant(rs_id=rs_id)
        return variants[0] if variants else None
    
    def get_variants_by_gene(self, gene: str) -> List[Dict]:
        """
        Get all variants associated with a gene.
        
        Args:
            gene: Gene symbol
            
        Returns:
            List of variants associated with the gene
        """
        return self.cache.lookup_variant(gene=gene)
