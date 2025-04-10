"""
Variant Enrichment Module for GASLIT-AF Variant Analysis.

This module provides advanced variant enrichment capabilities, integrating
multiple data sources to enhance variant annotations with a focus on
atrial fibrillation and related cardiac pathways.
"""

import os
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Optional, Any, Union, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

# Import GASLIT-AF modules
from .gene_lists import GASLIT_AF_GENES, KNOWN_SNPS
from .api_integration import VariantAPIIntegration
from .clinvar_integration import ClinVarIntegration
from .caching import AnalysisCache

# Configure logging
log = logging.getLogger("gaslit-af")

# Define AF-specific pathogenicity thresholds
AF_PATHOGENICITY_THRESHOLDS = {
    "cadd_phred": 15.0,  # CADD Phred score threshold for AF variants
    "sift_damaging": 0.05,  # SIFT score threshold (lower is more damaging)
    "polyphen_damaging": 0.85,  # PolyPhen score threshold (higher is more damaging)
    "gnomad_af_rare": 0.01  # gnomAD allele frequency threshold for rare variants
}

# Define AF-specific gene categories
AF_GENE_CATEGORIES = {
    "Ion Channels": [
        "KCNA5", "KCND3", "KCNE1", "KCNQ1", "HCN4", "SCN5A", "KCNJ5", "RYR2"
    ],
    "Structural Proteins": [
        "MYL4", "LMNA", "MYH7", "FLNA"
    ],
    "Transcription Factors": [
        "PITX2", "GATA4", "GATA5", "GATA6", "TBX3", "TBX5", "NKX2-5", "ZFHX3"
    ],
    "Signaling Molecules": [
        "PRKAA2", "CAMK2B", "SPEN", "KIAA1755", "GREM2", "NPPA", "SH3PXD2A"
    ],
    "Inflammatory Mediators": [
        "IL6R", "IL6", "IL1B", "NLRP3"
    ]
}

# Create reverse mapping from gene to category
AF_GENE_TO_CATEGORY = {}
for category, genes in AF_GENE_CATEGORIES.items():
    for gene in genes:
        AF_GENE_TO_CATEGORY[gene] = category


class VariantEnricher:
    """Enriches variant data with AF-specific annotations and insights."""
    
    def __init__(self, cache_dir: Optional[Path] = None, cache_ttl: int = 24):
        """
        Initialize the variant enricher.
        
        Args:
            cache_dir: Directory to cache enrichment data
            cache_ttl: Cache time-to-live in hours
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path("./cache/enrichment")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.cache_ttl = cache_ttl
        
        # Initialize API and ClinVar integrations
        self.api = VariantAPIIntegration(cache_dir=self.cache_dir, cache_ttl=cache_ttl)
        self.clinvar = ClinVarIntegration(cache_dir=self.cache_dir)
        
        # Initialize cache manager
        self.cache = AnalysisCache(cache_dir=self.cache_dir, max_age_hours=cache_ttl, enabled=True)
        
        # Track processed variants
        self.processed_variants = set()
        self.af_related_variants = set()
        
    def get_af_category(self, gene: str) -> str:
        """
        Get the AF-specific category for a gene.
        
        Args:
            gene: Gene symbol
            
        Returns:
            Category name or "Other" if not in AF categories
        """
        return AF_GENE_TO_CATEGORY.get(gene, "Other")
    
    def is_af_related_gene(self, gene: str) -> bool:
        """
        Check if a gene is related to atrial fibrillation.
        
        Args:
            gene: Gene symbol
            
        Returns:
            True if gene is AF-related, False otherwise
        """
        return gene in AF_GENE_TO_CATEGORY
    
    def is_pathogenic_for_af(self, variant_data: Dict) -> Tuple[bool, str]:
        """
        Determine if a variant is potentially pathogenic for AF.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Tuple of (is_pathogenic, reason)
        """
        # Check if gene is AF-related
        gene = variant_data.get("gene")
        if not gene or not self.is_af_related_gene(gene):
            return False, "Not in AF-related gene"
        
        # Check ClinVar significance
        clin_sig = variant_data.get("clinical_significance", "").lower()
        if "pathogenic" in clin_sig or "likely_pathogenic" in clin_sig:
            return True, "ClinVar pathogenic annotation"
        
        # Check pathogenicity scores
        scores = variant_data.get("pathogenicity_scores", {})
        
        # CADD score check
        cadd_phred = scores.get("cadd_phred")
        if cadd_phred and cadd_phred > AF_PATHOGENICITY_THRESHOLDS["cadd_phred"]:
            return True, f"High CADD score ({cadd_phred})"
        
        # SIFT score check (lower is more damaging)
        sift = scores.get("sift", "").lower()
        if sift and ("deleterious" in sift or "damaging" in sift):
            return True, f"Deleterious SIFT prediction ({sift})"
        
        # PolyPhen score check (higher is more damaging)
        polyphen = scores.get("polyphen", "").lower()
        if polyphen and ("probably_damaging" in polyphen or "possibly_damaging" in polyphen):
            return True, f"Damaging PolyPhen prediction ({polyphen})"
        
        # Check allele frequency (rare variants more likely to be pathogenic)
        af = variant_data.get("allele_frequency")
        if af and af < AF_PATHOGENICITY_THRESHOLDS["gnomad_af_rare"]:
            return True, f"Rare variant (AF={af})"
        
        return False, "No pathogenic evidence found"
    
    def enrich_variant(self, variant_id: str) -> Dict:
        """
        Enrich a single variant with AF-specific annotations.
        
        Args:
            variant_id: Variant identifier (e.g., rs123456)
            
        Returns:
            Dict: Enriched variant data
        """
        # Check cache first
        # Use variant_id as the vcf_path parameter (not actually a path but a unique identifier)
        # and 'variant_enrichment' as the analysis_type
        cached_data = self.cache.get(variant_id, "variant_enrichment")
        if cached_data:
            log.debug(f"Using cached enrichment for {variant_id}")
            return cached_data
        
        # Get basic variant details
        variant_data = self.api.get_variant_details(variant_id)
        
        # Add ClinVar details if available
        clinvar_details = self.clinvar.get_variant_details(variant_id)
        if "error" not in clinvar_details:
            variant_data["clinvar"] = clinvar_details
        
        # Add AF-specific annotations
        gene = variant_data.get("gene")
        if gene:
            variant_data["af_category"] = self.get_af_category(gene)
            variant_data["is_af_related"] = self.is_af_related_gene(gene)
            
            # Assess pathogenicity for AF
            is_pathogenic, reason = self.is_pathogenic_for_af(variant_data)
            variant_data["af_pathogenic"] = is_pathogenic
            variant_data["af_pathogenic_reason"] = reason
            
            # Track AF-related variants
            if variant_data["is_af_related"]:
                self.af_related_variants.add(variant_id)
        
        # Cache the enriched data
        self.cache.set(variant_data, variant_id, "variant_enrichment")
        self.processed_variants.add(variant_id)
        
        return variant_data
    
    def enrich_variants(self, variant_df: pd.DataFrame) -> pd.DataFrame:
        """
        Enrich a DataFrame of variants with AF-specific annotations.
        
        Args:
            variant_df: DataFrame containing variant data
            
        Returns:
            Enriched DataFrame
        """
        log.info(f"Enriching {len(variant_df)} variants with AF-specific annotations")
        
        # Check if rsid column exists
        if "rsid" not in variant_df.columns:
            log.warning("No rsid column found in variant DataFrame, cannot enrich")
            return variant_df
        
        # Create output columns if they don't exist
        for col in ["af_category", "is_af_related", "af_pathogenic", "af_pathogenic_reason"]:
            if col not in variant_df.columns:
                variant_df[col] = None
        
        # Process variants in parallel
        enriched_data = {}
        variant_ids = [vid for vid in variant_df["rsid"].unique() if vid and str(vid).startswith("rs")]
        
        with ThreadPoolExecutor(max_workers=10) as executor:
            future_to_variant = {
                executor.submit(self.enrich_variant, vid): vid 
                for vid in variant_ids
            }
            
            for future in as_completed(future_to_variant):
                variant_id = future_to_variant[future]
                try:
                    enriched_data[variant_id] = future.result()
                except Exception as e:
                    log.error(f"Error enriching variant {variant_id}: {e}")
        
        # Update DataFrame with enriched data
        for idx, row in variant_df.iterrows():
            rsid = row.get("rsid")
            if rsid in enriched_data:
                data = enriched_data[rsid]
                variant_df.at[idx, "af_category"] = data.get("af_category", "Other")
                variant_df.at[idx, "is_af_related"] = data.get("is_af_related", False)
                variant_df.at[idx, "af_pathogenic"] = data.get("af_pathogenic", False)
                variant_df.at[idx, "af_pathogenic_reason"] = data.get("af_pathogenic_reason", "")
        
        log.info(f"Enriched {len(enriched_data)} variants, found {len(self.af_related_variants)} AF-related variants")
        return variant_df
    
    def generate_af_enrichment_report(self, variant_df: pd.DataFrame, output_dir: str) -> str:
        """
        Generate a report focusing on AF-related variants.
        
        Args:
            variant_df: DataFrame of enriched variants
            output_dir: Output directory
            
        Returns:
            Path to the generated report
        """
        if variant_df.empty:
            log.warning("No variants to report")
            return ""
        
        # Filter for AF-related variants
        af_variants = variant_df[variant_df["is_af_related"] == True].copy()
        
        if af_variants.empty:
            log.warning("No AF-related variants found")
            return ""
        
        log.info(f"Generating AF enrichment report for {len(af_variants)} variants")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Group variants by category
        af_variants["af_category"] = af_variants["af_category"].fillna("Other")
        category_groups = af_variants.groupby("af_category")
        
        # Generate report
        report_path = os.path.join(output_dir, "af_enrichment_report.md")
        
        with open(report_path, 'w') as f:
            f.write("# Atrial Fibrillation Variant Enrichment Report\n\n")
            
            # Summary statistics
            f.write("## Summary\n\n")
            f.write(f"Total variants analyzed: {len(variant_df)}\n")
            f.write(f"AF-related variants identified: {len(af_variants)}\n")
            f.write(f"Potentially pathogenic AF variants: {len(af_variants[af_variants['af_pathogenic'] == True])}\n\n")
            
            # Category breakdown
            f.write("## Variant Distribution by AF Category\n\n")
            f.write("| Category | Variant Count | Pathogenic Count |\n")
            f.write("|----------|---------------|------------------|\n")
            
            for category, group in category_groups:
                pathogenic_count = len(group[group["af_pathogenic"] == True])
                f.write(f"| {category} | {len(group)} | {pathogenic_count} |\n")
            
            f.write("\n")
            
            # Detailed variant listing by category
            f.write("## Detailed Variant Analysis by Category\n\n")
            
            for category, group in category_groups:
                f.write(f"### {category}\n\n")
                f.write("| Gene | Variant | Pathogenic | Reason | Clinical Significance |\n")
                f.write("|------|---------|------------|--------|------------------------|\n")
                
                for _, row in group.iterrows():
                    gene = row.get("gene", "")
                    rsid = row.get("rsid", "")
                    pathogenic = "Yes" if row.get("af_pathogenic") else "No"
                    reason = row.get("af_pathogenic_reason", "")
                    clin_sig = row.get("clinvar_significance", "")
                    
                    f.write(f"| {gene} | {rsid} | {pathogenic} | {reason} | {clin_sig} |\n")
                
                f.write("\n")
            
            # Recommendations section
            f.write("## Clinical Implications and Recommendations\n\n")
            f.write("The variants identified in this report may have implications for atrial fibrillation risk ")
            f.write("and related cardiac conditions. Consider the following recommendations:\n\n")
            
            f.write("1. **Ion Channel Variants**: May affect cardiac conduction and repolarization. ")
            f.write("Consider ECG monitoring and antiarrhythmic medication response.\n\n")
            
            f.write("2. **Structural Protein Variants**: May affect cardiac remodeling. ")
            f.write("Consider echocardiographic assessment of atrial structure and function.\n\n")
            
            f.write("3. **Transcription Factor Variants**: May affect cardiac development. ")
            f.write("Consider family history assessment and genetic counseling.\n\n")
            
            f.write("4. **Signaling Molecule Variants**: May affect cellular signaling pathways. ")
            f.write("Consider response to rate control medications.\n\n")
            
            f.write("5. **Inflammatory Mediator Variants**: May affect inflammatory processes. ")
            f.write("Consider inflammatory biomarker assessment.\n\n")
        
        log.info(f"AF enrichment report generated: {report_path}")
        return report_path


def enrich_variants_with_af_data(variant_df: pd.DataFrame, 
                               cache_dir: Optional[Path] = None,
                               output_dir: Optional[str] = None) -> Tuple[pd.DataFrame, Optional[str]]:
    """
    Enrich variants with atrial fibrillation-specific annotations.
    
    Args:
        variant_df: DataFrame containing variant data
        cache_dir: Directory to cache enrichment data
        output_dir: Directory to save the enrichment report
        
    Returns:
        Tuple of (enriched_df, report_path)
    """
    enricher = VariantEnricher(cache_dir=cache_dir)
    
    # Enrich variants
    enriched_df = enricher.enrich_variants(variant_df)
    
    # Generate report if output directory is provided
    report_path = None
    if output_dir and not enriched_df.empty:
        report_path = enricher.generate_af_enrichment_report(enriched_df, output_dir)
    
    return enriched_df, report_path
