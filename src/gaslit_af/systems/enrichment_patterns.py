"""
Enrichment Patterns Module for GASLIT-AF Variant Analysis.

This module provides a modular, composable approach to variant enrichment patterns.
Each pattern can be applied independently or recursively combined for deeper analysis.
"""

import os
import json
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Optional, Any, Union, Tuple, Callable

from .gene_systems import get_gene_system_manager
from .variant_store import get_variant_store

# Configure logging
log = logging.getLogger("gaslit-af")

class EnrichmentPattern:
    """Base class for variant enrichment patterns."""
    
    def __init__(self, name: str, description: str):
        """
        Initialize the enrichment pattern.
        
        Args:
            name: Pattern name
            description: Pattern description
        """
        self.name = name
        self.description = description
    
    def enrich(self, variant_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Enrich a variant with additional information.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Enriched variant data
        """
        # Base implementation does nothing
        return variant_data
    
    def enrich_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Enrich a DataFrame of variants.
        
        Args:
            df: DataFrame of variants
            
        Returns:
            Enriched DataFrame
        """
        # Create a copy to avoid modifying the original
        result = df.copy()
        
        # Apply enrichment to each row
        for idx, row in result.iterrows():
            variant_data = row.to_dict()
            enriched_data = self.enrich(variant_data)
            
            # Update row with enriched data
            for key, value in enriched_data.items():
                if key not in result.columns:
                    result[key] = None
                result.at[idx, key] = value
        
        return result


class CompositeEnrichmentPattern(EnrichmentPattern):
    """Composite pattern that combines multiple enrichment patterns."""
    
    def __init__(self, name: str, description: str, patterns: List[EnrichmentPattern] = None):
        """
        Initialize the composite enrichment pattern.
        
        Args:
            name: Pattern name
            description: Pattern description
            patterns: List of enrichment patterns to apply
        """
        super().__init__(name, description)
        self.patterns = patterns or []
    
    def add_pattern(self, pattern: EnrichmentPattern):
        """
        Add an enrichment pattern to the composite.
        
        Args:
            pattern: Enrichment pattern to add
        """
        self.patterns.append(pattern)
    
    def enrich(self, variant_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply all enrichment patterns sequentially.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Enriched variant data
        """
        result = variant_data.copy()
        
        # Apply each pattern in sequence
        for pattern in self.patterns:
            result = pattern.enrich(result)
        
        return result


class PathogenicityEnrichmentPattern(EnrichmentPattern):
    """Enrichment pattern for assessing variant pathogenicity."""
    
    def __init__(self, name: str, description: str, thresholds: Dict[str, float] = None):
        """
        Initialize the pathogenicity enrichment pattern.
        
        Args:
            name: Pattern name
            description: Pattern description
            thresholds: Dictionary of pathogenicity thresholds
        """
        super().__init__(name, description)
        
        # Default thresholds
        self.thresholds = thresholds or {
            "cadd_phred": 15.0,  # CADD Phred score threshold
            "sift_damaging": 0.05,  # SIFT score threshold (lower is more damaging)
            "polyphen_damaging": 0.85,  # PolyPhen score threshold (higher is more damaging)
            "gnomad_af_rare": 0.01  # gnomAD allele frequency threshold for rare variants
        }
    
    def enrich(self, variant_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Assess pathogenicity of a variant.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Enriched variant data with pathogenicity assessment
        """
        result = variant_data.copy()
        
        # Initialize pathogenicity fields
        result["pathogenic"] = False
        result["pathogenic_reason"] = ""
        
        # Check ClinVar significance
        clin_sig = str(result.get("clinvar_significance", "")).lower()
        if "pathogenic" in clin_sig or "likely_pathogenic" in clin_sig:
            result["pathogenic"] = True
            result["pathogenic_reason"] = "ClinVar pathogenic annotation"
            return result
        
        # Check pathogenicity scores
        cadd_phred = result.get("cadd_phred")
        if cadd_phred and float(cadd_phred) > self.thresholds["cadd_phred"]:
            result["pathogenic"] = True
            result["pathogenic_reason"] = f"High CADD score ({cadd_phred})"
            return result
        
        # Check SIFT score (lower is more damaging)
        sift = str(result.get("sift", "")).lower()
        if sift and ("deleterious" in sift or "damaging" in sift):
            result["pathogenic"] = True
            result["pathogenic_reason"] = f"Deleterious SIFT prediction ({sift})"
            return result
        
        # Check PolyPhen score (higher is more damaging)
        polyphen = str(result.get("polyphen", "")).lower()
        if polyphen and ("probably_damaging" in polyphen or "possibly_damaging" in polyphen):
            result["pathogenic"] = True
            result["pathogenic_reason"] = f"Damaging PolyPhen prediction ({polyphen})"
            return result
        
        # Check allele frequency (rare variants more likely to be pathogenic)
        af = result.get("af") or result.get("allele_frequency")
        if af and float(af) < self.thresholds["gnomad_af_rare"]:
            result["pathogenic"] = True
            result["pathogenic_reason"] = f"Rare variant (AF={af})"
            return result
        
        return result


class SystemSpecificEnrichmentPattern(EnrichmentPattern):
    """Enrichment pattern specific to a biological system."""
    
    def __init__(self, system_id: str, name: Optional[str] = None, description: Optional[str] = None):
        """
        Initialize the system-specific enrichment pattern.
        
        Args:
            system_id: System identifier
            name: Pattern name (defaults to system name)
            description: Pattern description
        """
        self.system_id = system_id
        self.gene_systems = get_gene_system_manager()
        
        # Get system info
        system_info = self.gene_systems.systems.get(system_id, {})
        system_name = system_info.get("name", system_id)
        
        # Use provided name/description or defaults
        name = name or f"{system_name} Enrichment"
        description = description or f"Enrichment pattern for {system_name} variants"
        
        super().__init__(name, description)
        
        # Get genes for this system
        self.system_genes = set(self.gene_systems.get_genes_for_system(system_id))
    
    def enrich(self, variant_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Enrich a variant with system-specific information.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Enriched variant data
        """
        result = variant_data.copy()
        
        # Add system information
        gene = result.get("gene", "")
        result[f"{self.system_id}_related"] = gene in self.system_genes
        
        if gene in self.system_genes:
            result["system_id"] = self.system_id
            
            # System-specific enrichment logic can be added here
            # This is a placeholder for system-specific logic
            pass
        
        return result


class AFEnrichmentPattern(SystemSpecificEnrichmentPattern):
    """Enrichment pattern specific to atrial fibrillation."""
    
    def __init__(self):
        """Initialize the AF enrichment pattern."""
        super().__init__("cardiac_development", "AF Enrichment", "Enrichment pattern for atrial fibrillation variants")
        
        # AF-specific gene categories
        self.af_categories = {
            "ion_channels": ["KCNA5", "KCND3", "KCNE1", "KCNQ1", "HCN4", "SCN5A", "KCNJ5", "RYR2"],
            "structural": ["MYL4", "LMNA", "MYH7", "FLNA"],
            "transcription_factors": ["PITX2", "GATA4", "GATA5", "GATA6", "TBX3", "TBX5", "NKX2-5", "ZFHX3"],
            "signaling": ["PRKAA2", "CAMK2B", "SPEN", "KIAA1755", "GREM2", "NPPA", "SH3PXD2A"],
            "inflammatory": ["IL6R", "IL6", "IL1B", "NLRP3"]
        }
        
        # Create reverse mapping
        self.gene_to_category = {}
        for category, genes in self.af_categories.items():
            for gene in genes:
                self.gene_to_category[gene] = category
    
    def enrich(self, variant_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Enrich a variant with AF-specific information.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Enriched variant data
        """
        # Apply base system enrichment
        result = super().enrich(variant_data)
        
        # Add AF-specific fields
        gene = result.get("gene", "")
        result["af_related"] = gene in self.system_genes
        result["af_category"] = self.gene_to_category.get(gene, "other")
        
        # AF-specific pathogenicity assessment
        if result.get("af_related", False):
            # Check if already assessed as pathogenic
            if result.get("pathogenic", False):
                result["af_pathogenic"] = True
                result["af_pathogenic_reason"] = result.get("pathogenic_reason", "General pathogenicity")
            else:
                # AF-specific pathogenicity logic
                af_specific_pathogenic = False
                reason = ""
                
                # Ion channel variants with any impact are more likely to be pathogenic for AF
                if result.get("af_category") == "ion_channels" and result.get("impact") in ["MODERATE", "HIGH"]:
                    af_specific_pathogenic = True
                    reason = "Ion channel variant with moderate/high impact"
                
                # Transcription factor variants affecting DNA binding
                elif result.get("af_category") == "transcription_factors" and "dna binding" in str(result.get("consequence", "")).lower():
                    af_specific_pathogenic = True
                    reason = "Transcription factor variant affecting DNA binding"
                
                # Structural protein variants affecting protein structure
                elif result.get("af_category") == "structural" and "missense" in str(result.get("consequence", "")).lower():
                    af_specific_pathogenic = True
                    reason = "Structural protein missense variant"
                
                result["af_pathogenic"] = af_specific_pathogenic
                result["af_pathogenic_reason"] = reason
        
        return result


class MECFSEnrichmentPattern(SystemSpecificEnrichmentPattern):
    """Enrichment pattern specific to ME/CFS and post-viral syndromes."""
    
    def __init__(self):
        """Initialize the ME/CFS enrichment pattern."""
        super().__init__("mecfs_postviral", "ME/CFS Enrichment", "Enrichment pattern for ME/CFS and post-viral syndrome variants")
        
        # ME/CFS-specific gene categories
        self.mecfs_categories = {
            "mitochondrial": ["AKAP1", "S100PBP"],
            "immune_regulation": ["USP6NL"],
            "cell_signaling": ["CDON", "SULF2"]
        }
        
        # Create reverse mapping
        self.gene_to_category = {}
        for category, genes in self.mecfs_categories.items():
            for gene in genes:
                self.gene_to_category[gene] = category
    
    def enrich(self, variant_data: Dict[str, Any]) -> Dict[str, Any]:
        """
        Enrich a variant with ME/CFS-specific information.
        
        Args:
            variant_data: Variant data dictionary
            
        Returns:
            Enriched variant data
        """
        # Apply base system enrichment
        result = super().enrich(variant_data)
        
        # Add ME/CFS-specific fields
        gene = result.get("gene", "")
        result["mecfs_related"] = gene in self.system_genes
        result["mecfs_category"] = self.gene_to_category.get(gene, "other")
        
        # ME/CFS-specific pathogenicity assessment
        if result.get("mecfs_related", False):
            # Check if already assessed as pathogenic
            if result.get("pathogenic", False):
                result["mecfs_pathogenic"] = True
                result["mecfs_pathogenic_reason"] = result.get("pathogenic_reason", "General pathogenicity")
            else:
                # ME/CFS-specific pathogenicity logic
                mecfs_specific_pathogenic = False
                reason = ""
                
                # Mitochondrial function variants are key for ME/CFS
                if result.get("mecfs_category") == "mitochondrial" and result.get("impact") in ["MODERATE", "HIGH"]:
                    mecfs_specific_pathogenic = True
                    reason = "Mitochondrial function variant with moderate/high impact"
                
                # Immune regulation variants affecting protein function
                elif result.get("mecfs_category") == "immune_regulation" and "missense" in str(result.get("consequence", "")).lower():
                    mecfs_specific_pathogenic = True
                    reason = "Immune regulation variant affecting protein function"
                
                # Cell signaling variants in conserved domains
                elif result.get("mecfs_category") == "cell_signaling" and result.get("conserved_domain", False):
                    mecfs_specific_pathogenic = True
                    reason = "Cell signaling variant in conserved domain"
                
                result["mecfs_pathogenic"] = mecfs_specific_pathogenic
                result["mecfs_pathogenic_reason"] = reason
                
                # Check for viral susceptibility markers
                if gene in ["USP6NL", "CDON"]:
                    result["viral_susceptibility"] = True
                    result["viral_susceptibility_note"] = f"{gene} associated with differential response to viral infection"
        
        return result


# Factory function to create enrichment patterns
def create_enrichment_pattern(pattern_type: str, **kwargs) -> EnrichmentPattern:
    """
    Create an enrichment pattern of the specified type.
    
    Args:
        pattern_type: Type of enrichment pattern
        **kwargs: Additional arguments for the pattern
        
    Returns:
        EnrichmentPattern instance
    """
    if pattern_type == "pathogenicity":
        return PathogenicityEnrichmentPattern(
            name=kwargs.get("name", "Pathogenicity Enrichment"),
            description=kwargs.get("description", "Assesses variant pathogenicity"),
            thresholds=kwargs.get("thresholds")
        )
    
    elif pattern_type == "system":
        return SystemSpecificEnrichmentPattern(
            system_id=kwargs.get("system_id", "unknown"),
            name=kwargs.get("name"),
            description=kwargs.get("description")
        )
    
    elif pattern_type == "af":
        return AFEnrichmentPattern()
    
    elif pattern_type == "mecfs":
        return MECFSEnrichmentPattern()
    
    elif pattern_type == "composite":
        patterns = kwargs.get("patterns", [])
        return CompositeEnrichmentPattern(
            name=kwargs.get("name", "Composite Enrichment"),
            description=kwargs.get("description", "Combines multiple enrichment patterns"),
            patterns=patterns
        )
    
    else:
        log.warning(f"Unknown enrichment pattern type: {pattern_type}")
        return EnrichmentPattern(
            name=kwargs.get("name", "Unknown Enrichment"),
            description=kwargs.get("description", "Unknown enrichment pattern")
        )


# Create standard enrichment patterns
def create_standard_enrichment_patterns() -> Dict[str, EnrichmentPattern]:
    """
    Create standard enrichment patterns.
    
    Returns:
        Dictionary of enrichment patterns
    """
    patterns = {}
    
    # Basic pathogenicity pattern
    patterns["pathogenicity"] = create_enrichment_pattern("pathogenicity")
    
    # System-specific patterns
    gene_systems = get_gene_system_manager()
    for system_id in gene_systems.systems.keys():
        patterns[system_id] = create_enrichment_pattern("system", system_id=system_id)
    
    # Disease-specific patterns
    patterns["af"] = create_enrichment_pattern("af")
    patterns["mecfs"] = create_enrichment_pattern("mecfs")
    
    # Composite patterns
    patterns["comprehensive"] = create_enrichment_pattern(
        "composite",
        name="Comprehensive Enrichment",
        description="Comprehensive variant enrichment with all patterns",
        patterns=[
            patterns["pathogenicity"],
            patterns["af"],
            patterns["mecfs"]
        ]
    )
    
    return patterns
