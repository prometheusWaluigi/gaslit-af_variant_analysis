"""
Systems Module for GASLIT-AF Variant Analysis.

This module provides a modular, lightweight architecture for variant analysis
across different biological systems, with a focus on recursive pattern recognition
and quantum coherence across systems.
"""

from .gene_systems import GeneSystemManager, get_gene_system_manager
from .variant_store_simple import VariantStore, get_variant_store
from .enrichment_patterns import (
    EnrichmentPattern, 
    CompositeEnrichmentPattern,
    PathogenicityEnrichmentPattern,
    SystemSpecificEnrichmentPattern,
    AFEnrichmentPattern,
    MECFSEnrichmentPattern,
    create_enrichment_pattern,
    create_standard_enrichment_patterns
)

__all__ = [
    'GeneSystemManager',
    'get_gene_system_manager',
    'VariantStore',
    'get_variant_store',
    'EnrichmentPattern',
    'CompositeEnrichmentPattern',
    'PathogenicityEnrichmentPattern',
    'SystemSpecificEnrichmentPattern',
    'AFEnrichmentPattern',
    'MECFSEnrichmentPattern',
    'create_enrichment_pattern',
    'create_standard_enrichment_patterns'
]
