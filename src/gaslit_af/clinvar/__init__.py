"""
ClinVar Package for GASLIT-AF Variant Analysis.

This package is responsible for downloading, parsing, caching, and annotating
ClinVar data for use in the GASLIT-AF variant analysis pipeline.
"""

from src.gaslit_af.clinvar.downloader import ClinVarDownloader
from src.gaslit_af.clinvar.parser import ClinVarParser
from src.gaslit_af.clinvar.indexer import ClinVarIndexer
from src.gaslit_af.clinvar.cache_manager import ClinVarCache
from src.gaslit_af.clinvar.annotator import ClinVarAnnotator

# Main interface class
from src.gaslit_af.clinvar.annotator import ClinVarIntegration

__all__ = [
    'ClinVarDownloader',
    'ClinVarParser',
    'ClinVarIndexer',
    'ClinVarCache',
    'ClinVarAnnotator',
    'ClinVarIntegration',
]
