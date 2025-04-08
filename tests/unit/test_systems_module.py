"""
Unit tests for the Systems module in GASLIT-AF Variant Analysis.

Tests the modular gene systems, variant store, and enrichment patterns.
"""

import os
import json
import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import shutil

from src.gaslit_af.systems.gene_systems import GeneSystemManager, get_gene_system_manager
from src.gaslit_af.systems.variant_store_simple import VariantStore, get_variant_store
from src.gaslit_af.systems.enrichment_patterns import (
    EnrichmentPattern,
    CompositeEnrichmentPattern,
    PathogenicityEnrichmentPattern,
    SystemSpecificEnrichmentPattern,
    AFEnrichmentPattern,
    MECFSEnrichmentPattern,
    create_enrichment_pattern,
    create_standard_enrichment_patterns
)


class TestGeneSystemManager:
    """Tests for the GeneSystemManager class."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for testing."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def gene_system_manager(self, temp_dir):
        """Create a GeneSystemManager instance for testing."""
        systems_dir = Path(temp_dir) / "systems"
        return GeneSystemManager(systems_dir)
    
    def test_default_systems_creation(self, gene_system_manager):
        """Test that default systems are created correctly."""
        # Check that systems were created
        assert len(gene_system_manager.systems) > 0
        
        # Check that ME/CFS system exists
        assert "mecfs_postviral" in gene_system_manager.systems
        
        # Check that ME/CFS genes are included
        mecfs_genes = gene_system_manager.get_genes_for_system("mecfs_postviral")
        for gene in ["S100PBP", "AKAP1", "USP6NL", "CDON", "SULF2"]:
            assert gene in mecfs_genes
    
    def test_get_system_for_gene(self, gene_system_manager):
        """Test getting the system for a gene."""
        # Test ME/CFS genes
        assert gene_system_manager.get_system_for_gene("S100PBP") == "mecfs_postviral"
        assert gene_system_manager.get_system_for_gene("AKAP1") == "mecfs_postviral"
        
        # Test AF genes
        assert gene_system_manager.get_system_for_gene("PITX2") == "cardiac_development"
        
        # Test unknown gene
        assert gene_system_manager.get_system_for_gene("UNKNOWN_GENE") == "unknown"
    
    def test_add_gene_to_system(self, gene_system_manager):
        """Test adding a gene to a system."""
        # Add a new gene to the ME/CFS system
        result = gene_system_manager.add_gene_to_system("NEW_GENE", "mecfs_postviral")
        assert result is True
        
        # Check that the gene was added
        mecfs_genes = gene_system_manager.get_genes_for_system("mecfs_postviral")
        assert "NEW_GENE" in mecfs_genes
        
        # Check that the gene is mapped to the correct system
        assert gene_system_manager.get_system_for_gene("NEW_GENE") == "mecfs_postviral"


class TestVariantStore:
    """Tests for the VariantStore class."""
    
    @pytest.fixture
    def temp_dir(self):
        """Create a temporary directory for testing."""
        temp_dir = tempfile.mkdtemp()
        yield temp_dir
        shutil.rmtree(temp_dir)
    
    @pytest.fixture
    def variant_store(self, temp_dir):
        """Create a VariantStore instance for testing."""
        db_path = Path(temp_dir) / "variants.duckdb"
        systems_dir = Path(temp_dir) / "systems"
        
        # Create gene systems first
        gene_system_manager = GeneSystemManager(systems_dir)
        
        # Create variant store
        store = VariantStore(db_path, systems_dir)
        
        yield store
        
        # Clean up
        store.close()
    
    def test_store_variant(self, variant_store):
        """Test storing a variant."""
        # Create a test variant
        variant_data = {
            "variant_id": "test_variant",
            "chrom": "1",
            "pos": 1000,
            "ref": "A",
            "alt": "G",
            "gene": "S100PBP",  # ME/CFS gene
            "rsid": "rs123456",
            "clinvar_significance": "pathogenic",
            "af": 0.001
        }
        
        # Store the variant
        result = variant_store.store_variant(variant_data)
        assert result is True
        
        # Retrieve the variant
        stored_variant = variant_store.get_variant("test_variant")
        assert stored_variant is not None
        assert stored_variant["variant_id"] == "test_variant"
        assert stored_variant["gene"] == "S100PBP"
        assert stored_variant["system_id"] == "mecfs_postviral"
    
    def test_get_variants_by_system(self, variant_store):
        """Test getting variants by system."""
        # Store some test variants
        for i, gene in enumerate([
            "S100PBP",  # ME/CFS
            "AKAP1",    # ME/CFS
            "PITX2",    # AF
            "KCNA5"     # AF
        ]):
            variant_data = {
                "variant_id": f"variant_{i}",
                "chrom": "1",
                "pos": 1000 + i,
                "ref": "A",
                "alt": "G",
                "gene": gene,
                "rsid": f"rs{i}",
                "clinvar_significance": "pathogenic" if i % 2 == 0 else "benign",
                "af": 0.001
            }
            variant_store.store_variant(variant_data)
        
        # Get ME/CFS variants
        mecfs_variants = variant_store.get_variants_by_system("mecfs_postviral")
        assert len(mecfs_variants) == 2
        assert set(mecfs_variants["gene"].tolist()) == {"S100PBP", "AKAP1"}
    
        # Get AF variants (PITX2 in cardiac_development, KCNA5 in calcium_ion_channels)
        pitx2_variants = variant_store.get_variants_by_system("cardiac_development")
        kcna5_variants = variant_store.get_variants_by_system("calcium_ion_channels")
        assert len(pitx2_variants) == 1
        assert pitx2_variants["gene"].iloc[0] == "PITX2"
        assert len(kcna5_variants) == 1
        assert kcna5_variants["gene"].iloc[0] == "KCNA5"
    
    def test_get_pathogenic_variants(self, variant_store):
        """Test getting pathogenic variants."""
        # Store some test variants
        for i, (gene, pathogenic) in enumerate([
            ("S100PBP", True),   # ME/CFS, pathogenic
            ("AKAP1", False),    # ME/CFS, not pathogenic
            ("PITX2", True),     # AF, pathogenic
            ("KCNA5", False)     # AF, not pathogenic
        ]):
            variant_data = {
                "variant_id": f"variant_{i}",
                "chrom": "1",
                "pos": 1000 + i,
                "ref": "A",
                "alt": "G",
                "gene": gene,
                "rsid": f"rs{i}",
                "clinvar_significance": "pathogenic" if pathogenic else "benign",
                "pathogenic": pathogenic,
                "af": 0.001
            }
            variant_store.store_variant(variant_data)
        
        # Get all pathogenic variants
        pathogenic_variants = variant_store.get_pathogenic_variants()
        assert len(pathogenic_variants) == 2
        assert set(pathogenic_variants["gene"].tolist()) == {"S100PBP", "PITX2"}
        
        # Get ME/CFS pathogenic variants
        mecfs_pathogenic = variant_store.get_pathogenic_variants("mecfs_postviral")
        assert len(mecfs_pathogenic) == 1
        assert mecfs_pathogenic["gene"].iloc[0] == "S100PBP"


class TestEnrichmentPatterns:
    """Tests for the enrichment patterns."""
    
    def test_pathogenicity_enrichment(self):
        """Test pathogenicity enrichment pattern."""
        # Create the pattern
        pattern = PathogenicityEnrichmentPattern("Test Pathogenicity", "Test description")
        
        # Test pathogenic variant (ClinVar)
        variant_data = {
            "variant_id": "test_variant",
            "gene": "S100PBP",
            "clinvar_significance": "pathogenic"
        }
        enriched = pattern.enrich(variant_data)
        assert enriched["pathogenic"] is True
        assert "ClinVar" in enriched["pathogenic_reason"]
        
        # Test pathogenic variant (CADD)
        variant_data = {
            "variant_id": "test_variant",
            "gene": "S100PBP",
            "clinvar_significance": "uncertain",
            "cadd_phred": 25.0
        }
        enriched = pattern.enrich(variant_data)
        assert enriched["pathogenic"] is True
        assert "CADD" in enriched["pathogenic_reason"]
        
        # Test non-pathogenic variant
        variant_data = {
            "variant_id": "test_variant",
            "gene": "S100PBP",
            "clinvar_significance": "benign",
            "cadd_phred": 5.0
        }
        enriched = pattern.enrich(variant_data)
        assert enriched["pathogenic"] is False
    
    def test_mecfs_enrichment(self):
        """Test ME/CFS enrichment pattern."""
        # Create the pattern
        pattern = MECFSEnrichmentPattern()
        
        # Test ME/CFS-related variant
        variant_data = {
            "variant_id": "test_variant",
            "gene": "S100PBP",
            "impact": "HIGH",
            "consequence": "missense_variant"
        }
        enriched = pattern.enrich(variant_data)
        assert enriched["mecfs_related"] is True
        assert enriched["mecfs_category"] == "mitochondrial"
        assert enriched["mecfs_pathogenic"] is True
        
        # Test viral susceptibility
        variant_data = {
            "variant_id": "test_variant",
            "gene": "USP6NL",
            "impact": "MODERATE",
            "consequence": "missense_variant"
        }
        enriched = pattern.enrich(variant_data)
        assert enriched["mecfs_related"] is True
        assert enriched["viral_susceptibility"] is True
        
        # Test non-ME/CFS variant
        variant_data = {
            "variant_id": "test_variant",
            "gene": "PITX2"
        }
        enriched = pattern.enrich(variant_data)
        assert enriched["mecfs_related"] is False
    
    def test_composite_enrichment(self):
        """Test composite enrichment pattern."""
        # Create component patterns
        pathogenicity = PathogenicityEnrichmentPattern("Pathogenicity", "Test")
        mecfs = MECFSEnrichmentPattern()
        af = AFEnrichmentPattern()
        
        # Create composite pattern
        composite = CompositeEnrichmentPattern(
            "Composite", 
            "Test composite",
            [pathogenicity, mecfs, af]
        )
        
        # Test ME/CFS variant
        variant_data = {
            "variant_id": "test_variant",
            "gene": "S100PBP",
            "clinvar_significance": "pathogenic"
        }
        enriched = composite.enrich(variant_data)
        assert enriched["pathogenic"] is True
        assert enriched["mecfs_related"] is True
        assert enriched["af_related"] is False
        
        # Test AF variant
        variant_data = {
            "variant_id": "test_variant",
            "gene": "PITX2",
            "clinvar_significance": "pathogenic"
        }
        enriched = composite.enrich(variant_data)
        assert enriched["pathogenic"] is True
        assert enriched["mecfs_related"] is False
        assert enriched["af_related"] is True


def test_standard_patterns():
    """Test creating standard enrichment patterns."""
    patterns = create_standard_enrichment_patterns()
    
    # Check that all expected patterns exist
    assert "pathogenicity" in patterns
    assert "af" in patterns
    assert "mecfs" in patterns
    assert "comprehensive" in patterns
    
    # Check that the comprehensive pattern includes all components
    comprehensive = patterns["comprehensive"]
    assert isinstance(comprehensive, CompositeEnrichmentPattern)
    assert len(comprehensive.patterns) == 3
