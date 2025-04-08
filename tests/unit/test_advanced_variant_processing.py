"""
Unit tests for the advanced variant processing module.
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock
from src.gaslit_af.advanced_variant_processing import VariantProcessor, SNP_TO_GENE, KNOWN_SNPS

class TestAdvancedVariantProcessing:
    """Test the advanced variant processing module."""
    
    def test_variant_processor_initialization(self):
        """Test that VariantProcessor initializes correctly."""
        processor = VariantProcessor(threads=8)
        assert processor.threads == 8
        assert processor.queue is None
        assert isinstance(processor.gene_variants, dict)
        assert isinstance(processor.snp_variants, dict)
        assert isinstance(processor.rsid_map, dict)
    
    def test_known_snps_mapping(self):
        """Test that KNOWN_SNPS contains expected gene-SNP mappings."""
        # Check a few important genes and SNPs
        assert "CHRM2" in KNOWN_SNPS
        assert "rs8191992" in KNOWN_SNPS["CHRM2"]
        assert "rs2350780" in KNOWN_SNPS["CHRM2"]
        
        assert "DRD2" in KNOWN_SNPS
        assert "rs6277" in KNOWN_SNPS["DRD2"]
        
        assert "TFAM" in KNOWN_SNPS
        assert "rs1937" in KNOWN_SNPS["TFAM"]
        
        assert "ADGRV1" in KNOWN_SNPS
        assert "rs575602255" in KNOWN_SNPS["ADGRV1"]
        assert "rs555466095" in KNOWN_SNPS["ADGRV1"]
    
    def test_snp_to_gene_mapping(self):
        """Test that SNP_TO_GENE contains correct mappings."""
        # Check a few important SNPs
        assert "rs8191992" in SNP_TO_GENE
        assert SNP_TO_GENE["rs8191992"] == "CHRM2"
        
        assert "rs6277" in SNP_TO_GENE
        assert SNP_TO_GENE["rs6277"] == "DRD2"
        
        assert "rs1937" in SNP_TO_GENE
        assert SNP_TO_GENE["rs1937"] == "TFAM"
    
    def test_load_rsid_map(self):
        """Test loading rsID map from dbSNP file."""
        processor = VariantProcessor(threads=8)
        
        # Mock a dbSNP file that doesn't exist
        # The warning is logged, not raised as a UserWarning
        processor.load_rsid_map("nonexistent_file.vcf")
        
        # Assert that rsid_map is still empty
        assert len(processor.rsid_map) == 0
    
    @patch('pysam.VariantFile')
    def test_process_record(self, mock_variant_file):
        """Test processing a single VCF record."""
        processor = VariantProcessor(threads=8)
        
        # Create a mock record
        mock_record = MagicMock()
        mock_record.id = "rs8191992"  # A known SNP for CHRM2
        mock_record.chrom = "chr7"
        mock_record.pos = 136553416
        mock_record.ref = "A"
        mock_record.alts = ("G",)
        mock_record.qual = 100
        
        # Mock sample data
        mock_sample = MagicMock()
        mock_sample.gt_alleles = ["A", "G"]
        mock_record.samples = {"sample1": mock_sample}
        
        # Process the record
        target_genes = {"CHRM2"}
        gene, variant = processor._process_record(mock_record, target_genes)
        
        # Check the results
        assert gene == "CHRM2"
        assert variant["chrom"] == "chr7"
        assert variant["pos"] == 136553416
        assert variant["ref"] == "A"
        assert variant["alt"] == "G"
        assert variant["gene"] == "CHRM2"
        assert variant["rsid"] == "rs8191992"
        assert variant["genotype"] == "A/G"
        assert variant["quality"] == 100
    
    def test_annotate_variants(self):
        """Test annotating variants with additional information."""
        processor = VariantProcessor(threads=8)
        
        # Create a test DataFrame with variants
        variants = pd.DataFrame([
            {"gene": "CHRM2", "rsid": "rs8191992", "genotype": "A/G", "chrom": "chr7", "pos": 136553416},
            {"gene": "DRD2", "rsid": "rs6277", "genotype": "C/T", "chrom": "chr11", "pos": 113283459},
            {"gene": "UNKNOWN", "rsid": "rs12345", "genotype": "G/T", "chrom": "chr1", "pos": 12345}
        ])
        
        # Annotate the variants
        annotated = processor.annotate_variants(variants)
        
        # Check that annotations were added
        assert "impact" in annotated.columns
        assert "trait" in annotated.columns
        
        # Check specific annotations
        assert annotated.loc[0, "impact"] == "Executive function, memory, attention"
        assert annotated.loc[0, "trait"] == "Cognition & Brain Function"
        
        assert annotated.loc[1, "impact"] == "Dopamine modulation, cognitive flexibility"
        assert annotated.loc[1, "trait"] == "Cognition & Brain Function"
        
        # Unknown variants should have default annotations
        assert annotated.loc[2, "impact"] == "Unknown"
        assert annotated.loc[2, "trait"] == "Unknown"
    
    def test_generate_variant_report(self, output_dir):
        """Test generating a variant report."""
        processor = VariantProcessor(threads=8)
        
        # Create a test DataFrame with variants
        variants = pd.DataFrame([
            {"gene": "CHRM2", "rsid": "rs8191992", "genotype": "A/G", "chrom": "chr7", "pos": 136553416, 
             "impact": "Executive function, memory, attention", "trait": "Cognition & Brain Function"},
            {"gene": "DRD2", "rsid": "rs6277", "genotype": "C/T", "chrom": "chr11", "pos": 113283459,
             "impact": "Dopamine modulation, cognitive flexibility", "trait": "Cognition & Brain Function"},
            {"gene": "ADA", "rsid": "rs73598374", "genotype": "T/C", "chrom": "chr20", "pos": 44619661,
             "impact": "Deep sleep, longer delta wave cycles", "trait": "Sleep Traits"}
        ])
        
        # Generate the report
        report_path = processor.generate_variant_report(variants, output_dir)
        
        # Check that the report was created
        assert Path(report_path).exists()
        
        # Check the content of the report
        with open(report_path, 'r') as f:
            content = f.read()
            
            # Check that the report contains the expected sections
            assert "# Genomic Variant Analysis Report" in content
            assert "## Cognition & Brain Function" in content
            assert "## Sleep Traits" in content
            
            # Check that it contains the variant information
            assert "CHRM2" in content
            assert "rs8191992" in content
            assert "DRD2" in content
            assert "rs6277" in content
            assert "ADA" in content
            assert "rs73598374" in content
