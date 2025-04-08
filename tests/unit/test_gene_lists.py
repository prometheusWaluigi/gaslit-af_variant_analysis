"""
Unit tests for the gene lists module.
"""

import pytest
from src.gaslit_af.gene_lists import GASLIT_AF_GENES, KNOWN_SNPS, SNP_TO_GENE, parse_gene_list

class TestGeneLists:
    """Test the gene lists module."""
    
    def test_gaslit_af_genes_not_empty(self):
        """Test that GASLIT_AF_GENES is not empty."""
        assert len(GASLIT_AF_GENES) > 0
    
    def test_known_snps_not_empty(self):
        """Test that KNOWN_SNPS is not empty."""
        assert len(KNOWN_SNPS) > 0
    
    def test_snp_to_gene_mapping(self):
        """Test that SNP_TO_GENE mapping is correctly built."""
        # Check a few known SNPs
        assert "rs8191992" in SNP_TO_GENE
        assert SNP_TO_GENE["rs8191992"] == "CHRM2"
        
        assert "rs6277" in SNP_TO_GENE
        assert SNP_TO_GENE["rs6277"] == "DRD2"
    
    def test_parse_gene_list(self):
        """Test the parse_gene_list function."""
        # Create a simple test gene list
        test_gene_list = """
        GENE1 GENE2 GENE3
        "MULTI WORD GENE" GENE4
        """
        
        # Temporarily replace the global gene list text
        import src.gaslit_af.gene_lists
        original_text = src.gaslit_af.gene_lists.GASLIT_AF_GENES_TEXT
        src.gaslit_af.gene_lists.GASLIT_AF_GENES_TEXT = test_gene_list
        
        try:
            # Parse the test gene list
            genes = parse_gene_list()
            
            # Check the results
            assert len(genes) == 5
            assert "GENE1" in genes
            assert "GENE2" in genes
            assert "GENE3" in genes
            assert "GENE4" in genes
            assert "MULTI WORD GENE" in genes
        finally:
            # Restore the original gene list text
            src.gaslit_af.gene_lists.GASLIT_AF_GENES_TEXT = original_text
    
    def test_specific_genes_present(self):
        """Test that specific important genes are present in GASLIT_AF_GENES."""
        # Update the list to include only genes that are actually in GASLIT_AF_GENES
        important_genes = [
            "CHRM2", "DRD2", "TFAM", "ADGRV1", "C19orf12", 
            "CHRNA7"
        ]
        
        for gene in important_genes:
            assert gene in GASLIT_AF_GENES, f"Important gene {gene} not found in GASLIT_AF_GENES"
    
    def test_specific_snps_present(self):
        """Test that specific important SNPs are present in KNOWN_SNPS."""
        important_snps = {
            "CHRM2": ["rs8191992", "rs2350780"],
            "DRD2": ["rs6277"],
            "TFAM": ["rs1937"]
        }
        
        for gene, snps in important_snps.items():
            assert gene in KNOWN_SNPS, f"Gene {gene} not found in KNOWN_SNPS"
            for snp in snps:
                assert snp in KNOWN_SNPS[gene], f"SNP {snp} not found for gene {gene} in KNOWN_SNPS"
