"""
Unit tests for the API integration module.
"""

import os
import json
import unittest
import tempfile
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add project root to Python path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.gaslit_af.api_integration import (
    VariantAnnotator, 
    EnsemblAnnotator, 
    MyVariantAnnotator, 
    VariantAPIIntegration
)


class TestVariantAnnotator(unittest.TestCase):
    """Test cases for the VariantAnnotator base class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Initialize the annotator
        self.annotator = VariantAnnotator(cache_dir=self.temp_path, cache_ttl=24)
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_get_cache_path(self):
        """Test getting the cache file path."""
        cache_path = self.annotator._get_cache_path("rs123456", "test_source")
        expected_path = self.temp_path / "test_source_rs123456.json"
        self.assertEqual(cache_path, expected_path)
        
        # Test with special characters
        cache_path = self.annotator._get_cache_path("1:12345:A/G", "test_source")
        expected_path = self.temp_path / "test_source_1_12345_A_G.json"
        self.assertEqual(cache_path, expected_path)
    
    def test_is_cache_valid(self):
        """Test checking if cache is valid."""
        # Create a cache file
        cache_path = self.temp_path / "test_cache.json"
        with open(cache_path, 'w') as f:
            json.dump({"test": "data"}, f)
        
        # Cache should be valid
        self.assertTrue(self.annotator._is_cache_valid(cache_path))
        
        # Non-existent cache should be invalid
        non_existent_path = self.temp_path / "non_existent.json"
        self.assertFalse(self.annotator._is_cache_valid(non_existent_path))
    
    def test_load_from_cache(self):
        """Test loading data from cache."""
        # Create a cache file
        variant_id = "rs123456"
        source = "test_source"
        cache_data = {"test": "data"}
        
        cache_path = self.annotator._get_cache_path(variant_id, source)
        with open(cache_path, 'w') as f:
            json.dump(cache_data, f)
        
        # Load from cache
        loaded_data = self.annotator._load_from_cache(variant_id, source)
        self.assertEqual(loaded_data, cache_data)
        
        # Non-existent cache should return None
        loaded_data = self.annotator._load_from_cache("non_existent", source)
        self.assertIsNone(loaded_data)
    
    def test_save_to_cache(self):
        """Test saving data to cache."""
        variant_id = "rs123456"
        source = "test_source"
        cache_data = {"test": "data"}
        
        # Save to cache
        self.annotator._save_to_cache(variant_id, source, cache_data)
        
        # Check that the file was created
        cache_path = self.annotator._get_cache_path(variant_id, source)
        self.assertTrue(cache_path.exists())
        
        # Check the content
        with open(cache_path, 'r') as f:
            loaded_data = json.load(f)
        
        self.assertEqual(loaded_data, cache_data)


class TestEnsemblAnnotator(unittest.TestCase):
    """Test cases for the EnsemblAnnotator class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Initialize the annotator
        self.annotator = EnsemblAnnotator(cache_dir=self.temp_path, cache_ttl=24)
        
        # Mock data
        self.variant_id = "rs429358"
        self.vep_data = [
            {
                "input": "rs429358",
                "assembly_name": "GRCh38",
                "seq_region_name": "19",
                "start": 44908684,
                "end": 44908684,
                "strand": 1,
                "allele_string": "T/C",
                "transcript_consequences": [
                    {
                        "transcript_id": "ENST00000252486",
                        "gene_id": "ENSG00000130203",
                        "gene_symbol": "APOE",
                        "consequence_terms": ["missense_variant"],
                        "impact": "MODERATE",
                        "biotype": "protein_coding"
                    }
                ]
            }
        ]
        
        self.variant_info = {
            "name": "rs429358",
            "mappings": [
                {
                    "seq_region_name": "19",
                    "start": 44908684,
                    "end": 44908684,
                    "strand": 1,
                    "allele_string": "T/C"
                }
            ],
            "source": "dbSNP",
            "MAF": 0.2
        }
        
        self.gene_info = {
            "id": "ENSG00000130203",
            "display_name": "APOE",
            "description": "apolipoprotein E",
            "biotype": "protein_coding",
            "strand": 1,
            "seq_region_name": "19",
            "start": 44905796,
            "end": 44909393
        }
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    @patch('requests.get')
    def test_get_variant_consequences(self, mock_get):
        """Test getting variant consequences from Ensembl."""
        # Mock response
        mock_response = MagicMock()
        mock_response.json.return_value = self.vep_data
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response
        
        # Get variant consequences
        data = self.annotator.get_variant_consequences(self.variant_id)
        
        # Check that the API was called
        mock_get.assert_called_once_with(
            f"{self.annotator.base_url}/vep/human/id/{self.variant_id}",
            headers=self.annotator.headers
        )
        
        # Check the returned data
        self.assertEqual(data, self.vep_data)
        
        # Check that the data was cached
        cache_path = self.annotator._get_cache_path(self.variant_id, "ensembl_vep")
        self.assertTrue(cache_path.exists())
    
    @patch('requests.get')
    def test_get_variant_info(self, mock_get):
        """Test getting variant information from Ensembl."""
        # Mock response
        mock_response = MagicMock()
        mock_response.json.return_value = self.variant_info
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response
        
        # Get variant info
        data = self.annotator.get_variant_info(self.variant_id)
        
        # Check that the API was called
        mock_get.assert_called_once_with(
            f"{self.annotator.base_url}/variation/human/{self.variant_id}",
            headers=self.annotator.headers
        )
        
        # Check the returned data
        self.assertEqual(data, self.variant_info)
        
        # Check that the data was cached
        cache_path = self.annotator._get_cache_path(self.variant_id, "ensembl_variant")
        self.assertTrue(cache_path.exists())
    
    @patch('requests.get')
    def test_get_gene_info(self, mock_get):
        """Test getting gene information from Ensembl."""
        # Mock response
        mock_response = MagicMock()
        mock_response.json.return_value = self.gene_info
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response
        
        # Get gene info
        gene_id = "ENSG00000130203"
        data = self.annotator.get_gene_info(gene_id)
        
        # Check that the API was called
        mock_get.assert_called_once_with(
            f"{self.annotator.base_url}/lookup/id/{gene_id}",
            headers=self.annotator.headers
        )
        
        # Check the returned data
        self.assertEqual(data, self.gene_info)
        
        # Check that the data was cached
        cache_path = self.annotator._get_cache_path(gene_id, "ensembl_gene")
        self.assertTrue(cache_path.exists())


class TestMyVariantAnnotator(unittest.TestCase):
    """Test cases for the MyVariantAnnotator class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Initialize the annotator
        self.annotator = MyVariantAnnotator(cache_dir=self.temp_path, cache_ttl=24)
        
        # Mock data
        self.variant_id = "rs429358"
        self.variant_data = {
            "_id": "rs429358",
            "cadd": {
                "phred": 0.007
            },
            "clinvar": {
                "rcv": {
                    "clinical_significance": "Pathogenic"
                }
            },
            "gnomad_genome": {
                "af": 0.164436
            },
            "dbnsfp": {
                "sift": {
                    "pred": "T"
                },
                "polyphen2": {
                    "hdiv": {
                        "pred": "B"
                    }
                }
            }
        }
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    @patch('requests.get')
    def test_get_variant_info(self, mock_get):
        """Test getting variant information from MyVariant.info."""
        # Mock response
        mock_response = MagicMock()
        mock_response.json.return_value = self.variant_data
        mock_response.raise_for_status.return_value = None
        mock_get.return_value = mock_response
        
        # Get variant info
        data = self.annotator.get_variant_info(self.variant_id)
        
        # Check that the API was called
        mock_get.assert_called_once_with(
            f"{self.annotator.base_url}/variant/{self.variant_id}"
        )
        
        # Check the returned data
        self.assertEqual(data, self.variant_data)
        
        # Check that the data was cached
        cache_path = self.annotator._get_cache_path(self.variant_id, "myvariant")
        self.assertTrue(cache_path.exists())
    
    @patch('requests.post')
    def test_get_variants_info(self, mock_post):
        """Test getting information for multiple variants from MyVariant.info."""
        # Mock response
        mock_response = MagicMock()
        mock_response.json.return_value = [self.variant_data]
        mock_response.raise_for_status.return_value = None
        mock_post.return_value = mock_response
        
        # Get variants info
        variant_ids = [self.variant_id]
        data = self.annotator.get_variants_info(variant_ids)
        
        # Check that the API was called
        mock_post.assert_called_once_with(
            f"{self.annotator.base_url}/variant",
            data={"ids": self.variant_id}
        )
        
        # Check the returned data
        self.assertEqual(data, {self.variant_id: self.variant_data})
        
        # Check that the data was cached
        cache_path = self.annotator._get_cache_path(self.variant_id, "myvariant")
        self.assertTrue(cache_path.exists())


class TestVariantAPIIntegration(unittest.TestCase):
    """Test cases for the VariantAPIIntegration class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Initialize the API integration
        self.api = VariantAPIIntegration(cache_dir=self.temp_path, cache_ttl=24)
        
        # Mock data
        self.variant_id = "rs429358"
        self.ensembl_vep_data = [
            {
                "input": "rs429358",
                "assembly_name": "GRCh38",
                "seq_region_name": "19",
                "start": 44908684,
                "end": 44908684,
                "strand": 1,
                "allele_string": "T/C",
                "transcript_consequences": [
                    {
                        "transcript_id": "ENST00000252486",
                        "gene_id": "ENSG00000130203",
                        "gene_symbol": "APOE",
                        "consequence_terms": ["missense_variant"],
                        "impact": "MODERATE",
                        "biotype": "protein_coding"
                    }
                ]
            }
        ]
        
        self.ensembl_variant_data = {
            "name": "rs429358",
            "mappings": [
                {
                    "seq_region_name": "19",
                    "start": 44908684,
                    "end": 44908684,
                    "strand": 1,
                    "allele_string": "T/C"
                }
            ],
            "source": "dbSNP",
            "MAF": 0.2
        }
        
        self.myvariant_data = {
            "_id": "rs429358",
            "cadd": {
                "phred": 0.007
            },
            "clinvar": {
                "rcv": {
                    "clinical_significance": "Pathogenic"
                }
            },
            "gnomad_genome": {
                "af": 0.164436
            },
            "dbnsfp": {
                "sift": {
                    "pred": "T"
                },
                "polyphen2": {
                    "hdiv": {
                        "pred": "B"
                    }
                }
            }
        }
        
        # Create test variants DataFrame
        self.variants_df = pd.DataFrame({
            "rsid": [self.variant_id, "rs7412"],
            "gene": ["APOE", "APOE"],
            "chrom": ["19", "19"],
            "pos": [44908684, 44908822],
            "ref": ["T", "C"],
            "alt": ["C", "T"]
        })
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    @patch.object(EnsemblAnnotator, 'get_variant_consequences')
    @patch.object(EnsemblAnnotator, 'get_variant_info')
    @patch.object(MyVariantAnnotator, 'get_variant_info')
    def test_get_variant_details(self, mock_myvariant, mock_ensembl_variant, mock_ensembl_vep):
        """Test getting comprehensive details for a variant."""
        # Mock responses
        mock_ensembl_vep.return_value = self.ensembl_vep_data
        mock_ensembl_variant.return_value = self.ensembl_variant_data
        mock_myvariant.return_value = self.myvariant_data
        
        # Get variant details
        details = self.api.get_variant_details(self.variant_id)
        
        # Check that the APIs were called
        mock_ensembl_vep.assert_called_once_with(self.variant_id)
        mock_ensembl_variant.assert_called_once_with(self.variant_id)
        mock_myvariant.assert_called_once_with(self.variant_id)
        
        # Check the returned data
        self.assertEqual(details["variant_id"], self.variant_id)
        self.assertEqual(details["gene"], "APOE")
        self.assertEqual(details["consequence"], "missense_variant")
        self.assertEqual(details["impact"], "MODERATE")
        self.assertEqual(details["clinical_significance"], "Pathogenic")
        self.assertEqual(details["allele_frequency"], 0.164436)
        self.assertEqual(details["pathogenicity_scores"]["cadd_phred"], 0.007)
        self.assertEqual(details["pathogenicity_scores"]["sift"], "T")
        self.assertEqual(details["pathogenicity_scores"]["polyphen"], "B")
    
    @patch.object(EnsemblAnnotator, 'get_variant_consequences')
    @patch.object(MyVariantAnnotator, 'get_variants_info')
    def test_annotate_variants(self, mock_myvariant, mock_ensembl_vep):
        """Test annotating variants with data from external APIs."""
        # Mock responses
        mock_ensembl_vep.return_value = self.ensembl_vep_data
        mock_myvariant.return_value = {self.variant_id: self.myvariant_data}
        
        # Annotate variants
        annotated_df = self.api.annotate_variants(self.variants_df, sources=["ensembl", "myvariant"])
        
        # Check that the APIs were called
        self.assertEqual(mock_ensembl_vep.call_count, 2)  # Once for each variant
        mock_myvariant.assert_called_once()
        
        # Check the annotations
        self.assertIn("ensembl_most_severe", annotated_df.columns)
        self.assertIn("ensembl_consequence", annotated_df.columns)
        self.assertIn("ensembl_impact", annotated_df.columns)
        self.assertIn("clinvar_significance", annotated_df.columns)
        self.assertIn("cadd_phred", annotated_df.columns)
        self.assertIn("gnomad_af", annotated_df.columns)
        
        # Check values for the first variant
        self.assertEqual(annotated_df.loc[0, "ensembl_most_severe"], "missense_variant")
        self.assertEqual(annotated_df.loc[0, "ensembl_consequence"], "missense_variant")
        self.assertEqual(annotated_df.loc[0, "ensembl_impact"], "MODERATE")
        self.assertEqual(annotated_df.loc[0, "clinvar_significance"], "Pathogenic")
        self.assertEqual(annotated_df.loc[0, "cadd_phred"], 0.007)
        self.assertEqual(annotated_df.loc[0, "gnomad_af"], 0.164436)


if __name__ == "__main__":
    unittest.main()
