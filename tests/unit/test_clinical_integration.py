"""
Unit tests for the clinical integration module.
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

from src.gaslit_af.clinical_integration import ClinicalIntegration
from src.gaslit_af.clinical_variants import ClinicalVariantManager


class TestClinicalIntegration(unittest.TestCase):
    """Test cases for the ClinicalIntegration class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Create a test clinical data file
        self.clinical_data_path = self.temp_path / "test_clinical_data.json"
        self.clinical_data = {
            "conditions": [
                {
                    "name": "Test Condition 1",
                    "description": "A test condition for unit testing",
                    "symptoms": ["Symptom 1", "Symptom 2"],
                    "status": {
                        "risk": "Pathogenic",
                        "confidence": "High",
                        "classification": "D"
                    },
                    "genetic_data": {
                        "gene": "APOE",
                        "variant_id": "rs429358",
                        "rcv": "RCV000123456",
                        "genotype": "T/C"
                    },
                    "risk_assessment": {
                        "frequency": 0.01,
                        "version": "GRCh38"
                    }
                },
                {
                    "name": "Test Condition 2",
                    "description": "Another test condition",
                    "symptoms": ["Symptom 3", "Symptom 4"],
                    "status": {
                        "risk": "Benign",
                        "confidence": "Medium",
                        "classification": "S"
                    },
                    "genetic_data": {
                        "gene": "APOE",
                        "variant_id": "rs7412",
                        "rcv": "RCV000654321",
                        "genotype": "C/T"
                    },
                    "risk_assessment": {
                        "frequency": 0.05,
                        "version": "GRCh38"
                    }
                },
                {
                    "name": "Test Condition 3",
                    "description": "A VUS condition",
                    "symptoms": ["Symptom 5", "Symptom 6"],
                    "status": {
                        "risk": "Uncertain Significance",
                        "confidence": "Low",
                        "classification": "Other"
                    },
                    "genetic_data": {
                        "gene": "BRCA1",
                        "variant_id": "rs80357906",
                        "rcv": "RCV000789012",
                        "genotype": "G/A"
                    },
                    "risk_assessment": {
                        "frequency": 0.001,
                        "version": "GRCh38"
                    }
                }
            ]
        }
        
        with open(self.clinical_data_path, 'w') as f:
            json.dump(self.clinical_data, f)
        
        # Create test variants DataFrame
        self.variants_df = pd.DataFrame({
            "rsid": ["rs429358", "rs7412", "rs80357906", "rs999999"],
            "gene": ["APOE", "APOE", "BRCA1", "TEST"],
            "chrom": ["19", "19", "17", "1"],
            "pos": [44908684, 44908822, 43071077, 1000],
            "ref": ["T", "C", "G", "A"],
            "alt": ["C", "T", "A", "G"],
            "genotype": ["T/C", "C/T", "G/A", "A/G"]
        })
        
        # Initialize the integration
        self.integration = ClinicalIntegration(self.clinical_data_path)
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_load_clinical_data(self):
        """Test loading clinical data."""
        # Data should be loaded in setUp
        self.assertTrue(self.integration.clinical_data_loaded)
        
        # Test loading non-existent file
        integration = ClinicalIntegration()
        result = integration.load_clinical_data(self.temp_path / "non_existent.json")
        self.assertFalse(result)
        self.assertFalse(integration.clinical_data_loaded)
    
    def test_annotate_variants(self):
        """Test annotating variants with clinical information."""
        # Annotate variants
        annotated_df = self.integration.annotate_variants(self.variants_df)
        
        # Check annotations
        self.assertIn("clinical_significance", annotated_df.columns)
        self.assertIn("condition_name", annotated_df.columns)
        self.assertIn("condition_description", annotated_df.columns)
        self.assertIn("rcv", annotated_df.columns)
        
        # Check values
        self.assertEqual(annotated_df.loc[0, "clinical_significance"], "Pathogenic")
        self.assertEqual(annotated_df.loc[0, "condition_name"], "Test Condition 1")
        self.assertEqual(annotated_df.loc[0, "rcv"], "RCV000123456")
        
        self.assertEqual(annotated_df.loc[1, "clinical_significance"], "Benign")
        self.assertEqual(annotated_df.loc[1, "condition_name"], "Test Condition 2")
        
        self.assertEqual(annotated_df.loc[2, "clinical_significance"], "Uncertain Significance")
        self.assertEqual(annotated_df.loc[2, "condition_name"], "Test Condition 3")
        
        self.assertIsNone(annotated_df.loc[3, "clinical_significance"])
        self.assertIsNone(annotated_df.loc[3, "condition_name"])
    
    def test_generate_clinical_report(self):
        """Test generating a clinical report."""
        # Generate report
        report_path = self.integration.generate_clinical_report(self.variants_df, self.temp_path)
        
        # Check that the report was created
        self.assertIsNotNone(report_path)
        report_file = Path(report_path)
        self.assertTrue(report_file.exists())
        
        # Check report content
        with open(report_file, 'r') as f:
            report_content = f.read()
        
        self.assertIn("Clinical Variant Report", report_content)
        self.assertIn("Pathogenic/Likely Pathogenic Variants", report_content)
        self.assertIn("APOE", report_content)
        self.assertIn("rs429358", report_content)
        self.assertIn("Test Condition 1", report_content)
    
    def test_get_pathogenic_variants(self):
        """Test getting pathogenic variants."""
        # Get pathogenic variants
        pathogenic_df = self.integration.get_pathogenic_variants(self.variants_df)
        
        # Check results
        self.assertEqual(len(pathogenic_df), 1)
        self.assertEqual(pathogenic_df.iloc[0]["rsid"], "rs429358")
        self.assertEqual(pathogenic_df.iloc[0]["clinical_significance"], "Pathogenic")
    
    def test_get_variants_by_condition(self):
        """Test getting variants by condition."""
        # Get variants for Test Condition 1
        condition_df = self.integration.get_variants_by_condition(self.variants_df, "Test Condition 1")
        
        # Check results
        self.assertEqual(len(condition_df), 1)
        self.assertEqual(condition_df.iloc[0]["rsid"], "rs429358")
        self.assertEqual(condition_df.iloc[0]["condition_name"], "Test Condition 1")
        
        # Get variants for non-existent condition
        condition_df = self.integration.get_variants_by_condition(self.variants_df, "Non-existent Condition")
        self.assertEqual(len(condition_df), 0)
    
    def test_get_clinical_summary(self):
        """Test generating a clinical summary."""
        # Get clinical summary
        summary = self.integration.get_clinical_summary(self.variants_df)
        
        # Check summary
        self.assertEqual(summary["pathogenic_count"], 1)
        self.assertEqual(summary["benign_count"], 1)
        self.assertEqual(summary["vus_count"], 1)
        self.assertEqual(len(summary["conditions"]), 3)
        
        # Check condition details
        condition_names = [c["name"] for c in summary["conditions"]]
        self.assertIn("Test Condition 1", condition_names)
        self.assertIn("Test Condition 2", condition_names)
        self.assertIn("Test Condition 3", condition_names)
        
        # Check first condition
        condition = next(c for c in summary["conditions"] if c["name"] == "Test Condition 1")
        self.assertEqual(condition["variant_count"], 1)
        self.assertEqual(condition["genes"], ["APOE"])


if __name__ == "__main__":
    unittest.main()
