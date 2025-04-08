"""
Unit tests for the clinical variants module.
"""

import os
import json
import unittest
import tempfile
import pandas as pd
from pathlib import Path

# Add project root to Python path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.gaslit_af.clinical_variants import ClinicalVariantManager, create_example_conditions


class TestClinicalVariantManager(unittest.TestCase):
    """Test cases for the ClinicalVariantManager class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Create a test schema file
        self.schema_path = self.temp_path / "test_schema.json"
        self.schema = {
            "$schema": "http://json-schema.org/draft-07/schema#",
            "title": "Clinical Variants Schema",
            "type": "object",
            "required": ["conditions"],
            "properties": {
                "conditions": {
                    "type": "array",
                    "items": {"$ref": "#/definitions/Condition"}
                }
            },
            "definitions": {
                "Condition": {
                    "type": "object",
                    "required": ["name", "genetic_data", "status"],
                    "properties": {
                        "name": {"type": "string"},
                        "description": {"type": "string"},
                        "symptoms": {
                            "type": "array",
                            "items": {"type": "string"}
                        },
                        "status": {
                            "type": "object",
                            "required": ["risk", "confidence"],
                            "properties": {
                                "risk": {"type": "string"},
                                "confidence": {"type": "string"},
                                "classification": {"type": "string"}
                            }
                        },
                        "genetic_data": {
                            "type": "object",
                            "required": ["gene", "variant_id"],
                            "properties": {
                                "gene": {"type": "string"},
                                "variant_id": {"type": "string"},
                                "rcv": {"type": "string"},
                                "genotype": {"type": "string"}
                            }
                        },
                        "risk_assessment": {
                            "type": "object",
                            "properties": {
                                "frequency": {"type": "number"},
                                "version": {"type": "string"}
                            }
                        }
                    }
                }
            }
        }
        
        with open(self.schema_path, 'w') as f:
            json.dump(self.schema, f)
        
        # Create a test conditions file
        self.conditions_path = self.temp_path / "test_conditions.json"
        self.conditions = {
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
                        "gene": "TEST1",
                        "variant_id": "rs123456",
                        "rcv": "RCV000123456",
                        "genotype": "A/G"
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
                        "gene": "TEST2",
                        "variant_id": "rs654321",
                        "rcv": "RCV000654321",
                        "genotype": "C/T"
                    },
                    "risk_assessment": {
                        "frequency": 0.05,
                        "version": "GRCh38"
                    }
                }
            ]
        }
        
        with open(self.conditions_path, 'w') as f:
            json.dump(self.conditions, f)
        
        # Initialize the manager
        self.manager = ClinicalVariantManager(schema_path=self.schema_path)
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_load_schema(self):
        """Test loading the JSON schema."""
        schema = self.manager._load_schema()
        self.assertIsInstance(schema, dict)
        self.assertEqual(schema["title"], "Clinical Variants Schema")
    
    def test_load_conditions(self):
        """Test loading conditions from a file."""
        conditions = self.manager.load_conditions(self.conditions_path)
        self.assertIsInstance(conditions, dict)
        self.assertEqual(len(conditions["conditions"]), 2)
        self.assertEqual(conditions["conditions"][0]["name"], "Test Condition 1")
    
    def test_validate_conditions(self):
        """Test validating conditions against the schema."""
        # Valid conditions
        result = self.manager.validate_conditions(self.conditions)
        self.assertTrue(result)
        
        # Invalid conditions (missing required field)
        invalid_conditions = {
            "conditions": [
                {
                    "name": "Invalid Condition",
                    # Missing genetic_data
                    "status": {
                        "risk": "Pathogenic",
                        "confidence": "High"
                    }
                }
            ]
        }
        
        with self.assertRaises(Exception):
            self.manager.validate_conditions(invalid_conditions)
    
    def test_add_condition(self):
        """Test adding a new condition."""
        # Load initial conditions
        self.manager.load_conditions(self.conditions_path)
        initial_count = len(self.manager.conditions_data["conditions"])
        
        # Add a new condition
        new_condition = {
            "name": "New Test Condition",
            "description": "A new test condition",
            "symptoms": ["New Symptom 1", "New Symptom 2"],
            "status": {
                "risk": "Risk Factor",
                "confidence": "Low",
                "classification": "R"
            },
            "genetic_data": {
                "gene": "TEST3",
                "variant_id": "rs789012",
                "rcv": "RCV000789012",
                "genotype": "G/T"
            },
            "risk_assessment": {
                "frequency": 0.02,
                "version": "GRCh38"
            }
        }
        
        self.manager.add_condition(new_condition)
        
        # Check that the condition was added
        self.assertEqual(len(self.manager.conditions_data["conditions"]), initial_count + 1)
        self.assertEqual(self.manager.conditions_data["conditions"][-1]["name"], "New Test Condition")
    
    def test_save_conditions(self):
        """Test saving conditions to a file."""
        # Load initial conditions
        self.manager.load_conditions(self.conditions_path)
        
        # Add a new condition
        new_condition = {
            "name": "New Test Condition",
            "description": "A new test condition",
            "symptoms": ["New Symptom 1", "New Symptom 2"],
            "status": {
                "risk": "Risk Factor",
                "confidence": "Low",
                "classification": "R"
            },
            "genetic_data": {
                "gene": "TEST3",
                "variant_id": "rs789012",
                "rcv": "RCV000789012",
                "genotype": "G/T"
            },
            "risk_assessment": {
                "frequency": 0.02,
                "version": "GRCh38"
            }
        }
        
        self.manager.add_condition(new_condition)
        
        # Save conditions to a new file
        output_path = self.temp_path / "output_conditions.json"
        self.manager.save_conditions(output_path)
        
        # Check that the file was created
        self.assertTrue(output_path.exists())
        
        # Load the saved conditions and check content
        with open(output_path, 'r') as f:
            saved_conditions = json.load(f)
        
        self.assertEqual(len(saved_conditions["conditions"]), 3)
        self.assertEqual(saved_conditions["conditions"][2]["name"], "New Test Condition")
    
    def test_get_condition_by_gene(self):
        """Test getting conditions by gene."""
        # Load conditions
        self.manager.load_conditions(self.conditions_path)
        
        # Get conditions for TEST1
        conditions = self.manager.get_condition_by_gene("TEST1")
        self.assertEqual(len(conditions), 1)
        self.assertEqual(conditions[0]["name"], "Test Condition 1")
        
        # Get conditions for non-existent gene
        conditions = self.manager.get_condition_by_gene("NONEXISTENT")
        self.assertEqual(len(conditions), 0)
    
    def test_get_condition_by_variant_id(self):
        """Test getting conditions by variant ID."""
        # Load conditions
        self.manager.load_conditions(self.conditions_path)
        
        # Get conditions for rs123456
        conditions = self.manager.get_condition_by_variant_id("rs123456")
        self.assertEqual(len(conditions), 1)
        self.assertEqual(conditions[0]["name"], "Test Condition 1")
        
        # Get conditions for non-existent variant ID
        conditions = self.manager.get_condition_by_variant_id("rs999999")
        self.assertEqual(len(conditions), 0)
    
    def test_get_conditions_by_risk(self):
        """Test getting conditions by risk level."""
        # Load conditions
        self.manager.load_conditions(self.conditions_path)
        
        # Get pathogenic conditions
        conditions = self.manager.get_conditions_by_risk("Pathogenic")
        self.assertEqual(len(conditions), 1)
        self.assertEqual(conditions[0]["name"], "Test Condition 1")
        
        # Get benign conditions
        conditions = self.manager.get_conditions_by_risk("Benign")
        self.assertEqual(len(conditions), 1)
        self.assertEqual(conditions[0]["name"], "Test Condition 2")
    
    def test_get_pathogenic_conditions(self):
        """Test getting pathogenic conditions."""
        # Load conditions
        self.manager.load_conditions(self.conditions_path)
        
        # Get pathogenic conditions
        conditions = self.manager.get_pathogenic_conditions()
        self.assertEqual(len(conditions), 1)
        self.assertEqual(conditions[0]["name"], "Test Condition 1")
    
    def test_annotate_variants(self):
        """Test annotating variants with clinical information."""
        # Load conditions
        self.manager.load_conditions(self.conditions_path)
        
        # Create a test variants DataFrame
        variants_df = pd.DataFrame({
            "gene": ["TEST1", "TEST2", "TEST3"],
            "rsid": ["rs123456", "rs654321", "rs999999"],
            "chrom": ["1", "2", "3"],
            "pos": [1000, 2000, 3000],
            "ref": ["A", "C", "G"],
            "alt": ["G", "T", "A"]
        })
        
        # Annotate variants
        annotated_df = self.manager.annotate_variants(variants_df)
        
        # Check annotations
        self.assertIsNotNone(annotated_df["clinical_significance"][0])
        self.assertEqual(annotated_df["clinical_significance"][0], "Pathogenic")
        self.assertEqual(annotated_df["condition_name"][0], "Test Condition 1")
        
        self.assertIsNotNone(annotated_df["clinical_significance"][1])
        self.assertEqual(annotated_df["clinical_significance"][1], "Benign")
        self.assertEqual(annotated_df["condition_name"][1], "Test Condition 2")
        
        self.assertIsNone(annotated_df["clinical_significance"][2])
        self.assertIsNone(annotated_df["condition_name"][2])
    
    def test_generate_clinical_report(self):
        """Test generating a clinical report."""
        # Load conditions
        self.manager.load_conditions(self.conditions_path)
        
        # Create a test variants DataFrame
        variants_df = pd.DataFrame({
            "gene": ["TEST1", "TEST2", "TEST3"],
            "rsid": ["rs123456", "rs654321", "rs999999"],
            "chrom": ["1", "2", "3"],
            "pos": [1000, 2000, 3000],
            "ref": ["A", "C", "G"],
            "alt": ["G", "T", "A"],
            "genotype": ["A/G", "C/T", "G/A"]
        })
        
        # Generate report
        report_path = self.manager.generate_clinical_report(variants_df, self.temp_path)
        
        # Check that the report was created
        self.assertIsNotNone(report_path)
        report_file = Path(report_path)
        self.assertTrue(report_file.exists())
        
        # Check report content
        with open(report_file, 'r') as f:
            report_content = f.read()
        
        self.assertIn("Clinical Variant Report", report_content)
        self.assertIn("Pathogenic/Likely Pathogenic Variants", report_content)
        self.assertIn("TEST1", report_content)
        self.assertIn("rs123456", report_content)
    
    def test_create_example_conditions(self):
        """Test creating example conditions."""
        example_conditions = create_example_conditions()
        self.assertIsInstance(example_conditions, dict)
        self.assertIn("conditions", example_conditions)
        self.assertTrue(len(example_conditions["conditions"]) > 0)
        
        # Validate example conditions against the schema
        result = self.manager.validate_conditions(example_conditions)
        self.assertTrue(result)


if __name__ == "__main__":
    unittest.main()
