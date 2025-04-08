"""
Integration tests for the clinical variant and API integration modules.

This test verifies that the clinical variant data model and API integration
components work correctly with the main analysis pipeline.
"""

import os
import json
import unittest
import tempfile
import pandas as pd
import numpy as np
from pathlib import Path
from unittest.mock import patch, MagicMock

# Add project root to Python path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.gaslit_af.clinical_variants import ClinicalVariantManager, create_example_conditions
from src.gaslit_af.clinical_integration import ClinicalIntegration
from src.gaslit_af.api_integration import VariantAPIIntegration


class TestClinicalAPIIntegration(unittest.TestCase):
    """Integration tests for clinical variant and API integration."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Create test output directory
        self.output_dir = self.temp_path / "output"
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Create a test clinical data file
        self.clinical_data_path = self.temp_path / "clinical_data.json"
        self.clinical_data = create_example_conditions()
        
        with open(self.clinical_data_path, 'w') as f:
            json.dump(self.clinical_data, f)
        
        # Create a test variants DataFrame
        self.variants_df = pd.DataFrame({
            "rsid": ["rs429358", "rs7412", "rs80357906", "rs2350780", "rs6277"],
            "gene": ["APOE", "APOE", "BRCA1", "CHRM2", "DRD2"],
            "chrom": ["19", "19", "17", "7", "11"],
            "pos": [44908684, 44908822, 43071077, 136705002, 113412737],
            "ref": ["T", "C", "G", "G", "G"],
            "alt": ["C", "T", "A", "A", "A"],
            "genotype": ["T/C", "C/T", "G/A", "G/A", "G/A"],
            "quality": [100, 100, 100, 100, 100],
            "filter": ["PASS", "PASS", "PASS", "PASS", "PASS"],
            "info": [".", ".", ".", ".", "."]
        })
        
        # Save the variants DataFrame to a CSV file
        self.variants_csv = self.temp_path / "variants.csv"
        self.variants_df.to_csv(self.variants_csv, index=False)
        
        # Initialize components
        self.clinical_integration = ClinicalIntegration(self.clinical_data_path)
        self.api_integration = VariantAPIIntegration(cache_dir=self.temp_path / "api_cache")
        
        # Mock API responses
        self.mock_ensembl_data = {
            "rs429358": [
                {
                    "input": "rs429358",
                    "transcript_consequences": [
                        {
                            "gene_symbol": "APOE",
                            "consequence_terms": ["missense_variant"],
                            "impact": "MODERATE"
                        }
                    ]
                }
            ],
            "rs7412": [
                {
                    "input": "rs7412",
                    "transcript_consequences": [
                        {
                            "gene_symbol": "APOE",
                            "consequence_terms": ["missense_variant"],
                            "impact": "MODERATE"
                        }
                    ]
                }
            ],
            "rs2350780": [
                {
                    "input": "rs2350780",
                    "transcript_consequences": [
                        {
                            "gene_symbol": "CHRM2",
                            "consequence_terms": ["intron_variant"],
                            "impact": "MODIFIER"
                        }
                    ]
                }
            ],
            "rs6277": [
                {
                    "input": "rs6277",
                    "transcript_consequences": [
                        {
                            "gene_symbol": "DRD2",
                            "consequence_terms": ["synonymous_variant"],
                            "impact": "LOW"
                        }
                    ]
                }
            ]
        }
        
        self.mock_myvariant_data = {
            "rs429358": {
                "_id": "rs429358",
                "cadd": {"phred": 15.7},
                "clinvar": {"rcv": {"clinical_significance": "Pathogenic"}},
                "gnomad_genome": {"af": 0.164436},
                "dbnsfp": {
                    "sift": {"pred": "D"},
                    "polyphen2": {"hdiv": {"pred": "D"}}
                }
            },
            "rs7412": {
                "_id": "rs7412",
                "cadd": {"phred": 12.3},
                "clinvar": {"rcv": {"clinical_significance": "Benign"}},
                "gnomad_genome": {"af": 0.073},
                "dbnsfp": {
                    "sift": {"pred": "T"},
                    "polyphen2": {"hdiv": {"pred": "B"}}
                }
            },
            "rs2350780": {
                "_id": "rs2350780",
                "cadd": {"phred": 5.2},
                "gnomad_genome": {"af": 0.15},
                "dbnsfp": {
                    "sift": {"pred": "T"},
                    "polyphen2": {"hdiv": {"pred": "B"}}
                }
            },
            "rs6277": {
                "_id": "rs6277",
                "cadd": {"phred": 3.8},
                "gnomad_genome": {"af": 0.25},
                "dbnsfp": {
                    "sift": {"pred": "T"},
                    "polyphen2": {"hdiv": {"pred": "B"}}
                }
            }
        }
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    @patch.object(VariantAPIIntegration, 'get_variant_details')
    def test_end_to_end_pipeline(self, mock_get_variant_details):
        """Test the end-to-end pipeline with clinical and API integration."""
        # Mock API responses
        def mock_get_details(variant_id):
            if variant_id in self.mock_myvariant_data:
                mv_data = self.mock_myvariant_data[variant_id]
                ensembl_data = self.mock_ensembl_data.get(variant_id, [{}])[0]
                
                details = {
                    "variant_id": variant_id,
                    "gene": ensembl_data.get("transcript_consequences", [{}])[0].get("gene_symbol"),
                    "consequence": ensembl_data.get("transcript_consequences", [{}])[0].get("consequence_terms", [""])[0],
                    "impact": ensembl_data.get("transcript_consequences", [{}])[0].get("impact"),
                    "clinical_significance": mv_data.get("clinvar", {}).get("rcv", {}).get("clinical_significance"),
                    "allele_frequency": mv_data.get("gnomad_genome", {}).get("af"),
                    "pathogenicity_scores": {
                        "cadd_phred": mv_data.get("cadd", {}).get("phred"),
                        "sift": mv_data.get("dbnsfp", {}).get("sift", {}).get("pred"),
                        "polyphen": mv_data.get("dbnsfp", {}).get("polyphen2", {}).get("hdiv", {}).get("pred")
                    }
                }
                return details
            return {
                "variant_id": variant_id,
                "gene": None,
                "consequence": None,
                "impact": None,
                "clinical_significance": None,
                "allele_frequency": None,
                "pathogenicity_scores": {}
            }
        
        mock_get_variant_details.side_effect = mock_get_details
        
        # Step 1: Load the variants
        variants = pd.read_csv(self.variants_csv)
        self.assertEqual(len(variants), 5)
        
        # Step 2: Annotate with clinical data
        clinical_df = self.clinical_integration.annotate_variants(variants)
        
        # Verify clinical annotations
        self.assertIn("clinical_significance", clinical_df.columns)
        self.assertIn("condition_name", clinical_df.columns)
        
        # Check CHRM2 variant annotation
        chrm2_row = clinical_df[clinical_df["gene"] == "CHRM2"].iloc[0]
        self.assertEqual(chrm2_row["condition_name"], "CHRM2-Related Cognitive Function")
        self.assertEqual(chrm2_row["clinical_significance"], "Risk Factor")
        
        # Step 3: Generate clinical report
        report_path = self.clinical_integration.generate_clinical_report(clinical_df, self.output_dir)
        self.assertIsNotNone(report_path)
        self.assertTrue(Path(report_path).exists())
        
        # Step 4: Get clinical summary
        summary = self.clinical_integration.get_clinical_summary(clinical_df)
        self.assertIn("conditions", summary)
        self.assertTrue(len(summary["conditions"]) > 0)
        
        # Step 5: Annotate with API data
        for variant_id in variants["rsid"]:
            # Call get_variant_details for each variant to populate the mock
            self.api_integration.get_variant_details(variant_id)
        
        # Verify that the API was called for each variant
        self.assertEqual(mock_get_variant_details.call_count, 5)
        
        # Step 6: Combine clinical and API annotations
        # In a real scenario, this would be done in the workflow module
        # Here we'll simulate it by merging the data
        
        # First, get API annotations for each variant
        api_annotations = {}
        for variant_id in variants["rsid"]:
            api_annotations[variant_id] = self.api_integration.get_variant_details(variant_id)
        
        # Add API annotations to the DataFrame
        clinical_df["ensembl_consequence"] = clinical_df["rsid"].apply(
            lambda x: api_annotations.get(x, {}).get("consequence")
        )
        clinical_df["ensembl_impact"] = clinical_df["rsid"].apply(
            lambda x: api_annotations.get(x, {}).get("impact")
        )
        clinical_df["cadd_phred"] = clinical_df["rsid"].apply(
            lambda x: api_annotations.get(x, {}).get("pathogenicity_scores", {}).get("cadd_phred")
        )
        clinical_df["gnomad_af"] = clinical_df["rsid"].apply(
            lambda x: api_annotations.get(x, {}).get("allele_frequency")
        )
        
        # Verify API annotations
        self.assertIn("ensembl_consequence", clinical_df.columns)
        self.assertIn("cadd_phred", clinical_df.columns)
        
        # Check APOE variant annotation
        apoe_row = clinical_df[clinical_df["rsid"] == "rs429358"].iloc[0]
        self.assertEqual(apoe_row["ensembl_consequence"], "missense_variant")
        self.assertEqual(apoe_row["ensembl_impact"], "MODERATE")
        self.assertEqual(apoe_row["cadd_phred"], 15.7)
        
        # Step 7: Save the fully annotated data
        output_file = self.output_dir / "annotated_variants.csv"
        clinical_df.to_csv(output_file, index=False)
        self.assertTrue(output_file.exists())
        
        # Verify the saved file
        saved_df = pd.read_csv(output_file)
        self.assertEqual(len(saved_df), 5)
        self.assertIn("clinical_significance", saved_df.columns)
        self.assertIn("ensembl_consequence", saved_df.columns)
        self.assertIn("cadd_phred", saved_df.columns)


if __name__ == "__main__":
    unittest.main()
