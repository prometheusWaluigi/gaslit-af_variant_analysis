"""
Unit tests for the ClinVar integration module.
"""

import os
import json
import unittest
import tempfile
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock, mock_open
import io
import gzip

# Add project root to Python path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from src.gaslit_af.clinvar_integration import (
    ClinVarDownloader,
    ClinVarParser,
    ClinVarIntegration
)


class TestClinVarDownloader(unittest.TestCase):
    """Test cases for the ClinVarDownloader class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Initialize the downloader
        self.downloader = ClinVarDownloader(cache_dir=self.temp_path)
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    @patch('requests.get')
    def test_download_file(self, mock_get):
        """Test downloading a ClinVar file."""
        # Mock response
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.iter_content.return_value = [b'test data']
        mock_get.return_value = mock_response
        
        # Download file
        file_type = "variant_summary"
        file_path = self.downloader.download_file(file_type, force_download=True)
        
        # Check that the file was downloaded
        self.assertTrue(file_path.exists())
        self.assertEqual(file_path.parent, self.temp_path)
        self.assertEqual(file_path.name, "variant_summary.txt.gz")
        
        # Check that the request was made
        mock_get.assert_called_once()
        self.assertIn("variant_summary.txt.gz", mock_get.call_args[0][0])
    
    @patch('requests.head')
    def test_get_latest_release_date(self, mock_head):
        """Test getting the latest ClinVar release date."""
        # Mock response
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.headers = {'Last-Modified': 'Wed, 01 Jan 2023 12:00:00 GMT'}
        mock_head.return_value = mock_response
        
        # Get release date
        release_date = self.downloader.get_latest_release_date()
        
        # Check the release date
        self.assertEqual(release_date, "2023-01-01")
        
        # Check that the request was made
        mock_head.assert_called_once()
        self.assertIn("variant_summary.txt.gz", mock_head.call_args[0][0])


class TestClinVarParser(unittest.TestCase):
    """Test cases for the ClinVarParser class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Mock downloader
        self.mock_downloader = MagicMock()
        self.mock_downloader.cache_dir = self.temp_path
        # Add the file_types attribute needed by the parser
        self.mock_downloader.file_types = {
            "vcv_xml": "...",
            "rcv_xml": "...",
            "variant_summary": "...",
            "vcf_grch37": "...",
            "vcf_grch38": "...",
            "var_citations": "...",
            "cross_references": "..."
        }
        
        # Initialize the parser
        self.parser = ClinVarParser(self.mock_downloader)
        
        # Create test data
        self.variant_summary_data = """#AlleleID	Type	Name	GeneID	GeneSymbol	HGNC_ID	ClinicalSignificance	ClinSigSimple	LastEvaluated	RS# (dbSNP)	nsv/esv (dbVar)	RCVaccession	PhenotypeIDs	PhenotypeList	Origin	OriginSimple	Assembly	ChromosomeAccession	Chromosome	Start	Stop	ReferenceAllele	AlternateAllele	Cytogenetic	ReviewStatus	NumberSubmitters	Guidelines	TestedInGTR	OtherIDs	SubmitterCategories	VariationID	PositionVCF	ReferenceAlleleVCF	AlternateAlleleVCF
15	single nucleotide variant	NM_000059.3(BRCA2):c.7397T>C (p.Val2466Ala)	675	BRCA2	1101	Benign	1	2019-12-31	80359178		RCV000045034	MedGen:C0677776,OMIM:612555,Orphanet:ORPHA145	Breast-ovarian cancer, familial 2	germline	1	GRCh38	NC_000013.11	13	32356550	32356550	T	C	13q13.1	criteria provided, single submitter	1	1	N	LRG_293:g.30148T>C,LRG_293t1:c.7397T>C,NM_000059.3:c.7397T>C,NP_000050.2:p.Val2466Ala,NC_000013.10:g.32914766T>C,NC_000013.11:g.32356550T>C	research	12557	32356550	T	C
55	single nucleotide variant	NM_007294.3(BRCA1):c.5309C>T (p.Pro1770Leu)	672	BRCA1	1100	Uncertain significance	0	2019-12-31	80357284		RCV000045039	MedGen:C0677776,OMIM:604370,Orphanet:ORPHA145	Breast-ovarian cancer, familial 1	germline	1	GRCh38	NC_000017.11	17	43057051	43057051	G	A	17q21.31	criteria provided, single submitter	1	1	N	LRG_292:g.126456G>A,LRG_292t1:c.5309C>T,NM_007294.3:c.5309C>T,NP_009225.1:p.Pro1770Leu,NC_000017.10:g.41222975G>A,NC_000017.11:g.43057051G>A	research	14144	43057051	G	A
"""
        self.vcf_data = """##fileformat=VCFv4.1
##fileDate=2023-01-01
##source=ClinVar
##reference=GRCh38
##ID=<Description="ClinVar Variation ID">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance for this single variant">
##INFO=<ID=CLNREVSTAT,Number=.,Type=String,Description="ClinVar review status for the Variation ID">
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
13	32356550	15	T	C	.	.	CLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter
17	43057051	55	G	A	.	.	CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter
"""
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_parse_variant_summary(self):
        """Test parsing the variant summary file."""
        # Create a real temporary gzipped file
        mock_download_path = self.temp_path / "variant_summary.txt.gz"
        with gzip.open(mock_download_path, 'wb') as f_gz: # Write bytes
            f_gz.write(self.variant_summary_data.encode('utf-8'))
        
        # Mock downloader to return the path to the temp file
        self.mock_downloader.download_file.return_value = mock_download_path

        # Parse variant summary (actual function will call pd.read_csv)
        df = self.parser.parse_variant_summary(force_download=True)

        # Check that the downloader was called (using positional args as in the source)
        self.mock_downloader.download_file.assert_called_once_with("variant_summary", True)

        # Check the parsed data
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 2)
        self.assertIn("GeneSymbol", df.columns)
        self.assertEqual(df.loc[0, "GeneSymbol"], "BRCA2")

    @patch('builtins.open', new_callable=mock_open)
    def test_load_variant_summary(self, mock_open):
        """Test loading the processed variant summary data."""
        # Create a processed file
        processed_path = self.temp_path / "variant_summary_processed.parquet"
        
        # Mock parser method
        self.parser.parse_variant_summary = MagicMock()
        mock_df = pd.DataFrame({
            'GeneSymbol': ['BRCA2', 'BRCA1'],
            'ClinicalSignificance': ['Benign', 'Uncertain significance']
        })
        self.parser.parse_variant_summary.return_value = mock_df
        
        # Load variant summary
        df = self.parser.load_variant_summary(force_reprocess=True)
        
        # Check that the parser was called
        self.parser.parse_variant_summary.assert_called_once_with(force_download=True)
        
        # Check the loaded data
        self.assertEqual(len(df), 2)
        self.assertEqual(df.iloc[0]['GeneSymbol'], 'BRCA2')
        self.assertEqual(df.iloc[1]['GeneSymbol'], 'BRCA1')
    
    @patch('pandas.read_csv')
    @patch('gzip.open')
    def test_parse_vcf(self, mock_gzip_open, mock_read_csv):
        """Test parsing the ClinVar VCF file."""
        # Mock downloader
        mock_download_path = self.temp_path / "clinvar.vcf.gz"
        self.mock_downloader.download_file.return_value = mock_download_path

        # Mock gzip.open to return an iterator yielding STRINGS for header reading
        mock_gzip_open.return_value.__enter__.return_value = iter(self.vcf_data.splitlines(keepends=True))

        # Mock pandas.read_csv to return DataFrame INCLUDING the INFO column
        # The source code expects specific columns parsed by names=
        mock_read_csv.return_value = pd.DataFrame({
            'CHROM': ['13', '17'],
            'POS': [32356550, 43057051],
            'ID': ['15', '55'], # VCF IDs are typically strings
            'REF': ['T', 'G'],
            'ALT': ['C', 'A'],
            'QUAL': ['.', '.'], # Add dummy QUAL/FILTER/INFO
            'FILTER': ['PASS', 'PASS'],
            'INFO': [
                'ALLELEID=15;CLNDN=Breast-ovarian_cancer,_familial_2;CLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter',
                'ALLELEID=55;CLNDN=Breast-ovarian_cancer,_familial_1;CLNSIG=Uncertain_significance;CLNREVSTAT=criteria_provided,_single_submitter'
            ]
            # Note: The source code applies extract_clnsig, so the CLNSIG/CLNREVSTAT cols 
            #       in the *returned* df by parse_vcf will be derived from INFO.
            #       We don't need them in this *mocked* read_csv result.
        })

        # Parse VCF (will use mocked gzip.open for header, mocked read_csv for data)
        df = self.parser.parse_vcf(assembly="GRCh38", force_download=True)

        # Check downloader called
        self.mock_downloader.download_file.assert_called_once_with("vcf_grch38", True)

        # Check gzip.open was called correctly for the header read
        mock_gzip_open.assert_called_once_with(mock_download_path, 'rt')
        
        # Check that pandas.read_csv was called (by the source code)
        # We need to check the arguments passed by the source code to read_csv
        # Example: Check the filepath and the skiprows argument
        mock_read_csv.assert_called_once() 
        call_args, call_kwargs = mock_read_csv.call_args
        self.assertEqual(call_args[0], mock_download_path)
        self.assertEqual(call_kwargs.get('sep'), '\t')
        self.assertEqual(call_kwargs.get('comment'), '#')
        self.assertEqual(call_kwargs.get('compression'), 'gzip')
        # Check that the correct column names were passed to read_csv
        # Note: The source uses the header to define names, let's check that names were passed
        self.assertIn('names', call_kwargs)
        self.assertIsInstance(call_kwargs.get('names'), list) 

        # Check the final returned data (which is processed from our mock_read_csv result)
        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 2)
        # Verify columns created by the parse_vcf function after processing
        self.assertIn('CLNSIG', df.columns)
        self.assertEqual(df.loc[0, 'CLNSIG'], 'Benign')
        self.assertEqual(df.loc[1, 'CLNSIG'], 'Uncertain_significance')


class TestClinVarIntegration(unittest.TestCase):
    """Test cases for the ClinVarIntegration class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Initialize the integration
        self.integration = ClinVarIntegration(cache_dir=self.temp_path)
        
        # Mock data
        self.variant_summary_df = pd.DataFrame({
            '#AlleleID': [15, 55],
            'GeneSymbol': ['BRCA2', 'BRCA1'],
            'ClinicalSignificance': ['Benign', 'Uncertain significance'],
            'ReviewStatus': ['criteria provided, single submitter', 'criteria provided, single submitter'],
            'RS# (dbSNP)': [80359178, 80357284],
            'Chromosome': ['13', '17'],
            'Start': [32356550, 43057051],
            'Stop': [32356550, 43057051],
            'ReferenceAllele': ['T', 'G'],
            'AlternateAllele': ['C', 'A'],
            'PhenotypeList': ['Breast-ovarian cancer, familial 2', 'Breast-ovarian cancer, familial 1'],
            'LastEvaluated': ['2019-12-31', '2019-12-31'],
            'nsv/esv (dbVar)': ['', ''],
            'Origin': ['germline', 'germline'],
            'Assembly': ['GRCh38', 'GRCh38'],
            'Cytogenetic': ['13q13.1', '17q21.31']
        })
        
        self.vcf_df = pd.DataFrame({
            'CHROM': ['13', '17'],
            'POS': [32356550, 43057051],
            'ID': [15, 55],
            'REF': ['T', 'G'],
            'ALT': ['C', 'A'],
            'CLNSIG': ['Benign', 'Uncertain_significance'],
            'CLNREVSTAT': ['criteria_provided,_single_submitter', 'criteria_provided,_single_submitter']
        })
        
        # Create test variants DataFrame
        self.variants_df = pd.DataFrame({
            "rsid": ["rs80359178", "rs80357284"],
            "gene": ["BRCA2", "BRCA1"],
            "chrom": ["13", "17"],
            "pos": [32356550, 43057051],
            "ref": ["T", "G"],
            "alt": ["C", "A"]
        })
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    @patch.object(ClinVarIntegration, 'get_variant_summary')
    def test_get_variant_details(self, mock_get_variant_summary):
        """Test getting detailed information for a variant from ClinVar."""
        # Mock get_variant_summary
        mock_get_variant_summary.return_value = self.variant_summary_df
        
        # Get variant details
        variant_id = "rs80359178"
        details = self.integration.get_variant_details(variant_id)
        
        # Check that get_variant_summary was called
        mock_get_variant_summary.assert_called_once()
        
        # Check the returned details
        self.assertEqual(details["variant_id"], variant_id)
        self.assertEqual(details["clinvar_id"], 15)
        self.assertEqual(details["gene"], "BRCA2")
        self.assertEqual(details["clinical_significance"], "Benign")
        self.assertEqual(details["review_status"], "criteria provided, single submitter")
        self.assertEqual(details["chromosome"], "13")
        self.assertEqual(details["start"], 32356550)
        self.assertEqual(details["reference_allele"], "T")
        self.assertEqual(details["alternate_allele"], "C")
    
    @patch.object(ClinVarIntegration, 'get_variant_summary')
    def test_annotate_variants(self, mock_get_variant_summary):
        """Test annotating variants with ClinVar data."""
        # Mock get_variant_summary
        mock_get_variant_summary.return_value = self.variant_summary_df
        
        # Annotate variants
        annotated_df = self.integration.annotate_variants(self.variants_df)
        
        # Check that get_variant_summary was called
        mock_get_variant_summary.assert_called_once()
        
        # Check the annotations
        self.assertIn("clinvar_id", annotated_df.columns)
        self.assertIn("clinvar_significance", annotated_df.columns)
        self.assertIn("clinvar_review_status", annotated_df.columns)
        
        # Check values
        self.assertEqual(annotated_df.loc[0, "clinvar_id"], 15)
        self.assertEqual(annotated_df.loc[0, "clinvar_significance"], "Benign")
        self.assertEqual(annotated_df.loc[1, "clinvar_id"], 55)
        self.assertEqual(annotated_df.loc[1, "clinvar_significance"], "Uncertain significance")


if __name__ == "__main__":
    unittest.main()
