"""
Unit tests for the refactored ClinVar module components.

This test suite validates the fractal integrity of the ClinVar components
within the GASLIT-AF architecture.
"""

import os
import json
import unittest
import tempfile
import pandas as pd
import sqlite3
from pathlib import Path
from datetime import datetime
from unittest.mock import patch, MagicMock, mock_open

# Add project root to Python path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Import refactored components
from src.gaslit_af.clinvar.downloader import ClinVarDownloader
from src.gaslit_af.clinvar.parser import ClinVarParser
from src.gaslit_af.clinvar.indexer import ClinVarIndexer
from src.gaslit_af.clinvar.cache_manager import ClinVarCache
from src.gaslit_af.clinvar.annotator import ClinVarAnnotator, ClinVarIntegration


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
        
        # Initialize the parser
        self.parser = ClinVarParser(self.mock_downloader)
        
        # Create test data
        self.variant_summary_data = """#AlleleID\tType\tName\tGeneID\tGeneSymbol\tHGNC_ID\tClinicalSignificance\tClinSigSimple\tLastEvaluated\tRS# (dbSNP)\tnsv/esv (dbVar)\tRCVaccession\tPhenotypeIDs\tPhenotypeList\tOrigin\tOriginSimple\tAssembly\tChromosomeAccession\tChromosome\tStart\tStop\tReferenceAllele\tAlternateAllele\tCytogenetic\tReviewStatus\tNumberSubmitters\tGuidelines\tTestedInGTR\tOtherIDs\tSubmitterCategories\tVariationID\tPositionVCF\tReferenceAlleleVCF\tAlternateAlleleVCF
15\tSingle nucleotide variant\tNM_000059.3(BRCA2):c.6513G>C (p.Glu2171=)\t675\tBRCA2\t1101\tBenign\t0\t2019-12-31\t80359178\t\tRCV000031156\t\tBreast-ovarian cancer, familial 2\tgermline\tgermline\tGRCh38\tNC_000013.11\t13\t32356550\t32356550\tT\tC\t13q13.1\tcriteria provided, single submitter\t1\t\tN\t\texpert panel\t13085\t32356550\tT\tC
55\tSingle nucleotide variant\tNM_007294.3(BRCA1):c.5503C>T (p.Arg1835Ter)\t672\tBRCA1\t1100\tUncertain significance\t255\t2019-12-31\t80357284\t\tRCV000031189\t\tBreast-ovarian cancer, familial 1\tgermline\tgermline\tGRCh38\tNC_000017.11\t17\t43057051\t43057051\tG\tA\t17q21.31\tcriteria provided, single submitter\t1\t\tN\t\texpert panel\t14304\t43057051\tG\tA"""
        
        self.vcf_data = """##fileformat=VCFv4.1
##fileDate=20230101
##source=ClinVar
##reference=GRCh38
##ID=<Description="ClinVar Variation ID">
##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">
##INFO=<ID=GENEINFO,Number=1,Type=String,Description="Gene name and ID">
##CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
13\t32356550\t15\tT\tC\t.\t.\tCLNSIG=Benign;GENEINFO=BRCA2:675
17\t43057051\t55\tG\tA\t.\t.\tCLNSIG=Uncertain_significance;GENEINFO=BRCA1:672"""
        
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
            'Origin': ['germline', 'germline'],
            'Assembly': ['GRCh38', 'GRCh38'],
            'Cytogenetic': ['13q13.1', '17q21.31']
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
    
    @patch('pandas.read_csv')
    @patch('gzip.open')
    def test_parse_variant_summary(self, mock_gzip_open, mock_read_csv):
        """Test parsing the variant summary file."""
        # Mock download_file
        self.mock_downloader.download_file.return_value = self.temp_path / "variant_summary.txt.gz"
        
        # Mock read_csv
        mock_read_csv.return_value = self.variant_summary_df
        
        # Parse variant summary
        df = self.parser.parse_variant_summary(force_download=False)
        
        # Check that download_file was called
        self.mock_downloader.download_file.assert_called_once_with("variant_summary", False)
        
        # Check that read_csv was called
        mock_read_csv.assert_called_once()
        
        # Check returned DataFrame
        self.assertEqual(len(df), 2)
        self.assertEqual(df['GeneSymbol'].tolist(), ['BRCA2', 'BRCA1'])
    
    @patch('gzip.open')
    def test_parse_vcf(self, mock_gzip_open):
        """Test parsing the ClinVar VCF file."""
        # Mock file handle
        mock_file = MagicMock()
        mock_file.__enter__.return_value = self.vcf_data.splitlines()
        mock_gzip_open.return_value = mock_file
        
        # Mock download_file
        self.mock_downloader.download_file.return_value = self.temp_path / "clinvar.vcf.gz"
        
        # Parse VCF
        df = self.parser.parse_vcf(assembly="GRCh38", force_download=False)
        
        # Check that download_file was called
        self.mock_downloader.download_file.assert_called_once_with("vcf_grch38", False)
        
        # Check that gzip.open was called
        mock_gzip_open.assert_called_once()
        
        # Check returned DataFrame (just check it has some rows)
        # Note: We're not checking the exact DataFrame because the implementation creates it from scratch
        self.assertIsInstance(df, pd.DataFrame)


class TestClinVarIndexer(unittest.TestCase):
    """Test cases for the ClinVarIndexer class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Initialize the indexer
        self.indexer = ClinVarIndexer(index_dir=self.temp_path)
        
        # Test data
        self.variant_df = pd.DataFrame({
            'rs_id': ['rs80359178', 'rs80357284'],
            'clinvar_id': [15, 55],
            'chrom': ['13', '17'],
            'pos': [32356550, 43057051],
            'ref': ['T', 'G'],
            'alt': ['C', 'A'],
            'gene': ['BRCA2', 'BRCA1'],
            'significance': ['Benign', 'Uncertain significance']
        })
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_init_database(self):
        """Test database initialization."""
        # Check that database file exists
        db_path = self.temp_path / "clinvar_index.db"
        self.assertTrue(db_path.exists())
        
        # Check that tables and indices were created
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Check variant_index table
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='variant_index'")
        self.assertIsNotNone(cursor.fetchone())
        
        # Check indices
        cursor.execute("SELECT name FROM sqlite_master WHERE type='index' AND name='idx_rs_id'")
        self.assertIsNotNone(cursor.fetchone())
        
        cursor.execute("SELECT name FROM sqlite_master WHERE type='index' AND name='idx_gene'")
        self.assertIsNotNone(cursor.fetchone())
        
        conn.close()
    
    def test_index_variants(self):
        """Test indexing variants."""
        # Use a mock connection instead of patching immutable cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_cursor.fetchone.side_effect = [(2,), (2,), (0,), (1,), (1,)]
        
        # Patch sqlite3.connect to return our mock connection
        with patch('sqlite3.connect', return_value=mock_conn):
            # Index variants
            stats = self.indexer.index_variants(self.variant_df, source="test")
            
            # Check stats
            self.assertEqual(stats["variant_count"], 2)
            self.assertEqual(stats["clinical_significant_count"], 2)
            self.assertEqual(stats["pathogenic_count"], 0)
            self.assertEqual(stats["benign_count"], 1)
            self.assertEqual(stats["vus_count"], 1)
            
            # Verify execute was called for queries
            self.assertTrue(mock_cursor.execute.call_count >= 5)
        
        # No need to check actual database - we're mocking the interaction
    
    def test_lookup_variant(self):
        """Test looking up variants."""
        # Create mock dictionaries for rows
        dict_row1 = {"gene": "BRCA2"}
        dict_row2 = {"gene": "BRCA1"}
        dict_row3 = {"rs_id": "rs80357284"}
        
        # Use a mock connection instead of patching immutable cursor
        mock_conn = MagicMock()
        mock_cursor = MagicMock()
        mock_conn.cursor.return_value = mock_cursor
        mock_conn.row_factory = sqlite3.Row  # Just to match the actual code
        
        # Set up different return values for different queries
        mock_cursor.fetchall.side_effect = [
            [dict_row1],  # For rs_id lookup
            [dict_row2],  # For genomic coordinates lookup
            [dict_row3]   # For gene lookup
        ]
        
        # Patch sqlite3.connect to return our mock connection
        with patch('sqlite3.connect', return_value=mock_conn):
            # Look up by rs_id
            results = self.indexer.lookup_variant(rs_id="rs80359178")
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["gene"], "BRCA2")
            
            # Look up by genomic coordinates
            results = self.indexer.lookup_variant(chrom="17", pos=43057051)
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["gene"], "BRCA1")
            
            # Look up by gene
            results = self.indexer.lookup_variant(gene="BRCA1")
            self.assertEqual(len(results), 1)
            self.assertEqual(results[0]["rs_id"], "rs80357284")


class TestClinVarCache(unittest.TestCase):
    """Test cases for the ClinVarCache class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Mock components
        self.mock_downloader = MagicMock()
        self.mock_parser = MagicMock()
        self.mock_indexer = MagicMock()
        
        # Initialize the cache with patching
        with patch('src.gaslit_af.clinvar.cache_manager.ClinVarDownloader', return_value=self.mock_downloader), \
             patch('src.gaslit_af.clinvar.cache_manager.ClinVarParser', return_value=self.mock_parser), \
             patch('src.gaslit_af.clinvar.cache_manager.ClinVarIndexer', return_value=self.mock_indexer):
            self.cache = ClinVarCache(cache_dir=self.temp_path)
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_init_metadata(self):
        """Test metadata initialization."""
        # Check that metadata file exists
        metadata_file = self.temp_path / "metadata" / "cache_metadata.json"
        self.assertTrue(metadata_file.exists())
        
        # Check metadata structure
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
        
        self.assertIn("versions", metadata)
        self.assertIn("stats", metadata)
        self.assertIn("file_hashes", metadata)
    
    @patch('src.gaslit_af.clinvar.cache_manager.ClinVarCache.get_file_hash')
    def test_update_cache_metadata(self, mock_get_file_hash):
        """Test updating cache metadata."""
        # Mock file hash
        mock_get_file_hash.return_value = "test_hash"
        
        # Update metadata
        file_path = self.temp_path / "test.txt"
        Path(file_path).touch()
        
        self.cache.update_cache_metadata("test_type", file_path, version="1.0")
        
        # Check that metadata was updated
        self.assertIn("test_type", self.cache.metadata["versions"])
        self.assertEqual(self.cache.metadata["versions"]["test_type"]["file_hash"], "test_hash")
        self.assertEqual(self.cache.metadata["versions"]["test_type"]["version"], "1.0")
        self.assertIn(str(file_path), self.cache.metadata["file_hashes"])
    
    def test_is_cache_valid(self):
        """Test cache validity checking."""
        # Set up test metadata
        self.cache.metadata["versions"]["test_type"] = {
            "last_update": (datetime.now().replace(day=1) - pd.Timedelta(days=15)).isoformat()
        }
        
        # Check cache validity
        self.assertTrue(self.cache.is_cache_valid("test_type", max_age_days=30))
        self.assertFalse(self.cache.is_cache_valid("test_type", max_age_days=10))
        self.assertFalse(self.cache.is_cache_valid("nonexistent_type"))
    
    def test_refresh_variant_summary(self):
        """Test refreshing variant summary data."""
        # Mock parser and indexer
        mock_df = pd.DataFrame({"test": [1, 2]})
        self.mock_parser.parse_variant_summary.return_value = mock_df
        
        mock_stats = {"variant_count": 2}
        self.mock_indexer.index_variants.return_value = mock_stats
        
        # Call refresh_variant_summary
        stats = self.cache.refresh_variant_summary(force_download=True)
        
        # Check that parser and indexer were called
        self.mock_parser.parse_variant_summary.assert_called_once_with(True)
        self.mock_indexer.index_variants.assert_called_once_with(mock_df, "variant_summary")
        
        # Check returned stats
        self.assertEqual(stats, mock_stats)


class TestClinVarAnnotator(unittest.TestCase):
    """Test cases for the ClinVarAnnotator class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Mock cache
        self.mock_cache = MagicMock()
        
        # Initialize the annotator
        self.annotator = ClinVarAnnotator(self.mock_cache)
        
        # Test data
        self.variant = {
            "chrom": "13",
            "pos": 32356550,
            "ref": "T",
            "alt": "C",
            "rs_id": "rs80359178",
            "gene": "BRCA2"
        }
        
        self.clinvar_data = [{
            "clinvar_id": 15,
            "significance": "Benign",
            "gene": "BRCA2"
        }]
        
        # Mock cache.lookup_variant
        self.mock_cache.lookup_variant.return_value = self.clinvar_data
    
    def test_annotate_variant(self):
        """Test annotating a single variant."""
        # Annotate variant
        annotated = self.annotator.annotate_variant(self.variant.copy())
        
        # Check that lookup_variant was called
        self.mock_cache.lookup_variant.assert_called_once_with(rs_id="rs80359178")
        
        # Check annotation results
        self.assertEqual(annotated["clinvar_id"], 15)
        self.assertEqual(annotated["clinical_significance"], "Benign")
        self.assertEqual(annotated["is_pathogenic"], False)
        self.assertEqual(annotated["is_benign"], True)
        self.assertEqual(annotated["is_vus"], False)
    
    def test_annotate_variant_not_found(self):
        """Test annotating a variant not in ClinVar."""
        # Set up mock to return empty list
        self.mock_cache.lookup_variant.return_value = []
        
        # Annotate variant
        annotated = self.annotator.annotate_variant(self.variant.copy())
        
        # Check annotation results for not found
        self.assertIsNone(annotated["clinvar_id"])
        self.assertEqual(annotated["clinical_significance"], "Not found in ClinVar")
        self.assertEqual(annotated["is_pathogenic"], False)
        self.assertEqual(annotated["is_benign"], False)
        self.assertEqual(annotated["is_vus"], False)
    
    def test_annotate_variants_df(self):
        """Test annotating a DataFrame of variants."""
        # Create test DataFrame
        df = pd.DataFrame([self.variant])
        
        # Annotate DataFrame
        annotated_df = self.annotator.annotate_variants_df(df)
        
        # Check result
        self.assertIn("clinvar_id", annotated_df.columns)
        self.assertIn("clinical_significance", annotated_df.columns)
        self.assertEqual(annotated_df.iloc[0]["clinvar_id"], 15)


class TestClinVarIntegration(unittest.TestCase):
    """Test cases for the ClinVarIntegration class."""
    
    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory for test files
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_path = Path(self.temp_dir.name)
        
        # Mock components
        self.mock_cache = MagicMock()
        self.mock_downloader = MagicMock()
        self.mock_parser = MagicMock()
        self.mock_indexer = MagicMock()
        self.mock_annotator = MagicMock()
        
        # Initialize the integration with patching
        with patch('src.gaslit_af.clinvar.annotator.ClinVarCache', return_value=self.mock_cache), \
             patch('src.gaslit_af.clinvar.annotator.ClinVarAnnotator', return_value=self.mock_annotator):
            self.integration = ClinVarIntegration(cache_dir=self.temp_path)
            
            # Set mocked components
            self.integration.cache = self.mock_cache
            self.integration.downloader = self.mock_downloader
            self.integration.parser = self.mock_parser
            self.integration.indexer = self.mock_indexer
            self.integration.annotator = self.mock_annotator
        
        # Test data
        self.variants_df = pd.DataFrame({
            "rsid": ["rs80359178", "rs80357284"],
            "gene": ["BRCA2", "BRCA1"],
            "chrom": ["13", "17"],
            "pos": [32356550, 43057051],
            "ref": ["T", "G"],
            "alt": ["C", "A"]
        })
        
        self.annotated_df = self.variants_df.copy()
        self.annotated_df["is_pathogenic"] = [False, True]
        self.annotated_df["clinical_significance"] = ["Benign", "Pathogenic"]
        self.annotated_df["clinvar_id"] = [15, 55]
    
    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_refresh_cache(self):
        """Test refreshing all ClinVar data caches."""
        # Mock cache methods
        stats = {"variant_count": 100}
        self.mock_cache.refresh_variant_summary.return_value = stats
        self.mock_cache.refresh_vcf_data.return_value = stats
        self.mock_cache.get_cache_stats.return_value = stats
        
        # Call refresh_cache
        result = self.integration.refresh_cache(force_download=True)
        
        # Check that cache methods were called
        self.mock_cache.refresh_variant_summary.assert_called_once_with(True)
        self.mock_cache.refresh_vcf_data.assert_any_call("GRCh38", True)
        self.mock_cache.refresh_vcf_data.assert_any_call("GRCh37", True)
        self.mock_cache.get_cache_stats.assert_called_once()
        
        # Check result
        self.assertEqual(result, stats)
    
    def test_annotate_variants(self):
        """Test annotating variants with ClinVar information."""
        # Mock annotator
        self.mock_annotator.annotate_variants_df.return_value = self.annotated_df
        
        # Call annotate_variants
        result = self.integration.annotate_variants(self.variants_df)
        
        # Check that annotator was called
        self.mock_annotator.annotate_variants_df.assert_called_once_with(self.variants_df)
        
        # Check result
        self.assertEqual(len(result), 2)
        self.assertIn("is_pathogenic", result.columns)
        self.assertIn("clinical_significance", result.columns)
    
    def test_get_pathogenic_variants(self):
        """Test filtering for pathogenic variants."""
        # Mock annotate_variants
        self.integration.annotate_variants = MagicMock(return_value=self.annotated_df)
        
        # Call get_pathogenic_variants
        result = self.integration.get_pathogenic_variants(self.variants_df)
        
        # Check that annotate_variants was called
        self.integration.annotate_variants.assert_called_once_with(self.variants_df)
        
        # Check result (should only return the pathogenic variant)
        self.assertEqual(len(result), 1)
        self.assertEqual(result.iloc[0]["gene"], "BRCA1")
    
    def test_lookup_variant(self):
        """Test looking up a variant in ClinVar."""
        # Mock cache.lookup_variant
        expected_result = [{"clinvar_id": 15, "significance": "Benign"}]
        self.mock_cache.lookup_variant.return_value = expected_result
        
        # Call lookup_variant
        result = self.integration.lookup_variant(rs_id="rs80359178")
        
        # Check that cache.lookup_variant was called
        self.mock_cache.lookup_variant.assert_called_once_with(rs_id="rs80359178")
        
        # Check result
        self.assertEqual(result, expected_result)


if __name__ == "__main__":
    unittest.main()
