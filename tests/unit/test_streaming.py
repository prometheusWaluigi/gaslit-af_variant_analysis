"""
Unit tests for the streaming module.
"""

import os
import tempfile
import pytest
import pandas as pd
from pathlib import Path
from unittest.mock import patch, MagicMock

from src.gaslit_af.streaming import StreamingVCFProcessor, stream_process_vcf
from src.gaslit_af.exceptions import FileError, DataProcessingError


class TestStreamingVCFProcessor:
    """Test cases for the StreamingVCFProcessor class."""
    
    def setup_method(self):
        """Set up test fixtures."""
        # Create a temporary VCF file for testing
        self.temp_dir = tempfile.TemporaryDirectory()
        self.vcf_path = os.path.join(self.temp_dir.name, "test.vcf")
        
        # Create a minimal VCF file
        with open(self.vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("chr1\t100\t.\tA\tG\t100\tPASS\tANN=G|missense_variant|MODERATE|GENE1|GENE1|transcript|NM_001|protein_coding|1/5|c.100A>G|p.Thr34Ala\n")
            f.write("chr1\t200\t.\tC\tT\t100\tPASS\tANN=T|missense_variant|MODERATE|GENE2|GENE2|transcript|NM_002|protein_coding|2/5|c.200C>T|p.Arg67Trp\n")
            f.write("chr2\t300\t.\tG\tA\t100\tPASS\tANN=A|missense_variant|MODERATE|GENE3|GENE3|transcript|NM_003|protein_coding|3/5|c.300G>A|p.Met100Ile\n")
        
        # Set up target genes
        self.target_genes = {"GENE1", "GENE3"}
    
    def teardown_method(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()
    
    def test_init(self):
        """Test initialization of StreamingVCFProcessor."""
        processor = StreamingVCFProcessor(
            vcf_path=self.vcf_path,
            target_genes=self.target_genes
        )
        
        assert processor.vcf_path == self.vcf_path
        assert processor.target_genes == self.target_genes
        assert processor.max_ram_usage == 64
        assert processor.ram_buffer == 16
        assert processor.max_workers == 16
        assert processor.chunk_size == 5000000  # Updated for 16GB Arc A770 GPU
    
    def test_init_with_invalid_path(self):
        """Test initialization with invalid VCF path."""
        with pytest.raises(FileError):
            StreamingVCFProcessor(
                vcf_path="invalid_path.vcf",
                target_genes=self.target_genes
            )
    
    @patch("src.gaslit_af.streaming.process_batch_with_data")
    def test_process_chunk(self, mock_process_batch):
        """Test processing a chunk of VCF records."""
        # Mock the process_batch_with_data function
        mock_process_batch.return_value = (["GENE1"], [{"GENE": "GENE1"}])
        
        processor = StreamingVCFProcessor(
            vcf_path=self.vcf_path,
            target_genes=self.target_genes
        )
        
        # Create a mock chunk
        chunk = [MagicMock(), MagicMock()]
        
        # Process the chunk
        genes_found, variant_data = processor._process_chunk(chunk)
        
        # Check that process_batch_with_data was called with the correct arguments
        mock_process_batch.assert_called_once_with(chunk, self.target_genes)
        
        # Check the results
        assert genes_found == ["GENE1"]
        assert variant_data == [{"GENE": "GENE1"}]
    
    @patch("src.gaslit_af.streaming.update_gene_counts")
    @patch("src.gaslit_af.streaming.get_memory_usage")
    def test_stream_process(self, mock_get_memory_usage, mock_update_gene_counts):
        """Test streaming processing of a VCF file."""
        # Mock memory usage
        mock_get_memory_usage.return_value = 1.0
        
        processor = StreamingVCFProcessor(
            vcf_path=self.vcf_path,
            target_genes=self.target_genes,
            initial_chunk_size=1  # Process one record at a time for testing
        )
        
        # Process the VCF file
        gene_counts, variant_df = processor.stream_process()
        
        # Check that update_gene_counts was called
        assert mock_update_gene_counts.call_count > 0
        
        # Check the results
        assert isinstance(gene_counts, dict)
        assert isinstance(variant_df, pd.DataFrame)
        
        # We should have found GENE1 and GENE3, but not GENE2
        assert len(variant_df) == 2
        assert set(variant_df["GENE"].unique()) == {"GENE1", "GENE3"}
    
    def test_estimate_total_records(self):
        """Test estimating the total number of records in a VCF file."""
        processor = StreamingVCFProcessor(
            vcf_path=self.vcf_path,
            target_genes=self.target_genes
        )
        
        # Estimate total records
        total_records = processor._estimate_total_records()
        
        # We have 3 records in our test VCF
        assert total_records == 3
    
    def test_create_record_chunks(self):
        """Test creating chunks of VCF records."""
        processor = StreamingVCFProcessor(
            vcf_path=self.vcf_path,
            target_genes=self.target_genes
        )
        
        # Create a mock iterator
        mock_iterator = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        
        # Create chunks of size 3
        chunks = list(processor._create_record_chunks(mock_iterator, 3))
        
        # We should have 4 chunks: [1,2,3], [4,5,6], [7,8,9], [10]
        assert len(chunks) == 4
        assert chunks[0] == [1, 2, 3]
        assert chunks[1] == [4, 5, 6]
        assert chunks[2] == [7, 8, 9]
        assert chunks[3] == [10]


def test_stream_process_vcf():
    """Test the stream_process_vcf function."""
    # Create a temporary VCF file for testing
    with tempfile.TemporaryDirectory() as temp_dir:
        vcf_path = os.path.join(temp_dir, "test.vcf")
        
        # Create a minimal VCF file
        with open(vcf_path, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##INFO=<ID=ANN,Number=.,Type=String,Description=\"Functional annotations\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("chr1\t100\t.\tA\tG\t100\tPASS\tANN=G|missense_variant|MODERATE|GENE1|GENE1|transcript|NM_001|protein_coding|1/5|c.100A>G|p.Thr34Ala\n")
            f.write("chr1\t200\t.\tC\tT\t100\tPASS\tANN=T|missense_variant|MODERATE|GENE2|GENE2|transcript|NM_002|protein_coding|2/5|c.200C>T|p.Arg67Trp\n")
        
        # Set up target genes
        target_genes = {"GENE1"}
        
        # Define a mock progress callback
        def progress_callback(processed, total):
            assert processed <= total
        
        # Process the VCF file
        gene_counts, variant_df = stream_process_vcf(
            vcf_path=vcf_path,
            target_genes=target_genes,
            max_ram_usage=64,
            ram_buffer=16,
            max_workers=4,
            initial_chunk_size=1,
            progress_callback=progress_callback
        )
        
        # Check the results
        assert isinstance(gene_counts, dict)
        assert isinstance(variant_df, pd.DataFrame)
        
        # We should have found GENE1, but not GENE2
        assert len(variant_df) == 1
        assert set(variant_df["GENE"].unique()) == {"GENE1"}


def test_stream_process_vcf_with_invalid_path():
    """Test stream_process_vcf with an invalid VCF path."""
    with pytest.raises(FileError):
        stream_process_vcf(
            vcf_path="invalid_path.vcf",
            target_genes={"GENE1"}
        )
