"""
Unit tests for the CLI module.
"""

import pytest
import sys
from pathlib import Path
from unittest.mock import patch
from src.gaslit_af.cli import parse_args

class TestCLI:
    """Test the CLI module."""
    
    def test_parse_args_minimal(self):
        """Test parse_args with minimal arguments."""
        # Mock sys.argv with minimal arguments
        with patch('sys.argv', ['analyze.py', 'test.vcf']):
            args = parse_args()
            
            # Check that required arguments are set correctly
            assert args.vcf_path == 'test.vcf'
            
            # Check that optional arguments have default values
            assert args.output_dir == Path('output')
            assert args.batch_size == 2000000
            assert args.max_ram == 64
            assert args.ram_buffer == 16
            assert args.threads == 16
            assert args.cache_dir == Path('./cache')
            assert args.cache_max_age == 24
            assert not args.no_cache
            assert not args.system_analysis
            assert not args.use_pysam
            assert not args.known_variants_only
            assert not args.no_visualization
            assert not args.no_report
            assert not args.enhanced_report
            assert not args.open_browser
    
    def test_parse_args_full(self):
        """Test parse_args with all arguments."""
        # Mock sys.argv with all arguments
        with patch('sys.argv', [
            'analyze.py',
            'test.vcf',
            '--output-dir', 'custom_output',
            '--batch-size', '1000000',
            '--max-ram', '32',
            '--ram-buffer', '8',
            '--threads', '8',
            '--sample-limit', '50000',
            '--cache-dir', 'custom_cache',
            '--cache-max-age', '12',
            '--no-cache',
            '--system-analysis',
            '--use-pysam',
            '--dbsnp-path', 'dbsnp.vcf',
            '--known-variants-only',
            '--no-visualization',
            '--no-report',
            '--enhanced-report',
            '--open-browser'
        ]):
            args = parse_args()
            
            # Check that all arguments are set correctly
            assert args.vcf_path == 'test.vcf'
            assert args.output_dir == Path('custom_output')
            assert args.batch_size == 1000000
            assert args.max_ram == 32
            assert args.ram_buffer == 8
            assert args.threads == 8
            assert args.sample_limit == 50000
            assert args.cache_dir == Path('custom_cache')
            assert args.cache_max_age == 12
            assert args.no_cache
            assert args.system_analysis
            assert args.use_pysam
            assert args.dbsnp_path == 'dbsnp.vcf'
            assert args.known_variants_only
            assert args.no_visualization
            assert args.no_report
            assert args.enhanced_report
            assert args.open_browser
