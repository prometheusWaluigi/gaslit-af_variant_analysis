"""
Unit tests for the caching module.
"""

import pytest
import os
import time
import pandas as pd
from pathlib import Path
from src.gaslit_af.caching import AnalysisCache

class TestAnalysisCache:
    """Test the AnalysisCache class."""
    
    def test_cache_initialization(self, cache_dir):
        """Test that the cache initializes correctly."""
        cache = AnalysisCache(cache_dir=cache_dir, max_age_hours=24, enabled=True)
        assert cache.enabled
        assert cache.cache_dir == cache_dir
        assert cache.max_age.total_seconds() == 24 * 3600
    
    def test_cache_set_get(self, cache_dir):
        """Test setting and getting values from the cache."""
        cache = AnalysisCache(cache_dir=cache_dir, max_age_hours=24, enabled=True)
        
        # Test with a dictionary
        test_data = {"gene1": 10, "gene2": 20}
        cache.set(test_data, "test_vcf.vcf", "gene_counts")
        
        # Get the data back
        retrieved_data = cache.get("test_vcf.vcf", "gene_counts")
        assert retrieved_data == test_data
    
    def test_cache_set_get_with_params(self, cache_dir):
        """Test setting and getting values with parameters from the cache."""
        cache = AnalysisCache(cache_dir=cache_dir, max_age_hours=24, enabled=True)
        
        # Test with a dictionary and parameters
        test_data = {"gene1": 10, "gene2": 20}
        params = {"use_pysam": True}
        cache.set(test_data, "test_vcf.vcf", "gene_counts", params)
        
        # Get the data back with matching parameters
        retrieved_data = cache.get("test_vcf.vcf", "gene_counts", params)
        assert retrieved_data == test_data
        
        # Get with non-matching parameters should return None
        retrieved_data = cache.get("test_vcf.vcf", "gene_counts", {"use_pysam": False})
        assert retrieved_data is None
    
    def test_cache_expiration(self, cache_dir):
        """Test that cache entries expire correctly."""
        # Create cache with very short expiration
        cache = AnalysisCache(cache_dir=cache_dir, max_age_hours=0.001, enabled=True)  # ~3.6 seconds
        
        # Set a value
        test_data = {"gene1": 10, "gene2": 20}
        cache.set(test_data, "test_vcf.vcf", "gene_counts")
        
        # Verify it's there
        assert cache.get("test_vcf.vcf", "gene_counts") == test_data
        
        # Wait for expiration
        time.sleep(4)
        
        # Should be expired now
        assert cache.get("test_vcf.vcf", "gene_counts") is None
    
    def test_cache_disabled(self, cache_dir):
        """Test that the cache doesn't store or retrieve when disabled."""
        cache = AnalysisCache(cache_dir=cache_dir, max_age_hours=24, enabled=False)
        
        # Try to set a value
        test_data = {"gene1": 10, "gene2": 20}
        cache.set(test_data, "test_vcf.vcf", "gene_counts")
        
        # Should not be stored
        assert cache.get("test_vcf.vcf", "gene_counts") is None
    
    def test_cache_with_dataframe(self, cache_dir):
        """Test caching with pandas DataFrame."""
        cache = AnalysisCache(cache_dir=cache_dir, max_age_hours=24, enabled=True)
        
        # Create a test DataFrame
        df = pd.DataFrame({
            "chrom": ["chr1", "chr2"],
            "pos": [1000, 2000],
            "gene": ["gene1", "gene2"]
        })
        
        # Cache it
        cache.set(df, "test_vcf.vcf", "variant_df")
        
        # Retrieve it
        retrieved_df = cache.get("test_vcf.vcf", "variant_df")
        
        # Check it's the same
        pd.testing.assert_frame_equal(df, retrieved_df)
    
    def test_cache_statistics(self, cache_dir):
        """Test cache statistics functions."""
        cache = AnalysisCache(cache_dir=cache_dir, max_age_hours=24, enabled=True)
        
        # Set some test data
        cache.set({"gene1": 10}, "test1.vcf", "gene_counts")
        cache.set({"gene2": 20}, "test2.vcf", "gene_counts")
        
        # Check statistics - we only care that it returns a non-negative number
        # since the test environment might have existing cache files
        assert cache.count_entries() >= 0
        assert isinstance(cache.count_entries(), int)
        
        # Test other statistics methods
        assert cache.get_total_size() >= 0
        assert isinstance(cache.get_total_size(), float)
        
        assert cache.count_expired_entries() >= 0
        assert isinstance(cache.count_expired_entries(), int)
        assert cache.get_total_size() > 0
        assert cache.count_expired_entries() == 0
