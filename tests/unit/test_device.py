"""
Unit tests for the device module.
"""

import pytest
import os
import sys
from src.gaslit_af.device import get_memory_usage, check_memory_limits

class TestDeviceModule:
    """Test the device module functions."""
    
    def test_get_memory_usage(self):
        """Test that get_memory_usage returns a valid memory usage value."""
        memory_usage = get_memory_usage()
        assert isinstance(memory_usage, float)
        assert memory_usage > 0.0
    
    def test_check_memory_limits_within_limits(self):
        """Test check_memory_limits when memory usage is within limits."""
        # Set current usage to 10GB, max to 64GB, buffer to 16GB
        result = check_memory_limits(10.0, 64, 16)
        assert result is True
    
    def test_check_memory_limits_approaching_limits(self):
        """Test check_memory_limits when memory usage is approaching limits."""
        # Set current usage to 50GB, max to 64GB, buffer to 16GB
        # This should be approaching limits (50 > 64-16)
        result = check_memory_limits(50.0, 64, 16)
        assert result is False
    
    def test_check_memory_limits_exceeding_limits(self):
        """Test check_memory_limits when memory usage exceeds limits."""
        # Set current usage to 65GB, max to 64GB, buffer to 16GB
        result = check_memory_limits(65.0, 64, 16)
        assert result is False
