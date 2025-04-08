"""
Test configuration for GASLIT-AF Variant Analysis.
"""

import os
import sys
import pytest
from pathlib import Path

# Add project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Define fixtures that can be used across tests
@pytest.fixture(scope="session")
def test_data_dir():
    """Return the path to the test data directory."""
    return Path(__file__).parent / "data"

@pytest.fixture(scope="session")
def sample_vcf_path(test_data_dir):
    """Return the path to a sample VCF file for testing."""
    # For testing, prioritize using our sample test file
    vcf_path = test_data_dir / "sample.vcf"
    
    # If the sample test file doesn't exist, try to use the actual VCF file
    if not vcf_path.exists():
        vcf_path = Path("/home/k10/dev/windsage/gaslit-af_variant_analysis/data/KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.snp-indel.genome.vcf.gz")
        if not vcf_path.exists():
            pytest.skip(f"No test VCF file found. Neither the sample test file nor the real VCF exist.")
    
    return vcf_path

@pytest.fixture(scope="session")
def output_dir():
    """Return a temporary output directory for test results."""
    output_dir = Path("/home/k10/dev/windsage/gaslit-af_variant_analysis/tests/output")
    output_dir.mkdir(exist_ok=True, parents=True)
    return output_dir

@pytest.fixture(scope="session")
def cache_dir():
    """Return a temporary cache directory for testing."""
    cache_dir = Path("/home/k10/dev/windsage/gaslit-af_variant_analysis/tests/cache")
    cache_dir.mkdir(exist_ok=True, parents=True)
    return cache_dir

@pytest.fixture(scope="session")
def test_environment(test_data_dir, output_dir, cache_dir):
    """Set up a test environment with necessary directories and sample data."""
    # Create visualization directory
    viz_dir = output_dir / "visualizations"
    viz_dir.mkdir(exist_ok=True, parents=True)
    
    # Clean up any previous test files
    for file in output_dir.glob("*"):
        if file.is_file():
            file.unlink()
    
    for file in cache_dir.glob("*.cache"):
        file.unlink()
    
    # Return a dict with all the paths
    return {
        "test_data_dir": test_data_dir,
        "output_dir": output_dir,
        "cache_dir": cache_dir,
        "viz_dir": viz_dir
    }
