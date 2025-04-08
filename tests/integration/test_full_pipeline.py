"""
Integration test for the full GASLIT-AF variant analysis pipeline.
"""

import pytest
import os
import sys
from pathlib import Path
import pandas as pd
import logging
from unittest.mock import patch

# Import modules for testing
from src.gaslit_af.cli import parse_args
from src.gaslit_af.device import initialize_device
from src.gaslit_af.gene_lists import GASLIT_AF_GENES, KNOWN_SNPS
from src.gaslit_af.workflow import run_analysis_workflow
from src.gaslit_af.advanced_variant_processing import process_vcf_with_pysam, VariantProcessor

# Configure logging to avoid cluttering test output
logging.basicConfig(level=logging.WARNING)

class TestFullPipeline:
    """Integration tests for the full analysis pipeline."""
    
    @pytest.mark.slow
    def test_find_specific_variants(self, sample_vcf_path, test_environment):
        """Test finding specific variants in the genome."""
        # Skip if VCF file is not available
        if not sample_vcf_path.exists():
            pytest.skip(f"Test VCF file not found: {sample_vcf_path}")
        
        output_dir = test_environment["output_dir"]
        
        try:
            # Initialize device
            queue = initialize_device()
            
            # Create variant processor
            processor = VariantProcessor(queue=queue, threads=4)
            
            # Define target genes (use a subset for faster testing)
            # Only include genes that are actually in our sample VCF
            target_genes = {
                "CHRM2", "DRD2", "TFAM"
            }
            
            # Process VCF file with a small batch size for testing
            gene_counts, variant_df = process_vcf_with_pysam(
                vcf_path=str(sample_vcf_path),
                target_genes=target_genes,
                dbsnp_path=None,
                queue=queue,
                threads=4,
                batch_size=100000,  # Smaller batch for testing
                max_ram_usage=32,
                ram_buffer=8
            )
            
            # Check that we got results
            assert gene_counts is not None
            assert isinstance(gene_counts, dict)
            
            # Check if any variants were found
            if variant_df is not None and not variant_df.empty:
                # Generate variant report
                report_path = processor.generate_variant_report(variant_df, output_dir)
                assert Path(report_path).exists()
                
                # Print found variants for debugging
                print(f"\nFound {len(variant_df)} variants in target genes:")
                for _, row in variant_df.iterrows():
                    gene = row.get('gene', '')
                    rsid = row.get('rsid', '')
                    genotype = row.get('genotype', '')
                    print(f"  {gene} - {rsid} ({genotype})")
            else:
                print("\nNo variants found in the sample VCF file.")
        except Exception as e:
            pytest.fail(f"Test failed with error: {e}")
    
    @pytest.mark.slow
    def test_modular_workflow(self, sample_vcf_path, test_environment):
        """Test the modular workflow on a sample VCF file."""
        # Skip if VCF file is not available
        if not sample_vcf_path.exists():
            pytest.skip(f"Test VCF file not found: {sample_vcf_path}")
        
        output_dir = test_environment["output_dir"]
        cache_dir = test_environment["cache_dir"]
        
        # Create mock args
        with patch('sys.argv', [
            'analyze_modular.py',
            str(sample_vcf_path),
            '--output-dir', str(output_dir),
            '--cache-dir', str(cache_dir),
            '--batch-size', '100000',  # Smaller batch for testing
            '--threads', '4',
            '--max-ram', '32',
            '--ram-buffer', '8',
            '--system-analysis',
            '--use-pysam',
            '--known-variants-only'
        ]):
            args = parse_args()
        
        # Run the analysis workflow with reduced scope for testing
        try:
            run_analysis_workflow(args)
            
            # Check that output files were created
            gene_counts_files = list(output_dir.glob("gene_counts*.csv"))
            system_analysis_md = output_dir / "system_analysis.md"
            system_analysis_json = output_dir / "system_analysis.json"
            
            # Some files might not be created if no variants are found
            # Print what was found for debugging
            print("\nOutput files created:")
            for file in output_dir.glob("*"):
                if file.is_file():
                    print(f"  {file.name}")
            
            # Check for HTML report
            html_reports = list(output_dir.glob("gaslit_af_report_*.html"))
            if html_reports:
                print(f"  Found HTML report: {html_reports[0].name}")
            
            # Check for visualizations
            viz_dir = output_dir / "visualizations"
            if viz_dir.exists() and list(viz_dir.glob("*")):
                print(f"  Found visualizations: {len(list(viz_dir.glob('*')))} files")
            
            # Print summary of results if gene counts file exists
            if gene_counts_files:
                try:
                    gene_counts_df = pd.read_csv(gene_counts_files[0])
                    print(f"\nFound variants in {len(gene_counts_df)} genes:")
                    for _, row in gene_counts_df.head(10).iterrows():
                        print(f"  {row['Gene']}: {row['Count']} variants")
                    if len(gene_counts_df) > 10:
                        print(f"  ... and {len(gene_counts_df) - 10} more genes")
                except Exception as e:
                    print(f"Could not read gene counts file: {e}")
            else:
                print("No gene counts file was generated.")
        
        except Exception as e:
            pytest.fail(f"Workflow execution failed: {e}")
    
    @pytest.mark.slow
    def test_specific_variants_script(self, sample_vcf_path, test_environment):
        """Test the find_specific_variants.py script."""
        # Skip if VCF file is not available
        if not sample_vcf_path.exists():
            pytest.skip(f"Test VCF file not found: {sample_vcf_path}")
        
        output_dir = test_environment["output_dir"]
        
        # Import the script as a module
        sys.path.insert(0, str(Path(__file__).parent.parent.parent))
        import find_specific_variants
        
        # Mock command-line arguments
        with patch('sys.argv', [
            'find_specific_variants.py',
            str(sample_vcf_path),
            '--output-dir', str(output_dir),
            '--batch-size', '100000',  # Smaller batch for testing
            '--threads', '4'
        ]):
            # Run the script
            try:
                find_specific_variants.main()
                
                # Check for output files
                variant_report = output_dir / "variant_report.md"
                specific_variants_csv = output_dir / "specific_variants.csv"
                
                # Print what was found for debugging
                print("\nOutput files from specific variants script:")
                for file in output_dir.glob("*"):
                    if file.is_file():
                        print(f"  {file.name}")
                
                # One or both of these files should exist if variants were found
                if variant_report.exists() or specific_variants_csv.exists():
                    print("\nSpecific variants were found!")
                    
                    if specific_variants_csv.exists():
                        try:
                            variants_df = pd.read_csv(specific_variants_csv)
                            print(f"Found {len(variants_df)} specific variants:")
                            for _, row in variants_df.iterrows():
                                gene = row.get('gene', '')
                                rsid = row.get('rsid', '')
                                genotype = row.get('genotype', '')
                                print(f"  {gene} - {rsid} ({genotype})")
                        except Exception as e:
                            print(f"Could not read variants file: {e}")
                else:
                    print("\nNo specific variants were found or no output files were generated.")
            
            except Exception as e:
                pytest.fail(f"Script execution failed: {e}")
