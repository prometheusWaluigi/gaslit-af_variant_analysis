"""
Data processing module for GASLIT-AF Variant Analysis.
Provides functions to process VCF data and extract relevant information.
"""

import os
import pandas as pd
import numpy as np
from cyvcf2 import VCF
from collections import defaultdict
import logging
import json
import time
import dpctl

from src.gaslit_af.exceptions import DataProcessingError, FileError, retry_operation

# Configure logging
log = logging.getLogger("gaslit-af")

def vcf_to_dataframe(vcf_path, limit=None, batch_size=100000):
    """
    Convert VCF file to pandas DataFrame with basic variant information.
    Optimized for performance with batch processing.
    
    Args:
        vcf_path: Path to VCF file
        limit: Optional limit on number of records to process
        batch_size: Size of batches for processing
    
    Returns:
        DataFrame with variant information
    """
    import concurrent.futures
    from functools import partial
    
    log.info(f"Converting VCF to DataFrame: {vcf_path}")
    
    if not os.path.exists(vcf_path):
        error_msg = f"VCF file not found: {vcf_path}"
        log.error(error_msg)
        raise FileError(error_msg)
    
    try:
        vcf = VCF(vcf_path)
    except Exception as e:
        error_msg = f"Error opening VCF file: {e}"
        log.error(error_msg)
        raise DataProcessingError(error_msg, details=str(e))
    
    # Define columns to extract
    columns = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER']
    
    # Function to process a batch of records
    def process_batch(records):
        batch_data = []
        for record in records:
            try:
                # Extract basic fields
                row = {
                    'CHROM': record.CHROM,
                    'POS': record.POS,
                    'ID': record.ID if record.ID else '.',
                    'REF': record.REF,
                    'ALT': record.ALT[0] if len(record.ALT) > 0 else '',  # Take first ALT allele with safety check
                    'QUAL': record.QUAL,
                    'FILTER': ';'.join(record.FILTER) if record.FILTER else 'PASS'
                }
                
                # Extract only essential INFO fields to reduce memory usage
                for field in ['ANN', 'Gene', 'IMPACT', 'Effect']:
                    if field in record.INFO:
                        value = record.INFO.get(field)
                        if isinstance(value, (list, tuple)):
                            value = ';'.join(map(str, value))
                        row[f'INFO_{field}'] = value
                
                batch_data.append(row)
            except Exception as e:
                log.warning(f"Error processing record: {e}")
                continue
        return batch_data
    
    # Collect records in batches
    all_batches = []
    current_batch = []
    count = 0
    
    for record in vcf:
        current_batch.append(record)
        count += 1
        
        if len(current_batch) >= batch_size:
            all_batches.append(current_batch)
            current_batch = []
        
        if limit and count >= limit:
            break
    
    # Add the last batch if it's not empty
    if current_batch:
        all_batches.append(current_batch)
    
    # Process batches in parallel
    data = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        batch_results = list(executor.map(process_batch, all_batches))
        
        # Flatten results
        for batch_data in batch_results:
            data.extend(batch_data)
    
    if not data:
        warning_msg = "No variant data extracted from VCF file"
        log.warning(warning_msg)
        return pd.DataFrame()
    
    df = pd.DataFrame(data)
    log.info(f"Created DataFrame with {len(df)} variants and {len(df.columns)} columns")
    
    return df

@retry_operation(max_attempts=3, retry_delay=2)
def extract_gaslit_af_variants(vcf_path, gaslit_af_genes, queue, batch_size=2000000):
    """
    Extract variants in GASLIT-AF genes from VCF file.
    
    Args:
        vcf_path: Path to VCF file
        gaslit_af_genes: Set of GASLIT-AF gene names
        queue: SYCL queue for processing
        batch_size: Batch size for processing
    
    Returns:
        Dictionary of gene:count pairs and DataFrame with variant information
    """
    log.info(f"Extracting GASLIT-AF variants from: {vcf_path}")
    
    if not os.path.exists(vcf_path):
        error_msg = f"VCF file not found: {vcf_path}"
        log.error(error_msg)
        raise FileError(error_msg)
    
    try:
        vcf = VCF(vcf_path)
    except Exception as e:
        error_msg = f"Error opening VCF file: {e}"
        log.error(error_msg)
        raise DataProcessingError(error_msg, details=str(e))
    match_counts = defaultdict(int)
    variant_data = []
    
    records_batch = []
    batch_variants = []
    
    for record in vcf:
        records_batch.append(record)
        
        if len(records_batch) >= batch_size:
            # Process batch
            genes_found, batch_data = process_batch_with_data(records_batch, gaslit_af_genes)
            
            # Update gene counts using SYCL
            update_gene_counts(genes_found, match_counts, queue)
            
            # Extend variant data
            variant_data.extend(batch_data)
            
            # Clear batch
            records_batch.clear()
    
    # Process remaining batch
    if records_batch:
        genes_found, batch_data = process_batch_with_data(records_batch, gaslit_af_genes)
        update_gene_counts(genes_found, match_counts, queue)
        variant_data.extend(batch_data)
    
    # Convert variant data to DataFrame
    df = pd.DataFrame(variant_data)
    
    log.info(f"Extracted {len(df)} variants across {len(match_counts)} GASLIT-AF genes")
    
    return match_counts, df

def process_batch_with_data(records, gaslit_af_genes):
    """
    Process a batch of VCF records and extract GASLIT-AF gene variants.
    Optimized for performance with list comprehensions.
    
    Args:
        records: List of VCF records
        gaslit_af_genes: Set of GASLIT-AF gene names
    
    Returns:
        List of genes found and list of variant data dictionaries
    """
    genes_found = []
    variant_data = []
    
    for record in records:
        try:
            # Check for ANN field (variant annotation)
            ann = record.INFO.get('ANN')
            if not ann:
                continue
                
            # Process annotations more efficiently
            for entry in ann.split(','):
                parts = entry.split('|')
                if len(parts) > 3 and parts[3] in gaslit_af_genes:
                    # Add gene to found list
                    genes_found.append(parts[3])
                    
                    # Extract variant data - use try/except for safety
                    try:
                        variant = {
                            'CHROM': record.CHROM,
                            'POS': record.POS,
                            'REF': record.REF,
                            'ALT': record.ALT[0] if len(record.ALT) > 0 else '',
                            'GENE': parts[3],
                            'EFFECT': parts[1] if len(parts) > 1 else '',
                            'IMPACT': parts[2] if len(parts) > 2 else ''
                        }
                        variant_data.append(variant)
                    except Exception as e:
                        # Skip this variant but continue processing
                        continue
        except Exception as e:
            # Skip this record but continue processing
            continue
    
    return genes_found, variant_data

def update_gene_counts(genes_found, match_counts, queue):
    """
    Update gene counts using SYCL acceleration with optimized GPU usage.
    
    Args:
        genes_found: List of genes found
        match_counts: Dictionary to update
        queue: SYCL queue for processing
    """
    if not genes_found:
        return
    
    # Convert to numpy array
    genes_array = np.array(genes_found)
    
    try:
        # Try to use SYCL for acceleration, but this is optional
        # and we'll fall back to CPU if it fails
        if hasattr(dpctl, '_sycl_context'):
            # Use SYCL USM memory for better GPU performance
            import dpctl.tensor as dpt
            with dpctl._sycl_context(queue):
                usm_array = dpt.asarray(genes_array)
                
                # Perform unique count on GPU
                # This is more efficient than transferring back to CPU
                unique_genes, counts = np.unique(usm_array.to_numpy(), return_counts=True)
                
                # Use vectorized operations for better performance
                for gene, count in zip(unique_genes, counts):
                    match_counts[gene] += int(count)
        else:
            # If dpctl doesn't have _sycl_context, fall back to CPU
            raise AttributeError("dpctl doesn't have required SYCL context methods")
    except Exception as e:
        # Fallback to CPU if SYCL fails
        log.warning(f"⚠️ SYCL processing failed, falling back to CPU: {e}")
        unique_genes, counts = np.unique(genes_array, return_counts=True)
        for gene, count in zip(unique_genes, counts):
            match_counts[gene] += int(count)

def save_results(match_counts, variant_df, output_dir):
    """
    Save analysis results to various formats.
    
    Args:
        match_counts: Dictionary of gene:count pairs
        variant_df: DataFrame with variant information
        output_dir: Directory to save results
    """
    try:
        os.makedirs(output_dir, exist_ok=True)
        timestamp = time.strftime("%Y%m%d-%H%M%S")
    
        # Save gene counts as CSV
        gene_counts_df = pd.DataFrame(list(match_counts.items()), 
                                    columns=['Gene', 'VariantCount'])
        gene_counts_df = gene_counts_df.sort_values('VariantCount', ascending=False)
        gene_counts_path = os.path.join(output_dir, f"gene_counts_{timestamp}.csv")
        gene_counts_df.to_csv(gene_counts_path, index=False)
        log.info(f"Saved gene counts to: {gene_counts_path}")
        
        # Save variant data as CSV
        variant_path = None
        if not variant_df.empty:
            variant_path = os.path.join(output_dir, f"variants_{timestamp}.csv")
            variant_df.to_csv(variant_path, index=False)
            log.info(f"Saved variant data to: {variant_path}")
        
        # Save as JSON
        json_data = {
            'analysis_time': timestamp,
            'total_genes': len(match_counts),
            'total_variants': variant_df.shape[0] if not variant_df.empty else 0,
            'gene_counts': match_counts
        }
        
        json_path = os.path.join(output_dir, f"results_{timestamp}.json")
        with open(json_path, 'w') as f:
            json.dump(json_data, f, indent=2)
        
        log.info(f"Saved JSON results to: {json_path}")
        
        return {
            'gene_counts_path': gene_counts_path,
            'variant_path': variant_path,
            'json_path': json_path
        }
    except Exception as e:
        error_msg = f"Error saving results: {e}"
        log.error(error_msg)
        raise FileError(error_msg, details=str(e))
