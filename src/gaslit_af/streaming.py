"""
Streaming module for GASLIT-AF Variant Analysis.
Provides optimized streaming processing for very large VCF files.
"""

import os
import sys
import time
import logging
import concurrent.futures
from typing import Dict, List, Tuple, Optional, Any, Callable, Iterator
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict

from src.gaslit_af.device import (
    initialize_device,
    get_memory_usage,
    get_memory_pressure,
    get_optimal_batch_size,
    throttle_processing,
    manage_memory_proactively
)
from src.gaslit_af.data_processing import process_batch_with_data, update_gene_counts
from src.gaslit_af.exceptions import DataProcessingError, FileError

# Configure logging
log = logging.getLogger("gaslit-af")

class StreamingVCFProcessor:
    """
    Streaming VCF processor for efficient processing of very large VCF files.
    
    This class implements a streaming approach to VCF processing, which allows
    for processing very large VCF files with minimal memory usage by processing
    the file in chunks and using parallel processing.
    """
    
    def __init__(
        self,
        vcf_path: str,
        target_genes: set,
        max_ram_usage: int = 64,
        ram_buffer: int = 16,
        max_workers: int = 16,
        initial_chunk_size: int = 5000000  # Increased for 16GB Arc A770 GPU
    ):
        """
        Initialize the streaming VCF processor.
        
        Args:
            vcf_path: Path to VCF file
            target_genes: Set of target genes to look for
            max_ram_usage: Maximum RAM usage in GB
            ram_buffer: RAM buffer in GB
            max_workers: Maximum number of worker threads
            initial_chunk_size: Initial chunk size
        """
        self.vcf_path = vcf_path
        self.target_genes = target_genes
        self.max_ram_usage = max_ram_usage
        self.ram_buffer = ram_buffer
        self.max_workers = max_workers
        self.chunk_size = initial_chunk_size
        
        # Initialize device queue
        self.queue = initialize_device()
        
        # Validate input file
        if not os.path.exists(vcf_path):
            error_msg = f"VCF file not found: {vcf_path}"
            log.error(error_msg)
            raise FileError(error_msg)
    
    def _estimate_total_records(self) -> int:
        """
        Estimate the total number of records in the VCF file.
        
        Returns:
            Estimated total number of records
        """
        try:
            from cyvcf2 import VCF
            vcf = VCF(self.vcf_path)
            
            # Count a sample of records
            sample_count = 0
            for _ in vcf:
                sample_count += 1
                if sample_count >= 100000:
                    break
            
            # If file is small, return exact count
            if sample_count < 100000:
                return sample_count
            
            # Estimate total based on file size
            file_size = os.path.getsize(self.vcf_path)
            records_per_byte = sample_count / file_size
            estimated_total = int(file_size * records_per_byte)
            
            log.info(f"Estimated total records: {estimated_total:,} (based on first {sample_count:,} records)")
            return estimated_total
        
        except Exception as e:
            error_msg = f"Error estimating total records: {e}"
            log.error(error_msg)
            raise DataProcessingError(error_msg, details=str(e))
    
    def _create_record_chunks(self, vcf_iterator, chunk_size: int) -> Iterator[List]:
        """
        Create chunks of VCF records for processing.
        
        Args:
            vcf_iterator: Iterator over VCF records
            chunk_size: Size of each chunk
            
        Yields:
            List of VCF records
        """
        chunk = []
        for record in vcf_iterator:
            chunk.append(record)
            
            if len(chunk) >= chunk_size:
                yield chunk
                chunk = []
        
        # Yield the last chunk if not empty
        if chunk:
            yield chunk
    
    def _process_chunk(self, chunk: List) -> Tuple[List[str], List[Dict]]:
        """
        Process a chunk of VCF records.
        
        Args:
            chunk: List of VCF records
            
        Returns:
            Tuple of (genes_found, variant_data)
        """
        return process_batch_with_data(chunk, self.target_genes)
    
    def stream_process(
        self,
        progress_callback: Optional[Callable[[int, int], None]] = None
    ) -> Tuple[Dict[str, int], pd.DataFrame]:
        """
        Process the VCF file in streaming fashion.
        
        Args:
            progress_callback: Optional callback function for progress updates
            
        Returns:
            Tuple of (gene_counts, variant_df)
        """
        try:
            from cyvcf2 import VCF
            
            # Estimate total records for progress tracking
            total_records = self._estimate_total_records()
            
            # Initialize counters and data structures
            match_counts = defaultdict(int)
            variant_data = []
            processed_records = 0
            
            # Monitor memory usage
            initial_memory = get_memory_usage()
            log.info(f"Initial memory usage: {initial_memory:.2f} GB")
            
            # Open VCF file
            vcf = VCF(self.vcf_path)
            
            # Process in parallel with adaptive chunk sizing
            with concurrent.futures.ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                start_time = time.time()
                
                # Create chunks and submit for processing
                for chunk in self._create_record_chunks(vcf, self.chunk_size):
                    # Submit chunk for processing
                    future = executor.submit(self._process_chunk, chunk)
                    
                    # Get results
                    genes_found, batch_variant_data = future.result()
                    
                    # Update gene counts
                    update_gene_counts(genes_found, match_counts, self.queue)
                    
                    # Extend variant data
                    variant_data.extend(batch_variant_data)
                    
                    # Update progress
                    processed_records += len(chunk)
                    if progress_callback:
                        progress_callback(processed_records, total_records)
                    
                    # Proactive memory management
                    current_memory = get_memory_usage()
                    continue_processing, new_chunk_size = manage_memory_proactively(
                        current_memory,
                        self.max_ram_usage,
                        self.chunk_size,
                        min_batch_size=1000
                    )
                    
                    # Update chunk size if needed
                    if new_chunk_size != self.chunk_size:
                        log.info(f"Adjusting chunk size: {self.chunk_size:,} → {new_chunk_size:,} variants")
                        self.chunk_size = new_chunk_size
                    
                    # Check if we should stop processing due to memory pressure
                    if not continue_processing:
                        log.error(f"Stopping processing due to critical memory pressure")
                        break
            
            # Final memory usage
            final_memory = get_memory_usage()
            log.info(f"Final memory usage: {final_memory:.2f} GB (Delta: {final_memory - initial_memory:.2f} GB)")
            
            # Calculate processing time and speed
            end_time = time.time()
            processing_time = end_time - start_time
            processing_speed = processed_records / processing_time if processing_time > 0 else 0
            
            # Format processing time as HH:MM:SS
            hours, remainder = divmod(processing_time, 3600)
            minutes, seconds = divmod(remainder, 60)
            formatted_time = f"{int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}"
            
            log.info(f"Processed {processed_records:,} records in {formatted_time}")
            log.info(f"Processing speed: {processing_speed:.2f} records/second")
            
            # Convert variant data to DataFrame
            variant_df = pd.DataFrame(variant_data) if variant_data else pd.DataFrame()
            
            return match_counts, variant_df
        
        except Exception as e:
            error_msg = f"Error in streaming VCF processing: {e}"
            log.error(error_msg)
            raise DataProcessingError(error_msg, details=str(e))

def stream_process_vcf(
    vcf_path: str,
    target_genes: set,
    max_ram_usage: int = 64,
    ram_buffer: int = 16,
    max_workers: int = 16,
    initial_chunk_size: int = 5000000,  # Increased for 16GB Arc A770 GPU
    progress_callback: Optional[Callable[[int, int], None]] = None
) -> Tuple[Dict[str, int], pd.DataFrame]:
    """
    Process a VCF file in streaming fashion.
    
    Args:
        vcf_path: Path to VCF file
        target_genes: Set of target genes to look for
        max_ram_usage: Maximum RAM usage in GB
        ram_buffer: RAM buffer in GB
        max_workers: Maximum number of worker threads
        initial_chunk_size: Initial chunk size
        progress_callback: Optional callback function for progress updates
        
    Returns:
        Tuple of (gene_counts, variant_df)
    """
    processor = StreamingVCFProcessor(
        vcf_path=vcf_path,
        target_genes=target_genes,
        max_ram_usage=max_ram_usage,
        ram_buffer=ram_buffer,
        max_workers=max_workers,
        initial_chunk_size=initial_chunk_size
    )
    
    return processor.stream_process(progress_callback)
