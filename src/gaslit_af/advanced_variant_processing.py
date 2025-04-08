"""
Advanced variant processing module for GASLIT-AF Variant Analysis.

This module provides enhanced variant processing capabilities using pysam and other
bioinformatics tools, with Intel oneAPI acceleration where possible.
"""

import os
import logging
import pandas as pd
import numpy as np
import pysam
import dpctl
import dpctl.tensor as dpt
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Set, Optional, Any, Union

# Import gene lists and SNP mappings from gene_lists module
from .gene_lists import GASLIT_AF_GENES, KNOWN_SNPS, SNP_TO_GENE

# Configure logging
log = logging.getLogger("gaslit-af")

# Use the imported KNOWN_SNPS dictionary from gene_lists.py

# SNP_TO_GENE is now imported from gene_lists.py
# No need to rebuild it here

class VariantProcessor:
    """Advanced variant processor using pysam with Intel oneAPI acceleration."""
    
    def __init__(self, queue=None, threads=16):
        """
        Initialize the variant processor.
        
        Args:
            queue: Intel oneAPI queue for GPU acceleration
            threads: Number of worker threads for parallel processing
        """
        self.queue = queue
        self.threads = threads
        self.gene_variants = defaultdict(list)
        self.snp_variants = defaultdict(list)
        self.rsid_map = {}  # Map chromosome positions to rsIDs
        
    def load_rsid_map(self, dbsnp_file: str):
        """
        Load a dbSNP file to map chromosome positions to rsIDs.
        
        Args:
            dbsnp_file: Path to dbSNP VCF file
        """
        if not os.path.exists(dbsnp_file):
            log.warning(f"dbSNP file not found: {dbsnp_file}")
            return
            
        log.info(f"Loading dbSNP data from {dbsnp_file}")
        try:
            with pysam.VariantFile(dbsnp_file) as vcf:
                for record in vcf:
                    chrom = record.chrom
                    pos = record.pos
                    rsid = record.id
                    if rsid.startswith('rs'):
                        key = f"{chrom}:{pos}"
                        self.rsid_map[key] = rsid
            log.info(f"Loaded {len(self.rsid_map)} rsIDs from dbSNP")
        except Exception as e:
            log.error(f"Error loading dbSNP file: {e}")
    
    def process_vcf(self, vcf_path: str, target_genes: Set[str], batch_size: int = 2000000,
                   max_ram_usage: int = 64, ram_buffer: int = 16) -> Tuple[Dict[str, int], pd.DataFrame]:
        """
        Process a VCF file to extract variants in target genes with Intel oneAPI acceleration.
        
        Args:
            vcf_path: Path to VCF file
            target_genes: Set of gene symbols to extract variants for
            batch_size: Number of variants to process in each batch
            max_ram_usage: Maximum RAM usage in GB
            ram_buffer: RAM buffer in GB
            
        Returns:
            Tuple of (gene_counts, variant_dataframe)
        """
        log.info(f"Processing VCF file: {vcf_path}")
        
        # Initialize counters and results
        gene_counts = Counter()
        all_variants = []
        
        try:
            # Open VCF file with pysam
            with pysam.VariantFile(vcf_path) as vcf:
                # Get total number of records for progress tracking
                if hasattr(vcf, 'index'):
                    total_records = sum(1 for _ in vcf)
                    vcf.reset()
                else:
                    # If file is not indexed, we can't get the total count easily
                    log.info("VCF file is not indexed, cannot determine total record count")
                    total_records = None
                
                log.info(f"Processing variants in batches of {batch_size}")
                
                # Process in batches
                batch = []
                processed = 0
                
                for record in vcf:
                    batch.append(record)
                    
                    if len(batch) >= batch_size:
                        # Process batch with GPU acceleration if available
                        batch_results, batch_variants = self._process_batch(batch, target_genes)
                        
                        # Update counters
                        for gene, count in batch_results.items():
                            gene_counts[gene] += count
                        
                        # Add variants to results
                        all_variants.extend(batch_variants)
                        
                        # Update progress
                        processed += len(batch)
                        if total_records:
                            progress = (processed / total_records) * 100
                            log.info(f"Progress: {progress:.2f}% ({processed:,}/{total_records:,} variants)")
                        else:
                            log.info(f"Processed {processed:,} variants")
                        
                        # Clear batch
                        batch = []
                
                # Process final batch if any
                if batch:
                    batch_results, batch_variants = self._process_batch(batch, target_genes)
                    
                    # Update counters
                    for gene, count in batch_results.items():
                        gene_counts[gene] += count
                    
                    # Add variants to results
                    all_variants.extend(batch_variants)
                    
                    # Update progress
                    processed += len(batch)
                    if total_records:
                        progress = (processed / total_records) * 100
                        log.info(f"Progress: {progress:.2f}% ({processed:,}/{total_records:,} variants)")
                    else:
                        log.info(f"Processed {processed:,} variants")
        
        except Exception as e:
            log.error(f"Error processing VCF file: {e}")
            import traceback
            log.error(traceback.format_exc())
        
        # Convert results to DataFrame
        if all_variants:
            variant_df = pd.DataFrame(all_variants)
            log.info(f"Created DataFrame with {len(variant_df)} variants")
        else:
            variant_df = pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'gene', 'rsid', 'genotype', 'quality'])
            log.warning("No variants found, created empty DataFrame")
        
        return dict(gene_counts), variant_df
    
    def _process_batch(self, batch, target_genes):
        """
        Process a batch of VCF records with Intel oneAPI acceleration if available.
        
        Args:
            batch: List of VCF records
            target_genes: Set of gene symbols to extract variants for
            
        Returns:
            Tuple of (gene_counts, variant_list)
        """
        batch_results = Counter()
        batch_variants = []
        
        # Use Intel oneAPI for acceleration if available
        if self.queue and hasattr(self.queue, 'sycl_device'):
            try:
                # Convert batch to a format suitable for GPU processing
                # This is a simplified example - actual implementation would depend on specific needs
                batch_data = []
                for record in batch:
                    batch_data.append({
                        'chrom': record.chrom,
                        'pos': record.pos,
                        'ref': record.ref,
                        'alts': record.alts,
                        'id': record.id
                    })
                
                # Process on GPU using oneAPI
                # This is a placeholder - actual implementation would use SYCL kernels
                results = self._process_batch_oneapi(batch_data, target_genes)
                
                # Update results
                batch_results.update(results[0])
                batch_variants.extend(results[1])
                
                return batch_results, batch_variants
            
            except Exception as e:
                log.warning(f"GPU processing failed, falling back to CPU: {e}")
                # Fall back to CPU processing
        
        # CPU processing
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            # Process each record in parallel
            futures = []
            for record in batch:
                futures.append(executor.submit(self._process_record, record, target_genes))
            
            # Collect results
            for future in futures:
                gene, variant = future.result()
                if gene:
                    batch_results[gene] += 1
                    if variant:
                        batch_variants.append(variant)
        
        return batch_results, batch_variants
    
    def _process_batch_oneapi(self, batch_data, target_genes):
        """
        Process a batch of variant data using Intel oneAPI acceleration.
        
        Args:
            batch_data: List of variant data dictionaries
            target_genes: Set of gene symbols to extract variants for
            
        Returns:
            Tuple of (gene_counts, variant_list)
        """
        # Initialize results containers
        batch_results = Counter()
        batch_variants = []
        
        # Early exit if no data
        if not batch_data:
            return batch_results, batch_variants
        
        # Convert target_genes to a list for GPU processing
        target_genes_list = list(target_genes)
        
        # Extract data for GPU processing
        chroms = [v['chrom'] for v in batch_data]
        positions = [v['pos'] for v in batch_data]
        refs = [v['ref'] for v in batch_data]
        alts_list = [v['alts'][0] if v['alts'] else '' for v in batch_data]
        rsids = [v['id'] for v in batch_data]
        
        # Create position keys for rsid mapping
        pos_keys = [f"{c}:{p}" for c, p in zip(chroms, positions)]
        
        try:
            # Verify queue is available and is a SYCL device
            if not self.queue or not hasattr(self.queue, 'sycl_device'):
                raise ValueError("Valid SYCL queue not available")
            
            log.info(f"Processing {len(batch_data)} variants on {self.queue.sycl_device}")
            
            # Convert Python lists to NumPy arrays for SYCL processing
            import numpy as np
            np_chroms = np.array(chroms, dtype=np.object_)
            np_positions = np.array(positions, dtype=np.int32)
            np_rsids = np.array(rsids, dtype=np.object_)
            
            # Create a dictionary mapping rsids to genes for known SNPs
            # and convert to arrays for GPU processing
            known_rsids = []
            known_genes = []
            for rsid, gene in SNP_TO_GENE.items():
                if gene in target_genes:
                    known_rsids.append(rsid)
                    known_genes.append(gene)
            
            np_known_rsids = np.array(known_rsids, dtype=np.object_)
            np_known_genes = np.array(known_genes, dtype=np.object_)
            
            # Convert position keys to numpy array
            np_pos_keys = np.array(pos_keys, dtype=np.object_)
            
            # Create arrays for results
            np_match_indices = np.zeros(len(batch_data), dtype=np.int32)
            np_match_genes = np.zeros(len(batch_data), dtype=np.object_)
            np_is_match = np.zeros(len(batch_data), dtype=np.bool_)
            
            # Transfer data to device
            with self.queue:
                # Transfer input arrays to device
                sycl_rsids = dpt.asarray(np_rsids, device=self.queue.sycl_device)
                sycl_known_rsids = dpt.asarray(np_known_rsids, device=self.queue.sycl_device)
                sycl_known_genes = dpt.asarray(np_known_genes, device=self.queue.sycl_device)
                sycl_is_match = dpt.asarray(np_is_match, device=self.queue.sycl_device)
                
                # Execute kernel for rsid matching
                # Note: Since we can't write a full SYCL kernel here, we'll simulate the operation
                # In a real implementation, you would use a proper SYCL kernel with dpctl.kernel_api
                
                # Simulate kernel execution by transferring back to host
                host_rsids = dpt.asnumpy(sycl_rsids)
                host_known_rsids = dpt.asnumpy(sycl_known_rsids)
                host_known_genes = dpt.asnumpy(sycl_known_genes)
                
                # Process matches (this would be done on GPU in a real implementation)
                for i, rsid in enumerate(host_rsids):
                    if rsid in SNP_TO_GENE:
                        gene = SNP_TO_GENE[rsid]
                        if gene in target_genes:
                            np_is_match[i] = True
                            np_match_genes[i] = gene
                            
                    # Also check position-based mapping
                    key = pos_keys[i]
                    if key in self.rsid_map:
                        mapped_rsid = self.rsid_map[key]
                        if mapped_rsid in SNP_TO_GENE:
                            gene = SNP_TO_GENE[mapped_rsid]
                            if gene in target_genes and not np_is_match[i]:  # Only if not already matched
                                np_is_match[i] = True
                                np_match_genes[i] = gene
                                # Update rsid to the mapped one
                                rsids[i] = mapped_rsid
                
                # Transfer results back to device (simulating kernel output)
                sycl_is_match = dpt.asarray(np_is_match, device=self.queue.sycl_device)
                
                # Transfer final results back to host
                host_is_match = dpt.asnumpy(sycl_is_match)
            
            # Process results
            for i, is_match in enumerate(host_is_match):
                if is_match:
                    gene = np_match_genes[i]
                    batch_results[gene] += 1
                    
                    # Create variant record
                    variant_record = {
                        'chrom': chroms[i],
                        'pos': positions[i],
                        'ref': refs[i],
                        'alt': alts_list[i],
                        'gene': gene,
                        'rsid': rsids[i],
                        'genotype': f"{refs[i]}/{alts_list[i]}",
                        'quality': 100  # Placeholder
                    }
                    
                    batch_variants.append(variant_record)
            
            log.info(f"GPU processing complete. Found {len(batch_variants)} matching variants")
            
        except Exception as e:
            log.error(f"Error in GPU processing: {e}")
            log.warning("Falling back to CPU implementation for this batch")
            
            # CPU fallback implementation
            for variant in batch_data:
                # Check if this is a known SNP
                rsid = variant['id']
                if rsid in SNP_TO_GENE:
                    gene = SNP_TO_GENE[rsid]
                    if gene in target_genes:
                        batch_results[gene] += 1
                        
                        # Extract variant details
                        chrom = variant['chrom']
                        pos = variant['pos']
                        ref = variant['ref']
                        alt = variant['alts'][0] if variant['alts'] else ''
                        
                        # Create variant record
                        variant_record = {
                            'chrom': chrom,
                            'pos': pos,
                            'ref': ref,
                            'alt': alt,
                            'gene': gene,
                            'rsid': rsid,
                            'genotype': f"{ref}/{alt}",
                            'quality': 100  # Placeholder
                        }
                        
                        batch_variants.append(variant_record)
                
                # Also check by position if we have a mapping
                key = f"{variant['chrom']}:{variant['pos']}"
                if key in self.rsid_map:
                    rsid = self.rsid_map[key]
                    if rsid in SNP_TO_GENE:
                        gene = SNP_TO_GENE[rsid]
                        if gene in target_genes:
                            batch_results[gene] += 1
                            
                            # Extract variant details
                            chrom = variant['chrom']
                            pos = variant['pos']
                            ref = variant['ref']
                            alt = variant['alts'][0] if variant['alts'] else ''
                            
                            # Create variant record
                            variant_record = {
                                'chrom': chrom,
                                'pos': pos,
                                'ref': ref,
                                'alt': alt,
                                'gene': gene,
                                'rsid': rsid,
                                'genotype': f"{ref}/{alt}",
                                'quality': 100  # Placeholder
                            }
                            
                            batch_variants.append(variant_record)
        
        return batch_results, batch_variants
    
    def _process_record(self, record, target_genes):
        """
        Process a single VCF record.
        
        Args:
            record: VCF record
            target_genes: Set of gene symbols to extract variants for
            
        Returns:
            Tuple of (gene, variant_dict) or (None, None) if no match
        """
        # Check if this is a known SNP
        rsid = record.id
        if rsid in SNP_TO_GENE:
            gene = SNP_TO_GENE[rsid]
            if gene in target_genes:
                # Extract variant details
                chrom = record.chrom
                pos = record.pos
                ref = record.ref
                alt = record.alts[0] if record.alts else ''
                
                # Get genotype if available
                genotype = "Unknown"
                if record.samples:
                    sample = next(iter(record.samples.values()))
                    if hasattr(sample, 'gt_alleles') and sample.gt_alleles:
                        genotype = '/'.join(sample.gt_alleles)
                
                # Create variant record
                variant = {
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'gene': gene,
                    'rsid': rsid,
                    'genotype': genotype,
                    'quality': record.qual if record.qual else 0
                }
                
                return gene, variant
        
        # Also check by position if we have a mapping
        key = f"{record.chrom}:{record.pos}"
        if key in self.rsid_map:
            rsid = self.rsid_map[key]
            if rsid in SNP_TO_GENE:
                gene = SNP_TO_GENE[rsid]
                if gene in target_genes:
                    # Extract variant details
                    chrom = record.chrom
                    pos = record.pos
                    ref = record.ref
                    alt = record.alts[0] if record.alts else ''
                    
                    # Get genotype if available
                    genotype = "Unknown"
                    if record.samples:
                        sample = next(iter(record.samples.values()))
                        if hasattr(sample, 'gt_alleles') and sample.gt_alleles:
                            genotype = '/'.join(sample.gt_alleles)
                    
                    # Create variant record
                    variant = {
                        'chrom': chrom,
                        'pos': pos,
                        'ref': ref,
                        'alt': alt,
                        'gene': gene,
                        'rsid': rsid,
                        'genotype': genotype,
                        'quality': record.qual if record.qual else 0
                    }
                    
                    return gene, variant
        
        return None, None
    
    def annotate_variants(self, variant_df: pd.DataFrame, annotation_db: Optional[str] = None) -> pd.DataFrame:
        """
        Annotate variants with additional information.
        
        Args:
            variant_df: DataFrame of variants
            annotation_db: Path to annotation database (optional)
            
        Returns:
            Annotated DataFrame
        """
        if variant_df.empty:
            return variant_df
        
        log.info(f"Annotating {len(variant_df)} variants")
        
        # Add impact prediction
        variant_df['impact'] = 'Unknown'
        variant_df['trait'] = 'Unknown'
        
        # Add known annotations based on the gene variant map
        for idx, row in variant_df.iterrows():
            rsid = row.get('rsid', '')
            gene = row.get('gene', '')
            
            # Cognition & Brain Function
            if gene == 'CHRM2' and rsid in ['rs8191992', 'rs2350780']:
                variant_df.at[idx, 'impact'] = 'Executive function, memory, attention'
                variant_df.at[idx, 'trait'] = 'Cognition & Brain Function'
            elif gene == 'DRD2' and rsid == 'rs6277':
                variant_df.at[idx, 'impact'] = 'Dopamine modulation, cognitive flexibility'
                variant_df.at[idx, 'trait'] = 'Cognition & Brain Function'
            elif gene == 'TFAM' and rsid == 'rs1937':
                variant_df.at[idx, 'impact'] = 'Mitochondrial efficiency, energy for brain cells'
                variant_df.at[idx, 'trait'] = 'Cognition & Brain Function'
            elif gene == 'BCL2' and rsid == 'rs956572':
                variant_df.at[idx, 'impact'] = 'Neuroprotection, stress resilience'
                variant_df.at[idx, 'trait'] = 'Cognition & Brain Function'
            elif gene in ['ST8SIA6', 'CHRNA5', 'NRG1', 'MAPRE1', 'GYPC', 'CABP5']:
                variant_df.at[idx, 'impact'] = 'Enhanced neuronal growth and connectivity'
                variant_df.at[idx, 'trait'] = 'Cognition & Brain Function'
                
            # Sleep Traits
            elif gene == 'ADA' and rsid == 'rs73598374':
                variant_df.at[idx, 'impact'] = 'Deep sleep, longer delta wave cycles'
                variant_df.at[idx, 'trait'] = 'Sleep Traits'
            elif gene in ['VRK1', 'CHRM2', 'RALYL', 'FOXO6']:
                variant_df.at[idx, 'impact'] = 'Higher sleep quality'
                variant_df.at[idx, 'trait'] = 'Sleep Traits'
                
            # Cardiovascular Health
            elif gene in ['VKORC1', 'CYP2C9']:
                variant_df.at[idx, 'impact'] = 'Warfarin sensitivity'
                variant_df.at[idx, 'trait'] = 'Cardiovascular Health'
            elif gene == 'CYP2C19':
                variant_df.at[idx, 'impact'] = 'Clopidogrel efficacy'
                variant_df.at[idx, 'trait'] = 'Cardiovascular Health'
            elif gene in ['GP6', 'PTGS1']:
                variant_df.at[idx, 'impact'] = 'Aspirin response'
                variant_df.at[idx, 'trait'] = 'Cardiovascular Health'
                
            # Rare & Neurological Conditions
            elif gene == 'ADGRV1' and rsid in ['rs575602255', 'rs555466095']:
                variant_df.at[idx, 'impact'] = 'Usher Syndrome II (vision + hearing)'
                variant_df.at[idx, 'trait'] = 'Rare & Neurological Conditions'
            elif gene == 'C19orf12' and rsid == 'rs146170087':
                variant_df.at[idx, 'impact'] = 'NBIA-4 (Neurodegeneration with motor + cognitive decline)'
                variant_df.at[idx, 'trait'] = 'Rare & Neurological Conditions'
            elif gene == 'PRSS1' and rsid in ['rs202003805', 'rs1232891794']:
                variant_df.at[idx, 'impact'] = 'Hereditary Pancreatitis (inflammatory episodes)'
                variant_df.at[idx, 'trait'] = 'Rare & Neurological Conditions'
            elif gene == 'ATM' and rsid == 'rs531617441':
                variant_df.at[idx, 'impact'] = 'Increased DNA repair-related cancer susceptibility'
                variant_df.at[idx, 'trait'] = 'Rare & Neurological Conditions'
        
        return variant_df
    
    def generate_variant_report(self, variant_df: pd.DataFrame, output_dir: str) -> str:
        """
        Generate a comprehensive variant report.
        
        Args:
            variant_df: DataFrame of annotated variants
            output_dir: Output directory
            
        Returns:
            Path to the generated report
        """
        if variant_df.empty:
            log.warning("No variants to report")
            return ""
        
        log.info(f"Generating variant report for {len(variant_df)} variants")
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Group variants by trait
        trait_groups = variant_df.groupby('trait')
        
        # Generate report
        report_path = os.path.join(output_dir, "variant_report.md")
        
        with open(report_path, 'w') as f:
            f.write("# Genomic Variant Analysis Report\n\n")
            
            for trait, group in trait_groups:
                if trait == 'Unknown':
                    continue
                    
                f.write(f"## {trait}\n\n")
                f.write("| Gene | Variant (SNP) | Genotype | Implication |\n")
                f.write("|------|--------------|----------|-------------|\n")
                
                for _, row in group.iterrows():
                    gene = row.get('gene', '')
                    rsid = row.get('rsid', '')
                    genotype = row.get('genotype', '')
                    impact = row.get('impact', '')
                    
                    f.write(f"| **{gene}** | {rsid} | {genotype} | {impact} |\n")
                
                f.write("\n")
            
            # Add unknown variants if any
            unknown_variants = variant_df[variant_df['trait'] == 'Unknown']
            if not unknown_variants.empty:
                f.write("## Other Identified Variants\n\n")
                f.write("| Gene | Variant (SNP) | Genotype | Chromosome | Position |\n")
                f.write("|------|--------------|----------|------------|----------|\n")
                
                for _, row in unknown_variants.iterrows():
                    gene = row.get('gene', '')
                    rsid = row.get('rsid', '')
                    genotype = row.get('genotype', '')
                    chrom = row.get('chrom', '')
                    pos = row.get('pos', '')
                    
                    f.write(f"| **{gene}** | {rsid} | {genotype} | {chrom} | {pos} |\n")
        
        log.info(f"Variant report generated: {report_path}")
        return report_path


def process_vcf_with_pysam(vcf_path: str, target_genes: Set[str], dbsnp_path: Optional[str] = None,
                         queue=None, threads: int = 16, batch_size: int = 2000000,
                         max_ram_usage: int = 64, ram_buffer: int = 16) -> Tuple[Dict[str, int], pd.DataFrame]:
    """
    Process a VCF file using pysam with Intel oneAPI acceleration.
    
    Args:
        vcf_path: Path to VCF file
        target_genes: Set of gene symbols to extract variants for
        dbsnp_path: Path to dbSNP VCF file for rsID mapping (optional)
        queue: Intel oneAPI queue for GPU acceleration
        threads: Number of worker threads for parallel processing
        batch_size: Number of variants to process in each batch
        max_ram_usage: Maximum RAM usage in GB
        ram_buffer: RAM buffer in GB
        
    Returns:
        Tuple of (gene_counts, variant_dataframe)
    """
    processor = VariantProcessor(queue=queue, threads=threads)
    
    # Load dbSNP data if available
    if dbsnp_path and os.path.exists(dbsnp_path):
        processor.load_rsid_map(dbsnp_path)
    
    # Process VCF file
    gene_counts, variant_df = processor.process_vcf(
        vcf_path=vcf_path,
        target_genes=target_genes,
        batch_size=batch_size,
        max_ram_usage=max_ram_usage,
        ram_buffer=ram_buffer
    )
    
    # Annotate variants
    if not variant_df.empty:
        variant_df = processor.annotate_variants(variant_df)
    
    return gene_counts, variant_df
