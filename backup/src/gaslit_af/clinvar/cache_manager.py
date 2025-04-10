"""
ClinVar Cache Manager Module for GASLIT-AF Variant Analysis.

This module manages the caching of ClinVar data including metadata management,
versioning, integrity checks, and coordinating the downloader, parser, and indexer.
"""

import os
import json
import gzip
import hashlib
import logging
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, Any, List, Union

from src.gaslit_af.clinvar.downloader import ClinVarDownloader
from src.gaslit_af.clinvar.parser import ClinVarParser
from src.gaslit_af.clinvar.indexer import ClinVarIndexer

# Configure logging
log = logging.getLogger("gaslit-af")

class ClinVarCache:
    """Enhanced caching system for ClinVar data."""
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the ClinVar cache.
        
        Args:
            cache_dir: Directory to store cached data
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path("./cache/clinvar")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Cache subdirectories
        self.raw_dir = self.cache_dir / "raw"
        self.processed_dir = self.cache_dir / "processed"
        self.index_dir = self.cache_dir / "index"
        self.metadata_dir = self.cache_dir / "metadata"
        
        # Create subdirectories
        self.raw_dir.mkdir(exist_ok=True)
        self.processed_dir.mkdir(exist_ok=True)
        self.index_dir.mkdir(exist_ok=True)
        self.metadata_dir.mkdir(exist_ok=True)
        
        # Initialize components
        self.downloader = ClinVarDownloader(self.raw_dir)
        self.parser = ClinVarParser(self.downloader)
        self.indexer = ClinVarIndexer(self.index_dir)
        
        # Metadata file
        self.metadata_file = self.metadata_dir / "cache_metadata.json"
        self._init_metadata()
        
    def _init_metadata(self):
        """Initialize cache metadata."""
        if self.metadata_file.exists():
            try:
                with open(self.metadata_file, 'r') as f:
                    self.metadata = json.load(f)
                log.info(f"Loaded cache metadata from {self.metadata_file}")
            except Exception as e:
                log.error(f"Error loading metadata: {e}")
                self.metadata = self._create_default_metadata()
                self._save_metadata()
        else:
            self.metadata = self._create_default_metadata()
            self._save_metadata()
            
    def _create_default_metadata(self) -> Dict:
        """
        Create default metadata structure.
        
        Returns:
            Default metadata dictionary
        """
        return {
            "versions": {},
            "stats": {
                "variant_count": 0,
                "clinical_significant_count": 0,
                "pathogenic_count": 0,
                "benign_count": 0,
                "vus_count": 0,
                "last_updated": datetime.now().isoformat()
            },
            "file_hashes": {}
        }
        
    def _save_metadata(self):
        """Save cache metadata to file."""
        try:
            with open(self.metadata_file, 'w') as f:
                json.dump(self.metadata, f, indent=2)
            log.debug(f"Saved cache metadata to {self.metadata_file}")
        except Exception as e:
            log.error(f"Error saving metadata: {e}")
            
    def get_file_hash(self, file_path: Path) -> str:
        """
        Calculate MD5 hash of a file.
        
        Args:
            file_path: Path to the file
            
        Returns:
            MD5 hash as a string
        """
        md5 = hashlib.md5()
        
        # Handle gzipped files
        if str(file_path).endswith('.gz'):
            with gzip.open(file_path, 'rb') as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    md5.update(chunk)
        else:
            with open(file_path, 'rb') as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    md5.update(chunk)
                    
        return md5.hexdigest()
        
    def is_cache_valid(self, cache_type: str, max_age_days: int = 30) -> bool:
        """
        Check if a specific cache type is valid.
        
        Args:
            cache_type: Type of cache to check (e.g., 'variant_summary', 'vcf_grch38')
            max_age_days: Maximum age in days for the cache to be considered valid
            
        Returns:
            True if cache is valid, False otherwise
        """
        if cache_type not in self.metadata["versions"]:
            return False
        
        last_update = self.metadata["versions"].get(cache_type, {}).get("last_update")
        if not last_update:
            return False
        
        # Convert to datetime
        last_update_dt = datetime.fromisoformat(last_update)
        age = datetime.now() - last_update_dt
        
        return age.days < max_age_days
        
    def update_cache_metadata(self, cache_type: str, file_path: Path, version: str = None):
        """
        Update cache metadata after processing a file.
        
        Args:
            cache_type: Type of cache being updated
            file_path: Path to the processed file
            version: Version of the data (optional)
        """
        file_hash = self.get_file_hash(file_path)
        
        if cache_type not in self.metadata["versions"]:
            self.metadata["versions"][cache_type] = {}
        
        self.metadata["versions"][cache_type].update({
            "last_update": datetime.now().isoformat(),
            "file_path": str(file_path),
            "file_hash": file_hash,
            "version": version or datetime.now().strftime("%Y%m%d")
        })
        
        self.metadata["file_hashes"][str(file_path)] = file_hash
        self._save_metadata()
        
    def get_cache_stats(self) -> Dict:
        """
        Get statistics about the cache.
        
        Returns:
            Dictionary with cache statistics
        """
        return self.metadata["stats"]
        
    def refresh_variant_summary(self, force_download: bool = False) -> Dict:
        """
        Download and process the latest variant summary data.
        
        Args:
            force_download: Whether to force download even if cache is valid
            
        Returns:
            Dictionary with cache statistics
        """
        log.info("Refreshing variant summary data...")
        
        # Download and parse
        variant_df = self.parser.parse_variant_summary(force_download)
        
        # Index the data
        stats = self.indexer.index_variants(variant_df, "variant_summary")
        
        # Update metadata
        processed_path = self.processed_dir / "variant_summary_processed.parquet"
        if processed_path.exists():
            self.update_cache_metadata("variant_summary", processed_path)
            
        # Update stats
        self.metadata["stats"] = stats
        self._save_metadata()
        
        return stats
        
    def refresh_vcf_data(self, assembly: str = "GRCh38", force_download: bool = False) -> Dict:
        """
        Download and process the latest VCF data.
        
        Args:
            assembly: Genome assembly (GRCh37 or GRCh38)
            force_download: Whether to force download even if cache is valid
            
        Returns:
            Dictionary with cache statistics
        """
        log.info(f"Refreshing VCF data for {assembly}...")
        
        # Download and parse
        vcf_df = self.parser.parse_vcf(assembly, force_download)
        
        # Index the data
        stats = self.indexer.index_variants(vcf_df, f"vcf_{assembly.lower()}")
        
        # Update metadata
        processed_path = self.processed_dir / f"clinvar_vcf_{assembly.lower()}_processed.parquet"
        if processed_path.exists():
            self.update_cache_metadata(f"vcf_{assembly.lower()}", processed_path)
            
        # Update stats
        self.metadata["stats"].update(stats)
        self._save_metadata()
        
        return stats
        
    def lookup_variant(self, **kwargs) -> List[Dict]:
        """
        Look up variants in the index.
        
        Args:
            **kwargs: Query parameters (rs_id, clinvar_id, chrom, pos, ref, alt, gene)
            
        Returns:
            List of matching variants
        """
        return self.indexer.lookup_variant(**kwargs)
        
    def clear_cache(self, cache_type: Optional[str] = None):
        """
        Clear the cache.
        
        Args:
            cache_type: Type of cache to clear (None for all)
        """
        if cache_type:
            # Clear specific cache type
            if cache_type in self.metadata["versions"]:
                file_path = self.metadata["versions"][cache_type].get("file_path")
                if file_path and Path(file_path).exists():
                    Path(file_path).unlink()
                
                del self.metadata["versions"][cache_type]
                self._save_metadata()
                
                log.info(f"Cleared cache for {cache_type}")
        else:
            # Clear all caches
            for subdir in [self.raw_dir, self.processed_dir]:
                for file in subdir.glob("*"):
                    if file.is_file():
                        file.unlink()
            
            # Reset database by reinitializing the indexer
            self.indexer = ClinVarIndexer(self.index_dir)
            
            # Reset metadata
            self.metadata = self._create_default_metadata()
            self._save_metadata()
            
            log.info("Cleared all ClinVar caches")
