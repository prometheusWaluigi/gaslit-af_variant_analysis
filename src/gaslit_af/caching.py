"""
Caching module for GASLIT-AF Variant Analysis.
Provides functions to cache intermediate analysis results for improved performance.
"""

import os
import json
import pickle
import hashlib
import time
import logging
from pathlib import Path
from datetime import datetime, timedelta

# Configure logging
log = logging.getLogger("gaslit-af")

class AnalysisCache:
    """Cache manager for analysis results to improve performance on repeated runs."""
    
    def __init__(self, cache_dir="./cache", max_age_hours=24, enabled=True):
        """
        Initialize the cache manager.
        
        Args:
            cache_dir: Directory to store cache files
            max_age_hours: Maximum age of cache files in hours before invalidation
            enabled: Whether caching is enabled
        """
        self.cache_dir = Path(cache_dir)
        self.max_age = timedelta(hours=max_age_hours)
        self.enabled = enabled
        
        if self.enabled:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
            log.info(f"Cache initialized at {self.cache_dir} (max age: {max_age_hours} hours)")
    
    def count_entries(self):
        """Count the number of cache entries.
        
        Returns:
            int: Number of cache entries
        """
        if not self.enabled or not self.cache_dir.exists():
            return 0
        
        return len(list(self.cache_dir.glob('*.cache')))
    
    def get_total_size(self):
        """Get the total size of all cache entries in MB.
        
        Returns:
            float: Total size in MB
        """
        if not self.enabled or not self.cache_dir.exists():
            return 0.0
        
        total_bytes = sum(f.stat().st_size for f in self.cache_dir.glob('*.cache'))
        return total_bytes / (1024 * 1024)  # Convert to MB
    
    def count_expired_entries(self):
        """Count the number of expired cache entries.
        
        Returns:
            int: Number of expired cache entries
        """
        if not self.enabled or not self.cache_dir.exists():
            return 0
        
        now = datetime.now()
        expired_count = 0
        
        for cache_file in self.cache_dir.glob('*.cache'):
            try:
                # Get file modification time
                mtime = datetime.fromtimestamp(cache_file.stat().st_mtime)
                if now - mtime > self.max_age:
                    expired_count += 1
            except Exception:
                # If we can't determine the age, consider it expired
                expired_count += 1
        
        return expired_count
    
    def _get_cache_key(self, vcf_path, analysis_type, params=None):
        """
        Generate a unique cache key based on input parameters.
        
        Args:
            vcf_path: Path to VCF file
            analysis_type: Type of analysis (e.g., 'gene_counts', 'variant_data')
            params: Additional parameters that affect the analysis
        
        Returns:
            String cache key
        """
        # Get VCF file metadata
        try:
            file_stat = os.stat(vcf_path)
            file_size = file_stat.st_size
            file_mtime = file_stat.st_mtime
        except (FileNotFoundError, OSError):
            file_size = 0
            file_mtime = 0
        
        # Create a unique string based on inputs
        key_parts = [
            os.path.abspath(vcf_path),
            str(file_size),
            str(file_mtime),
            analysis_type
        ]
        
        # Add additional parameters if provided
        if params:
            if isinstance(params, dict):
                for k, v in sorted(params.items()):
                    key_parts.append(f"{k}:{v}")
            else:
                key_parts.append(str(params))
        
        # Create a hash of the combined string
        key_string = "|".join(key_parts)
        return hashlib.md5(key_string.encode()).hexdigest()
    
    def _get_cache_path(self, cache_key, analysis_type):
        """Get the file path for a cache entry."""
        return self.cache_dir / f"{analysis_type}_{cache_key}.cache"
    
    def get(self, vcf_path, analysis_type, params=None):
        """
        Retrieve cached results if available and valid.
        
        Args:
            vcf_path: Path to VCF file
            analysis_type: Type of analysis
            params: Additional parameters
        
        Returns:
            Cached data if available, None otherwise
        """
        if not self.enabled:
            return None
        
        cache_key = self._get_cache_key(vcf_path, analysis_type, params)
        cache_path = self._get_cache_path(cache_key, analysis_type)
        
        if not cache_path.exists():
            return None
        
        # Check if cache is expired
        cache_time = datetime.fromtimestamp(cache_path.stat().st_mtime)
        if datetime.now() - cache_time > self.max_age:
            log.info(f"Cache expired for {analysis_type} (age: {datetime.now() - cache_time})")
            return None
        
        try:
            with open(cache_path, 'rb') as f:
                data = pickle.load(f)
                log.info(f"Cache hit for {analysis_type} (key: {cache_key[:8]}...)")
                return data
        except (pickle.PickleError, EOFError, Exception) as e:
            log.warning(f"Error reading cache: {e}")
            return None
    
    def set(self, data, vcf_path, analysis_type, params=None):
        """
        Store results in cache.
        
        Args:
            data: Data to cache
            vcf_path: Path to VCF file
            analysis_type: Type of analysis
            params: Additional parameters
        
        Returns:
            True if successful, False otherwise
        """
        if not self.enabled:
            return False
        
        cache_key = self._get_cache_key(vcf_path, analysis_type, params)
        cache_path = self._get_cache_path(cache_key, analysis_type)
        
        try:
            with open(cache_path, 'wb') as f:
                pickle.dump(data, f)
            log.info(f"Cached {analysis_type} results (key: {cache_key[:8]}...)")
            return True
        except Exception as e:
            log.warning(f"Error writing to cache: {e}")
            return False
    
    def invalidate(self, vcf_path=None, analysis_type=None):
        """
        Invalidate cache entries.
        
        Args:
            vcf_path: Path to VCF file (if None, all files for the analysis_type)
            analysis_type: Type of analysis (if None, all types for the vcf_path)
        
        Returns:
            Number of cache entries invalidated
        """
        if not self.enabled:
            return 0
        
        count = 0
        
        if vcf_path is None and analysis_type is None:
            # Clear all cache
            for cache_file in self.cache_dir.glob("*.cache"):
                cache_file.unlink()
                count += 1
            log.info(f"Cleared entire cache ({count} entries)")
        
        elif vcf_path is not None and analysis_type is None:
            # Clear all cache for a specific VCF file
            for cache_file in self.cache_dir.glob("*.cache"):
                cache_key = cache_file.stem.split('_', 1)[1]
                if vcf_path in cache_key:
                    cache_file.unlink()
                    count += 1
            log.info(f"Cleared cache for {vcf_path} ({count} entries)")
        
        elif vcf_path is None and analysis_type is not None:
            # Clear all cache for a specific analysis type
            for cache_file in self.cache_dir.glob(f"{analysis_type}_*.cache"):
                cache_file.unlink()
                count += 1
            log.info(f"Cleared cache for {analysis_type} ({count} entries)")
        
        else:
            # Clear specific cache
            cache_key = self._get_cache_key(vcf_path, analysis_type)
            cache_path = self._get_cache_path(cache_key, analysis_type)
            if cache_path.exists():
                cache_path.unlink()
                count = 1
                log.info(f"Cleared cache for {analysis_type} on {vcf_path}")
        
        return count
    
    def clean_expired(self):
        """
        Remove all expired cache entries.
        
        Returns:
            Number of expired entries removed
        """
        if not self.enabled:
            return 0
        
        count = 0
        now = datetime.now()
        
        for cache_file in self.cache_dir.glob("*.cache"):
            cache_time = datetime.fromtimestamp(cache_file.stat().st_mtime)
            if now - cache_time > self.max_age:
                cache_file.unlink()
                count += 1
        
        if count > 0:
            log.info(f"Cleaned {count} expired cache entries")
        
        return count
    
    def get_stats(self):
        """
        Get cache statistics.
        
        Returns:
            Dictionary with cache statistics
        """
        if not self.enabled:
            return {"enabled": False}
        
        stats = {
            "enabled": True,
            "cache_dir": str(self.cache_dir),
            "max_age_hours": self.max_age.total_seconds() / 3600,
            "total_entries": 0,
            "total_size_mb": 0,
            "expired_entries": 0,
            "analysis_types": {}
        }
        
        now = datetime.now()
        
        for cache_file in self.cache_dir.glob("*.cache"):
            stats["total_entries"] += 1
            stats["total_size_mb"] += cache_file.stat().st_size / (1024 * 1024)
            
            # Check if expired
            cache_time = datetime.fromtimestamp(cache_file.stat().st_mtime)
            if now - cache_time > self.max_age:
                stats["expired_entries"] += 1
            
            # Count by analysis type
            try:
                analysis_type = cache_file.stem.split('_', 1)[0]
                if analysis_type not in stats["analysis_types"]:
                    stats["analysis_types"][analysis_type] = 0
                stats["analysis_types"][analysis_type] += 1
            except:
                pass
        
        stats["total_size_mb"] = round(stats["total_size_mb"], 2)
        
        return stats
