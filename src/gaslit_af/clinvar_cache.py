"""
ClinVar Cache Module for GASLIT-AF Variant Analysis.

This module provides enhanced caching capabilities for ClinVar data,
including versioning, incremental updates, and indexing for faster lookups.
"""

import os
import json
import gzip
import hashlib
import logging
import shutil
import sqlite3
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from datetime import datetime, timedelta

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
        
        # Initialize SQLite database for indexing
        self.db_path = self.index_dir / "clinvar_index.db"
        self._init_database()
        
        # Metadata file
        self.metadata_file = self.metadata_dir / "cache_metadata.json"
        self._init_metadata()
    
    def _init_database(self):
        """Initialize the SQLite database for indexing."""
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Create tables if they don't exist
        cursor.execute('''
        CREATE TABLE IF NOT EXISTS variant_index (
            id INTEGER PRIMARY KEY,
            rs_id TEXT,
            clinvar_id INTEGER,
            chrom TEXT,
            pos INTEGER,
            ref TEXT,
            alt TEXT,
            gene TEXT,
            significance TEXT,
            last_updated TEXT
        )
        ''')
        
        # Create indices for faster lookups
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_rs_id ON variant_index(rs_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_clinvar_id ON variant_index(clinvar_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_chrom_pos ON variant_index(chrom, pos)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_gene ON variant_index(gene)')
        
        conn.commit()
        conn.close()
    
    def _init_metadata(self):
        """Initialize or load cache metadata."""
        if self.metadata_file.exists():
            try:
                with open(self.metadata_file, 'r') as f:
                    self.metadata = json.load(f)
            except Exception as e:
                log.warning(f"Error loading metadata: {e}. Creating new metadata.")
                self.metadata = self._create_default_metadata()
        else:
            self.metadata = self._create_default_metadata()
            self._save_metadata()
    
    def _create_default_metadata(self) -> Dict:
        """Create default metadata structure."""
        return {
            "last_update": None,
            "versions": {},
            "file_hashes": {},
            "stats": {
                "variant_count": 0,
                "clinical_significant_count": 0,
                "pathogenic_count": 0,
                "benign_count": 0,
                "vus_count": 0
            }
        }
    
    def _save_metadata(self):
        """Save metadata to file."""
        with open(self.metadata_file, 'w') as f:
            json.dump(self.metadata, f, indent=2)
    
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
    
    def index_variants(self, df: pd.DataFrame, source: str):
        """
        Index variants in the SQLite database.
        
        Args:
            df: DataFrame containing variant data
            source: Source of the data (e.g., 'variant_summary', 'vcf')
        """
        conn = sqlite3.connect(self.db_path)
        cursor = conn.cursor()
        
        # Clear existing data if this is a full update
        if source in ['variant_summary', 'vcf_grch37', 'vcf_grch38']:
            cursor.execute('DELETE FROM variant_index WHERE 1=1')
        
        # Prepare data for insertion
        if source == 'variant_summary':
            # Map column names to database fields
            data = []
            for _, row in df.iterrows():
                rs_id = f"rs{int(row['RS# (dbSNP)'])}" if not pd.isna(row['RS# (dbSNP)']) else None
                clinvar_id = int(row['#AlleleID']) if not pd.isna(row['#AlleleID']) else None
                chrom = str(row['Chromosome']) if not pd.isna(row['Chromosome']) else None
                pos = int(row['Start']) if not pd.isna(row['Start']) else None
                ref = str(row['ReferenceAllele']) if not pd.isna(row['ReferenceAllele']) else None
                alt = str(row['AlternateAllele']) if not pd.isna(row['AlternateAllele']) else None
                gene = str(row['GeneSymbol']) if not pd.isna(row['GeneSymbol']) else None
                significance = str(row['ClinicalSignificance']) if not pd.isna(row['ClinicalSignificance']) else None
                last_updated = datetime.now().isoformat()
                
                data.append((rs_id, clinvar_id, chrom, pos, ref, alt, gene, significance, last_updated))
        
        elif source.startswith('vcf_'):
            # Map VCF columns to database fields
            data = []
            for _, row in df.iterrows():
                rs_id = None  # VCF doesn't typically include rsIDs
                clinvar_id = int(row['ID']) if not pd.isna(row['ID']) else None
                chrom = str(row['CHROM']) if not pd.isna(row['CHROM']) else None
                pos = int(row['POS']) if not pd.isna(row['POS']) else None
                ref = str(row['REF']) if not pd.isna(row['REF']) else None
                alt = str(row['ALT']) if not pd.isna(row['ALT']) else None
                gene = None  # VCF doesn't typically include gene symbols
                significance = str(row['CLNSIG']) if not pd.isna(row['CLNSIG']) else None
                last_updated = datetime.now().isoformat()
                
                data.append((rs_id, clinvar_id, chrom, pos, ref, alt, gene, significance, last_updated))
        
        # Insert data in batches
        cursor.executemany(
            'INSERT INTO variant_index (rs_id, clinvar_id, chrom, pos, ref, alt, gene, significance, last_updated) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)',
            data
        )
        
        # Update statistics
        cursor.execute('SELECT COUNT(*) FROM variant_index')
        variant_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM variant_index WHERE significance != 'Uncertain significance' AND significance != 'not provided' AND significance IS NOT NULL")
        clinical_significant_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM variant_index WHERE significance LIKE '%athogenic%'")
        pathogenic_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM variant_index WHERE significance LIKE '%enign%'")
        benign_count = cursor.fetchone()[0]
        
        cursor.execute("SELECT COUNT(*) FROM variant_index WHERE significance LIKE '%ncertain%'")
        vus_count = cursor.fetchone()[0]
        
        # Update metadata
        self.metadata["stats"] = {
            "variant_count": variant_count,
            "clinical_significant_count": clinical_significant_count,
            "pathogenic_count": pathogenic_count,
            "benign_count": benign_count,
            "vus_count": vus_count,
            "last_updated": datetime.now().isoformat()
        }
        
        conn.commit()
        conn.close()
        self._save_metadata()
        
        log.info(f"Indexed {len(data)} variants from {source}")
        log.info(f"Total variants in database: {variant_count}")
    
    def lookup_variant(self, **kwargs) -> List[Dict]:
        """
        Look up variants in the index.
        
        Args:
            **kwargs: Query parameters (rs_id, clinvar_id, chrom, pos, ref, alt, gene)
            
        Returns:
            List of matching variants
        """
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        
        # Build query
        query = "SELECT * FROM variant_index WHERE 1=1"
        params = []
        
        if 'rs_id' in kwargs and kwargs['rs_id']:
            query += " AND rs_id = ?"
            params.append(kwargs['rs_id'])
        
        if 'clinvar_id' in kwargs and kwargs['clinvar_id']:
            query += " AND clinvar_id = ?"
            params.append(int(kwargs['clinvar_id']))
        
        if 'chrom' in kwargs and kwargs['chrom'] and 'pos' in kwargs and kwargs['pos']:
            query += " AND chrom = ? AND pos = ?"
            params.append(str(kwargs['chrom']))
            params.append(int(kwargs['pos']))
            
            if 'ref' in kwargs and kwargs['ref'] and 'alt' in kwargs and kwargs['alt']:
                query += " AND ref = ? AND alt = ?"
                params.append(str(kwargs['ref']))
                params.append(str(kwargs['alt']))
        
        if 'gene' in kwargs and kwargs['gene']:
            query += " AND gene = ?"
            params.append(str(kwargs['gene']))
        
        # Execute query
        cursor.execute(query, params)
        results = [dict(row) for row in cursor.fetchall()]
        
        conn.close()
        return results
    
    def get_cache_stats(self) -> Dict:
        """
        Get statistics about the cache.
        
        Returns:
            Dictionary with cache statistics
        """
        return self.metadata["stats"]
    
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
            
            # Reset database
            if self.db_path.exists():
                self.db_path.unlink()
                self._init_database()
            
            # Reset metadata
            self.metadata = self._create_default_metadata()
            self._save_metadata()
            
            log.info("Cleared all ClinVar caches")


# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Initialize cache
    cache = ClinVarCache()
    
    # Print cache statistics
    stats = cache.get_cache_stats()
    print("Cache Statistics:")
    for key, value in stats.items():
        print(f"  {key}: {value}")
