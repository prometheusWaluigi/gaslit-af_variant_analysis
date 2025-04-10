"""
ClinVar Indexer Module for GASLIT-AF Variant Analysis.

This module manages the SQLite database for indexing variants for fast lookups.
"""

import logging
import sqlite3
import pandas as pd
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict, List, Any, Union

# Configure logging
log = logging.getLogger("gaslit-af")

class ClinVarIndexer:
    """Manages the SQLite database for ClinVar data."""
    
    def __init__(self, index_dir: Optional[Path] = None):
        """
        Initialize the ClinVar indexer.
        
        Args:
            index_dir: Directory to store the index database
        """
        self.index_dir = Path(index_dir) if index_dir else Path("./cache/clinvar/index")
        self.index_dir.mkdir(parents=True, exist_ok=True)
        
        # Database path
        self.db_path = self.index_dir / "clinvar_index.db"
        
        # Initialize the database
        self._init_database()
    
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
        data = []
        if source == 'variant_summary':
            # Map column names to database fields
            for _, row in df.iterrows():
                try:
                    data.append((
                        row.get('RS# (dbSNP)', ''),
                        row.get('VariationID', 0),
                        row.get('Chromosome', ''),
                        row.get('Start', 0),
                        row.get('ReferenceAllele', ''),
                        row.get('AlternateAllele', ''),
                        row.get('GeneSymbol', ''),
                        row.get('ClinicalSignificance', ''),
                        datetime.now().isoformat()
                    ))
                except Exception as e:
                    log.error(f"Error preparing row for indexing: {e}")
                    continue
        elif source.startswith('vcf_'):
            # Map VCF columns to database fields
            for _, row in df.iterrows():
                try:
                    data.append((
                        row.get('id', ''),
                        row.get('clinvar_id', 0),
                        row.get('chrom', ''),
                        row.get('pos', 0),
                        row.get('ref', ''),
                        row.get('alt', ''),
                        row.get('gene', ''),
                        row.get('significance', ''),
                        datetime.now().isoformat()
                    ))
                except Exception as e:
                    log.error(f"Error preparing row for indexing: {e}")
                    continue
        
        # Insert data in batches
        cursor.executemany(
            'INSERT INTO variant_index (rs_id, clinvar_id, chrom, pos, ref, alt, gene, significance, last_updated) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)', 
            data
        )
        
        # Get stats
        cursor.execute('SELECT COUNT(*) FROM variant_index')
        variant_count = cursor.fetchone()[0]
        
        cursor.execute('SELECT COUNT(*) FROM variant_index WHERE significance NOT IN ("", "not provided", "Uncertain significance", "not reported", NULL)')
        clinical_significant_count = cursor.fetchone()[0]
        
        cursor.execute('SELECT COUNT(*) FROM variant_index WHERE significance LIKE "%pathogenic%" OR significance LIKE "%Pathogenic%"')
        pathogenic_count = cursor.fetchone()[0]
        
        cursor.execute('SELECT COUNT(*) FROM variant_index WHERE significance LIKE "%benign%" OR significance LIKE "%Benign%"')
        benign_count = cursor.fetchone()[0]
        
        cursor.execute('SELECT COUNT(*) FROM variant_index WHERE significance LIKE "%uncertain%" OR significance LIKE "%Uncertain%"')
        vus_count = cursor.fetchone()[0]
        
        conn.commit()
        conn.close()
        
        log.info(f"Indexed {len(data)} variants from {source}")
        log.info(f"Total variants in database: {variant_count}")
        
        # Return stats for cache metadata
        return {
            "variant_count": variant_count,
            "clinical_significant_count": clinical_significant_count,
            "pathogenic_count": pathogenic_count,
            "benign_count": benign_count,
            "vus_count": vus_count,
            "last_updated": datetime.now().isoformat()
        }
    
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
