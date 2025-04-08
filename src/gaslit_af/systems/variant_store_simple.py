"""
Simplified Variant Store Module for GASLIT-AF Variant Analysis.

This module provides a lightweight, high-performance variant storage system using DuckDB.
Variants are stored by biological system, allowing for efficient querying and recursive pattern analysis.
"""

import os
import json
import logging
import pandas as pd
import duckdb
from pathlib import Path
from typing import Dict, List, Set, Optional, Any, Union, Tuple

from .gene_systems import get_gene_system_manager

# Configure logging
log = logging.getLogger("gaslit-af")

class VariantStore:
    """Lightweight variant storage using DuckDB for high-performance queries."""
    
    def __init__(self, db_path: Optional[Path] = None, systems_dir: Optional[Path] = None):
        """
        Initialize the variant store.
        
        Args:
            db_path: Path to DuckDB database file
            systems_dir: Directory containing gene system JSON files
        """
        self.db_path = Path(db_path) if db_path else Path("./data/variants.duckdb")
        self.db_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Initialize gene system manager
        self.gene_systems = get_gene_system_manager(systems_dir)
        
        # Initialize DuckDB connection
        self._initialize_db()
    
    def _initialize_db(self):
        """Initialize the DuckDB database and tables."""
        log.info(f"Initializing variant store at {self.db_path}")
        
        # Connect to DuckDB
        self.conn = duckdb.connect(str(self.db_path))
        
        # Create tables if they don't exist
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS variants (
                variant_id VARCHAR PRIMARY KEY,
                chrom VARCHAR,
                pos INTEGER,
                ref VARCHAR,
                alt VARCHAR,
                gene VARCHAR,
                system_id VARCHAR,
                rsid VARCHAR,
                clinvar_id VARCHAR,
                clinvar_significance VARCHAR,
                af DOUBLE,
                impact VARCHAR,
                pathogenic BOOLEAN,
                pathogenic_reason VARCHAR,
                metadata JSON
            )
        """)
        
        self.conn.execute("""
            CREATE TABLE IF NOT EXISTS system_variants (
                system_id VARCHAR,
                variant_id VARCHAR,
                gene VARCHAR,
                pathogenic BOOLEAN,
                PRIMARY KEY (system_id, variant_id)
            )
        """)
        
        self.conn.execute("""
            CREATE INDEX IF NOT EXISTS idx_variants_gene ON variants(gene);
            CREATE INDEX IF NOT EXISTS idx_variants_system ON variants(system_id);
            CREATE INDEX IF NOT EXISTS idx_variants_rsid ON variants(rsid);
            CREATE INDEX IF NOT EXISTS idx_system_variants_system ON system_variants(system_id);
        """)
    
    def store_variant(self, variant_data: Dict[str, Any]) -> bool:
        """
        Store a variant in the database.
        
        Args:
            variant_data: Dictionary containing variant data
            
        Returns:
            True if successful, False otherwise
        """
        try:
            # Extract required fields
            variant_id = variant_data.get("variant_id")
            if not variant_id:
                # Generate a variant ID if not provided
                chrom = variant_data.get("chrom", "")
                pos = variant_data.get("pos", 0)
                ref = variant_data.get("ref", "")
                alt = variant_data.get("alt", "")
                variant_id = f"{chrom}:{pos}:{ref}:{alt}"
                variant_data["variant_id"] = variant_id
            
            # Get gene and system
            gene = variant_data.get("gene", "")
            system_id = self.gene_systems.get_system_for_gene(gene)
            
            # Store remaining data as JSON metadata
            metadata = {k: v for k, v in variant_data.items() if k not in [
                "variant_id", "chrom", "pos", "ref", "alt", "gene", "rsid", 
                "clinvar_id", "clinvar_significance", "af", "impact", 
                "pathogenic", "pathogenic_reason", "system_id"
            ]}
            
            # Direct SQL insert for variants table
            self.conn.execute("""
                INSERT INTO variants 
                (variant_id, chrom, pos, ref, alt, gene, system_id, rsid, 
                clinvar_id, clinvar_significance, af, impact, pathogenic, 
                pathogenic_reason, metadata)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT (variant_id) DO UPDATE SET
                chrom = excluded.chrom,
                pos = excluded.pos,
                ref = excluded.ref,
                alt = excluded.alt,
                gene = excluded.gene,
                system_id = excluded.system_id,
                rsid = excluded.rsid,
                clinvar_id = excluded.clinvar_id,
                clinvar_significance = excluded.clinvar_significance,
                af = excluded.af,
                impact = excluded.impact,
                pathogenic = excluded.pathogenic,
                pathogenic_reason = excluded.pathogenic_reason,
                metadata = excluded.metadata
            """, [
                variant_id,
                variant_data.get("chrom", ""),
                variant_data.get("pos", 0),
                variant_data.get("ref", ""),
                variant_data.get("alt", ""),
                gene,
                system_id,
                variant_data.get("rsid", ""),
                variant_data.get("clinvar_id", ""),
                variant_data.get("clinvar_significance", ""),
                variant_data.get("af", 0.0),
                variant_data.get("impact", ""),
                variant_data.get("pathogenic", False),
                variant_data.get("pathogenic_reason", ""),
                json.dumps(metadata)
            ])
            
            # Direct SQL insert for system_variants table
            self.conn.execute("""
                INSERT INTO system_variants 
                (system_id, variant_id, gene, pathogenic)
                VALUES (?, ?, ?, ?)
                ON CONFLICT (system_id, variant_id) DO UPDATE SET
                gene = excluded.gene,
                pathogenic = excluded.pathogenic
            """, [
                system_id,
                variant_id,
                gene,
                variant_data.get("pathogenic", False)
            ])
            
            return True
        
        except Exception as e:
            log.error(f"Error storing variant: {e}")
            return False
    
    def store_variants(self, variant_df: pd.DataFrame) -> int:
        """
        Store multiple variants from a DataFrame.
        
        Args:
            variant_df: DataFrame containing variant data
            
        Returns:
            Number of variants successfully stored
        """
        success_count = 0
        
        # Process each row
        for _, row in variant_df.iterrows():
            variant_data = row.to_dict()
            if self.store_variant(variant_data):
                success_count += 1
        
        log.info(f"Stored {success_count} out of {len(variant_df)} variants")
        return success_count
    
    def get_variant(self, variant_id: str) -> Optional[Dict[str, Any]]:
        """
        Get a variant by ID.
        
        Args:
            variant_id: Variant identifier
            
        Returns:
            Dictionary containing variant data or None if not found
        """
        try:
            # Use fetch_df to get a DataFrame and then convert to dict
            result_df = self.conn.execute(
                "SELECT * FROM variants WHERE variant_id = ?", [variant_id]
            ).fetch_df()
            
            if not result_df.empty:
                # Convert first row to dictionary
                variant_data = result_df.iloc[0].to_dict()
                
                # Parse metadata
                if "metadata" in variant_data and variant_data["metadata"]:
                    try:
                        metadata = json.loads(variant_data["metadata"])
                        variant_data.update(metadata)
                        del variant_data["metadata"]
                    except:
                        pass
                
                return variant_data
            
            return None
        
        except Exception as e:
            log.error(f"Error retrieving variant {variant_id}: {e}")
            return None
    
    def get_variants_by_gene(self, gene: str) -> pd.DataFrame:
        """
        Get all variants for a gene.
        
        Args:
            gene: Gene symbol
            
        Returns:
            DataFrame containing variant data
        """
        try:
            query = "SELECT * FROM variants WHERE gene = ?"
            result = self.conn.execute(query, [gene]).fetch_df()
            
            # Parse metadata
            if not result.empty and "metadata" in result.columns:
                result = self._expand_metadata(result)
            
            return result
        
        except Exception as e:
            log.error(f"Error retrieving variants for gene {gene}: {e}")
            return pd.DataFrame()
    
    def get_variants_by_system(self, system_id: str) -> pd.DataFrame:
        """
        Get all variants for a biological system.
        
        Args:
            system_id: System identifier
            
        Returns:
            DataFrame containing variant data
        """
        try:
            # Directly query the variants table by system_id
            query = "SELECT * FROM variants WHERE system_id = ?"
            result = self.conn.execute(query, [system_id]).fetch_df()
            
            # Parse metadata
            if not result.empty and "metadata" in result.columns:
                result = self._expand_metadata(result)
            
            return result
        
        except Exception as e:
            log.error(f"Error retrieving variants for system {system_id}: {e}")
            return pd.DataFrame()
    
    def get_pathogenic_variants(self, system_id: Optional[str] = None) -> pd.DataFrame:
        """
        Get all pathogenic variants, optionally filtered by system.
        
        Args:
            system_id: System identifier (optional)
            
        Returns:
            DataFrame containing variant data
        """
        try:
            if system_id:
                query = "SELECT * FROM variants WHERE pathogenic = TRUE AND system_id = ?"
                result = self.conn.execute(query, [system_id]).fetch_df()
            else:
                query = "SELECT * FROM variants WHERE pathogenic = TRUE"
                result = self.conn.execute(query).fetch_df()
            
            # Parse metadata
            if not result.empty and "metadata" in result.columns:
                result = self._expand_metadata(result)
            
            return result
        
        except Exception as e:
            log.error(f"Error retrieving pathogenic variants: {e}")
            return pd.DataFrame()
    
    def _expand_metadata(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Expand JSON metadata column into separate columns.
        
        Args:
            df: DataFrame with metadata column
            
        Returns:
            DataFrame with expanded metadata
        """
        # Create a copy to avoid modifying the original
        result = df.copy()
        
        # Parse metadata for each row
        for idx, row in result.iterrows():
            if pd.notna(row.get("metadata")) and row["metadata"]:
                try:
                    metadata = json.loads(row["metadata"])
                    for key, value in metadata.items():
                        result.at[idx, key] = value
                except:
                    pass
        
        # Drop metadata column
        if "metadata" in result.columns:
            result = result.drop("metadata", axis=1)
        
        return result
    
    def get_system_summary(self, system_id: str) -> Dict[str, Any]:
        """
        Get a summary of variants for a biological system.
        
        Args:
            system_id: System identifier
            
        Returns:
            Dictionary containing summary statistics
        """
        try:
            # Get system info
            system_info = self.gene_systems.systems.get(system_id, {})
            system_name = system_info.get("name", system_id)
            
            # Get variant counts
            total_query = "SELECT COUNT(*) as count FROM system_variants WHERE system_id = ?"
            total_count = self.conn.execute(total_query, [system_id]).fetchone()[0]
            
            pathogenic_query = "SELECT COUNT(*) as count FROM system_variants WHERE system_id = ? AND pathogenic = TRUE"
            pathogenic_count = self.conn.execute(pathogenic_query, [system_id]).fetchone()[0]
            
            gene_query = """
                SELECT gene, COUNT(*) as count 
                FROM system_variants 
                WHERE system_id = ? 
                GROUP BY gene 
                ORDER BY count DESC
            """
            gene_counts = self.conn.execute(gene_query, [system_id]).fetch_df()
            
            return {
                "system_id": system_id,
                "system_name": system_name,
                "total_variants": int(total_count),
                "pathogenic_variants": int(pathogenic_count),
                "pathogenic_percentage": round(pathogenic_count / total_count * 100 if total_count > 0 else 0, 2),
                "gene_counts": gene_counts.to_dict(orient="records") if not gene_counts.empty else []
            }
        
        except Exception as e:
            log.error(f"Error generating summary for system {system_id}: {e}")
            return {
                "system_id": system_id,
                "error": str(e)
            }
    
    def close(self):
        """Close the database connection."""
        if hasattr(self, 'conn'):
            self.conn.close()


# Singleton instance
_variant_store = None

def get_variant_store(db_path: Optional[Path] = None, systems_dir: Optional[Path] = None) -> VariantStore:
    """
    Get the singleton variant store instance.
    
    Args:
        db_path: Path to DuckDB database file
        systems_dir: Directory containing gene system JSON files
        
    Returns:
        VariantStore instance
    """
    global _variant_store
    if _variant_store is None:
        _variant_store = VariantStore(db_path, systems_dir)
    return _variant_store
