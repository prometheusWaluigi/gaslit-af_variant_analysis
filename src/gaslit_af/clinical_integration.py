"""
Clinical Integration Module for GASLIT-AF Variant Analysis.

This module integrates clinical variant data with the variant analysis workflow.
"""

import os
import json
import logging
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional, Union

from src.gaslit_af.clinical_variants import ClinicalVariantManager

# Configure logging
log = logging.getLogger("gaslit-af")

class ClinicalIntegration:
    """Integrates clinical variant data with the variant analysis workflow."""
    
    def __init__(self, clinical_data_path: Optional[Union[str, Path]] = None):
        """
        Initialize the clinical integration module.
        
        Args:
            clinical_data_path: Path to the clinical data JSON file. If None, no clinical data is loaded.
        """
        self.manager = ClinicalVariantManager()
        self.clinical_data_loaded = False
        
        if clinical_data_path:
            self.load_clinical_data(clinical_data_path)
    
    def load_clinical_data(self, clinical_data_path: Union[str, Path]) -> bool:
        """
        Load clinical data from a JSON file.
        
        Args:
            clinical_data_path: Path to the clinical data JSON file
            
        Returns:
            bool: True if data was loaded successfully, False otherwise
        """
        try:
            self.manager.load_conditions(clinical_data_path)
            self.clinical_data_loaded = True
            log.info(f"Loaded clinical data from {clinical_data_path}")
            return True
        except Exception as e:
            log.error(f"Error loading clinical data from {clinical_data_path}: {e}")
            self.clinical_data_loaded = False
            return False
    
    def annotate_variants(self, variants_df) -> pd.DataFrame:
        """
        Annotate variants with clinical information.
        
        Args:
            variants_df: DataFrame of variants
            
        Returns:
            pd.DataFrame: Annotated variants DataFrame
        """
        if not self.clinical_data_loaded:
            log.warning("No clinical data loaded, cannot annotate variants")
            return variants_df
        
        return self.manager.annotate_variants(variants_df)
    
    def generate_clinical_report(self, variants_df, output_dir: Union[str, Path]) -> Optional[str]:
        """
        Generate a clinical report based on the variants and clinical data.
        
        Args:
            variants_df: DataFrame of variants
            output_dir: Directory to save the report
            
        Returns:
            Optional[str]: Path to the generated report, or None if no report was generated
        """
        if not self.clinical_data_loaded:
            log.warning("No clinical data loaded, cannot generate clinical report")
            return None
        
        return self.manager.generate_clinical_report(variants_df, output_dir)
    
    def get_pathogenic_variants(self, variants_df) -> pd.DataFrame:
        """
        Filter variants to only include those that are pathogenic or likely pathogenic.
        
        Args:
            variants_df: DataFrame of variants
            
        Returns:
            pd.DataFrame: DataFrame containing only pathogenic variants
        """
        if not self.clinical_data_loaded or variants_df is None or variants_df.empty:
            return pd.DataFrame()
        
        # Annotate variants first
        annotated_df = self.annotate_variants(variants_df)
        
        # Filter for pathogenic variants
        return annotated_df[annotated_df['clinical_significance'].isin(['Pathogenic', 'Likely Pathogenic'])]
    
    def get_variants_by_condition(self, variants_df, condition_name: str) -> pd.DataFrame:
        """
        Filter variants to only include those associated with a specific condition.
        
        Args:
            variants_df: DataFrame of variants
            condition_name: Name of the condition to filter for
            
        Returns:
            pd.DataFrame: DataFrame containing only variants associated with the condition
        """
        if not self.clinical_data_loaded or variants_df is None or variants_df.empty:
            return pd.DataFrame()
        
        # Annotate variants first
        annotated_df = self.annotate_variants(variants_df)
        
        # Filter for condition
        return annotated_df[annotated_df['condition_name'] == condition_name]
    
    def get_clinical_summary(self, variants_df) -> Dict:
        """
        Generate a summary of clinical findings from the variants.
        
        Args:
            variants_df: DataFrame of variants
            
        Returns:
            Dict: Summary of clinical findings
        """
        if not self.clinical_data_loaded or variants_df is None or variants_df.empty:
            return {"pathogenic_count": 0, "benign_count": 0, "vus_count": 0, "conditions": []}
        
        # Annotate variants first
        annotated_df = self.annotate_variants(variants_df)
        
        # Count variants by clinical significance
        pathogenic_count = annotated_df[annotated_df['clinical_significance'].isin(
            ['Pathogenic', 'Likely Pathogenic'])].shape[0]
        
        benign_count = annotated_df[annotated_df['clinical_significance'].isin(
            ['Benign', 'Likely Benign'])].shape[0]
        
        vus_count = annotated_df[annotated_df['clinical_significance'] == 'Uncertain Significance'].shape[0]
        
        # Get unique conditions
        conditions = []
        for condition_name in annotated_df['condition_name'].dropna().unique():
            condition_variants = annotated_df[annotated_df['condition_name'] == condition_name]
            conditions.append({
                "name": condition_name,
                "description": condition_variants['condition_description'].iloc[0],
                "variant_count": condition_variants.shape[0],
                "genes": condition_variants['gene'].unique().tolist()
            })
        
        return {
            "pathogenic_count": pathogenic_count,
            "benign_count": benign_count,
            "vus_count": vus_count,
            "conditions": conditions
        }
