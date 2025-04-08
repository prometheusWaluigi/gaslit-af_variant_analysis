"""
API Integration Module for GASLIT-AF Variant Analysis.

This module provides integration with external variant annotation APIs:
- Ensembl REST API
- MyVariant.info API
- ClinVar data
- gnomAD data (via BigQuery if configured)
"""

import os
import json
import time
import logging
import requests
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional, Union, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
log = logging.getLogger("gaslit-af")

class VariantAnnotator:
    """Base class for variant annotation services."""
    
    def __init__(self, cache_dir: Optional[Path] = None, cache_ttl: int = 24):
        """
        Initialize the variant annotator.
        
        Args:
            cache_dir: Directory to cache API responses
            cache_ttl: Cache time-to-live in hours
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path("./cache/api")
        self.cache_ttl = cache_ttl
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
    def _get_cache_path(self, variant_id: str, source: str) -> Path:
        """
        Get the cache file path for a variant.
        
        Args:
            variant_id: Variant identifier (e.g., rs123456)
            source: Source of the annotation (e.g., 'ensembl', 'myvariant')
            
        Returns:
            Path: Path to the cache file
        """
        # Create a safe filename from the variant ID
        safe_id = variant_id.replace(":", "_").replace("/", "_")
        return self.cache_dir / f"{source}_{safe_id}.json"
    
    def _is_cache_valid(self, cache_path: Path) -> bool:
        """
        Check if a cache file is valid (exists and not expired).
        
        Args:
            cache_path: Path to the cache file
            
        Returns:
            bool: True if cache is valid, False otherwise
        """
        if not cache_path.exists():
            return False
        
        # Check if cache has expired
        cache_age = time.time() - cache_path.stat().st_mtime
        cache_age_hours = cache_age / 3600
        
        return cache_age_hours < self.cache_ttl
    
    def _load_from_cache(self, variant_id: str, source: str) -> Optional[Dict]:
        """
        Load annotation data from cache if available.
        
        Args:
            variant_id: Variant identifier
            source: Source of the annotation
            
        Returns:
            Optional[Dict]: Cached data if available, None otherwise
        """
        cache_path = self._get_cache_path(variant_id, source)
        
        if self._is_cache_valid(cache_path):
            try:
                with open(cache_path, 'r') as f:
                    return json.load(f)
            except Exception as e:
                log.warning(f"Error loading cache for {variant_id} from {source}: {e}")
                return None
        
        return None
    
    def _save_to_cache(self, variant_id: str, source: str, data: Dict) -> None:
        """
        Save annotation data to cache.
        
        Args:
            variant_id: Variant identifier
            source: Source of the annotation
            data: Annotation data to cache
        """
        cache_path = self._get_cache_path(variant_id, source)
        
        try:
            with open(cache_path, 'w') as f:
                json.dump(data, f, indent=2)
        except Exception as e:
            log.warning(f"Error saving cache for {variant_id} from {source}: {e}")


class EnsemblAnnotator(VariantAnnotator):
    """Annotator for Ensembl REST API."""
    
    def __init__(self, cache_dir: Optional[Path] = None, cache_ttl: int = 24):
        """Initialize the Ensembl annotator."""
        super().__init__(cache_dir, cache_ttl)
        self.base_url = "https://rest.ensembl.org"
        self.headers = {"Content-Type": "application/json"}
        
    def get_variant_consequences(self, variant_id: str) -> Dict:
        """
        Get variant consequences from Ensembl VEP.
        
        Args:
            variant_id: Variant identifier (e.g., rs123456)
            
        Returns:
            Dict: Variant consequences data
        """
        # Check cache first
        cached_data = self._load_from_cache(variant_id, "ensembl_vep")
        if cached_data:
            log.info(f"Using cached Ensembl VEP data for {variant_id}")
            return cached_data
        
        # Make API request
        url = f"{self.base_url}/vep/human/id/{variant_id}"
        
        try:
            log.info(f"Fetching Ensembl VEP data for {variant_id}")
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()
            data = response.json()
            
            # Cache the response
            self._save_to_cache(variant_id, "ensembl_vep", data)
            
            return data
        except Exception as e:
            log.error(f"Error fetching Ensembl VEP data for {variant_id}: {e}")
            return {}
    
    def get_variant_info(self, variant_id: str) -> Dict:
        """
        Get variant information from Ensembl.
        
        Args:
            variant_id: Variant identifier (e.g., rs123456)
            
        Returns:
            Dict: Variant information
        """
        # Check cache first
        cached_data = self._load_from_cache(variant_id, "ensembl_variant")
        if cached_data:
            log.info(f"Using cached Ensembl variant data for {variant_id}")
            return cached_data
        
        # Make API request
        url = f"{self.base_url}/variation/human/{variant_id}"
        
        try:
            log.info(f"Fetching Ensembl variant data for {variant_id}")
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()
            data = response.json()
            
            # Cache the response
            self._save_to_cache(variant_id, "ensembl_variant", data)
            
            return data
        except Exception as e:
            log.error(f"Error fetching Ensembl variant data for {variant_id}: {e}")
            return {}
    
    def get_gene_info(self, gene_id: str) -> Dict:
        """
        Get gene information from Ensembl.
        
        Args:
            gene_id: Gene identifier (e.g., ENSG00000139618)
            
        Returns:
            Dict: Gene information
        """
        # Check cache first
        cached_data = self._load_from_cache(gene_id, "ensembl_gene")
        if cached_data:
            log.info(f"Using cached Ensembl gene data for {gene_id}")
            return cached_data
        
        # Make API request
        url = f"{self.base_url}/lookup/id/{gene_id}"
        
        try:
            log.info(f"Fetching Ensembl gene data for {gene_id}")
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()
            data = response.json()
            
            # Cache the response
            self._save_to_cache(gene_id, "ensembl_gene", data)
            
            return data
        except Exception as e:
            log.error(f"Error fetching Ensembl gene data for {gene_id}: {e}")
            return {}


class MyVariantAnnotator(VariantAnnotator):
    """Annotator for MyVariant.info API."""
    
    def __init__(self, cache_dir: Optional[Path] = None, cache_ttl: int = 24):
        """Initialize the MyVariant annotator."""
        super().__init__(cache_dir, cache_ttl)
        self.base_url = "https://myvariant.info/v1"
        
    def get_variant_info(self, variant_id: str) -> Dict:
        """
        Get variant information from MyVariant.info.
        
        Args:
            variant_id: Variant identifier (e.g., rs123456)
            
        Returns:
            Dict: Variant information
        """
        # Check cache first
        cached_data = self._load_from_cache(variant_id, "myvariant")
        if cached_data:
            log.info(f"Using cached MyVariant data for {variant_id}")
            return cached_data
        
        # Make API request
        url = f"{self.base_url}/variant/{variant_id}"
        
        try:
            log.info(f"Fetching MyVariant data for {variant_id}")
            response = requests.get(url)
            response.raise_for_status()
            data = response.json()
            
            # Cache the response
            self._save_to_cache(variant_id, "myvariant", data)
            
            return data
        except Exception as e:
            log.error(f"Error fetching MyVariant data for {variant_id}: {e}")
            return {}
    
    def get_variants_info(self, variant_ids: List[str]) -> Dict[str, Dict]:
        """
        Get information for multiple variants from MyVariant.info.
        
        Args:
            variant_ids: List of variant identifiers
            
        Returns:
            Dict[str, Dict]: Dictionary mapping variant IDs to their information
        """
        results = {}
        
        # Check which variants need to be fetched
        to_fetch = []
        for variant_id in variant_ids:
            cached_data = self._load_from_cache(variant_id, "myvariant")
            if cached_data:
                log.info(f"Using cached MyVariant data for {variant_id}")
                results[variant_id] = cached_data
            else:
                to_fetch.append(variant_id)
        
        if not to_fetch:
            return results
        
        # Fetch variants in batches of 10
        batch_size = 10
        for i in range(0, len(to_fetch), batch_size):
            batch = to_fetch[i:i+batch_size]
            ids_str = ",".join(batch)
            url = f"{self.base_url}/variant"
            
            try:
                log.info(f"Fetching MyVariant data for batch of {len(batch)} variants")
                response = requests.post(url, data={"ids": ids_str})
                response.raise_for_status()
                data = response.json()
                
                # Process and cache each variant
                for variant_data in data:
                    if "_id" in variant_data:
                        variant_id = variant_data["_id"]
                        results[variant_id] = variant_data
                        self._save_to_cache(variant_id, "myvariant", variant_data)
            except Exception as e:
                log.error(f"Error fetching MyVariant data for batch: {e}")
        
        return results


class VariantAPIIntegration:
    """Integration with multiple variant annotation APIs."""
    
    def __init__(self, cache_dir: Optional[Path] = None, cache_ttl: int = 24):
        """
        Initialize the variant API integration.
        
        Args:
            cache_dir: Directory to cache API responses
            cache_ttl: Cache time-to-live in hours
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path("./cache/api")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize annotators
        self.ensembl = EnsemblAnnotator(self.cache_dir / "ensembl", cache_ttl)
        self.myvariant = MyVariantAnnotator(self.cache_dir / "myvariant", cache_ttl)
        
        # Maximum number of concurrent API requests
        self.max_workers = 5
    
    def annotate_variants(self, variants_df: pd.DataFrame, sources: List[str] = None) -> pd.DataFrame:
        """
        Annotate variants with data from external APIs.
        
        Args:
            variants_df: DataFrame of variants
            sources: List of sources to use for annotation (default: all)
            
        Returns:
            pd.DataFrame: Annotated variants DataFrame
        """
        if variants_df is None or variants_df.empty:
            return variants_df
        
        # Create a copy to avoid modifying the original
        annotated_df = variants_df.copy()
        
        # Default sources
        if sources is None:
            sources = ["ensembl", "myvariant"]
        
        # Get list of variant IDs (rsIDs)
        rsid_col = next((col for col in ["rsid", "variant_id", "id"] if col in annotated_df.columns), None)
        
        if rsid_col is None:
            log.warning("No variant ID column found in DataFrame, cannot annotate")
            return annotated_df
        
        # Filter out missing or invalid rsIDs
        valid_variants = annotated_df[annotated_df[rsid_col].notna() & 
                                    annotated_df[rsid_col].str.startswith('rs', na=False)]
        
        if valid_variants.empty:
            log.warning("No valid rsIDs found in DataFrame, cannot annotate")
            return annotated_df
        
        # Get unique variant IDs
        variant_ids = valid_variants[rsid_col].unique().tolist()
        log.info(f"Annotating {len(variant_ids)} unique variants")
        
        # Annotate with selected sources
        if "ensembl" in sources:
            self._annotate_with_ensembl(annotated_df, variant_ids, rsid_col)
        
        if "myvariant" in sources:
            self._annotate_with_myvariant(annotated_df, variant_ids, rsid_col)
        
        return annotated_df
    
    def _annotate_with_ensembl(self, df: pd.DataFrame, variant_ids: List[str], rsid_col: str) -> None:
        """
        Annotate variants with Ensembl data.
        
        Args:
            df: DataFrame to annotate
            variant_ids: List of variant IDs
            rsid_col: Column name for variant IDs
        """
        log.info("Annotating variants with Ensembl data")
        
        # Add Ensembl annotation columns
        df['ensembl_most_severe'] = None
        df['ensembl_consequence'] = None
        df['ensembl_impact'] = None
        df['ensembl_gene_id'] = None
        df['ensembl_transcript_id'] = None
        df['ensembl_biotype'] = None
        
        # Fetch data for each variant
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit tasks
            future_to_variant = {
                executor.submit(self.ensembl.get_variant_consequences, variant_id): variant_id
                for variant_id in variant_ids
            }
            
            # Process results as they complete
            for future in as_completed(future_to_variant):
                variant_id = future_to_variant[future]
                try:
                    data = future.result()
                    if data and isinstance(data, list) and len(data) > 0:
                        # Find the most severe consequence
                        most_severe = None
                        highest_impact = 0
                        impact_scores = {
                            "HIGH": 4,
                            "MODERATE": 3,
                            "LOW": 2,
                            "MODIFIER": 1
                        }
                        
                        for transcript in data[0].get("transcript_consequences", []):
                            consequences = transcript.get("consequence_terms", [])
                            impact = transcript.get("impact")
                            
                            if impact in impact_scores:
                                impact_score = impact_scores[impact]
                                if impact_score > highest_impact:
                                    highest_impact = impact_score
                                    most_severe = consequences[0] if consequences else None
                        
                        # Update DataFrame for this variant
                        mask = df[rsid_col] == variant_id
                        if most_severe:
                            df.loc[mask, 'ensembl_most_severe'] = most_severe
                        
                        # Get the first transcript consequence (if any)
                        if data[0].get("transcript_consequences"):
                            tc = data[0]["transcript_consequences"][0]
                            df.loc[mask, 'ensembl_consequence'] = ", ".join(tc.get("consequence_terms", []))
                            df.loc[mask, 'ensembl_impact'] = tc.get("impact")
                            df.loc[mask, 'ensembl_gene_id'] = tc.get("gene_id")
                            df.loc[mask, 'ensembl_transcript_id'] = tc.get("transcript_id")
                            df.loc[mask, 'ensembl_biotype'] = tc.get("biotype")
                
                except Exception as e:
                    log.error(f"Error processing Ensembl data for {variant_id}: {e}")
    
    def _annotate_with_myvariant(self, df: pd.DataFrame, variant_ids: List[str], rsid_col: str) -> None:
        """
        Annotate variants with MyVariant.info data.
        
        Args:
            df: DataFrame to annotate
            variant_ids: List of variant IDs
            rsid_col: Column name for variant IDs
        """
        log.info("Annotating variants with MyVariant.info data")
        
        # Add MyVariant annotation columns
        df['clinvar_significance'] = None
        df['cadd_phred'] = None
        df['gnomad_af'] = None
        df['sift_pred'] = None
        df['polyphen_pred'] = None
        
        # Fetch data for all variants
        variants_data = self.myvariant.get_variants_info(variant_ids)
        
        # Process results
        for variant_id, data in variants_data.items():
            if not data:
                continue
            
            # Update DataFrame for this variant
            mask = df[rsid_col] == variant_id
            
            # Extract ClinVar data
            if 'clinvar' in data:
                clinvar = data['clinvar']
                if 'rcv' in clinvar and 'clinical_significance' in clinvar['rcv']:
                    df.loc[mask, 'clinvar_significance'] = clinvar['rcv']['clinical_significance']
            
            # Extract CADD score
            if 'cadd' in data and 'phred' in data['cadd']:
                df.loc[mask, 'cadd_phred'] = data['cadd']['phred']
            
            # Extract gnomAD allele frequency
            if 'gnomad_genome' in data and 'af' in data['gnomad_genome']:
                df.loc[mask, 'gnomad_af'] = data['gnomad_genome']['af']
            
            # Extract SIFT prediction
            if 'dbnsfp' in data and 'sift' in data['dbnsfp'] and 'pred' in data['dbnsfp']['sift']:
                sift_pred = data['dbnsfp']['sift']['pred']
                if isinstance(sift_pred, list) and sift_pred:
                    df.loc[mask, 'sift_pred'] = sift_pred[0]
                else:
                    df.loc[mask, 'sift_pred'] = sift_pred
            
            # Extract PolyPhen prediction
            if 'dbnsfp' in data and 'polyphen2' in data['dbnsfp'] and 'hdiv' in data['dbnsfp']['polyphen2'] and 'pred' in data['dbnsfp']['polyphen2']['hdiv']:
                polyphen_pred = data['dbnsfp']['polyphen2']['hdiv']['pred']
                if isinstance(polyphen_pred, list) and polyphen_pred:
                    df.loc[mask, 'polyphen_pred'] = polyphen_pred[0]
                else:
                    df.loc[mask, 'polyphen_pred'] = polyphen_pred
    
    def get_variant_details(self, variant_id: str) -> Dict:
        """
        Get comprehensive details for a variant from multiple sources.
        
        Args:
            variant_id: Variant identifier (e.g., rs123456)
            
        Returns:
            Dict: Comprehensive variant details
        """
        details = {
            "variant_id": variant_id,
            "sources": {}
        }
        
        # Get Ensembl data
        ensembl_vep = self.ensembl.get_variant_consequences(variant_id)
        ensembl_variant = self.ensembl.get_variant_info(variant_id)
        
        if ensembl_vep:
            details["sources"]["ensembl_vep"] = ensembl_vep
        
        if ensembl_variant:
            details["sources"]["ensembl_variant"] = ensembl_variant
        
        # Get MyVariant data
        myvariant_data = self.myvariant.get_variant_info(variant_id)
        
        if myvariant_data:
            details["sources"]["myvariant"] = myvariant_data
        
        # Extract key information
        self._extract_key_details(details)
        
        return details
    
    def _extract_key_details(self, details: Dict) -> None:
        """
        Extract key details from the raw API responses.
        
        Args:
            details: Variant details dictionary to update
        """
        # Initialize key fields
        details["gene"] = None
        details["consequence"] = None
        details["impact"] = None
        details["clinical_significance"] = None
        details["allele_frequency"] = None
        details["pathogenicity_scores"] = {}
        
        # Extract from Ensembl VEP
        if "ensembl_vep" in details["sources"] and details["sources"]["ensembl_vep"]:
            vep_data = details["sources"]["ensembl_vep"]
            if isinstance(vep_data, list) and vep_data:
                vep_data = vep_data[0]
                
                # Get gene
                if "transcript_consequences" in vep_data and vep_data["transcript_consequences"]:
                    tc = vep_data["transcript_consequences"][0]
                    details["gene"] = tc.get("gene_symbol") or tc.get("gene_id")
                    details["consequence"] = ", ".join(tc.get("consequence_terms", []))
                    details["impact"] = tc.get("impact")
        
        # Extract from MyVariant
        if "myvariant" in details["sources"] and details["sources"]["myvariant"]:
            mv_data = details["sources"]["myvariant"]
            
            # Get clinical significance
            if "clinvar" in mv_data and "rcv" in mv_data["clinvar"] and "clinical_significance" in mv_data["clinvar"]["rcv"]:
                details["clinical_significance"] = mv_data["clinvar"]["rcv"]["clinical_significance"]
            
            # Get allele frequency
            if "gnomad_genome" in mv_data and "af" in mv_data["gnomad_genome"]:
                details["allele_frequency"] = mv_data["gnomad_genome"]["af"]
            
            # Get pathogenicity scores
            if "cadd" in mv_data and "phred" in mv_data["cadd"]:
                details["pathogenicity_scores"]["cadd_phred"] = mv_data["cadd"]["phred"]
            
            if "dbnsfp" in mv_data:
                dbnsfp = mv_data["dbnsfp"]
                
                if "sift" in dbnsfp and "pred" in dbnsfp["sift"]:
                    sift_pred = dbnsfp["sift"]["pred"]
                    if isinstance(sift_pred, list) and sift_pred:
                        details["pathogenicity_scores"]["sift"] = sift_pred[0]
                    else:
                        details["pathogenicity_scores"]["sift"] = sift_pred
                
                if "polyphen2" in dbnsfp and "hdiv" in dbnsfp["polyphen2"] and "pred" in dbnsfp["polyphen2"]["hdiv"]:
                    polyphen_pred = dbnsfp["polyphen2"]["hdiv"]["pred"]
                    if isinstance(polyphen_pred, list) and polyphen_pred:
                        details["pathogenicity_scores"]["polyphen"] = polyphen_pred[0]
                    else:
                        details["pathogenicity_scores"]["polyphen"] = polyphen_pred


# Example usage
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    # Initialize API integration
    api = VariantAPIIntegration()
    
    # Example: Get details for a variant
    variant_id = "rs429358"  # APOE variant
    details = api.get_variant_details(variant_id)
    
    print(f"Variant: {variant_id}")
    print(f"Gene: {details['gene']}")
    print(f"Consequence: {details['consequence']}")
    print(f"Impact: {details['impact']}")
    print(f"Clinical Significance: {details['clinical_significance']}")
    print(f"Allele Frequency: {details['allele_frequency']}")
    print("Pathogenicity Scores:")
    for score, value in details['pathogenicity_scores'].items():
        print(f"  {score}: {value}")
