"""
Gene Mapping Enhancement Module for GASLIT-AF Variant Analysis.

This module provides enhanced gene mapping capabilities to create a deeper
quantum coherence bridge between personal genomic architecture and the
GASLIT-AF theoretical model parameters.
"""

import os
import json
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Set, Tuple, Any, Optional
from collections import defaultdict

# Import GASLIT-AF modules
from src.gaslit_af.gene_lists import GASLIT_AF_GENES, KNOWN_SNPS
from src.gaslit_af.biological_systems import BIOLOGICAL_SYSTEMS, get_system_for_gene

# Configure logging
log = logging.getLogger("gaslit-af")

# Define gene mapping databases
GENE_DATABASES = [
    "NCBI Gene",
    "Ensembl",
    "HGNC",
    "UniProt",
    "OMIM",
    "GeneCards"
]

# Define gene mapping resources
GENE_RESOURCES = [
    "Gene Ontology",
    "KEGG Pathways",
    "Reactome",
    "BioCyc",
    "STRING",
    "BioGRID"
]

# Define GASLIT-AF system to pathway mappings
SYSTEM_PATHWAY_MAP = {
    "Immune & Inflammatory Pathways": [
        "Cytokine signaling",
        "Innate immune response",
        "Adaptive immune response",
        "Inflammatory response",
        "Complement cascade",
        "NF-kB signaling",
        "JAK-STAT signaling"
    ],
    "Autonomic & Neurotransmitter Pathways": [
        "Cholinergic signaling",
        "Adrenergic signaling",
        "Dopaminergic signaling",
        "Serotonergic signaling",
        "GABA signaling",
        "Glutamatergic signaling",
        "Neuropeptide signaling"
    ],
    "Structural & Connective Tissue Integrity": [
        "Collagen biosynthesis",
        "Extracellular matrix organization",
        "Cell adhesion",
        "Cytoskeleton organization",
        "Focal adhesion",
        "Integrin signaling",
        "Elastin metabolism"
    ],
    "Metabolic, Mitochondrial & Oxidative Stress": [
        "Oxidative phosphorylation",
        "TCA cycle",
        "Fatty acid metabolism",
        "Glucose metabolism",
        "Mitochondrial dynamics",
        "Redox homeostasis",
        "Antioxidant response"
    ],
    "Endocannabinoid System (ECS)": [
        "Endocannabinoid synthesis",
        "Endocannabinoid degradation",
        "Cannabinoid receptor signaling",
        "Retrograde endocannabinoid signaling",
        "Anandamide metabolism",
        "2-AG metabolism"
    ],
    "Calcium & Ion Channels": [
        "Calcium signaling",
        "Potassium channels",
        "Sodium channels",
        "Chloride channels",
        "TRP channels",
        "Voltage-gated ion channels",
        "Ligand-gated ion channels"
    ],
    "Mast Cell Activation & Histamine Metabolism": [
        "Histamine metabolism",
        "Mast cell degranulation",
        "IgE signaling",
        "Tryptase activity",
        "Histamine receptor signaling",
        "Leukotriene biosynthesis"
    ],
    "Kynurenine Pathway": [
        "Tryptophan metabolism",
        "Kynurenine synthesis",
        "Quinolinic acid synthesis",
        "Kynurenic acid synthesis",
        "IDO activity",
        "TDO activity",
        "NAD+ biosynthesis"
    ]
}

# Define GASLIT-AF parameter to biological function mappings
PARAMETER_FUNCTION_MAP = {
    "γ": [  # Genetic fragility
        "DNA repair",
        "Genome stability",
        "Telomere maintenance",
        "Chromatin remodeling",
        "Epigenetic regulation",
        "Transposon silencing"
    ],
    "Λ": [  # Allostatic load
        "Stress response",
        "HPA axis regulation",
        "Cortisol metabolism",
        "Inflammatory regulation",
        "Autonomic balance",
        "Energy homeostasis"
    ],
    "Ω": [  # Endocannabinoid buffering
        "Endocannabinoid synthesis",
        "Endocannabinoid degradation",
        "Cannabinoid receptor signaling",
        "Retrograde signaling",
        "Synaptic plasticity",
        "Neuroimmune modulation"
    ],
    "Χ": [  # Physiological coherence
        "Autonomic regulation",
        "Heart rate variability",
        "Circadian rhythm",
        "Neural synchrony",
        "Metabolic flexibility",
        "Homeostatic resilience"
    ],
    "σ": [  # Entropy production
        "Mitochondrial efficiency",
        "Redox balance",
        "ATP synthesis",
        "Proton leak",
        "Metabolic rate",
        "Thermogenesis"
    ]
}

class GeneMapper:
    """
    Enhanced gene mapper for GASLIT-AF Variant Analysis.
    
    This class provides methods to map genes to biological systems,
    pathways, and GASLIT-AF parameters with enhanced precision.
    """
    
    def __init__(self, gene_db_path: Optional[str] = None,
                external_mappings_path: Optional[str] = None):
        """
        Initialize gene mapper.
        
        Args:
            gene_db_path: Path to gene database file (JSON)
            external_mappings_path: Path to external mappings file (JSON)
        """
        self.gene_db = {}
        self.external_mappings = {}
        
        # Load gene database if provided
        if gene_db_path and os.path.exists(gene_db_path):
            try:
                with open(gene_db_path, 'r', encoding='utf-8') as f:
                    self.gene_db = json.load(f)
                log.info(f"Loaded gene database from {gene_db_path}")
            except Exception as e:
                log.error(f"Error loading gene database: {e}")
        
        # Load external mappings if provided
        if external_mappings_path and os.path.exists(external_mappings_path):
            try:
                with open(external_mappings_path, 'r', encoding='utf-8') as f:
                    self.external_mappings = json.load(f)
                log.info(f"Loaded external mappings from {external_mappings_path}")
            except Exception as e:
                log.error(f"Error loading external mappings: {e}")
        
        # Initialize system gene sets
        self.system_genes = {system: set(genes) for system, genes in BIOLOGICAL_SYSTEMS.items()}
        
        # Initialize reverse mappings
        self.rsid_to_gene = {}
        for gene, rsids in KNOWN_SNPS.items():
            for rsid in rsids:
                self.rsid_to_gene[rsid] = gene
    
    def map_rsid_to_gene(self, rsid: str) -> Optional[str]:
        """
        Map rsID to gene symbol.
        
        Args:
            rsid: rsID (e.g., rs429358)
            
        Returns:
            Gene symbol or None if not found
        """
        # Check GASLIT-AF known SNPs
        if rsid in self.rsid_to_gene:
            return self.rsid_to_gene[rsid]
        
        # Check external mappings
        if rsid in self.external_mappings.get('rsid_to_gene', {}):
            return self.external_mappings['rsid_to_gene'][rsid]
        
        # Check gene database
        if rsid in self.gene_db.get('rsid_to_gene', {}):
            return self.gene_db['rsid_to_gene'][rsid]
        
        return None
    
    def map_gene_to_system(self, gene: str) -> str:
        """
        Map gene to biological system with enhanced precision.
        
        Args:
            gene: Gene symbol
            
        Returns:
            Biological system name
        """
        # Check GASLIT-AF biological systems
        for system, genes in self.system_genes.items():
            if gene in genes:
                return system
        
        # Check external mappings
        if gene in self.external_mappings.get('gene_to_system', {}):
            return self.external_mappings['gene_to_system'][gene]
        
        # Check gene database
        if gene in self.gene_db.get('gene_to_system', {}):
            return self.gene_db['gene_to_system'][gene]
        
        # Default to "Other"
        return "Other"
    
    def map_gene_to_pathways(self, gene: str) -> List[str]:
        """
        Map gene to biological pathways.
        
        Args:
            gene: Gene symbol
            
        Returns:
            List of pathway names
        """
        pathways = []
        
        # Check external mappings
        if gene in self.external_mappings.get('gene_to_pathway', {}):
            pathways.extend(self.external_mappings['gene_to_pathway'][gene])
        
        # Check gene database
        if gene in self.gene_db.get('gene_to_pathway', {}):
            pathways.extend(self.gene_db['gene_to_pathway'][gene])
        
        # If no pathways found, use system pathways
        if not pathways:
            system = self.map_gene_to_system(gene)
            if system in SYSTEM_PATHWAY_MAP:
                # Add first two pathways as default
                pathways.extend(SYSTEM_PATHWAY_MAP[system][:2])
        
        return pathways
    
    def map_gene_to_parameters(self, gene: str) -> Dict[str, float]:
        """
        Map gene to GASLIT-AF parameters.
        
        Args:
            gene: Gene symbol
            
        Returns:
            Dictionary of parameter scores
        """
        # Initialize parameter scores
        params = {
            "γ": 0.0,  # Genetic fragility
            "Λ": 0.0,  # Allostatic load
            "Ω": 0.0,  # Endocannabinoid buffering
            "Χ": 0.0,  # Physiological coherence
            "σ": 0.0   # Entropy production
        }
        
        # Check external mappings
        if gene in self.external_mappings.get('gene_to_parameter', {}):
            for param, score in self.external_mappings['gene_to_parameter'][gene].items():
                if param in params:
                    params[param] = score
        
        # Check gene database
        if gene in self.gene_db.get('gene_to_parameter', {}):
            for param, score in self.gene_db['gene_to_parameter'][gene].items():
                if param in params:
                    params[param] = score
        
        # If no parameter scores found, estimate based on system
        if all(score == 0.0 for score in params.values()):
            system = self.map_gene_to_system(gene)
            
            # Set default parameter scores based on system
            if system == "Immune & Inflammatory Pathways":
                params["Λ"] = 0.7  # High allostatic load
                params["σ"] = 0.6  # High entropy production
            elif system == "Autonomic & Neurotransmitter Pathways":
                params["Χ"] = 0.7  # High physiological coherence
                params["Λ"] = 0.6  # Moderate allostatic load
            elif system == "Structural & Connective Tissue Integrity":
                params["γ"] = 0.7  # High genetic fragility
                params["Χ"] = 0.5  # Moderate physiological coherence
            elif system == "Metabolic, Mitochondrial & Oxidative Stress":
                params["σ"] = 0.8  # Very high entropy production
                params["Λ"] = 0.5  # Moderate allostatic load
            elif system == "Endocannabinoid System (ECS)":
                params["Ω"] = 0.9  # Very high endocannabinoid buffering
                params["Χ"] = 0.7  # High physiological coherence
            elif system == "Calcium & Ion Channels":
                params["Χ"] = 0.8  # Very high physiological coherence
                params["σ"] = 0.5  # Moderate entropy production
            elif system == "Mast Cell Activation & Histamine Metabolism":
                params["Λ"] = 0.8  # Very high allostatic load
                params["σ"] = 0.7  # High entropy production
            elif system == "Kynurenine Pathway":
                params["Ω"] = 0.6  # Moderate endocannabinoid buffering
                params["Λ"] = 0.7  # High allostatic load
        
        return params
    
    def map_variant_to_impact(self, variant: Dict[str, Any]) -> Dict[str, Any]:
        """
        Map variant to functional impact.
        
        Args:
            variant: Variant dictionary
            
        Returns:
            Variant dictionary with enhanced impact information
        """
        # Create a copy of the variant
        enhanced_variant = variant.copy()
        
        # Extract key information
        rsid = variant.get('variant_id', '')
        gene = variant.get('gene', '')
        classification = variant.get('classification', '')
        
        # Map rsID to gene if gene is missing
        if not gene and rsid:
            gene = self.map_rsid_to_gene(rsid)
            if gene:
                enhanced_variant['gene'] = gene
        
        # Map gene to system
        if gene:
            system = self.map_gene_to_system(gene)
            enhanced_variant['biological_system'] = system
            
            # Map gene to pathways
            pathways = self.map_gene_to_pathways(gene)
            enhanced_variant['pathways'] = pathways
            
            # Map gene to parameters
            param_scores = self.map_gene_to_parameters(gene)
            
            # Adjust parameter scores based on variant classification
            if classification == "Pathogenic":
                # Increase genetic fragility for pathogenic variants
                param_scores["γ"] = min(1.0, param_scores["γ"] + 0.3)
                # Increase entropy production
                param_scores["σ"] = min(1.0, param_scores["σ"] + 0.2)
            elif classification == "Likely pathogenic":
                # Moderately increase genetic fragility
                param_scores["γ"] = min(1.0, param_scores["γ"] + 0.2)
                # Moderately increase entropy production
                param_scores["σ"] = min(1.0, param_scores["σ"] + 0.1)
            
            enhanced_variant['parameter_scores'] = param_scores
        
        return enhanced_variant
    
    def enhance_variants(self, variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """
        Enhance variant data with improved gene mapping.
        
        Args:
            variants: List of variant dictionaries
            
        Returns:
            List of enhanced variant dictionaries
        """
        enhanced_variants = []
        
        log.info(f"Enhancing gene mapping for {len(variants)} variants")
        
        for variant in variants:
            enhanced_variant = self.map_variant_to_impact(variant)
            enhanced_variants.append(enhanced_variant)
        
        # Count variants by system
        system_counts = defaultdict(int)
        for variant in enhanced_variants:
            system = variant.get('biological_system', 'Other')
            system_counts[system] += 1
        
        # Log system distribution
        log.info("Enhanced variant distribution by biological system:")
        for system, count in sorted(system_counts.items(), key=lambda x: -x[1]):
            log.info(f"  - {system}: {count} variants")
        
        return enhanced_variants
    
    def create_gene_mapping_database(self, output_path: str,
                                   variant_data: Optional[pd.DataFrame] = None) -> None:
        """
        Create a gene mapping database from variant data and external resources.
        
        Args:
            output_path: Path to save the gene mapping database
            variant_data: DataFrame with variant data (optional)
        """
        # Initialize database structure
        db = {
            "rsid_to_gene": {},
            "gene_to_system": {},
            "gene_to_pathway": {},
            "gene_to_parameter": {},
            "gene_metadata": {}
        }
        
        # Add GASLIT-AF gene mappings
        for system, genes in BIOLOGICAL_SYSTEMS.items():
            for gene in genes:
                db["gene_to_system"][gene] = system
        
        # Add GASLIT-AF rsID mappings
        for gene, rsids in KNOWN_SNPS.items():
            for rsid in rsids:
                db["rsid_to_gene"][rsid] = gene
        
        # Add system pathway mappings
        for system, pathways in SYSTEM_PATHWAY_MAP.items():
            for gene in BIOLOGICAL_SYSTEMS.get(system, []):
                db["gene_to_pathway"][gene] = pathways[:3]  # Add top 3 pathways
        
        # Add parameter function mappings
        for param, functions in PARAMETER_FUNCTION_MAP.items():
            for system, genes in BIOLOGICAL_SYSTEMS.items():
                for gene in genes:
                    if gene not in db["gene_to_parameter"]:
                        db["gene_to_parameter"][gene] = {}
                    
                    # Set default parameter score based on system
                    if param == "γ" and system == "Structural & Connective Tissue Integrity":
                        db["gene_to_parameter"][gene][param] = 0.7
                    elif param == "Λ" and system == "Immune & Inflammatory Pathways":
                        db["gene_to_parameter"][gene][param] = 0.7
                    elif param == "Ω" and system == "Endocannabinoid System (ECS)":
                        db["gene_to_parameter"][gene][param] = 0.9
                    elif param == "Χ" and system == "Autonomic & Neurotransmitter Pathways":
                        db["gene_to_parameter"][gene][param] = 0.7
                    elif param == "σ" and system == "Metabolic, Mitochondrial & Oxidative Stress":
                        db["gene_to_parameter"][gene][param] = 0.8
                    else:
                        db["gene_to_parameter"][gene][param] = 0.3
        
        # Add variant data if provided
        if variant_data is not None and not variant_data.empty:
            # Extract unique genes and rsIDs
            genes = variant_data['gene'].dropna().unique()
            rsids = variant_data['rsid'].dropna().unique()
            
            # Add gene metadata
            for gene in genes:
                if gene and gene not in db["gene_metadata"]:
                    db["gene_metadata"][gene] = {
                        "symbol": gene,
                        "description": "",
                        "aliases": [],
                        "chromosome": "",
                        "position": "",
                        "strand": "",
                        "type": ""
                    }
            
            # Add rsID to gene mappings
            for _, row in variant_data.iterrows():
                rsid = row.get('rsid', '')
                gene = row.get('gene', '')
                if rsid and gene and rsid not in db["rsid_to_gene"]:
                    db["rsid_to_gene"][rsid] = gene
        
        # Save database to file
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(db, f, indent=2)
        
        log.info(f"Created gene mapping database at {output_path}")
        log.info(f"  - rsID to gene mappings: {len(db['rsid_to_gene'])}")
        log.info(f"  - Gene to system mappings: {len(db['gene_to_system'])}")
        log.info(f"  - Gene to pathway mappings: {len(db['gene_to_pathway'])}")
        log.info(f"  - Gene to parameter mappings: {len(db['gene_to_parameter'])}")
        log.info(f"  - Gene metadata entries: {len(db['gene_metadata'])}")

# Function to enhance gene mapping for variant data
def enhance_gene_mapping(variants: List[Dict[str, Any]],
                       gene_db_path: Optional[str] = None,
                       external_mappings_path: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    Enhance gene mapping for variant data.
    
    Args:
        variants: List of variant dictionaries
        gene_db_path: Path to gene database file (JSON)
        external_mappings_path: Path to external mappings file (JSON)
        
    Returns:
        List of enhanced variant dictionaries
    """
    mapper = GeneMapper(gene_db_path, external_mappings_path)
    return mapper.enhance_variants(variants)

# Function to create a gene mapping database
def create_gene_mapping_database(output_path: str,
                               variant_data: Optional[pd.DataFrame] = None) -> None:
    """
    Create a gene mapping database.
    
    Args:
        output_path: Path to save the gene mapping database
        variant_data: DataFrame with variant data (optional)
    """
    mapper = GeneMapper()
    mapper.create_gene_mapping_database(output_path, variant_data)
