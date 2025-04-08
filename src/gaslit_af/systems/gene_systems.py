"""
Gene Systems Module for GASLIT-AF Variant Analysis.

This module provides a modular, JSON-based approach to organizing genes by biological systems.
Each system can be loaded independently, allowing for recursive composition of analysis patterns.
"""

import os
import json
import logging
from pathlib import Path
from typing import Dict, List, Set, Optional, Any, Union

# Configure logging
log = logging.getLogger("gaslit-af")

class GeneSystemManager:
    """Manages gene systems with a lightweight, modular approach."""
    
    def __init__(self, systems_dir: Optional[Path] = None):
        """
        Initialize the gene system manager.
        
        Args:
            systems_dir: Directory containing gene system JSON files
        """
        self.systems_dir = Path(systems_dir) if systems_dir else Path("./data/systems")
        self.systems_dir.mkdir(parents=True, exist_ok=True)
        self.systems = {}
        self.gene_to_system = {}
        self.load_systems()
    
    def load_systems(self):
        """Load all gene systems from JSON files."""
        log.info(f"Loading gene systems from {self.systems_dir}")
        
        # Default systems if no files exist
        if not any(self.systems_dir.glob("*.json")):
            self._create_default_systems()
        
        # Load all JSON files
        for system_file in self.systems_dir.glob("*.json"):
            try:
                with open(system_file, 'r') as f:
                    system_data = json.load(f)
                    system_name = system_file.stem
                    self.systems[system_name] = system_data
                    
                    # Create reverse mapping
                    if "genes" in system_data:
                        for gene in system_data["genes"]:
                            self.gene_to_system[gene] = system_name
                    
                log.info(f"Loaded system {system_name} with {len(system_data.get('genes', []))} genes")
            except Exception as e:
                log.error(f"Error loading system file {system_file}: {e}")
    
    def _create_default_systems(self):
        """Create default gene systems if none exist."""
        default_systems = {
            "immune_inflammatory": {
                "name": "Immune & Inflammatory System",
                "description": "Genes involved in immune response and inflammation",
                "genes": ["IDO2", "AHR", "AHRR", "IL36RN", "CFH", "MBL2", "NLRP3", 
                         "IL1B", "IL6", "IL17", "IL13", "IL4", "HLA-DQB1", "PTPN22", 
                         "CTLA4", "ASXL1", "CBL", "DNMT3B", "ETV6", "IDH1", "IL6R"]
            },
            "autonomic_neurotransmitter": {
                "name": "Autonomic & Neurotransmitter System",
                "description": "Genes involved in autonomic function and neurotransmission",
                "genes": ["COMT", "CHRM2", "DRD2", "GABRA1", "CHRNA7", "ADRB1", "ADRB2", 
                         "NOS3", "GNB3", "SLC6A2", "NET", "EZH2", "SLC6A4", "HTR2A", 
                         "TAAR1", "OPRM1", "GCH1", "TRPV2", "MYT1L", "NRXN3"]
            },
            "structural_connective": {
                "name": "Structural & Connective Tissue",
                "description": "Genes involved in structural integrity and connective tissue",
                "genes": ["TNXB", "ADAMTS10", "SELENON", "NEB", "MYH7", "MAPRE1", "ADGRV1", 
                         "PLXNA2", "COL3A1", "FBN1", "FLNA", "COL5A1", "FKBP14", "PLOD1", 
                         "CDON", "SULF2"]
            },
            "metabolic": {
                "name": "Metabolic System",
                "description": "Genes involved in metabolic processes",
                "genes": ["APOE", "PCSK9", "UGT1A1", "HNF1A", "ABCC8", "TFAM", "C19orf12", 
                         "MT-ATP6", "MT-ATP8", "PDHA1", "SDHB", "NAMPT", "NMRK1", "PGC1A", 
                         "PRKAA2"]
            },
            "endocannabinoid": {
                "name": "Endocannabinoid System",
                "description": "Genes involved in endocannabinoid signaling",
                "genes": ["CNR1", "CNR2", "FAAH", "MGLL"]
            },
            "calcium_ion_channels": {
                "name": "Calcium & Ion Channels",
                "description": "Genes involved in calcium signaling and ion channel function",
                "genes": ["ITPR1", "KCNJ5", "RYR2", "KCNA5", "KCND3", "KCNE1", "KCNQ1", 
                         "HCN4", "CAMK2B"]
            },
            "mast_cell": {
                "name": "Mast Cell Activation",
                "description": "Genes involved in mast cell activation and histamine metabolism",
                "genes": ["TPSAB1", "KIT", "HNMT", "TET2"]
            },
            "kynurenine": {
                "name": "Kynurenine Pathway",
                "description": "Genes involved in kynurenine pathway metabolism",
                "genes": ["IDO1", "KMO", "KYNU", "TDO2", "HAAO", "ARNT", "BECN1", "ATG5"]
            },
            "vascular_ras": {
                "name": "Vascular & RAS System",
                "description": "Genes involved in vascular function and renin-angiotensin system",
                "genes": ["ROCK1", "ROCK2", "ARG1", "ACE", "ACE2", "TGFβ1", "TGFβ2", "TGFβ3", 
                         "GDF-15", "NPPA"]
            },
            "mitochondrial_stress": {
                "name": "Mitochondrial & Cellular Stress",
                "description": "Genes involved in mitochondrial function and cellular stress response",
                "genes": ["DRP1", "SIRT1", "IFNL1", "PGE2", "ATG13", "NEFL", "S100B", "TWEAK", 
                         "S100PBP", "AKAP1", "USP6NL"]
            },
            "cardiac_development": {
                "name": "Cardiac Development & Conduction",
                "description": "Genes involved in cardiac development and electrical conduction",
                "genes": ["PITX2", "SPEN", "KIAA1755", "GATA4", "GATA5", "GATA6", "TBX3", "TBX5", 
                         "NKX2-5", "ZFHX3", "GREM2", "NPPA", "SCN5A", "SH3PXD2A", "MYL4", "LMNA"]
            },
            "mecfs_postviral": {
                "name": "ME/CFS & Post-Viral Syndromes",
                "description": "Genes associated with ME/CFS and post-viral syndromes",
                "genes": ["S100PBP", "AKAP1", "USP6NL", "CDON", "SULF2"]
            }
        }
        
        # Save default systems
        for system_id, system_data in default_systems.items():
            system_file = self.systems_dir / f"{system_id}.json"
            with open(system_file, 'w') as f:
                json.dump(system_data, f, indent=2)
            
            # Add to in-memory systems
            self.systems[system_id] = system_data
            
            # Create reverse mapping
            for gene in system_data["genes"]:
                self.gene_to_system[gene] = system_id
        
        log.info(f"Created {len(default_systems)} default gene systems")
    
    def get_system_for_gene(self, gene: str) -> str:
        """
        Get the system for a gene.
        
        Args:
            gene: Gene symbol
            
        Returns:
            System ID or "unknown" if not found
        """
        return self.gene_to_system.get(gene, "unknown")
    
    def get_genes_for_system(self, system_id: str) -> List[str]:
        """
        Get all genes for a system.
        
        Args:
            system_id: System identifier
            
        Returns:
            List of gene symbols
        """
        if system_id in self.systems and "genes" in self.systems[system_id]:
            return self.systems[system_id]["genes"]
        return []
    
    def get_all_genes(self) -> Set[str]:
        """
        Get all genes across all systems.
        
        Returns:
            Set of all gene symbols
        """
        all_genes = set()
        for system_data in self.systems.values():
            if "genes" in system_data:
                all_genes.update(system_data["genes"])
        return all_genes
    
    def add_gene_to_system(self, gene: str, system_id: str, save: bool = True) -> bool:
        """
        Add a gene to a system.
        
        Args:
            gene: Gene symbol
            system_id: System identifier
            save: Whether to save changes to file
            
        Returns:
            True if successful, False otherwise
        """
        if system_id not in self.systems:
            log.error(f"System {system_id} not found")
            return False
        
        if "genes" not in self.systems[system_id]:
            self.systems[system_id]["genes"] = []
        
        if gene not in self.systems[system_id]["genes"]:
            self.systems[system_id]["genes"].append(gene)
            self.gene_to_system[gene] = system_id
            
            if save:
                self._save_system(system_id)
            
            log.info(f"Added gene {gene} to system {system_id}")
            return True
        
        return False
    
    def _save_system(self, system_id: str) -> bool:
        """
        Save a system to its JSON file.
        
        Args:
            system_id: System identifier
            
        Returns:
            True if successful, False otherwise
        """
        if system_id not in self.systems:
            return False
        
        system_file = self.systems_dir / f"{system_id}.json"
        try:
            with open(system_file, 'w') as f:
                json.dump(self.systems[system_id], f, indent=2)
            return True
        except Exception as e:
            log.error(f"Error saving system {system_id}: {e}")
            return False
    
    def create_system(self, system_id: str, name: str, description: str, genes: List[str] = None) -> bool:
        """
        Create a new gene system.
        
        Args:
            system_id: System identifier
            name: Display name
            description: System description
            genes: List of gene symbols
            
        Returns:
            True if successful, False otherwise
        """
        if system_id in self.systems:
            log.error(f"System {system_id} already exists")
            return False
        
        self.systems[system_id] = {
            "name": name,
            "description": description,
            "genes": genes or []
        }
        
        # Update reverse mapping
        for gene in self.systems[system_id]["genes"]:
            self.gene_to_system[gene] = system_id
        
        # Save to file
        return self._save_system(system_id)


# Singleton instance
_gene_system_manager = None

def get_gene_system_manager(systems_dir: Optional[Path] = None) -> GeneSystemManager:
    """
    Get the singleton gene system manager instance.
    
    Args:
        systems_dir: Directory containing gene system JSON files
        
    Returns:
        GeneSystemManager instance
    """
    global _gene_system_manager
    if _gene_system_manager is None:
        _gene_system_manager = GeneSystemManager(systems_dir)
    return _gene_system_manager
