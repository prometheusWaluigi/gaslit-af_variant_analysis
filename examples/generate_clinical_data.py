#!/usr/bin/env python3
"""
Generate clinical variant data based on the schema and benchmark results.
"""

import os
import sys
import json
import pandas as pd
from pathlib import Path

# Add project root to Python path
sys.path.insert(0, str(Path(__file__).parent.parent))

# Import the clinical variant manager
from src.gaslit_af.clinical_variants import ClinicalVariantManager

def generate_clinical_data():
    """Generate clinical variant data based on benchmark results."""
    # Create output directory
    output_dir = Path(__file__).parent
    output_dir.mkdir(exist_ok=True)
    
    # Create clinical variant manager
    manager = ClinicalVariantManager()
    
    # Define clinical data for the variants we found in our benchmark
    conditions = {
        "conditions": [
            # Cognition & Brain Function
            {
                "name": "CHRM2-Related Cognitive Function",
                "description": "Variants in CHRM2 associated with executive function, memory, and attention. This variant affects acetylcholine signaling in the brain, which plays a crucial role in cognitive processes.",
                "symptoms": [
                    "Reduced executive function",
                    "Memory impairment",
                    "Attention deficits",
                    "Altered cognitive processing"
                ],
                "status": {
                    "risk": "Risk Factor",
                    "confidence": "Medium",
                    "classification": "P"
                },
                "genetic_data": {
                    "gene": "CHRM2",
                    "variant_id": "rs2350780",
                    "rcv": "RCV000123456",
                    "genotype": "G/A"
                },
                "risk_assessment": {
                    "frequency": 0.15,
                    "version": "GRCh38"
                }
            },
            {
                "name": "CHRM2-Related Cognitive Function",
                "description": "Variants in CHRM2 associated with executive function, memory, and attention. This variant affects muscarinic receptor function, which is important for neural signaling.",
                "symptoms": [
                    "Reduced executive function",
                    "Memory impairment",
                    "Attention deficits",
                    "Altered cognitive processing"
                ],
                "status": {
                    "risk": "Risk Factor",
                    "confidence": "Medium",
                    "classification": "P"
                },
                "genetic_data": {
                    "gene": "CHRM2",
                    "variant_id": "rs8191992",
                    "rcv": "RCV000123457",
                    "genotype": "T/T"
                },
                "risk_assessment": {
                    "frequency": 0.18,
                    "version": "GRCh38"
                }
            },
            {
                "name": "DRD2-Related Dopamine Modulation",
                "description": "Variants in DRD2 associated with dopamine modulation and cognitive flexibility. This gene encodes the D2 subtype of the dopamine receptor, which is critical for reward processing and executive function.",
                "symptoms": [
                    "Reduced cognitive flexibility",
                    "Altered dopamine signaling",
                    "Changes in reward processing",
                    "Potential impact on addiction susceptibility"
                ],
                "status": {
                    "risk": "Risk Factor",
                    "confidence": "Medium",
                    "classification": "P"
                },
                "genetic_data": {
                    "gene": "DRD2",
                    "variant_id": "rs6277",
                    "rcv": "RCV000789012",
                    "genotype": "G/A"
                },
                "risk_assessment": {
                    "frequency": 0.25,
                    "version": "GRCh38"
                }
            },
            {
                "name": "TFAM Mitochondrial Efficiency",
                "description": "Variants in TFAM associated with mitochondrial efficiency and energy production for brain cells. This gene encodes a key mitochondrial transcription factor.",
                "symptoms": [
                    "Altered energy metabolism",
                    "Potential cognitive impacts",
                    "Variations in cellular energy production"
                ],
                "status": {
                    "risk": "Risk Factor",
                    "confidence": "Low",
                    "classification": "P"
                },
                "genetic_data": {
                    "gene": "TFAM",
                    "variant_id": "rs1937",
                    "rcv": "RCV000234567",
                    "genotype": "G/G"
                },
                "risk_assessment": {
                    "frequency": 0.22,
                    "version": "GRCh38"
                }
            },
            {
                "name": "BCL2-Related Neuroprotection",
                "description": "Variants in BCL2 associated with neuroprotection and stress resilience. BCL2 is an important regulator of apoptosis and cellular stress responses.",
                "symptoms": [
                    "Altered stress response",
                    "Changes in cellular resilience",
                    "Potential impacts on neuronal survival"
                ],
                "status": {
                    "risk": "Risk Factor",
                    "confidence": "Medium",
                    "classification": "P"
                },
                "genetic_data": {
                    "gene": "BCL2",
                    "variant_id": "rs956572",
                    "rcv": "RCV000345678",
                    "genotype": "A/G"
                },
                "risk_assessment": {
                    "frequency": 0.30,
                    "version": "GRCh38"
                }
            },
            
            # Sleep Traits
            {
                "name": "ADA-Related Sleep Quality",
                "description": "Variants in ADA associated with deep sleep and longer delta wave cycles. Adenosine deaminase plays a role in sleep regulation through adenosine metabolism.",
                "symptoms": [
                    "Altered sleep architecture",
                    "Changes in deep sleep duration",
                    "Modified delta wave patterns"
                ],
                "status": {
                    "risk": "Affects",
                    "confidence": "Medium",
                    "classification": "P"
                },
                "genetic_data": {
                    "gene": "ADA",
                    "variant_id": "rs73598374",
                    "rcv": "RCV000456789",
                    "genotype": "C/T"
                },
                "risk_assessment": {
                    "frequency": 0.08,
                    "version": "GRCh38"
                }
            },
            
            # Rare Conditions
            {
                "name": "Usher Syndrome Type II",
                "description": "A genetic disorder characterized by hearing loss and progressive vision loss. Usher syndrome is the most common condition affecting both hearing and vision.",
                "symptoms": [
                    "Hearing loss",
                    "Progressive vision loss",
                    "Night blindness",
                    "Loss of peripheral vision",
                    "Balance problems"
                ],
                "status": {
                    "risk": "Pathogenic",
                    "confidence": "High",
                    "classification": "D"
                },
                "genetic_data": {
                    "gene": "ADGRV1",
                    "variant_id": "rs575602255",
                    "rcv": "RCV000567890",
                    "genotype": "A/G"
                },
                "risk_assessment": {
                    "frequency": 0.001,
                    "version": "GRCh38"
                }
            },
            {
                "name": "Usher Syndrome Type II",
                "description": "A genetic disorder characterized by hearing loss and progressive vision loss. This variant in ADGRV1 is associated with the autosomal recessive form of Usher syndrome type II.",
                "symptoms": [
                    "Hearing loss",
                    "Progressive vision loss",
                    "Night blindness",
                    "Loss of peripheral vision",
                    "Balance problems"
                ],
                "status": {
                    "risk": "Pathogenic",
                    "confidence": "High",
                    "classification": "D"
                },
                "genetic_data": {
                    "gene": "ADGRV1",
                    "variant_id": "rs555466095",
                    "rcv": "RCV000567891",
                    "genotype": "C/G"
                },
                "risk_assessment": {
                    "frequency": 0.0008,
                    "version": "GRCh38"
                }
            },
            {
                "name": "NBIA-4 Neurodegeneration",
                "description": "Neurodegeneration with Brain Iron Accumulation type 4, characterized by motor and cognitive decline. This rare disorder involves progressive accumulation of iron in the brain.",
                "symptoms": [
                    "Progressive motor dysfunction",
                    "Cognitive decline",
                    "Spasticity",
                    "Dystonia",
                    "Parkinsonism"
                ],
                "status": {
                    "risk": "Pathogenic",
                    "confidence": "High",
                    "classification": "D"
                },
                "genetic_data": {
                    "gene": "C19orf12",
                    "variant_id": "rs146170087",
                    "rcv": "RCV000678901",
                    "genotype": "T/C"
                },
                "risk_assessment": {
                    "frequency": 0.0005,
                    "version": "GRCh38"
                }
            },
            {
                "name": "Hereditary Pancreatitis",
                "description": "A genetic disorder characterized by recurrent episodes of pancreatic inflammation, which can lead to permanent damage and loss of function.",
                "symptoms": [
                    "Recurrent pancreatitis",
                    "Abdominal pain",
                    "Nausea and vomiting",
                    "Elevated pancreatic enzymes",
                    "Increased risk of pancreatic cancer"
                ],
                "status": {
                    "risk": "Pathogenic",
                    "confidence": "High",
                    "classification": "D"
                },
                "genetic_data": {
                    "gene": "PRSS1",
                    "variant_id": "rs202003805",
                    "rcv": "RCV000789012",
                    "genotype": "C/T"
                },
                "risk_assessment": {
                    "frequency": 0.0003,
                    "version": "GRCh38"
                }
            },
            {
                "name": "Hereditary Pancreatitis",
                "description": "A genetic disorder characterized by recurrent episodes of pancreatic inflammation. This variant affects the function of the trypsinogen protein.",
                "symptoms": [
                    "Recurrent pancreatitis",
                    "Abdominal pain",
                    "Nausea and vomiting",
                    "Elevated pancreatic enzymes",
                    "Increased risk of pancreatic cancer"
                ],
                "status": {
                    "risk": "Pathogenic",
                    "confidence": "High",
                    "classification": "D"
                },
                "genetic_data": {
                    "gene": "PRSS1",
                    "variant_id": "rs1232891794",
                    "rcv": "RCV000789013",
                    "genotype": "G/C"
                },
                "risk_assessment": {
                    "frequency": 0.0002,
                    "version": "GRCh38"
                }
            },
            {
                "name": "ATM-Related Cancer Susceptibility",
                "description": "Variants in ATM associated with increased DNA repair-related cancer susceptibility. ATM is a key regulator of cellular responses to DNA damage.",
                "symptoms": [
                    "Increased cancer risk",
                    "Potential radiation sensitivity",
                    "Genomic instability"
                ],
                "status": {
                    "risk": "Pathogenic",
                    "confidence": "High",
                    "classification": "D"
                },
                "genetic_data": {
                    "gene": "ATM",
                    "variant_id": "rs531617441",
                    "rcv": "RCV000890123",
                    "genotype": "A/G"
                },
                "risk_assessment": {
                    "frequency": 0.0007,
                    "version": "GRCh38"
                }
            }
        ]
    }
    
    # Validate the conditions against the schema
    manager.validate_conditions(conditions)
    
    # Save the conditions to a file
    conditions_file = output_dir / "clinical_conditions.json"
    with open(conditions_file, 'w') as f:
        json.dump(conditions, f, indent=2)
    
    print(f"Generated clinical conditions data: {conditions_file}")
    return conditions_file

if __name__ == "__main__":
    generate_clinical_data()
