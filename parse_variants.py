#!/usr/bin/env python3
"""
Variant Parser for GASLIT-AF Framework

This script parses a variants.txt file and converts it to a structured JSON format
that creates a quantum coherence bridge between personal genomic architecture
and the GASLIT-AF theoretical model parameters.

The JSON output is stored in a gitignored location to maintain privacy while
enabling recursive integration with the GASLIT-AF analysis framework.
"""

import os
import re
import json
import logging
import argparse
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af")

# Define GASLIT-AF parameter mappings
GASLIT_PARAMETERS = {
    "γ": ["genetic_fragility", "pathogenic", "likely_pathogenic", "risk_allele"],
    "Λ": ["allostatic_load", "penetrance", "expressivity", "confidence"],
    "Ω": ["endocannabinoid_buffering", "protective_allele", "resilience"],
    "Χ": ["physiological_coherence", "system_impact", "pathway_disruption"],
    "σ": ["entropy_production", "variant_impact", "functional_consequence"]
}

def parse_args():
    """Parse command-line arguments for variant parsing."""
    parser = argparse.ArgumentParser(description="Variant Parser for GASLIT-AF Framework")
    
    parser.add_argument("--input", type=str, default="variants.txt",
                      help="Input variants file (default: variants.txt)")
    parser.add_argument("--output", type=str, default="data/personal_variants.json",
                      help="Output JSON file (default: data/personal_variants.json)")
    parser.add_argument("--format", type=str, choices=["standard", "recursive"], default="recursive",
                      help="Output format (default: recursive)")
    parser.add_argument("--include-gaslit-params", action="store_true",
                      help="Include GASLIT-AF parameter mappings in output")
    
    return parser.parse_args()

def extract_variant_blocks(content: str) -> List[str]:
    """
    Extract individual variant blocks from the content.
    
    Args:
        content: Raw file content
        
    Returns:
        List of variant blocks
    """
    # Split content by condition headers or empty lines
    # This is a heuristic approach and may need adjustment based on exact file format
    blocks = []
    current_block = []
    in_block = False
    
    lines = content.split('\n')
    for i, line in enumerate(lines):
        # Skip header lines
        if i < 10:
            continue
            
        # New condition starts a new block
        if re.match(r'^[A-Z][a-zA-Z\s]+$', line.strip()) and len(line.strip()) > 0:
            if in_block and current_block:
                blocks.append('\n'.join(current_block))
                current_block = []
            in_block = True
            current_block.append(line.strip())
        elif in_block:
            current_block.append(line.strip())
    
    # Add the last block
    if current_block:
        blocks.append('\n'.join(current_block))
    
    return blocks

def parse_variant_block(block: str) -> Dict[str, Any]:
    """
    Parse a variant block into a structured dictionary.
    
    Args:
        block: Variant block text
        
    Returns:
        Dictionary with parsed variant information
    """
    lines = block.split('\n')
    variant_data = {
        "condition": "",
        "status": "",
        "classification": "",
        "confidence": "",
        "gene": "",
        "variant_id": "",
        "rcv": "",
        "genotype": "",
        "frequency": "",
        "risk_allele": "",
        "description": "",
        "symptoms": []
    }
    
    # Extract condition from first line
    if lines and lines[0].strip():
        variant_data["condition"] = lines[0].strip()
    
    # Extract other fields
    for i, line in enumerate(lines):
        line = line.strip()
        
        # Skip empty lines
        if not line:
            continue
            
        # Extract status
        if line in ["Risk (D)", "Carrier (R)", "Not a carrier", "Typical"]:
            variant_data["status"] = line
            
        # Extract classification
        elif line in ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]:
            variant_data["classification"] = line
            
        # Extract confidence
        elif line.startswith(("High", "Medium", "Low")) and "(" in line:
            variant_data["confidence"] = line
            
        # Extract gene and variant ID
        elif re.match(r'^[A-Z0-9]+$', line) and i+1 < len(lines) and lines[i+1].strip().startswith("rs"):
            variant_data["gene"] = line
            variant_data["variant_id"] = lines[i+1].strip()
            
        # Extract RCV
        elif line.startswith("RCV"):
            variant_data["rcv"] = line
            
        # Extract genotype
        elif re.match(r'^[ACGT]{1,2}$', line):
            variant_data["genotype"] = line
            
        # Extract frequency
        elif re.match(r'^0\.\d+$', line):
            variant_data["frequency"] = line
            
        # Extract risk allele
        elif re.match(r'^[ACGT]$', line) and variant_data["frequency"]:
            variant_data["risk_allele"] = line
            
        # Extract description hint
        elif "description" in line.lower():
            variant_data["description"] = "Available"
            
        # Extract symptoms hint
        elif "symptoms" in line.lower():
            variant_data["symptoms"] = ["Available"]
    
    return variant_data

def map_to_gaslit_parameters(variant: Dict[str, Any]) -> Dict[str, Any]:
    """
    Map variant data to GASLIT-AF parameters.
    
    Args:
        variant: Parsed variant data
        
    Returns:
        Variant data with GASLIT-AF parameter mappings
    """
    # Create a copy of the variant data
    enriched_variant = variant.copy()
    
    # Add GASLIT-AF parameter mappings
    gaslit_params = {}
    
    # γ (genetic fragility)
    gamma_score = 0.0
    if variant["classification"] == "Pathogenic":
        gamma_score = 1.0
    elif variant["classification"] == "Likely pathogenic":
        gamma_score = 0.8
    elif variant["classification"] == "Uncertain significance":
        gamma_score = 0.5
    elif variant["classification"] == "Likely benign":
        gamma_score = 0.2
    elif variant["classification"] == "Benign":
        gamma_score = 0.0
        
    # Adjust by confidence
    confidence_multiplier = 1.0
    if "High" in variant["confidence"]:
        confidence_multiplier = 1.0
    elif "Medium" in variant["confidence"]:
        confidence_multiplier = 0.8
    elif "Low" in variant["confidence"]:
        confidence_multiplier = 0.6
        
    gamma_score *= confidence_multiplier
    gaslit_params["γ"] = round(gamma_score, 2)
    
    # Λ (allostatic load)
    lambda_score = 0.0
    if variant["status"] == "Risk (D)":
        lambda_score = 0.8
    elif variant["status"] == "Carrier (R)":
        lambda_score = 0.4
    elif variant["status"] == "Not a carrier":
        lambda_score = 0.0
    elif variant["status"] == "Typical":
        lambda_score = 0.2
        
    # Adjust by frequency
    if variant["frequency"]:
        try:
            freq = float(variant["frequency"])
            if freq < 0.01:
                lambda_score *= 1.2  # Rare variants may have higher impact
            elif freq > 0.1:
                lambda_score *= 0.8  # Common variants may have lower impact
        except ValueError:
            pass
            
    gaslit_params["Λ"] = round(lambda_score, 2)
    
    # Ω (endocannabinoid buffering)
    # This is more speculative and would require additional data
    omega_score = 0.5  # Default neutral value
    gaslit_params["Ω"] = round(omega_score, 2)
    
    # Χ (physiological coherence)
    # Based on gene function and pathway involvement
    chi_score = 0.5  # Default neutral value
    gaslit_params["Χ"] = round(chi_score, 2)
    
    # σ (entropy production)
    # Based on variant impact and functional consequence
    sigma_score = gamma_score * lambda_score
    gaslit_params["σ"] = round(sigma_score, 2)
    
    enriched_variant["gaslit_parameters"] = gaslit_params
    return enriched_variant

def create_standard_output(variants: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Create standard JSON output format.
    
    Args:
        variants: List of parsed variants
        
    Returns:
        Standard format JSON object
    """
    return {
        "format_version": "1.0",
        "variant_count": len(variants),
        "variants": variants
    }

def create_recursive_output(variants: List[Dict[str, Any]]) -> Dict[str, Any]:
    """
    Create recursive JSON output format aligned with GASLIT-AF model.
    
    Args:
        variants: List of parsed variants
        
    Returns:
        Recursive format JSON object
    """
    # Group variants by biological system
    systems = {
        "immune_inflammatory": [],
        "autonomic_neurotransmitter": [],
        "structural_connective": [],
        "metabolic": [],
        "endocannabinoid": [],
        "calcium_ion_channels": [],
        "mast_cell_activation": [],
        "kynurenine_pathway": [],
        "other": []
    }
    
    # Map genes to systems (simplified mapping)
    gene_to_system = {
        "APOE": "metabolic",
        "PCSK9": "metabolic",
        "CHRM2": "autonomic_neurotransmitter",
        "DRD2": "autonomic_neurotransmitter",
        "ADGRV1": "structural_connective",
        "C19orf12": "metabolic"
        # Add more mappings as needed
    }
    
    # Group variants by system
    for variant in variants:
        gene = variant.get("gene", "")
        system = gene_to_system.get(gene, "other")
        systems[system].append(variant)
    
    # Create recursive structure
    return {
        "format_version": "1.0",
        "variant_count": len(variants),
        "gaslit_af_model": {
            "parameters": {
                "γ": {
                    "description": "Genetic fragility",
                    "variants": sorted([v for v in variants if v.get("gaslit_parameters", {}).get("γ", 0) > 0.5], 
                                     key=lambda x: x.get("gaslit_parameters", {}).get("γ", 0), reverse=True)
                },
                "Λ": {
                    "description": "Allostatic load",
                    "variants": sorted([v for v in variants if v.get("gaslit_parameters", {}).get("Λ", 0) > 0.5],
                                     key=lambda x: x.get("gaslit_parameters", {}).get("Λ", 0), reverse=True)
                },
                "Ω": {
                    "description": "Endocannabinoid buffering capacity",
                    "variants": sorted([v for v in variants if v.get("gaslit_parameters", {}).get("Ω", 0) > 0.5],
                                     key=lambda x: x.get("gaslit_parameters", {}).get("Ω", 0), reverse=True)
                },
                "Χ": {
                    "description": "Physiological coherence",
                    "variants": sorted([v for v in variants if v.get("gaslit_parameters", {}).get("Χ", 0) > 0.5],
                                     key=lambda x: x.get("gaslit_parameters", {}).get("Χ", 0), reverse=True)
                },
                "σ": {
                    "description": "Entropy production",
                    "variants": sorted([v for v in variants if v.get("gaslit_parameters", {}).get("σ", 0) > 0.5],
                                     key=lambda x: x.get("gaslit_parameters", {}).get("σ", 0), reverse=True)
                }
            },
            "biological_systems": {
                "immune_inflammatory": {
                    "description": "Immune & Inflammatory System",
                    "variants": systems["immune_inflammatory"]
                },
                "autonomic_neurotransmitter": {
                    "description": "Autonomic & Neurotransmitter System",
                    "variants": systems["autonomic_neurotransmitter"]
                },
                "structural_connective": {
                    "description": "Structural & Connective Tissue",
                    "variants": systems["structural_connective"]
                },
                "metabolic": {
                    "description": "Metabolic System",
                    "variants": systems["metabolic"]
                },
                "endocannabinoid": {
                    "description": "Endocannabinoid System",
                    "variants": systems["endocannabinoid"]
                },
                "calcium_ion_channels": {
                    "description": "Calcium & Ion Channels",
                    "variants": systems["calcium_ion_channels"]
                },
                "mast_cell_activation": {
                    "description": "Mast Cell Activation",
                    "variants": systems["mast_cell_activation"]
                },
                "kynurenine_pathway": {
                    "description": "Kynurenine Pathway",
                    "variants": systems["kynurenine_pathway"]
                },
                "other": {
                    "description": "Other Systems",
                    "variants": systems["other"]
                }
            }
        }
    }

def main():
    """Main entry point for variant parsing."""
    args = parse_args()
    
    # Ensure input file exists
    input_path = Path(args.input)
    if not input_path.exists():
        log.error(f"Input file not found: {input_path}")
        return
    
    # Create output directory if it doesn't exist
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Read input file
    log.info(f"Reading variant data from {input_path}")
    with open(input_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Extract variant blocks
    log.info("Extracting variant blocks")
    blocks = extract_variant_blocks(content)
    log.info(f"Found {len(blocks)} variant blocks")
    
    # Parse variant blocks
    variants = []
    with Progress(
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        TimeElapsedColumn()
    ) as progress:
        task = progress.add_task("[cyan]Parsing variants...", total=len(blocks))
        
        for block in blocks:
            variant_data = parse_variant_block(block)
            
            # Map to GASLIT-AF parameters if requested
            if args.include_gaslit_params:
                variant_data = map_to_gaslit_parameters(variant_data)
                
            variants.append(variant_data)
            progress.update(task, advance=1)
    
    # Create output based on format
    if args.format == "standard":
        output_data = create_standard_output(variants)
    else:  # recursive
        output_data = create_recursive_output(variants)
    
    # Write output file
    log.info(f"Writing {len(variants)} variants to {output_path}")
    with open(output_path, 'w', encoding='utf-8') as f:
        json.dump(output_data, f, indent=2)
    
    # Update .gitignore to ensure personal data is not committed
    gitignore_path = Path(".gitignore")
    gitignore_entry = str(output_path)
    
    if gitignore_path.exists():
        with open(gitignore_path, 'r', encoding='utf-8') as f:
            gitignore_content = f.read()
            
        if gitignore_entry not in gitignore_content:
            log.info(f"Adding {gitignore_entry} to .gitignore")
            with open(gitignore_path, 'a', encoding='utf-8') as f:
                f.write(f"\n# Personal variant data\n{gitignore_entry}\n")
    else:
        log.info(f"Creating .gitignore with {gitignore_entry}")
        with open(gitignore_path, 'w', encoding='utf-8') as f:
            f.write(f"# Personal variant data\n{gitignore_entry}\n")
    
    log.info(f"[bold green]Variant parsing completed successfully![/]")
    log.info(f"Parsed {len(variants)} variants into {args.format} format")
    log.info(f"Output saved to {output_path} (added to .gitignore)")
    
    # Provide next steps
    log.info("\nTo analyze your variants with GASLIT-AF, run:")
    log.info(f"python analyze_personal_variants.py --variant-file {output_path}")

if __name__ == "__main__":
    main()
