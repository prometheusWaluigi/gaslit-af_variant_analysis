#!/usr/bin/env python3
"""
GASLIT-AF Parameter Mapping

This module establishes a quantum coherence bridge between genomic variant
distributions and the core theoretical parameters of the GASLIT-AF model:

γ = genetic fragility
Λ = allostatic load
Ω = endocannabinoid buffering capacity
Χ = physiological coherence
σ = entropy production

The module creates a recursive mapping that quantifies how each biological
system's variant load influences specific model parameters, establishing
a fractal connection between genomic architecture and physiological expression.
"""

import os
import sys
import json
import logging
import argparse
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from rich.console import Console
from rich.logging import RichHandler
from rich.progress import Progress

# Configure rich console and logging
console = Console()
logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(rich_tracebacks=True, console=console)]
)
log = logging.getLogger("gaslit-af-parameter-mapping")

# Define GASLIT-AF parameters
PARAMETERS = {
    'γ': 'genetic fragility',
    'Λ': 'allostatic load',
    'Ω': 'endocannabinoid buffering capacity',
    'Χ': 'physiological coherence',
    'σ': 'entropy production'
}

# Define parameter influence weights for each biological system
# These weights represent the recursive mapping between biological systems
# and GASLIT-AF theoretical parameters
SYSTEM_PARAMETER_WEIGHTS = {
    "Immune & Inflammatory Pathways": {
        'γ': 0.15,  # Moderate influence on genetic fragility
        'Λ': 0.40,  # Strong influence on allostatic load
        'Ω': 0.10,  # Minor influence on endocannabinoid buffering
        'Χ': 0.20,  # Moderate influence on physiological coherence
        'σ': 0.15   # Moderate influence on entropy production
    },
    "Autonomic & Neurotransmitter Pathways": {
        'γ': 0.10,  # Minor influence on genetic fragility
        'Λ': 0.25,  # Moderate influence on allostatic load
        'Ω': 0.15,  # Moderate influence on endocannabinoid buffering
        'Χ': 0.40,  # Strong influence on physiological coherence
        'σ': 0.10   # Minor influence on entropy production
    },
    "Structural & Connective Tissue Integrity": {
        'γ': 0.50,  # Strong influence on genetic fragility
        'Λ': 0.15,  # Moderate influence on allostatic load
        'Ω': 0.05,  # Minor influence on endocannabinoid buffering
        'Χ': 0.20,  # Moderate influence on physiological coherence
        'σ': 0.10   # Minor influence on entropy production
    },
    "Metabolic, Mitochondrial & Oxidative Stress": {
        'γ': 0.10,  # Minor influence on genetic fragility
        'Λ': 0.15,  # Moderate influence on allostatic load
        'Ω': 0.15,  # Moderate influence on endocannabinoid buffering
        'Χ': 0.10,  # Minor influence on physiological coherence
        'σ': 0.50   # Strong influence on entropy production
    },
    "Endocannabinoid System (ECS)": {
        'γ': 0.05,  # Minor influence on genetic fragility
        'Λ': 0.15,  # Moderate influence on allostatic load
        'Ω': 0.60,  # Strong influence on endocannabinoid buffering
        'Χ': 0.15,  # Moderate influence on physiological coherence
        'σ': 0.05   # Minor influence on entropy production
    },
    "Calcium & Ion Channels": {
        'γ': 0.10,  # Minor influence on genetic fragility
        'Λ': 0.15,  # Moderate influence on allostatic load
        'Ω': 0.10,  # Minor influence on endocannabinoid buffering
        'Χ': 0.50,  # Strong influence on physiological coherence
        'σ': 0.15   # Moderate influence on entropy production
    },
    "Mast Cell Activation & Histamine Metabolism": {
        'γ': 0.10,  # Minor influence on genetic fragility
        'Λ': 0.40,  # Strong influence on allostatic load
        'Ω': 0.20,  # Moderate influence on endocannabinoid buffering
        'Χ': 0.20,  # Moderate influence on physiological coherence
        'σ': 0.10   # Minor influence on entropy production
    },
    "Kynurenine Pathway": {
        'γ': 0.10,  # Minor influence on genetic fragility
        'Λ': 0.25,  # Moderate influence on allostatic load
        'Ω': 0.15,  # Moderate influence on endocannabinoid buffering
        'Χ': 0.15,  # Moderate influence on physiological coherence
        'σ': 0.35   # Strong influence on entropy production
    }
}

class ParameterMapper:
    """
    GASLIT-AF Parameter Mapper
    
    This class establishes a quantum coherence bridge between genomic variant
    distributions and the core theoretical parameters of the GASLIT-AF model.
    """
    
    def __init__(self, system_analysis_path: str):
        """
        Initialize the Parameter Mapper.
        
        Args:
            system_analysis_path: Path to system analysis JSON file
        """
        self.system_analysis_path = system_analysis_path
        self.system_analysis = None
        self.parameter_values = {}
        self.parameter_contributions = {}
        self.normalized_parameter_values = {}
        
        # Load system analysis
        self._load_system_analysis()
        
        log.info("Initialized GASLIT-AF Parameter Mapper")
    
    def _load_system_analysis(self):
        """
        Load system analysis from JSON file.
        """
        try:
            with open(self.system_analysis_path, 'r') as f:
                self.system_analysis = json.load(f)
            
            log.info(f"Loaded system analysis from {self.system_analysis_path}")
            log.info(f"Total variants: {self.system_analysis['total_variants']}")
            
            # Log system counts
            for system, count in self.system_analysis['system_counts'].items():
                log.info(f"  - {system}: {count} variants ({self.system_analysis['system_percentages'][system]:.2f}%)")
        
        except Exception as e:
            log.error(f"Error loading system analysis: {e}")
            sys.exit(1)
    
    def map_parameters(self):
        """
        Map system variant distributions to GASLIT-AF parameters.
        """
        log.info("Mapping system variant distributions to GASLIT-AF parameters")
        
        # Initialize parameter values and contributions
        for param in PARAMETERS.keys():
            self.parameter_values[param] = 0.0
            self.parameter_contributions[param] = {}
        
        # Calculate parameter values based on system variant distributions
        total_variants = self.system_analysis['total_variants']
        
        for system, count in self.system_analysis['system_counts'].items():
            # Skip systems not in weights
            if system not in SYSTEM_PARAMETER_WEIGHTS:
                log.warning(f"System not found in weights: {system}")
                continue
            
            # Calculate system weight based on variant count
            system_weight = count / total_variants
            
            # Calculate parameter contributions
            for param, weight in SYSTEM_PARAMETER_WEIGHTS[system].items():
                # Calculate contribution
                contribution = system_weight * weight
                
                # Update parameter value
                self.parameter_values[param] += contribution
                
                # Store contribution
                self.parameter_contributions[param][system] = contribution
        
        # Normalize parameter values to 0-10 scale
        self._normalize_parameters()
        
        # Log parameter values
        for param, value in self.normalized_parameter_values.items():
            log.info(f"{param} ({PARAMETERS[param]}): {value:.2f}")
    
    def _normalize_parameters(self):
        """
        Normalize parameter values to 0-10 scale.
        """
        # Calculate normalization factor
        # We want to scale the parameters so that the average is 5.0
        total = sum(self.parameter_values.values())
        avg = total / len(self.parameter_values)
        norm_factor = 5.0 / avg
        
        # Normalize parameters
        for param, value in self.parameter_values.items():
            normalized_value = value * norm_factor
            
            # Ensure value is between 0 and 10
            normalized_value = max(0.0, min(10.0, normalized_value))
            
            self.normalized_parameter_values[param] = normalized_value
    
    def generate_parameter_report(self, output_dir: str):
        """
        Generate parameter mapping report.
        
        Args:
            output_dir: Directory to save report
        """
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Save parameter values to JSON
        parameter_values_path = os.path.join(output_dir, 'parameter_values.json')
        with open(parameter_values_path, 'w') as f:
            json.dump({
                'parameter_values': self.parameter_values,
                'normalized_parameter_values': self.normalized_parameter_values,
                'parameter_contributions': self.parameter_contributions
            }, f, indent=2)
        
        log.info(f"Saved parameter values to {parameter_values_path}")
        
        # Generate parameter radar chart
        self._generate_parameter_radar_chart(output_dir)
        
        # Generate parameter contribution heatmap
        self._generate_parameter_contribution_heatmap(output_dir)
        
        # Generate parameter contribution bar chart
        self._generate_parameter_contribution_bar_chart(output_dir)
    
    def _generate_parameter_radar_chart(self, output_dir: str):
        """
        Generate parameter radar chart.
        
        Args:
            output_dir: Directory to save chart
        """
        try:
            # Set up radar chart
            params = list(PARAMETERS.keys())
            param_names = [f"{p} ({PARAMETERS[p]})" for p in params]
            values = [self.normalized_parameter_values[p] for p in params]
            
            # Create figure
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot(111, polar=True)
            
            # Set up radar chart
            angles = np.linspace(0, 2*np.pi, len(params), endpoint=False).tolist()
            values += values[:1]  # Close the loop
            angles += angles[:1]  # Close the loop
            
            # Plot radar chart
            ax.plot(angles, values, 'o-', linewidth=2)
            ax.fill(angles, values, alpha=0.25)
            
            # Set up chart labels
            ax.set_thetagrids(np.degrees(angles[:-1]), param_names)
            ax.set_ylim(0, 10)
            ax.set_yticks(np.arange(0, 11, 2))
            ax.set_yticklabels([str(x) for x in np.arange(0, 11, 2)])
            ax.grid(True)
            
            # Set title
            plt.title('GASLIT-AF Parameter Values', size=15, y=1.1)
            
            # Save chart
            radar_chart_path = os.path.join(output_dir, 'parameter_radar_chart.png')
            plt.tight_layout()
            plt.savefig(radar_chart_path, dpi=300)
            plt.close()
            
            log.info(f"Generated parameter radar chart: {radar_chart_path}")
        
        except Exception as e:
            log.error(f"Error generating parameter radar chart: {e}")
    
    def _generate_parameter_contribution_heatmap(self, output_dir: str):
        """
        Generate parameter contribution heatmap.
        
        Args:
            output_dir: Directory to save heatmap
        """
        try:
            # Create contribution matrix
            systems = list(SYSTEM_PARAMETER_WEIGHTS.keys())
            params = list(PARAMETERS.keys())
            
            # Initialize matrix
            matrix = np.zeros((len(systems), len(params)))
            
            # Fill matrix with contributions
            for i, system in enumerate(systems):
                for j, param in enumerate(params):
                    if system in self.parameter_contributions[param]:
                        matrix[i, j] = self.parameter_contributions[param][system]
            
            # Create figure
            plt.figure(figsize=(12, 10))
            
            # Create custom colormap
            colors = [(0.95, 0.95, 0.95), (0.8, 0.8, 1), (0, 0.5, 1)]
            cmap = LinearSegmentedColormap.from_list('custom_blue', colors, N=100)
            
            # Plot heatmap
            sns.heatmap(
                matrix,
                annot=True,
                fmt='.3f',
                cmap=cmap,
                xticklabels=[f"{p} ({PARAMETERS[p]})" for p in params],
                yticklabels=systems,
                cbar_kws={'label': 'Contribution'}
            )
            
            # Set title and labels
            plt.title('System Contributions to GASLIT-AF Parameters', size=15)
            plt.xlabel('Parameters')
            plt.ylabel('Biological Systems')
            
            # Adjust layout
            plt.tight_layout()
            
            # Save heatmap
            heatmap_path = os.path.join(output_dir, 'parameter_contribution_heatmap.png')
            plt.savefig(heatmap_path, dpi=300)
            plt.close()
            
            log.info(f"Generated parameter contribution heatmap: {heatmap_path}")
        
        except Exception as e:
            log.error(f"Error generating parameter contribution heatmap: {e}")
    
    def _generate_parameter_contribution_bar_chart(self, output_dir: str):
        """
        Generate parameter contribution bar chart.
        
        Args:
            output_dir: Directory to save chart
        """
        try:
            # Create figure
            fig, axes = plt.subplots(len(PARAMETERS), 1, figsize=(12, 15))
            
            # Plot bar chart for each parameter
            for i, (param, full_name) in enumerate(PARAMETERS.items()):
                # Get contributions
                contributions = self.parameter_contributions[param]
                
                # Sort contributions
                sorted_contributions = sorted(
                    contributions.items(),
                    key=lambda x: x[1],
                    reverse=True
                )
                
                # Extract systems and values
                systems = [s[0] for s in sorted_contributions]
                values = [s[1] for s in sorted_contributions]
                
                # Plot bar chart
                axes[i].barh(systems, values)
                
                # Set title and labels
                axes[i].set_title(f"{param} ({full_name})")
                axes[i].set_xlabel('Contribution')
                
                # Add value labels
                for j, v in enumerate(values):
                    axes[i].text(v + 0.001, j, f"{v:.3f}", va='center')
            
            # Adjust layout
            plt.tight_layout()
            
            # Save chart
            bar_chart_path = os.path.join(output_dir, 'parameter_contribution_bar_chart.png')
            plt.savefig(bar_chart_path, dpi=300)
            plt.close()
            
            log.info(f"Generated parameter contribution bar chart: {bar_chart_path}")
        
        except Exception as e:
            log.error(f"Error generating parameter contribution bar chart: {e}")

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="GASLIT-AF Parameter Mapping",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--system-analysis",
        default="analysis_results/direct_analysis/system_analysis.json",
        help="Path to system analysis JSON file"
    )
    
    parser.add_argument(
        "--output-dir",
        default="analysis_results/parameter_mapping",
        help="Directory to save parameter mapping results"
    )
    
    return parser.parse_args()

def main():
    """
    Main entry point for GASLIT-AF Parameter Mapping.
    """
    args = parse_args()
    
    # Create parameter mapper
    mapper = ParameterMapper(args.system_analysis)
    
    # Map parameters
    mapper.map_parameters()
    
    # Generate parameter report
    mapper.generate_parameter_report(args.output_dir)
    
    log.info(f"Parameter mapping complete. Results saved to {args.output_dir}")

if __name__ == "__main__":
    main()
