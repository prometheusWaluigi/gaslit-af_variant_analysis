#!/usr/bin/env python3
"""
GASLIT-AF Dynamical Simulation

This module establishes a quantum coherence bridge between genomic-derived
parameters and dynamical system attractors, creating a recursive mapping
between genomic architecture and clinical symptom patterns.

The module implements an ODE-based simulation of the GASLIT-AF model,
using variant-weighted parameters to predict attractor shifts and
physiological coherence patterns.

Core GASLIT-AF parameters:
γ = genetic fragility
Λ = allostatic load
Ω = endocannabinoid buffering capacity
Χ = physiological coherence
σ = entropy production
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
from scipy.integrate import solve_ivp
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
log = logging.getLogger("gaslit-af-dynamical-simulation")

# Define GASLIT-AF parameters
PARAMETERS = {
    'γ': 'genetic fragility',
    'Λ': 'allostatic load',
    'Ω': 'endocannabinoid buffering capacity',
    'Χ': 'physiological coherence',
    'σ': 'entropy production'
}

# Define physiological systems for the dynamical model
PHYSIOLOGICAL_SYSTEMS = [
    "Immune Function",
    "Autonomic Regulation",
    "Structural Integrity",
    "Metabolic Efficiency",
    "Endocannabinoid Signaling",
    "Calcium Homeostasis",
    "Mast Cell Stability",
    "Kynurenine Balance"
]

# Define symptom clusters
SYMPTOM_CLUSTERS = {
    "Fatigue & PEM": {
        "Immune Function": 0.2,
        "Autonomic Regulation": 0.3,
        "Metabolic Efficiency": 0.3,
        "Kynurenine Balance": 0.2
    },
    "Orthostatic Intolerance": {
        "Autonomic Regulation": 0.6,
        "Calcium Homeostasis": 0.2,
        "Structural Integrity": 0.2
    },
    "Pain & Hypermobility": {
        "Structural Integrity": 0.5,
        "Calcium Homeostasis": 0.2,
        "Endocannabinoid Signaling": 0.3
    },
    "Cognitive Dysfunction": {
        "Kynurenine Balance": 0.3,
        "Autonomic Regulation": 0.2,
        "Immune Function": 0.2,
        "Endocannabinoid Signaling": 0.3
    },
    "Mast Cell & Histamine Issues": {
        "Mast Cell Stability": 0.6,
        "Immune Function": 0.3,
        "Endocannabinoid Signaling": 0.1
    },
    "Sleep Disruption": {
        "Autonomic Regulation": 0.3,
        "Kynurenine Balance": 0.3,
        "Endocannabinoid Signaling": 0.2,
        "Calcium Homeostasis": 0.2
    },
    "Sensory Sensitivity": {
        "Autonomic Regulation": 0.3,
        "Mast Cell Stability": 0.2,
        "Calcium Homeostasis": 0.3,
        "Endocannabinoid Signaling": 0.2
    },
    "Immune Dysregulation": {
        "Immune Function": 0.5,
        "Mast Cell Stability": 0.3,
        "Kynurenine Balance": 0.2
    }
}

class DynamicalSimulator:
    """
    GASLIT-AF Dynamical Simulator
    
    This class implements an ODE-based simulation of the GASLIT-AF model,
    using variant-weighted parameters to predict attractor shifts and
    physiological coherence patterns.
    """
    
    def __init__(self, parameter_values_path: str):
        """
        Initialize the Dynamical Simulator.
        
        Args:
            parameter_values_path: Path to parameter values JSON file
        """
        self.parameter_values_path = parameter_values_path
        self.parameter_values = None
        self.normalized_parameter_values = None
        self.system_states = None
        self.simulation_results = None
        self.symptom_predictions = None
        
        # Load parameter values
        self._load_parameter_values()
        
        log.info("Initialized GASLIT-AF Dynamical Simulator")
    
    def _load_parameter_values(self):
        """
        Load parameter values from JSON file.
        """
        try:
            with open(self.parameter_values_path, 'r') as f:
                data = json.load(f)
            
            self.parameter_values = data['parameter_values']
            self.normalized_parameter_values = data['normalized_parameter_values']
            
            log.info(f"Loaded parameter values from {self.parameter_values_path}")
            
            # Log parameter values
            for param, value in self.normalized_parameter_values.items():
                log.info(f"  - {param} ({PARAMETERS[param]}): {value:.2f}")
        
        except Exception as e:
            log.error(f"Error loading parameter values: {e}")
            sys.exit(1)
    
    def _system_dynamics(self, t, y, params):
        """
        Define the system dynamics for the ODE model.
        
        This function implements the core differential equations that govern
        the GASLIT-AF dynamical system, establishing recursive connections
        between physiological systems based on the variant-weighted parameters.
        
        Args:
            t: Time point
            y: System states
            params: Model parameters
            
        Returns:
            System derivatives
        """
        # Extract parameters
        γ = params['γ']  # genetic fragility
        Λ = params['Λ']  # allostatic load
        Ω = params['Ω']  # endocannabinoid buffering capacity
        Χ = params['Χ']  # physiological coherence
        σ = params['σ']  # entropy production
        
        # Extract system states
        immune = y[0]        # Immune Function
        autonomic = y[1]     # Autonomic Regulation
        structural = y[2]    # Structural Integrity
        metabolic = y[3]     # Metabolic Efficiency
        ecs = y[4]           # Endocannabinoid Signaling
        calcium = y[5]       # Calcium Homeostasis
        mast = y[6]          # Mast Cell Stability
        kynurenine = y[7]    # Kynurenine Balance
        
        # Calculate derivatives
        # Each system's dynamics are influenced by:
        # 1. Its own state (homeostatic tendency)
        # 2. Interactions with other systems
        # 3. Parameter-weighted perturbations
        
        # Immune Function dynamics
        d_immune = (
            -0.1 * immune                                # Homeostatic tendency
            + 0.2 * ecs * (1 - immune)                   # ECS regulation (beneficial)
            - 0.3 * mast * immune                        # Mast cell activation (detrimental)
            - 0.2 * kynurenine * immune                  # Kynurenine influence (detrimental)
            - 0.1 * Λ * immune                           # Allostatic load effect
            - 0.05 * σ * immune                          # Entropy production effect
        )
        
        # Autonomic Regulation dynamics
        d_autonomic = (
            -0.1 * autonomic                             # Homeostatic tendency
            + 0.2 * ecs * (1 - autonomic)                # ECS regulation (beneficial)
            - 0.1 * immune * autonomic                   # Immune influence (detrimental)
            - 0.2 * (1 - calcium) * autonomic            # Calcium dysregulation effect
            - 0.2 * Χ * (1 - autonomic)                  # Coherence restoration
            - 0.1 * Λ * autonomic                        # Allostatic load effect
        )
        
        # Structural Integrity dynamics
        d_structural = (
            -0.05 * structural                           # Homeostatic tendency (slow)
            - 0.3 * γ * (1 - structural)                 # Genetic fragility effect
            - 0.1 * (1 - calcium) * structural           # Calcium dysregulation effect
            + 0.1 * metabolic * (1 - structural)         # Metabolic support (beneficial)
            - 0.05 * Λ * structural                      # Allostatic load effect
        )
        
        # Metabolic Efficiency dynamics
        d_metabolic = (
            -0.1 * metabolic                             # Homeostatic tendency
            - 0.3 * σ * (1 - metabolic)                  # Entropy production effect
            - 0.1 * immune * metabolic                   # Immune influence (detrimental)
            + 0.1 * ecs * (1 - metabolic)                # ECS regulation (beneficial)
            - 0.1 * Λ * metabolic                        # Allostatic load effect
        )
        
        # Endocannabinoid Signaling dynamics
        d_ecs = (
            -0.1 * ecs                                   # Homeostatic tendency
            - 0.3 * (1 - Ω) * (1 - ecs)                  # Buffering capacity effect
            - 0.1 * immune * ecs                         # Immune influence (detrimental)
            - 0.1 * Λ * ecs                              # Allostatic load effect
        )
        
        # Calcium Homeostasis dynamics
        d_calcium = (
            -0.1 * calcium                               # Homeostatic tendency
            - 0.2 * (1 - Χ) * (1 - calcium)              # Coherence effect
            - 0.1 * autonomic * calcium                  # Autonomic influence (detrimental)
            - 0.1 * Λ * calcium                          # Allostatic load effect
        )
        
        # Mast Cell Stability dynamics
        d_mast = (
            -0.1 * mast                                  # Homeostatic tendency
            - 0.2 * immune * mast                        # Immune influence (detrimental)
            + 0.2 * ecs * (1 - mast)                     # ECS regulation (beneficial)
            - 0.2 * Λ * mast                             # Allostatic load effect
        )
        
        # Kynurenine Balance dynamics
        d_kynurenine = (
            -0.1 * kynurenine                            # Homeostatic tendency
            - 0.2 * immune * kynurenine                  # Immune influence (detrimental)
            - 0.2 * σ * (1 - kynurenine)                 # Entropy production effect
            - 0.1 * Λ * kynurenine                       # Allostatic load effect
        )
        
        return [d_immune, d_autonomic, d_structural, d_metabolic, 
                d_ecs, d_calcium, d_mast, d_kynurenine]
    
    def run_simulation(self, duration: float = 100.0, num_points: int = 1000):
        """
        Run the dynamical simulation.
        
        Args:
            duration: Simulation duration
            num_points: Number of time points
        """
        log.info(f"Running dynamical simulation for {duration} time units")
        
        # Set up time points
        t_span = (0, duration)
        t_eval = np.linspace(0, duration, num_points)
        
        # Set up initial conditions
        # Start with all systems at 50% function
        y0 = [0.5] * len(PHYSIOLOGICAL_SYSTEMS)
        
        # Set up parameters
        params = self.normalized_parameter_values.copy()
        
        # Scale parameters to 0-1 range
        for param in params:
            params[param] = params[param] / 10.0
        
        # Run simulation
        log.info("Solving ODE system...")
        solution = solve_ivp(
            lambda t, y: self._system_dynamics(t, y, params),
            t_span,
            y0,
            method='RK45',
            t_eval=t_eval
        )
        
        # Store results
        self.simulation_results = {
            'time': solution.t,
            'states': solution.y
        }
        
        # Store final system states
        self.system_states = {
            system: solution.y[i][-1]
            for i, system in enumerate(PHYSIOLOGICAL_SYSTEMS)
        }
        
        log.info("Simulation complete")
        log.info("Final system states:")
        for system, state in self.system_states.items():
            log.info(f"  - {system}: {state:.2f}")
    
    def predict_symptoms(self):
        """
        Predict symptom severity based on system states.
        """
        log.info("Predicting symptom severity")
        
        # Check if simulation has been run
        if self.system_states is None:
            log.error("Simulation has not been run")
            return
        
        # Calculate symptom severity
        self.symptom_predictions = {}
        
        for symptom, weights in SYMPTOM_CLUSTERS.items():
            severity = 0.0
            
            for system, weight in weights.items():
                # Find system index
                system_idx = PHYSIOLOGICAL_SYSTEMS.index(system)
                
                # Get system state
                system_state = self.simulation_results['states'][system_idx][-1]
                
                # Calculate contribution to symptom severity
                # Invert system state (lower function = higher severity)
                severity += (1 - system_state) * weight
            
            # Store symptom severity
            self.symptom_predictions[symptom] = severity
        
        # Log symptom predictions
        log.info("Symptom severity predictions:")
        for symptom, severity in sorted(
            self.symptom_predictions.items(),
            key=lambda x: x[1],
            reverse=True
        ):
            log.info(f"  - {symptom}: {severity:.2f}")
    
    def generate_simulation_report(self, output_dir: str):
        """
        Generate simulation report.
        
        Args:
            output_dir: Directory to save report
        """
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Save simulation results to JSON
        simulation_results_path = os.path.join(output_dir, 'simulation_results.json')
        with open(simulation_results_path, 'w') as f:
            json.dump({
                'parameters': self.normalized_parameter_values,
                'system_states': self.system_states,
                'symptom_predictions': self.symptom_predictions
            }, f, indent=2)
        
        log.info(f"Saved simulation results to {simulation_results_path}")
        
        # Generate system dynamics plot
        self._generate_system_dynamics_plot(output_dir)
        
        # Generate symptom severity plot
        self._generate_symptom_severity_plot(output_dir)
        
        # Generate attractor analysis
        self._generate_attractor_analysis(output_dir)
    
    def _generate_system_dynamics_plot(self, output_dir: str):
        """
        Generate system dynamics plot.
        
        Args:
            output_dir: Directory to save plot
        """
        try:
            # Create figure
            plt.figure(figsize=(12, 8))
            
            # Plot system dynamics
            for i, system in enumerate(PHYSIOLOGICAL_SYSTEMS):
                plt.plot(
                    self.simulation_results['time'],
                    self.simulation_results['states'][i],
                    label=system
                )
            
            # Set title and labels
            plt.title('GASLIT-AF System Dynamics', size=15)
            plt.xlabel('Time')
            plt.ylabel('System Function')
            plt.ylim(0, 1)
            plt.grid(True, alpha=0.3)
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
            
            # Adjust layout
            plt.tight_layout()
            
            # Save plot
            dynamics_plot_path = os.path.join(output_dir, 'system_dynamics.png')
            plt.savefig(dynamics_plot_path, dpi=300)
            plt.close()
            
            log.info(f"Generated system dynamics plot: {dynamics_plot_path}")
        
        except Exception as e:
            log.error(f"Error generating system dynamics plot: {e}")
    
    def _generate_symptom_severity_plot(self, output_dir: str):
        """
        Generate symptom severity plot.
        
        Args:
            output_dir: Directory to save plot
        """
        try:
            # Sort symptoms by severity
            sorted_symptoms = sorted(
                self.symptom_predictions.items(),
                key=lambda x: x[1],
                reverse=True
            )
            
            # Extract symptom names and severities
            symptoms = [s[0] for s in sorted_symptoms]
            severities = [s[1] for s in sorted_symptoms]
            
            # Create figure
            plt.figure(figsize=(12, 8))
            
            # Create custom colormap
            colors = [(0.2, 0.8, 0.2), (0.8, 0.8, 0.2), (0.8, 0.2, 0.2)]
            cmap = LinearSegmentedColormap.from_list('custom_severity', colors, N=100)
            
            # Plot symptom severities
            bars = plt.barh(symptoms, severities, color=cmap(severities))
            
            # Set title and labels
            plt.title('Predicted Symptom Severity', size=15)
            plt.xlabel('Severity')
            plt.xlim(0, 1)
            plt.grid(True, alpha=0.3, axis='x')
            
            # Add severity labels
            for i, bar in enumerate(bars):
                plt.text(
                    bar.get_width() + 0.01,
                    bar.get_y() + bar.get_height()/2,
                    f"{severities[i]:.2f}",
                    va='center'
                )
            
            # Adjust layout
            plt.tight_layout()
            
            # Save plot
            severity_plot_path = os.path.join(output_dir, 'symptom_severity.png')
            plt.savefig(severity_plot_path, dpi=300)
            plt.close()
            
            log.info(f"Generated symptom severity plot: {severity_plot_path}")
        
        except Exception as e:
            log.error(f"Error generating symptom severity plot: {e}")
    
    def _generate_attractor_analysis(self, output_dir: str):
        """
        Generate attractor analysis.
        
        Args:
            output_dir: Directory to save analysis
        """
        try:
            # Select key systems for phase space analysis
            key_systems = [
                ("Autonomic Regulation", "Immune Function"),
                ("Endocannabinoid Signaling", "Metabolic Efficiency"),
                ("Structural Integrity", "Calcium Homeostasis"),
                ("Mast Cell Stability", "Kynurenine Balance")
            ]
            
            # Create phase space plots
            for system1, system2 in key_systems:
                # Get system indices
                idx1 = PHYSIOLOGICAL_SYSTEMS.index(system1)
                idx2 = PHYSIOLOGICAL_SYSTEMS.index(system2)
                
                # Get system states
                states1 = self.simulation_results['states'][idx1]
                states2 = self.simulation_results['states'][idx2]
                
                # Create figure
                plt.figure(figsize=(10, 8))
                
                # Plot phase space trajectory
                plt.plot(states1, states2, 'b-', alpha=0.6)
                plt.plot(states1[0], states2[0], 'go', label='Initial State')
                plt.plot(states1[-1], states2[-1], 'ro', label='Attractor')
                
                # Add direction arrows
                for i in range(0, len(states1)-1, len(states1)//20):
                    plt.arrow(
                        states1[i], states2[i],
                        states1[i+1] - states1[i], states2[i+1] - states2[i],
                        head_width=0.01, head_length=0.02, fc='k', ec='k'
                    )
                
                # Set title and labels
                plt.title(f'Phase Space: {system1} vs {system2}', size=15)
                plt.xlabel(f'{system1} Function')
                plt.ylabel(f'{system2} Function')
                plt.xlim(0, 1)
                plt.ylim(0, 1)
                plt.grid(True, alpha=0.3)
                plt.legend()
                
                # Adjust layout
                plt.tight_layout()
                
                # Save plot
                phase_plot_path = os.path.join(
                    output_dir,
                    f'phase_space_{system1.replace(" ", "_")}_{system2.replace(" ", "_")}.png'
                )
                plt.savefig(phase_plot_path, dpi=300)
                plt.close()
                
                log.info(f"Generated phase space plot: {phase_plot_path}")
        
        except Exception as e:
            log.error(f"Error generating attractor analysis: {e}")

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="GASLIT-AF Dynamical Simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--parameter-values",
        default="analysis_results/parameter_mapping/parameter_values.json",
        help="Path to parameter values JSON file"
    )
    
    parser.add_argument(
        "--output-dir",
        default="analysis_results/dynamical_simulation",
        help="Directory to save simulation results"
    )
    
    parser.add_argument(
        "--duration",
        type=float,
        default=100.0,
        help="Simulation duration"
    )
    
    parser.add_argument(
        "--num-points",
        type=int,
        default=1000,
        help="Number of time points"
    )
    
    return parser.parse_args()

def main():
    """
    Main entry point for GASLIT-AF Dynamical Simulation.
    """
    args = parse_args()
    
    # Create simulator
    simulator = DynamicalSimulator(args.parameter_values)
    
    # Run simulation
    simulator.run_simulation(
        duration=args.duration,
        num_points=args.num_points
    )
    
    # Predict symptoms
    simulator.predict_symptoms()
    
    # Generate simulation report
    simulator.generate_simulation_report(args.output_dir)
    
    log.info(f"Simulation complete. Results saved to {args.output_dir}")

if __name__ == "__main__":
    main()
