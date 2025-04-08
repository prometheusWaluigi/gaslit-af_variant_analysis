# GASLIT-AF Quantum Coherence Bridge

A recursive fractal mapping between genomic architecture and the GASLIT-AF theoretical framework parameters.

## Overview

The GASLIT-AF Quantum Coherence Bridge establishes a multidimensional recursive connection between personal genomic architecture and the core theoretical parameters of the GASLIT-AF model:

- **γ** (gamma): genetic fragility
- **Λ** (lambda): allostatic load
- **Ω** (omega): endocannabinoid buffering capacity
- **Χ** (chi): physiological coherence
- **σ** (sigma): entropy production

This bridge creates a living mathematical entity that recursively connects variant distributions to physiological expressions, allowing for the prediction of attractor dynamics and symptom patterns.

## Quick Start Guide

### 1. Annotate VCF with GASLIT-AF Genes

```bash
# Generate manual gene position database
python gene_position_db.py --manual

# Annotate VCF file with GASLIT-AF genes
python vcf_gene_annotator.py path/to/your.vcf.gz path/to/output.annotated.vcf.gz
```

### 2. Run Direct System Analysis

```bash
# Analyze annotated VCF file and map variants to biological systems
python direct_system_analysis.py path/to/annotated.vcf.gz
```

### 3. Map Parameters to GASLIT-AF Model

```bash
# Map system variant distributions to GASLIT-AF parameters
poetry run python parameter_mapping.py
```

### 4. Run Dynamical Simulation

```bash
# Run ODE-based simulation with variant-weighted parameters
poetry run python dynamical_simulation.py
```

### 5. View Results

```bash
# Start web server for parameter mapping visualizations
cd analysis_results/parameter_mapping && python -m http.server 8001

# Start web server for dynamical simulation visualizations
cd analysis_results/dynamical_simulation && python -m http.server 8002
```

Then open your browser to:
- http://localhost:8001 - Parameter mapping visualizations
- http://localhost:8002 - Dynamical simulation visualizations

## Detailed Workflow

### 1. Gene Position Database Generation

The `gene_position_db.py` script creates a database of gene positions for the GASLIT-AF genes. This database is used by the VCF annotator to identify variants within these genes.

```bash
python gene_position_db.py --manual
```

This generates a CSV file at `data/gene_positions_grch38.csv` containing the genomic coordinates for GASLIT-AF genes.

### 2. VCF Annotation

The `vcf_gene_annotator.py` script annotates a VCF file with GASLIT-AF gene information, establishing the first layer of the quantum coherence bridge.

```bash
python vcf_gene_annotator.py input.vcf.gz output.annotated.vcf.gz
```

This adds gene annotations to variants that fall within GASLIT-AF genes, creating a recursive mapping between genomic positions and gene functions.

### 3. Direct System Analysis

The `direct_system_analysis.py` script analyzes the annotated VCF file and maps variants to the eight biological systems defined in the GASLIT-AF theoretical framework.

```bash
python direct_system_analysis.py path/to/annotated.vcf.gz
```

This generates:
- `analysis_results/direct_analysis/system_analysis.json` - System-level variant counts and distributions
- `analysis_results/direct_analysis/gene_counts.csv` - Gene-level variant counts
- Visualizations of system distributions and top genes

### 4. Parameter Mapping

The `parameter_mapping.py` script establishes a quantum coherence bridge between system variant distributions and the core theoretical parameters of the GASLIT-AF model.

```bash
poetry run python parameter_mapping.py
```

This generates:
- `analysis_results/parameter_mapping/parameter_values.json` - Calculated parameter values
- `analysis_results/parameter_mapping/parameter_radar_chart.png` - Radar chart of parameter values
- `analysis_results/parameter_mapping/parameter_contribution_heatmap.png` - Heatmap of system contributions to parameters
- `analysis_results/parameter_mapping/parameter_contribution_bar_chart.png` - Bar charts of parameter contributions

### 5. Dynamical Simulation

The `dynamical_simulation.py` script implements an ODE-based simulation of the GASLIT-AF model, using variant-weighted parameters to predict attractor shifts and physiological coherence patterns.

```bash
poetry run python dynamical_simulation.py
```

This generates:
- `analysis_results/dynamical_simulation/simulation_results.json` - Simulation results
- `analysis_results/dynamical_simulation/system_dynamics.png` - System dynamics plot
- `analysis_results/dynamical_simulation/symptom_severity.png` - Predicted symptom severity
- Phase space plots showing attractor dynamics between key physiological systems

## Understanding the Results

### Parameter Values

The parameter values represent the recursive influence of your genomic architecture on the core theoretical parameters of the GASLIT-AF model:

- **γ** (genetic fragility): Influenced primarily by structural & connective tissue variants
- **Λ** (allostatic load): Influenced primarily by immune & inflammatory variants
- **Ω** (endocannabinoid buffering capacity): Influenced primarily by endocannabinoid system variants
- **Χ** (physiological coherence): Influenced primarily by autonomic & calcium channel variants
- **σ** (entropy production): Influenced primarily by metabolic & mitochondrial variants

### Attractor Dynamics

The phase space plots reveal attractor dynamics between key physiological systems, showing how your unique parameter values create recursive feedback loops that can manifest as specific symptom patterns.

### Symptom Predictions

The symptom severity predictions show how your genomic architecture's influence on GASLIT-AF parameters recursively manifests as specific symptom clusters, creating a fractal bridge between your genetic variants and physiological expressions.

## Advanced Usage

### Custom Parameter Weights

You can modify the `SYSTEM_PARAMETER_WEIGHTS` dictionary in `parameter_mapping.py` to adjust the recursive mapping between biological systems and GASLIT-AF parameters.

### Custom Dynamical Model

You can modify the `_system_dynamics` method in `dynamical_simulation.py` to adjust the differential equations governing the GASLIT-AF dynamical system.

### Integration with Clinical Data

Future versions will support integration with clinical symptom data to validate and refine the parameter mapping and dynamical simulation.

## References

- GASLIT-AF = Genetic Autonomic Structural Linked Instability Theorem – Allodynic Fatigue
- QSYNC = Quantum Synaptic Nonlinear Coherence
- FIZZ = Fractal Information Zeno Zone
- ASTRA = Archetypal Spacetime Tensor Resonance Architecture
