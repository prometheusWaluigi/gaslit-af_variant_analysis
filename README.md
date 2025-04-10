# GASLIT-AF Variant Analysis & Genomic Reconnaissance Pipeline

[![Poetry](https://img.shields.io/endpoint?url=https://python-poetry.org/badge/v0.json)](https://python-poetry.org/)
[![Python Version](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) <!-- Add appropriate license -->

**A high-performance pipeline for analyzing genomic variants (VCF) within the context of GASLIT-AF gene clusters, leveraging Intel oneAPI for acceleration. This project serves as the foundation for a multi-phase genomic reconnaissance effort aimed at decoding complex systemic dysregulation patterns.**

This project initially focuses on identifying genetic variants (primarily SNPs) within specific biological pathways relevant to the **Genetic Autonomic Structural Linked Instability Theorem – Allodynic Fatigue (GASLIT-AF)** model. However, its ultimate goal is to evolve into a comprehensive **GASLIT-AF Genomic Reconnaissance** platform, integrating multi-omic data and advanced analysis techniques to map the recursive attractor states underlying chronic conditions.

**Current Foundational Capabilities:**

*   Identify genetic variants within predefined GASLIT-AF gene lists from VCF files.
*   Integrate external databases like ClinVar for variant annotation.
*   Perform biological system-level analysis to understand pathway impacts.
*   Generate comprehensive reports and visualizations.
*   Utilize GPU acceleration via Intel oneAPI/SYCL (with CPU fallback) for efficient processing of large datasets.

## Project Vision: The GASLIT-AF Genomic Recon Protocol (v0.1)

While the current pipeline provides valuable insights into baseline genetic predispositions (γ), understanding the full GASLIT-AF dynamics requires a deeper, multi-dimensional approach. The long-term vision is guided by the **GASLIT-AF Genomic Recon Protocol**, encompassing:

**I. Deep Genomic Archeology (RCCX Region Forensics):** Detecting structural variants (CNVs, HERVs, chimeras) in the highly complex 6p21.3 RCCX locus (C4A/B, TNXB, CYP21A2).
**II. Epigenetic Tripwire Mapping:** Profiling methylation patterns in key feedback modulators (e.g., COMT, AHR/AHRR, NR3C1, IL6, TNF-α, FAAH, CNR1).
**III. Functional ECS Variant Sweep:** Characterizing endocannabinoid system (Ω) buffering capacity via genetic variants (FAAH, CNR1/2, MGLL, etc.).
**IV. Neuroimmune-IDO2-Kynurenine Nexus:** Investigating the IDO/TDO-Kynurenine-AHR pathway for metabolic trapping signatures (IDO1/2, KMO, KYNU, AHR, AHRR).
**V. Multi-Omic Loop Mapping:** Integrating transcriptomics, proteomics, and metabolomics to identify systemic attractor state signatures.
**VI. AI-Powered Variant Constellation Search:** Employing combinatorial analysis to find non-obvious SNP constellations driving recursive collapse.
**VII. Personalized GASLIT-AF Collapse Simulator:** Building a digital twin to model individual state-space trajectories Ψ(t) and simulate intervention effects.

**Current Status:** The pipeline currently implements the foundational layer, focusing on efficient variant detection within specified gene lists (related to steps III, IV, and baseline γ assessment). Future development will progressively integrate capabilities outlined in the Recon Protocol.

## Core Concepts (Foundational Pipeline)
*   **GASLIT-AF**: This pipeline is designed to explore the genetic underpinnings hypothesized by the GASLIT-AF framework, focusing on the comprehensive gene clusters (as defined in `src/gaslit_af/gene_lists.py`) associated with chronic multisystem conditions involving immune/inflammatory, autonomic/neurotransmitter, structural/connective tissue, metabolic, endocannabinoid, ion channel, mast cell, and kynurenine pathways. **Notably, the pipeline leverages a comprehensive gene list, enabling a broad analysis of genetic variants within the GASLIT-AF context.**
*   **VCF Analysis**: Processes standard Variant Call Format (VCF) files to extract relevant single nucleotide polymorphisms (SNPs) and other variants.
*   **Biological Systems**: Groups genes into functional pathways to analyze the systemic impact of identified variants (current implementation).

## Key Features (Foundational Pipeline)
*   **Intel oneAPI Acceleration**: Leverages `dpctl` and potentially other oneAPI libraries for SYCL-based GPU acceleration, **successfully tested and operational on Intel Arc A770 GPUs**, significantly speeding up variant processing on compatible Intel hardware. Includes automatic CPU fallback for broader compatibility.
*   **Memory-Optimized Chunking**: Processes large VCF files efficiently by reading and analyzing data in memory-bounded chunks.
*   **Parallel Processing**: Utilizes multi-threading for I/O and CPU-bound tasks.
*   **Modular Architecture**: Code is organized into logical modules (`workflow`, `data_processing`, `clinvar_integration`, `visualization`, etc.) for clarity and maintainability.
*   **ClinVar Integration**: Downloads, caches, and integrates ClinVar data for variant annotation (pathogenicity, associated conditions).
*   **Caching**: Caches intermediate results (like parsed ClinVar data, API responses) to speed up subsequent runs.
*   **Biological System Analysis**: Aggregates variant findings based on predefined biological systems/pathways.
*   **Comprehensive Reporting**: Generates static HTML reports and optional enhanced interactive reports (using Plotly) summarizing findings, variant counts, system distributions, and potentially clinical correlations.
*   **Direct VCF Analysis**: Provides specialized scripts for efficient analysis of pre-annotated VCF files, extracting gene information directly from the ANN field.
*   **Batch Processing**: Supports processing multiple VCF files in a single run, with automatic generation of combined reports for cross-sample comparison.
*   **Interactive Visualizations**: Creates interactive HTML visualizations using Plotly for system distributions, gene contributions, and variant patterns.
*   **Extensible**: Designed for future integration with modules addressing the broader GASLIT-AF Recon Protocol (e.g., CNV analysis, methylation data, multi-omic integration).

## Architecture Overview

### Core Library

The core logic resides in the `src/gaslit_af` directory:

*   `cli.py`: Handles command-line argument parsing.
*   `workflow.py`: Orchestrates the main analysis pipeline.
*   `device.py`: Manages device selection (GPU/CPU) for oneAPI.
*   `gene_lists.py` / `biological_systems.py`: Define the target genes and pathways for GASLIT-AF.
*   `data_processing.py`: Primary VCF parsing and variant extraction logic.
*   `streaming.py`: Provides memory-efficient streaming processing for large VCF files.
*   `caching.py`: Implements caching mechanisms.
*   `visualization.py`: Generates plots and figures (using Matplotlib, Seaborn, Plotly).
*   `reporting.py`: Creates HTML output reports.
*   `annovar_integration.py`: Optional integration with ANNOVAR for advanced variant annotation.
*   `rccx_analysis.py`: Specialized analysis of the RCCX region (6p21.3).
*   **(Future Modules):** Placeholder for methylation integration, advanced AI modeling, etc.

### Analysis Scripts

The repository includes several high-level scripts for different analysis workflows:

*   `direct_system_analysis.py`: Standalone script for analyzing pre-annotated VCF files, extracting gene information directly from the ANN field.
*   `run_direct_analysis.py`: Batch processor for running direct analysis on multiple VCF files in a directory.
*   `generate_combined_report.py`: Creates a comprehensive HTML report combining results from multiple analyses.
*   `run_full_analysis.py`: End-to-end pipeline that runs direct analysis on all VCF files and generates a combined report.

## Setup & Installation

This project uses [Poetry](https://python-poetry.org/) for dependency management.

1.  **Clone the repository:**
    ```bash
    git clone <repository-url>
    cd gaslit-af_variant_analysis
    ```
2.  **Ensure Poetry is installed:**
    ```bash
    pip install poetry
    ```
    *(Refer to official Poetry documentation for detailed installation instructions)*
3.  **Install dependencies:**
    ```bash
    poetry install
    ```
    *(This will create a virtual environment and install all required packages specified in `pyproject.toml`)*
4.  **Ensure Intel oneAPI Base Toolkit is installed** (for GPU acceleration): Follow Intel's installation guide for your operating system. Make sure the environment variables (`ONEAPI_ROOT`, etc.) are set correctly, or source the `setvars.sh` script.

## Usage

### Direct VCF Analysis Pipeline

The recommended workflow is to use the direct analysis pipeline, which is optimized for pre-annotated VCF files and provides a more streamlined experience.

#### Full Analysis Pipeline (Recommended)

The simplest way to analyze multiple VCF files is to use the `run_full_analysis.py` script:

```bash
python run_full_analysis.py --data-dir data --file-pattern "*.vcf.gz" --open-browser
```

This script:
1. Processes all VCF files in the specified directory
2. Performs biological system analysis on each file
3. Generates visualizations for system distributions and top genes
4. Creates a combined HTML report with interactive tabs for each sample

**Common Options:**
* `--data-dir <dir>`: Directory containing VCF files (default: "data")
* `--file-pattern <pattern>`: Pattern to match VCF files (default: "*.vcf.gz")
* `--output-dir <dir>`: Base directory for analysis results (default: "analysis_results")
* `--chunk-size <size>`: Number of variants to process at once (default: 10000000)
* `--open-browser`: Automatically open the report in a browser when complete

### Direct VCF Analysis

For more efficient analysis of annotated VCF files, the pipeline provides several specialized scripts:

#### 1. Direct Analysis of a Single VCF File

```bash
python direct_system_analysis.py <path/to/annotated.vcf.gz> --output-dir <output_directory>
```

This script extracts gene information directly from the ANN field in annotated VCF files and maps variants to GASLIT-AF biological systems.

#### 2. Batch Analysis of Multiple VCF Files

```bash
python run_direct_analysis.py --data-dir <directory_with_vcf_files> --file-pattern "*.vcf.gz"
```

This script processes all VCF files in the specified directory that match the given pattern, performing direct system analysis on each file.

#### 3. Generate Combined Report

```bash
python generate_combined_report.py --results-dir <analysis_results_directory>
```

This script generates a comprehensive HTML report combining the results from multiple direct analyses, including system distributions, top genes, and variant counts.

#### 4. Full Analysis Pipeline

```bash
python run_full_analysis.py --data-dir <directory_with_vcf_files> --open-browser
```

This script runs the complete pipeline:
1. Performs direct analysis on all VCF files in the specified directory
2. Generates a combined HTML report of all results
3. Optionally opens the report in a browser

The full analysis pipeline is particularly useful for comparing multiple samples or analyzing different types of variants (SNPs, CNVs, SVs) from the same sample.

## Outputs

### Analysis Outputs

### Direct Analysis Outputs

The direct analysis pipeline generates a more structured set of outputs:

*   **Per-VCF Analysis Directory**: For each analyzed VCF file, a directory is created (e.g., `analysis_results/direct_sample1.vcf`) containing:
    *   `system_analysis.json`: Detailed JSON with system counts, percentages, and gene mappings
    *   `gene_counts.csv`: CSV file with gene-level variant counts and system assignments
    *   `system_distribution.html/.png`: Interactive and static visualizations of variant distribution across biological systems
    *   `system_distribution_pie.html/.png`: Pie chart representation of system distribution
    *   `top_genes.html/.png`: Bar chart of the top genes by variant count

*   **Combined Report**: When using the full pipeline or `generate_combined_report.py`:
    *   `gaslit_af_combined_report.html`: A comprehensive HTML report that integrates results from all analyzed VCF files
    *   Interactive tabs for each sample
    *   Summary tables comparing variant counts and distributions across samples
    *   Embedded visualizations for each sample
    *   System-level comparisons and gene-level details

## Development & Testing

*   Tests are located in the `tests/` directory and use `pytest`.
*   Run tests using:
    ```bash
    poetry run pytest
    ```
    *(Or use the `run_tests.py` script)*
*   Code follows functional programming principles where practical and aims for modularity.
*   Contributions are welcome! Please follow standard fork-and-pull-request workflows, especially for developing modules aligned with the GASLIT-AF Recon Protocol.

## VCF File Requirements

### Standard Analysis Mode

For the standard analysis mode using `analyze_modular.py`, the pipeline accepts:

* Standard VCF files (version 4.0+)
* Gzipped VCF files (.vcf.gz)
* Both annotated and unannotated VCF files
* Whole genome, exome, or targeted panel VCF files

The pipeline will extract variants that match genes in the GASLIT-AF gene lists and perform annotation as needed.

### Direct Analysis Mode

For the direct analysis mode using `direct_system_analysis.py`, `run_direct_analysis.py`, or `run_full_analysis.py`, the pipeline requires:

* Pre-annotated VCF files with gene information in the ANN field
* Annotation should follow the standard format: `ANN=Gene|GeneID|...`
* Compatible with annotation tools like SnpEff, VEP, or ANNOVAR

The direct analysis mode is significantly faster as it extracts gene information directly from the annotation field without requiring additional lookups or API calls.

### Supported Variant Types

The pipeline currently focuses on:

* SNPs (Single Nucleotide Polymorphisms)
* Small indels (insertions/deletions)
* CNVs (Copy Number Variations) - basic support
* SVs (Structural Variants) - basic support

Future versions will expand support for more comprehensive analysis of structural variants, particularly in the RCCX region.

## Roadmap / Future Directions

1.  **Phase I Integration:** Develop modules for RCCX structural variant analysis (CNV detection specific to 6p21.3).
2.  **Phase II Integration:** Incorporate tools and workflows for analyzing targeted methylation data (WGBS/RRBS/Array).
3.  **Expand Gene Lists & Pathways:** Refine and expand gene sets based on ongoing research, particularly for ECS and Kynurenine pathways.
4.  **Multi-Omic Data Input:** Design interfaces for integrating transcriptomic, proteomic, and metabolomic datasets.
5.  **Advanced AI Models:** Implement combinatorial variant analysis and attractor state modeling (Phase VI/VII).
6.  **Refine oneAPI Usage:** Explore further optimization opportunities with SYCL/DPC++ for complex algorithms.
7.  **Enhanced VCF Support:** Improve handling of complex structural variants and multi-sample VCF files.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. <!-- Create a LICENSE file if one doesn't exist -->
