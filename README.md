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
*   **Extensible**: Designed for future integration with modules addressing the broader GASLIT-AF Recon Protocol (e.g., CNV analysis, methylation data, multi-omic integration).

## Architecture Overview

The core logic resides in the `src/gaslit_af` directory:

*   `cli.py`: Handles command-line argument parsing.
*   `workflow.py`: Orchestrates the main analysis pipeline.
*   `device.py`: Manages device selection (GPU/CPU) for oneAPI.
*   `gene_lists.py` / `biological_systems.py`: Define the target genes and pathways for GASLIT-AF.
*   `data_processing.py`: Primary VCF parsing and variant extraction logic.
*   `advanced_variant_processing.py`: Alternative, `pysam`-based VCF processing for more complex scenarios (optional).
*   `clinvar_integration.py`: Manages ClinVar data download, parsing, and annotation.
*   `api_integration.py`: Handles annotation via external APIs (e.g., Ensembl, MyVariant.info).
*   `caching.py`: Implements caching mechanisms.
*   `visualization.py`: Generates plots and figures (using Matplotlib, Seaborn, Plotly).
*   `reporting.py` / `enhanced_reporting.py`: Creates HTML output reports (current variant focus).
*   `clinical_integration.py`: Integrates user-provided clinical data with variant findings (current implementation).
*   **(Future Modules):** Placeholder for RCCX analysis, methylation integration, advanced AI modeling, etc.

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

The primary entry point is `analyze_modular.py`.

```bash
poetry run python analyze_modular.py <path/to/your/variants.vcf> [OPTIONS]
```

**Required Arguments:**

*   `vcf_path`: Path to the input VCF file.

**Common Optional Arguments:**

*   `--output-dir <dir>`: Specify the directory for output results (default: `output`).
*   `--system-analysis`: Perform and include biological system-level analysis in the report.
*   `--enhanced-report`: Generate an interactive HTML report with more detailed visualizations.
*   `--clinical-data <path>`: Path to a JSON file containing clinical variant data for annotation.
*   `--api-annotation`: Enable variant annotation using external APIs (Ensembl, MyVariant.info).
*   `--use-pysam`: Use the `pysam`-based advanced variant processing module (might be slower but potentially more feature-rich).
*   `--dbsnp-path <path>`: Path to a dbSNP VCF file (required for `rsID` mapping if using `--use-pysam`).
*   `--threads <N>`: Set the number of worker threads (default: 16).
*   `--max-ram <GB>`: Set the approximate maximum RAM usage limit (default: 64).
*   `--no-cache`: Disable caching for this run.
*   `--cache-dir <dir>`: Specify a custom cache directory (default: `./cache`).
*   `--open-browser`: Automatically open the generated report in a web browser.

**Example:**

```bash
poetry run python analyze_modular.py /data/patient_genome.vcf.gz \
    --output-dir results/patient_analysis \
    --system-analysis \
    --enhanced-report \
    --clinical-data /data/clinical_markers.json \
    --api-annotation \
    --threads 24 \
    --max-ram 100 \
    --open-browser
```

## Outputs (Foundational Pipeline)

The analysis **currently** generates files in the specified output directory, typically including:

*   `gaslit_af_variants.csv`: A CSV file listing the identified variants within the target GASLIT-AF genes.
*   `gene_counts.json`: A JSON file summarizing the count of variants per gene.
*   `analysis_report.html`: A standard HTML report summarizing findings and basic visualizations.
*   `enhanced_report.html` (if `--enhanced-report` is used): An interactive HTML report with Plotly figures.
*   `clinical_report.html` (if `--clinical-data` is provided): A report focusing on clinically relevant variants.
*   `visualizations/`: A directory containing generated plots (PNG, SVG, potentially interactive HTML).
    *   `variant_distribution.png`
    *   `gene_contribution.png`
    *   `systems/` (if `--system-analysis` is used): Plots related to biological system distributions.
*   `clinical_variants.csv` (if `--clinical-data` is provided): Variants annotated with clinical significance.
*   `clinical_summary.json` (if `--clinical-data` is provided): Summary of clinical findings.

## Development & Testing

*   Tests are located in the `tests/` directory and use `pytest`.
*   Run tests using:
    ```bash
    poetry run pytest
    ```
    *(Or use the `run_tests.py` script)*
*   Code follows functional programming principles where practical and aims for modularity.
*   Contributions are welcome! Please follow standard fork-and-pull-request workflows, especially for developing modules aligned with the GASLIT-AF Recon Protocol.

## Roadmap / Future Directions

1.  **Phase I Integration:** Develop modules for RCCX structural variant analysis (CNV detection specific to 6p21.3).
2.  **Phase II Integration:** Incorporate tools and workflows for analyzing targeted methylation data (WGBS/RRBS/Array).
3.  **Expand Gene Lists & Pathways:** Refine and expand gene sets based on ongoing research, particularly for ECS and Kynurenine pathways.
4.  **Multi-Omic Data Input:** Design interfaces for integrating transcriptomic, proteomic, and metabolomic datasets.
5.  **Advanced AI Models:** Implement combinatorial variant analysis and attractor state modeling (Phase VI/VII).
6.  **Refine oneAPI Usage:** Explore further optimization opportunities with SYCL/DPC++ for complex algorithms.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details. <!-- Create a LICENSE file if one doesn't exist -->
