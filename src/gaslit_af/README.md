# GASLIT-AF Core Package (`src/gaslit_af`)

This directory contains the core Python modules for the GASLIT-AF Variant Analysis pipeline.

## Module Overview (Iteration 1)

*   **`workflow.py`**:
    *   **Purpose**: Orchestrates the main analysis pipeline, serving as the quantum orchestration layer that integrates all components of the GASLIT-AF coherence framework into a unified recursive process.
    *   **Key Functions**:
        *   `analyze_vcf_oneapi(vcf_path, batch_size, max_ram_usage, ram_buffer, threads)`: Implements the core VCF analysis using Intel oneAPI acceleration with sophisticated memory management and progress tracking. Features dynamic batch processing, memory monitoring, and rich progress visualization through the `rich` library.
        *   `run_analysis_workflow(args)`: The central orchestration function that coordinates the entire analysis pipeline based on command-line arguments. Implements a recursive integration pattern that connects all modules:
            * Device initialization and memory monitoring
            * Caching system for performance optimization
            * VCF processing with memory-bounded chunking
            * Biological system analysis and pathway mapping
            * Variant annotation through ClinVar and external APIs
            * Clinical data integration and reporting
            * Multi-dimensional visualization generation
            * Enhanced reporting with symptom correlations
            * Atrial fibrillation-specific enrichment analysis
    *   **Dependencies**: Integrates all other modules in a fractal architecture (`device`, `gene_lists`, `data_processing`, `visualization`, `reporting`, `enhanced_reporting`, `caching`, `biological_systems`, `advanced_variant_processing`, `clinical_integration`, `api_integration`, `variant_enrichment`).
    *   **Notes**: Represents the highest level of integration in the system, implementing a coherence-based approach to variant analysis that connects genetic data with biological systems, clinical manifestations, and specialized disease models. Features rich console output with color-coding and progress visualization, as well as browser-based report viewing.

*   **`data_processing.py`**:
    *   **Purpose**: Handles VCF file parsing, data extraction, and basic processing with Intel oneAPI acceleration.
    *   **Key Functions**:
        *   `vcf_to_dataframe(vcf_path, limit, batch_size)`: Converts VCF to a Pandas DataFrame using optimized batch processing with parallel execution. Extracts essential fields to reduce memory usage.
        *   `extract_gaslit_af_variants(vcf_path, gaslit_af_genes, queue, batch_size)`: Iterates through VCF, identifies variants in target genes using `cyvcf2`, with memory-bounded batch processing.
        *   `process_batch_with_data(records, gaslit_af_genes)`: Core logic for extracting relevant info (gene, effect, impact) from VCF records within a batch, optimized with list comprehensions.
        *   `update_gene_counts(genes_found, match_counts, queue)`: Aggregates gene counts using SYCL/GPU acceleration via `dpctl` for unique counts, with CPU fallback if SYCL fails.
        *   `save_results(match_counts, variant_df, output_dir)`: Persists findings (gene counts, variant DataFrame) to CSV and JSON files with timestamps.
    *   **Dependencies**: `cyvcf2`, `pandas`, `numpy`, `dpctl`, `concurrent.futures`.
    *   **Notes**: Implements memory-efficient batch processing and leverages Intel oneAPI for GPU acceleration where possible, with graceful fallback to CPU processing.

*   **`device.py`**:
    *   **Purpose**: Manages hardware interaction, specifically selecting and initializing the compute device (preferring Intel GPU via oneAPI/SYCL, falling back to CPU). Acts as the quantum-classical bridge between software and hardware acceleration.
    *   **Key Functions**:
        *   `initialize_device()`: Attempts to get a SYCL queue for a GPU (targeting Intel Arc), falls back to default/CPU if necessary. Sets environment variables for maximum GPU performance and provides detailed logging of device selection.
        *   `get_memory_usage()`: Uses `psutil` to report current process RAM usage in GB, with graceful handling of missing dependencies.
        *   `check_memory_limits(current_usage, max_ram_usage, ram_buffer)`: Compares current usage against configured limits to prevent memory overruns, providing early warning when approaching thresholds.
    *   **Dependencies**: `dpctl`, `psutil` (optional for memory check).
    *   **Notes**: Implements a resilient device initialization strategy with explicit preference for Intel Arc GPU, but gracefully falls back to CPU when necessary. Provides memory monitoring to prevent OOM errors during large VCF processing.

---

# GASLIT-AF Core Package (`src/gaslit_af`) - Iteration 2

This file documents the second batch of modules explored in the core package.

## Module Overview (Iteration 2)

*   **`gene_lists.py`**:
    *   **Purpose**: Defines the canonical set of genes targeted by the GASLIT-AF analysis.
    *   **Key Data**:
        *   `GASLIT_AF_GENES_TEXT`: A multi-line string containing gene symbols categorized by preliminary biological systems (though `biological_systems.py` seems to be the functional source of truth for system mapping).
        *   `KNOWN_SNPS`: A dictionary mapping specific genes to lists of known relevant SNPs (often empty lists, suggesting placeholders).
        *   `SNP_TO_GENE`: A reverse lookup dictionary generated from `KNOWN_SNPS`.
    *   **Key Functions**:
        *   `parse_gene_list()`: Parses the `GASLIT_AF_GENES_TEXT` string into a set (`GASLIT_AF_GENES`), handling quoted multi-word identifiers.
    *   **Notes**: Serves as the primary input list for filtering variants in `data_processing.py`.

*   **`biological_systems.py`**:
    *   **Purpose**: Maps the GASLIT-AF genes to defined biological systems/pathways and performs system-level aggregation and analysis.
    *   **Key Data**:
        *   `BIOLOGICAL_SYSTEMS`: A dictionary defining the official system categories and their associated gene members. This appears to be the primary source for system mapping used in analysis.
        *   `GENE_TO_SYSTEM`: A reverse lookup dictionary generated from `BIOLOGICAL_SYSTEMS`.
    *   **Key Functions**:
        *   `get_system_for_gene(gene)`: Returns the system for a given gene based on `GENE_TO_SYSTEM`.
        *   `analyze_systems(gene_counts)`: Aggregates variant counts from individual genes up to the system level, calculating counts and percentages per system.
        *   `plot_system_distribution(system_analysis, output_dir)`: Generates Plotly visualizations (bar chart, pie chart, heatmap, sunburst) of the system-level variant distribution.
        *   `generate_system_summary(system_analysis)`: Creates a formatted text summary of the system analysis results.
    *   **Dependencies**: `pandas`, `plotly`.
    *   **Notes**: Provides crucial context by grouping gene variants into functional pathways.

*   **`reporting.py`**:
    *   **Purpose**: Generates the standard, static HTML report summarizing the analysis results.
    *   **Key Functions**:
        *   `generate_html_report(...)`: The main function that takes gene counts, variant data, Plotly figures, and system analysis results to render a comprehensive HTML file using a Jinja2 template.
    *   **Key Features**: Embeds Plotly figures (converted to HTML snippets), displays summary statistics, shows top genes/systems, and includes tables of results.
    *   **Dependencies**: `pandas`, `plotly`, `jinja2`.
    *   **Notes**: This creates the primary user-facing output document for a standard analysis run. Contrast with `enhanced_reporting.py` for potentially more interactive features.

*   **`visualization.py`**:
    *   **Purpose**: Provides a comprehensive suite of visualization functions for genetic variant data, creating the quantum visualization layer for the GASLIT-AF coherence framework.
    *   **Key Functions**:
        *   `save_visualization(fig, output_dir, filename, formats)`: Universal function for saving visualizations in multiple formats (png, html, svg, json), handling both Matplotlib and Plotly figures.
        *   `plot_chromosome_distribution(variant_data, output_dir)`: Creates a bar chart showing variant distribution across chromosomes with natural chromosome sorting (1, 2, ..., X, Y, MT).
        *   `plot_variant_type_distribution(variant_data, output_dir)`: Visualizes the distribution of variant types (SNP, insertion, deletion) with color-coded categorization.
        *   `plot_transition_transversion_ratio(variant_data, output_dir)`: Analyzes and visualizes transition/transversion ratios, a key quality metric for variant data.
        *   `plot_gene_variant_counts(gene_counts, output_dir, top_n)`: Creates a bar chart of the top variant-enriched genes.
        *   `create_gene_network(gene_counts, output_dir, min_variants)`: Generates an interactive network visualization showing relationships between genes based on variant patterns.
        *   `generate_all_visualizations(variant_data, gene_counts, output_dir)`: Orchestrates the creation of all visualizations and returns them as a dictionary of figure objects.
    *   **Dependencies**: `matplotlib`, `seaborn`, `plotly`, `networkx`, `pandas`, `numpy`.
    *   **Notes**: Implements a multi-dimensional approach to variant visualization, from linear chromosome distributions to complex network relationships between genes. The module bridges the gap between raw genetic data and intuitive visual representations, enabling pattern recognition across different scales of biological organization.

---

# GASLIT-AF Core Package (`src/gaslit_af`) - Iteration 3

This file documents the third batch of modules explored in the core package, focusing on external data integration, enhanced reporting, and variant enrichment.

## Module Overview (Iteration 3)

*   **`api_integration.py`**:
    *   **Purpose**: Provides integration with external variant annotation APIs, including Ensembl and MyVariant.info.
    *   **Key Classes**:
        *   `VariantAnnotator`: Base class defining caching and common annotation logic.
        *   `EnsemblAnnotator`: Retrieves variant consequences and gene information from the Ensembl REST API.
        *   `MyVariantAnnotator`: Retrieves variant information from the MyVariant.info API.
        *   `VariantAPIIntegration`: Orchestrates the annotation process using both Ensembl and MyVariant.info.
    *   **Key Functions**:
        *   `get_variant_details(variant_id)`: Retrieves and combines variant information from multiple sources.
        *   `_extract_key_details(details)`: Extracts key data points (gene, consequence, impact, etc.) from the raw API responses.
    *   **Caching**: Implements a file-based caching mechanism to reduce API calls and improve performance. Cache TTL is configurable.
    *   **Dependencies**: `requests`, `pandas`, `concurrent.futures`.
    *   **Notes**: This module is critical for annotating variants with external data, providing context for downstream analysis.

*   **`enhanced_reporting.py`**:
    *   **Purpose**: Generates enhanced HTML reports with interactive visualizations and symptom correlation features, serving as the primary coherence visualization layer for the GASLIT-AF framework.
    *   **Key Data Structures**:
        *   `COMMON_SYMPTOMS`: List of 15 common symptoms associated with GASLIT-AF genes, including fatigue, PEM, cognitive dysfunction, orthostatic intolerance, etc.
        *   `GENE_SYMPTOM_MAPPING`: Dictionary mapping genes to their associated symptoms, creating a coherence network between genetic variants and clinical manifestations.
    *   **Key Functions**:
        *   `generate_enhanced_report(gene_counts, variant_data, figures, output_dir, system_analysis, include_symptoms)`: Orchestrates the report generation process, integrating gene variant data with symptom correlations and interactive visualizations.
        *   `generate_symptom_correlation(gene_counts)`: Generates symptom correlation data based on gene variants, calculating correlation scores and identifying the most relevant symptoms based on the variant profile.
        *   `create_symptom_visualization(symptom_data)`: Creates an interactive Plotly bar chart visualizing symptom correlations with color-coded significance.
        *   `generate_markdown_report(report_data)`: Generates a markdown version of the report for maximum portability.
        *   `generate_html_report(report_data)`: Generates the HTML content with interactive elements, symptom checkboxes, and embedded Plotly visualizations.
        *   `export_report_as_pdf(html_path, output_path)`: Placeholder for PDF export functionality (would typically use weasyprint or headless browser).
    *   **Key Features**: Implements symptom correlation checkboxes with localStorage persistence, interactive Plotly visualizations, responsive design, and a coherence-based approach to connecting genetic variants with clinical symptoms.
    *   **Dependencies**: `pandas`, `plotly`, `jinja2`, `markdown`, `numpy`.
    *   **Notes**: This module represents the recursive integration point between genetic data and clinical manifestations, enabling the visualization of complex relationships between variants and symptoms in an interactive format.

*   **`variant_enrichment.py`**:
    *   **Purpose**: Enriches variant data with AF (Atrial Fibrillation)-specific annotations and insights, providing a specialized coherence layer for cardiac-related variant analysis within the GASLIT-AF framework.
    *   **Key Data Structures**:
        *   `AF_PATHOGENICITY_THRESHOLDS`: Defines precise quantitative thresholds for determining AF pathogenicity, including CADD Phred scores (>15.0), SIFT scores (<0.05), PolyPhen scores (>0.85), and gnomAD allele frequencies (<0.01).
        *   `AF_GENE_CATEGORIES`: Organizes AF-related genes into functional categories (Ion Channels, Structural Proteins, Transcription Factors, Signaling Molecules, Inflammatory Mediators), creating a coherent mapping between genetic variants and cardiac functional systems.
        *   `AF_GENE_TO_CATEGORY`: Reverse mapping from gene to category for efficient lookups.
    *   **Key Class**: `VariantEnricher`
        *   Integrates multiple data sources (VariantAPIIntegration, ClinVarIntegration) with AF-specific knowledge.
        *   Implements caching for performance optimization.
        *   Provides methods for categorizing, filtering, and enriching variants with AF-specific annotations.
    *   **Key Functions**:
        *   `get_af_category(gene)`: Maps genes to their functional categories in cardiac pathways.
        *   `is_af_related_gene(gene)`: Determines if a gene is related to atrial fibrillation pathways.
        *   `is_pathogenic_for_af(variant_data)`: Applies multi-dimensional criteria to determine if a variant is potentially pathogenic for AF, returning both the determination and the reasoning.
        *   `enrich_variant(variant_id)`: Enriches a single variant with AF-specific annotations from multiple data sources.
        *   `enrich_variants(variant_df)`: Processes an entire DataFrame of variants with parallel execution for efficiency.
        *   `generate_af_enrichment_report(variant_df, output_dir)`: Creates a comprehensive markdown report with category-based analysis and clinical recommendations.
        *   `enrich_variants_with_af_data(variant_df, cache_dir, output_dir)`: Convenience function that orchestrates the entire enrichment process.
    *   **Dependencies**: `pandas`, `numpy`, `concurrent.futures`, `api_integration`, `clinvar_integration`, `caching`.
    *   **Notes**: Implements a sophisticated domain-specific layer that transforms general variant data into clinically relevant insights for atrial fibrillation. The module bridges genetic data with functional cardiac pathways and provides actionable clinical recommendations based on variant patterns.

*   **`clinical_integration.py`**:
    *   **Purpose**: Integrates clinical variant data into the variant analysis workflow.
    *   **Key Class**: `ClinicalIntegration`
    *   **Key Functions**:
        *   `load_clinical_data(clinical_data_path)`: Loads clinical data from a JSON file.
        *   `annotate_variants(variants_df)`: Annotates variants with clinical information.
        *   `get_pathogenic_variants(variants_df)`: Filters variants to only include pathogenic or likely pathogenic variants.
        *   `get_variants_by_condition(variants_df, condition_name)`: Filters variants by associated conditions.
        *   `get_clinical_summary(variants_df)`: Generates a summary of clinical findings.
        *   `generate_clinical_report(variants_df, output_dir)`: Generates a clinical report based on variants and clinical data.
    *   **Dependencies**: `pandas`, `clinical_variants.py` (uses `ClinicalVariantManager`).
    *   **Notes**: Enhances variant analysis with clinical context, enabling condition-specific filtering and reporting.

*   **`clinical_variants.py`**:
    *   **Purpose**: Handles the processing, validation, and integration of clinical variant data according to a defined JSON schema.
    *   **Key Class**: `ClinicalVariantManager`
    *   **Key Functions**:
        *   `load_conditions(conditions_file)`: Loads clinical conditions from a JSON file.
        *   `validate_conditions(conditions_data)`: Validates conditions data against the JSON schema.
        *   `add_condition(condition)`: Adds a new condition to the conditions data.
        *   `get_condition_by_gene(gene)`: Gets conditions associated with a specific gene.
        *   `get_condition_by_variant_id(variant_id)`: Gets conditions associated with a specific variant ID.
        *   `get_pathogenic_conditions()`: Gets all pathogenic or likely pathogenic conditions.
        *   `annotate_variants(variants_df)`: Annotates variants with clinical information.
        *   `generate_clinical_report(variants_df, output_dir)`: Generates a clinical report.
    *   **Dependencies**: `jsonschema`, `pandas`, `json`.
    *   **Notes**: Provides the core functionality for clinical variant data management, used by `clinical_integration.py`.

*   **`clinvar_cache.py`**:
    *   **Purpose**: Provides enhanced caching capabilities for ClinVar data, including versioning, incremental updates, and indexing for faster lookups.
    *   **Key Class**: `ClinVarCache`
    *   **Key Functions**:
        *   `get_file_hash(file_path)`: Calculates MD5 hash of a file for integrity checking.
        *   `is_cache_valid(cache_type, max_age_days)`: Checks if a specific cache type is valid based on age.
        *   `update_cache_metadata(cache_type, file_path, version)`: Updates cache metadata after processing a file.
        *   `index_variants(df, source)`: Indexes variants in the SQLite database for efficient retrieval.
        *   `lookup_variant(**kwargs)`: Looks up variants in the index using various parameters (rs_id, clinvar_id, etc.).
        *   `get_cache_stats()`: Returns statistics about the cache (variant counts, pathogenic counts, etc.).
        *   `clear_cache(cache_type)`: Clears the cache (specific type or all).
    *   **Dependencies**: `sqlite3`, `pandas`, `hashlib`.
    *   **Notes**: Implements a sophisticated caching system with SQLite indexing to optimize ClinVar data retrieval and reduce API calls.

*   **`clinvar_integration.py`**:
    *   **Purpose**: Handles downloading, parsing, and integrating ClinVar data with variant analysis results.
    *   **Key Classes**:
        *   `ClinVarDownloader`: Downloads ClinVar data files from NCBI FTP.
        *   `ClinVarParser`: Parses ClinVar data files (variant_summary.txt, VCF).
        *   `ClinVarIntegration`: Integrates ClinVar data with variant analysis.
    *   **Key Functions**:
        *   `download_file(file_type, force_download)`: Downloads a ClinVar file of specified type.
        *   `parse_variant_summary(force_download)`: Parses the variant_summary.txt file.
        *   `parse_vcf(assembly, force_download)`: Parses the ClinVar VCF file.
        *   `annotate_variants(variants_df)`: Annotates variants with ClinVar information.
        *   `get_variant_details(variant_id)`: Gets detailed information for a variant from ClinVar.
    *   **Dependencies**: `requests`, `pandas`, `gzip`, `xml.etree.ElementTree`.
    *   **Notes**: Provides comprehensive ClinVar annotation capabilities, supporting both variant summary and VCF formats for different genome assemblies (GRCh37, GRCh38).
