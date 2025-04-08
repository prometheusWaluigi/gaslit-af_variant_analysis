# GASLIT-AF Systems Submodule

This directory encapsulates the logic for analyzing variants within the context of the broader biological systems defined by the GASLIT-AF framework. It provides the foundational components for mapping genes to their respective systems, efficiently storing and querying system-aware variant data using DuckDB, and applying modular enrichment patterns to annotate variants based on pathogenicity, system relevance, and specific disease contexts (like AF or ME/CFS).

## Modules

*   **`__init__.py`**:
    *   Acts as the public interface for the `systems` submodule.
    *   Exports key components: `GeneSystemManager` (from `gene_systems`), `VariantStore` (from `variant_store_simple`), various `EnrichmentPattern` classes, and factory functions (`create_enrichment_pattern`, `create_standard_enrichment_patterns`) from `enrichment_patterns`.

*   **`gene_systems.py`**:
    *   Defines the `GeneSystemManager` class (implemented as a singleton via `get_gene_system_manager`).
    *   Manages the mapping between individual genes and the GASLIT-AF biological systems (e.g., Immune & Inflammatory, Autonomic & Neurotransmitter, Structural & Connective Tissue, etc.).
    *   Loads system definitions from JSON files (defaulting to `./data/systems/`) or uses hardcoded defaults if files are absent.
    *   Provides methods for querying gene-system relationships (`get_system_for_gene`, `get_genes_for_system`, `get_all_genes`).

*   **`variant_store_simple.py`**:
    *   Provides the primary, active implementation of the variant storage mechanism via the `VariantStore` class (implemented as a singleton via `get_variant_store`).
    *   Uses **DuckDB** as a backend (defaulting to `./data/variants.duckdb`) for efficient, file-based storage and querying.
    *   Defines database schemas (`variants`, `system_variants`) to store variant details (position, alleles, gene, annotations) and link them to their biological system (obtained via `GeneSystemManager`).
    *   Includes methods to store variants (`store_variant`, `store_variants`) and query them based on various criteria (`get_variant`, `get_variants_by_gene`, `get_variants_by_system`, `get_pathogenic_variants`).
    *   Offers functions to generate summary statistics per system (`get_system_summary`) and across all systems (`get_all_systems_summary`).

*   **`enrichment_patterns.py`**:
    *   Defines a composable framework for annotating or enriching variant data based on specific rules or contexts.
    *   Includes a base `EnrichmentPattern` class and specialized subclasses:
        *   `CompositeEnrichmentPattern`: Applies multiple patterns sequentially.
        *   `PathogenicityEnrichmentPattern`: Assesses pathogenicity based on configurable thresholds (e.g., CADD score, allele frequency, ClinVar significance).
        *   `SystemSpecificEnrichmentPattern`: Adds context specific to a variant's biological system.
        *   `AFEnrichmentPattern` / `MECFSEnrichmentPattern`: Apply rules specific to Atrial Fibrillation or ME/CFS context (e.g., flagging variants in known risk genes or pathways).
    *   Provides factory functions (`create_enrichment_pattern`, `create_standard_enrichment_patterns`) to easily instantiate and combine these patterns.

*   **`variant_store.py`**:
    *   _(Note: Appears to be an alternative or older DuckDB-based VariantStore implementation. Based on `__init__.py`, it is **not** the currently active store used by the submodule. `variant_store_simple.py` provides the active implementation.)_

## Purpose

The modules herein facilitate the crucial analytical step of moving beyond individual variant calls to understanding their potential collective impact within the interconnected biological systems relevant to the GASLIT-AF model. The use of DuckDB provides performance, while the enrichment patterns allow for flexible, context-aware annotation.
