# FASTQ to VCF Pipeline for GASLIT-AF Variant Analysis

This pipeline extends the GASLIT-AF Variant Analysis project to process raw FASTQ files through alignment, variant calling, and integration with the existing GASLIT-AF analysis framework. The pipeline leverages Intel oneAPI for optimization where possible.

## Overview

The pipeline consists of the following steps:

1. **Preprocessing (FASTQ → aligned BAM)** using BWA-MEM2
   - Optimized via Intel oneAPI HPC Toolkit
   - Leverages AVX-512 instructions on compatible CPUs

2. **Alignment QC & BAM processing** using Samtools
   - Sort, index, and prepare BAM files for variant calling

3. **Variant calling (BAM → VCF)** using DeepVariant
   - GPU-accelerated with Intel Arc GPU via OpenVINO backend (if available)
   - Produces high-quality variant calls with confidence estimates

4. **Integration with GASLIT-AF analysis**
   - Automatically runs the existing GASLIT-AF analysis on the generated VCF file
   - Produces visualizations and reports for biological system analysis

## Prerequisites

- Python 3.9+
- Intel oneAPI Base Toolkit (optional, for optimization)
- BWA-MEM2 (for alignment)
- Samtools (for BAM processing)
- DeepVariant (for variant calling) or Docker (to run DeepVariant container)
- Reference genome (hg38 or T2T-CHM13)

## Installation

### 1. Clone the repository (if you haven't already)

```bash
git clone <repository-url>
cd gaslit-af_variant_analysis
```

### 2. Install dependencies

The project uses Poetry for dependency management:

```bash
pip install poetry
poetry install
```

### 3. Download reference genome and install tools

We provide a script to download the reference genome and optionally install the necessary tools:

```bash
# Download reference genome only
python download_reference.py --reference hg38 --output-dir reference

# Download reference genome and install tools using Conda
python download_reference.py --reference hg38 --output-dir reference --install-tools --conda-env bioinformatics
```

This will:
- Download the specified reference genome (hg38 or T2T-CHM13)
- Optionally install BWA-MEM2, Samtools, and other tools using Conda or apt-get
- Provide instructions for the next steps

### 4. Index the reference genome

After downloading the reference genome, you need to index it for BWA-MEM2:

```bash
# If using Conda
conda activate bioinformatics

# Index the reference genome
bwa-mem2 index reference/hg38.fa
```

## Usage

### Basic Usage

```bash
python fastq_to_vcf_pipeline.py \
  --fastq1 data/KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.1.fq.gz \
  --fastq2 data/KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.2.fq.gz \
  --reference reference/hg38.fa \
  --output-dir pipeline_output \
  --run-gaslit-analysis
```

### Advanced Usage

```bash
python fastq_to_vcf_pipeline.py \
  --fastq1 data/KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.1.fq.gz \
  --fastq2 data/KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.2.fq.gz \
  --reference reference/hg38.fa \
  --output-dir pipeline_output \
  --sample-name KetanRaturi-WGS \
  --threads 32 \
  --memory 128 \
  --use-gpu \
  --run-gaslit-analysis
```

### Command-line Arguments

#### Required Arguments
- `--fastq1`: Path to first FASTQ file (R1)
- `--fastq2`: Path to second FASTQ file (R2)
- `--reference`: Path to reference genome FASTA file

#### Optional Arguments
- `--output-dir`: Output directory for intermediate and final files (default: "pipeline_output")
- `--sample-name`: Sample name (default: derived from FASTQ filename)
- `--bwa-path`: Path to BWA-MEM2 executable (default: "bwa-mem2")
- `--samtools-path`: Path to Samtools executable (default: "samtools")
- `--deepvariant-path`: Path to DeepVariant executable or Docker image (default: "deepvariant")
- `--threads`: Number of threads to use for parallel processing (default: 16)
- `--memory`: Maximum memory usage in GB (default: 64)
- `--use-gpu`: Use GPU acceleration for DeepVariant (requires OpenVINO)
- `--skip-alignment`: Skip alignment step (use existing BAM file)
- `--skip-variant-calling`: Skip variant calling step (use existing VCF file)
- `--run-gaslit-analysis`: Run GASLIT-AF analysis on the generated VCF file

## Pipeline Output

The pipeline generates the following output files in the specified output directory:

- `{sample_name}.sorted.bam`: Sorted and indexed BAM file
- `{sample_name}.sorted.bam.bai`: BAM index file
- `{sample_name}.vcf.gz`: Compressed VCF file with variant calls
- `{sample_name}.vcf.gz.tbi`: VCF index file
- `gaslit_analysis_{sample_name}/`: Directory containing GASLIT-AF analysis results
  - `system_analysis.json`: JSON file with system-level analysis results
  - `gene_counts.csv`: CSV file with gene-level variant counts
  - `system_distribution.html/.png`: Visualization of variant distribution across biological systems
  - `system_distribution_pie.html/.png`: Pie chart of system distribution
  - `top_genes.html/.png`: Bar chart of top genes by variant count

## Integration with GASLIT-AF Analysis

After generating the VCF file, the pipeline can automatically run the GASLIT-AF analysis using the `direct_system_analysis.py` script. This will:

1. Extract gene information from the VCF file
2. Map variants to GASLIT-AF biological systems
3. Generate visualizations and reports
4. Save the results in the `gaslit_analysis_{sample_name}/` directory

You can then use the existing GASLIT-AF tools to generate combined reports and perform further analysis:

```bash
# Generate combined report
python generate_combined_report.py --results-dir pipeline_output
```

## Performance Optimization

The pipeline is designed to leverage Intel oneAPI for optimization where possible:

- BWA-MEM2 alignment is optimized with Intel oneAPI HPC Toolkit
- DeepVariant can use GPU acceleration via OpenVINO backend on Intel Arc GPUs
- Memory usage is carefully managed to avoid out-of-memory errors
- Parallel processing is used throughout the pipeline

## Troubleshooting

### Common Issues

1. **Out of memory errors during alignment**
   - Reduce the number of threads or increase the memory limit
   - Use a smaller batch size for processing

2. **DeepVariant Docker issues**
   - Ensure Docker is installed and running
   - Make sure the current user has permission to run Docker commands
   - Check that the reference genome and BAM file are accessible to Docker

3. **Missing tools**
   - Use the `download_reference.py` script with the `--install-tools` option to install the necessary tools
   - Alternatively, install the tools manually using Conda or apt-get

### Logs and Debugging

The pipeline generates detailed logs that can help diagnose issues. Look for error messages in the console output and check the following:

- BWA-MEM2 alignment logs
- Samtools processing logs
- DeepVariant logs (in the output directory)
- GASLIT-AF analysis logs

## References

- [BWA-MEM2](https://github.com/bwa-mem2/bwa-mem2)
- [Samtools](http://www.htslib.org/)
- [DeepVariant](https://github.com/google/deepvariant)
- [Intel oneAPI](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
- [GASLIT-AF Variant Analysis](https://github.com/yourusername/gaslit-af_variant_analysis)
