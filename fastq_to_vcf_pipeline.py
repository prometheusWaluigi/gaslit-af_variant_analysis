#!/usr/bin/env python3
"""
FASTQ to VCF Pipeline for GASLIT-AF Variant Analysis

This script implements a complete genomic analysis pipeline:
1. Preprocessing (FASTQ → aligned BAM) using BWA-MEM2
2. Alignment QC & BAM processing using Samtools
3. Variant calling (BAM → VCF) using DeepVariant
4. Integration with GASLIT-AF analysis

The pipeline leverages Intel oneAPI for optimization where possible.
"""

import os
import sys
import subprocess
import argparse
import logging
from pathlib import Path
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
log = logging.getLogger("gaslit-af-fastq-pipeline")

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="FASTQ to VCF Pipeline for GASLIT-AF Variant Analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Input/output arguments
    parser.add_argument(
        "--fastq1",
        required=True,
        help="Path to first FASTQ file (R1)"
    )
    parser.add_argument(
        "--fastq2",
        required=True,
        help="Path to second FASTQ file (R2)"
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Path to reference genome FASTA file"
    )
    parser.add_argument(
        "--output-dir",
        default="pipeline_output",
        help="Output directory for intermediate and final files"
    )
    parser.add_argument(
        "--sample-name",
        help="Sample name (default: derived from FASTQ filename)"
    )
    
    # Tool paths
    parser.add_argument(
        "--bwa-path",
        default="bwa-mem2",
        help="Path to BWA-MEM2 executable"
    )
    parser.add_argument(
        "--samtools-path",
        default="samtools",
        help="Path to Samtools executable"
    )
    parser.add_argument(
        "--deepvariant-path",
        default="deepvariant",
        help="Path to DeepVariant executable or Docker image"
    )
    
    # Performance arguments
    parser.add_argument(
        "--threads",
        type=int,
        default=16,
        help="Number of threads to use for parallel processing"
    )
    parser.add_argument(
        "--memory",
        type=int,
        default=64,
        help="Maximum memory usage in GB"
    )
    parser.add_argument(
        "--use-gpu",
        action="store_true",
        help="Use GPU acceleration for DeepVariant (requires OpenVINO)"
    )
    
    # Pipeline control
    parser.add_argument(
        "--skip-alignment",
        action="store_true",
        help="Skip alignment step (use existing BAM file)"
    )
    parser.add_argument(
        "--skip-variant-calling",
        action="store_true",
        help="Skip variant calling step (use existing VCF file)"
    )
    parser.add_argument(
        "--run-gaslit-analysis",
        action="store_true",
        help="Run GASLIT-AF analysis on the generated VCF file"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert string paths to Path objects
    args.fastq1 = Path(args.fastq1)
    args.fastq2 = Path(args.fastq2)
    args.reference = Path(args.reference)
    args.output_dir = Path(args.output_dir)
    
    # Derive sample name from FASTQ filename if not provided
    if not args.sample_name:
        args.sample_name = args.fastq1.stem.split('.')[0]
        if args.sample_name.endswith('_1') or args.sample_name.endswith('.1'):
            args.sample_name = args.sample_name[:-2]
    
    return args

def check_dependencies(args):
    """Check if required tools are installed."""
    log.info("Checking dependencies...")
    
    # Check if BWA-MEM2 is installed
    try:
        subprocess.run([args.bwa_path, "version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        log.info(f"✓ BWA-MEM2 found at {args.bwa_path}")
    except (subprocess.SubprocessError, FileNotFoundError):
        log.error(f"✗ BWA-MEM2 not found at {args.bwa_path}")
        log.error("Please install BWA-MEM2 or provide the correct path using --bwa-path")
        return False
    
    # Check if Samtools is installed
    try:
        subprocess.run([args.samtools_path, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        log.info(f"✓ Samtools found at {args.samtools_path}")
    except (subprocess.SubprocessError, FileNotFoundError):
        log.error(f"✗ Samtools not found at {args.samtools_path}")
        log.error("Please install Samtools or provide the correct path using --samtools-path")
        return False
    
    # Check if DeepVariant is installed (this is more complex as it might be a Docker image)
    # Skip this check if variant calling is skipped
    if not args.skip_variant_calling:
        # For simplicity, we'll just check if the command exists
        if args.deepvariant_path == "deepvariant":
            try:
                # Check if run_deepvariant exists
                subprocess.run(["which", "run_deepvariant"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                log.info(f"✓ DeepVariant found (run_deepvariant)")
            except (subprocess.SubprocessError, FileNotFoundError):
                # Check if Docker is installed (for DeepVariant Docker image)
                try:
                    subprocess.run(["sudo", "docker", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                    log.info(f"✓ Docker found (will use DeepVariant Docker image with sudo)")
                except (subprocess.SubprocessError, FileNotFoundError):
                    log.error(f"✗ DeepVariant not found and Docker not installed")
                    log.error("Please install DeepVariant or Docker, or provide the correct path using --deepvariant-path")
                    return False
        else:
            # Custom DeepVariant path provided, check if it exists
            if not Path(args.deepvariant_path).exists():
                log.error(f"✗ DeepVariant not found at {args.deepvariant_path}")
                log.error("Please install DeepVariant or provide the correct path using --deepvariant-path")
                return False
            log.info(f"✓ Using custom DeepVariant path: {args.deepvariant_path}")
    else:
        log.info("Skipping DeepVariant check as variant calling is disabled")
    
    # Check if reference genome exists
    if not args.reference.exists():
        log.error(f"✗ Reference genome not found at {args.reference}")
        log.error("Please provide the correct path to the reference genome using --reference")
        return False
    log.info(f"✓ Reference genome found at {args.reference}")
    
    # Check if reference genome is indexed
    if not args.skip_alignment:
        reference_index = Path(f"{args.reference}.bwt")
        if not reference_index.exists():
            log.warning(f"! Reference genome index not found at {reference_index}")
            log.warning("Will need to index reference genome before alignment")
        else:
            log.info(f"✓ Reference genome index found at {reference_index}")
    
    # Check if FASTQ files exist
    if not args.skip_alignment:
        if not args.fastq1.exists():
            log.error(f"✗ FASTQ file 1 not found at {args.fastq1}")
            log.error("Please provide the correct path to the FASTQ file using --fastq1")
            return False
        log.info(f"✓ FASTQ file 1 found at {args.fastq1}")
        
        if not args.fastq2.exists():
            log.error(f"✗ FASTQ file 2 not found at {args.fastq2}")
            log.error("Please provide the correct path to the FASTQ file using --fastq2")
            return False
        log.info(f"✓ FASTQ file 2 found at {args.fastq2}")
    
    # Check if Intel oneAPI is installed
    oneapi_path = Path("/opt/intel/oneapi")
    if oneapi_path.exists():
        log.info(f"✓ Intel oneAPI found at {oneapi_path}")
        
        # Check if OpenVINO is installed
        openvino_path = oneapi_path / "openvino"
        if openvino_path.exists():
            log.info(f"✓ OpenVINO found at {openvino_path}")
        else:
            log.warning(f"! OpenVINO not found at {openvino_path}")
            if args.use_gpu:
                log.warning("GPU acceleration may not be available without OpenVINO")
    else:
        log.warning(f"! Intel oneAPI not found at {oneapi_path}")
        log.warning("Performance optimizations may not be available without Intel oneAPI")
    
    return True

def index_reference_genome(args):
    """Index the reference genome using BWA-MEM2."""
    log.info(f"Indexing reference genome: {args.reference}")
    
    # Check if reference genome is already indexed
    reference_index = Path(f"{args.reference}.bwt")
    if reference_index.exists():
        log.info(f"Reference genome already indexed at {reference_index}")
        return True
    
    # Index reference genome
    try:
        cmd = [args.bwa_path, "index", str(args.reference)]
        log.info(f"Running command: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True)
        
        log.info(f"Reference genome indexed successfully")
        return True
    except subprocess.SubprocessError as e:
        log.error(f"Error indexing reference genome: {e}")
        return False

def align_fastq_to_bam(args):
    """Align FASTQ files to reference genome using BWA-MEM2."""
    if args.skip_alignment:
        log.info("Skipping alignment step as requested")
        return True
    
    log.info(f"Aligning FASTQ files to reference genome")
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Define output files
    sam_file = args.output_dir / f"{args.sample_name}.sam"
    bam_file = args.output_dir / f"{args.sample_name}.bam"
    sorted_bam_file = args.output_dir / f"{args.sample_name}.sorted.bam"
    
    # Step 1: Align FASTQ to SAM using BWA-MEM2
    try:
        # Set up Intel oneAPI environment if available
        env = os.environ.copy()
        oneapi_vars = Path("/opt/intel/oneapi/setvars.sh")
        if oneapi_vars.exists():
            log.info(f"Setting up Intel oneAPI environment")
            # We can't source the script directly in Python, so we'll use a workaround
            # This is a simplified approach and may not capture all environment variables
            cmd = [f"source {oneapi_vars} && env"]
            proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            for line in proc.stdout:
                (key, _, value) = line.decode().partition("=")
                env[key] = value.rstrip()
        
        # Align FASTQ to SAM
        cmd = [
            args.bwa_path, "mem",
            "-t", str(args.threads),
            "-R", f"@RG\\tID:{args.sample_name}\\tSM:{args.sample_name}\\tPL:ILLUMINA",
            str(args.reference),
            str(args.fastq1),
            str(args.fastq2)
        ]
        log.info(f"Running alignment: {' '.join(cmd)}")
        
        with open(sam_file, "w") as f:
            subprocess.run(cmd, stdout=f, env=env, check=True)
        
        log.info(f"Alignment completed successfully: {sam_file}")
        
        # Step 2: Convert SAM to BAM
        cmd = [
            args.samtools_path, "view",
            "-@", str(args.threads),
            "-b", "-o", str(bam_file),
            str(sam_file)
        ]
        log.info(f"Converting SAM to BAM: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True)
        
        log.info(f"SAM to BAM conversion completed successfully: {bam_file}")
        
        # Step 3: Sort BAM file
        cmd = [
            args.samtools_path, "sort",
            "-@", str(args.threads),
            "-m", f"{args.memory // args.threads}G",
            "-o", str(sorted_bam_file),
            str(bam_file)
        ]
        log.info(f"Sorting BAM file: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True)
        
        log.info(f"BAM sorting completed successfully: {sorted_bam_file}")
        
        # Step 4: Index sorted BAM file
        cmd = [
            args.samtools_path, "index",
            "-@", str(args.threads),
            str(sorted_bam_file)
        ]
        log.info(f"Indexing sorted BAM file: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True)
        
        log.info(f"BAM indexing completed successfully: {sorted_bam_file}.bai")
        
        # Step 5: Remove intermediate files
        if sam_file.exists():
            log.info(f"Removing intermediate SAM file: {sam_file}")
            sam_file.unlink()
        
        if bam_file.exists():
            log.info(f"Removing intermediate BAM file: {bam_file}")
            bam_file.unlink()
        
        return True
    except subprocess.SubprocessError as e:
        log.error(f"Error in alignment pipeline: {e}")
        return False

def call_variants(args):
    """Call variants using DeepVariant."""
    if args.skip_variant_calling:
        log.info("Skipping variant calling step as requested")
        return True
    
    log.info(f"Calling variants using DeepVariant")
    
    # Define input and output files
    sorted_bam_file = args.output_dir / f"{args.sample_name}.sorted.bam"
    vcf_output = args.output_dir / f"{args.sample_name}.vcf.gz"
    
    # Check if input BAM file exists
    if not sorted_bam_file.exists():
        log.error(f"Input BAM file not found: {sorted_bam_file}")
        return False
    
    try:
        # Determine if we're using native DeepVariant or Docker
        if args.deepvariant_path == "deepvariant":
            # Check if run_deepvariant exists
            try:
                subprocess.run(["which", "run_deepvariant"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
                # Use native DeepVariant
                cmd = [
                    "run_deepvariant",
                    "--model_type=WGS",  # Whole Genome Sequencing
                    f"--ref={args.reference}",
                    f"--reads={sorted_bam_file}",
                    f"--output_vcf={vcf_output}",
                    f"--num_shards={args.threads}"
                ]
                
                # Add GPU acceleration if requested
                if args.use_gpu:
                    cmd.append("--use_gpu")
                
                log.info(f"Running DeepVariant: {' '.join(cmd)}")
                subprocess.run(cmd, check=True)
            except (subprocess.SubprocessError, FileNotFoundError):
                # Use Docker
                log.info("Using DeepVariant Docker image")
                
                # Create output directory for DeepVariant
                deepvariant_output = args.output_dir / "deepvariant_output"
                deepvariant_output.mkdir(parents=True, exist_ok=True)
                
                # Get absolute paths for Docker volume mapping
                ref_abs = args.reference.absolute()
                bam_abs = sorted_bam_file.absolute()
                output_abs = deepvariant_output.absolute()
                
                # Determine Docker image based on GPU usage
                if args.use_gpu:
                    docker_image = "google/deepvariant:1.5.0-gpu"
                    runtime_arg = "--gpus all"
                else:
                    docker_image = "google/deepvariant:1.5.0"
                    runtime_arg = ""
                
                # Run DeepVariant using Docker with sudo
                cmd = [
                    "sudo", "docker", "run",
                    runtime_arg,
                    "-v", f"{ref_abs.parent}:/input-ref",
                    "-v", f"{bam_abs.parent}:/input-bam",
                    "-v", f"{output_abs}:/output",
                    docker_image,
                    "/opt/deepvariant/bin/run_deepvariant",
                    "--model_type=WGS",
                    f"--ref=/input-ref/{ref_abs.name}",
                    f"--reads=/input-bam/{bam_abs.name}",
                    "--output_vcf=/output/output.vcf.gz",
                    f"--num_shards={args.threads}"
                ]
                
                # Remove empty strings from command
                cmd = [x for x in cmd if x]
                
                log.info(f"Running DeepVariant Docker: {' '.join(cmd)}")
                subprocess.run(cmd, check=True)
                
                # Copy output files to expected location
                docker_output_vcf = deepvariant_output / "output.vcf.gz"
                if docker_output_vcf.exists():
                    import shutil
                    shutil.copy(docker_output_vcf, vcf_output)
                    log.info(f"Copied DeepVariant output to {vcf_output}")
                    
                    # Also copy the index file if it exists
                    docker_output_vcf_index = deepvariant_output / "output.vcf.gz.tbi"
                    if docker_output_vcf_index.exists():
                        shutil.copy(docker_output_vcf_index, f"{vcf_output}.tbi")
                        log.info(f"Copied DeepVariant output index to {vcf_output}.tbi")
                else:
                    log.error(f"DeepVariant output not found: {docker_output_vcf}")
                    return False
        else:
            # Use custom DeepVariant path
            cmd = [
                args.deepvariant_path,
                "--model_type=WGS",
                f"--ref={args.reference}",
                f"--reads={sorted_bam_file}",
                f"--output_vcf={vcf_output}",
                f"--num_shards={args.threads}"
            ]
            
            # Add GPU acceleration if requested
            if args.use_gpu:
                cmd.append("--use_gpu")
            
            log.info(f"Running DeepVariant: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
        
        log.info(f"Variant calling completed successfully: {vcf_output}")
        return True
    except subprocess.SubprocessError as e:
        log.error(f"Error calling variants: {e}")
        return False

def run_gaslit_analysis(args):
    """Run GASLIT-AF analysis on the generated VCF file."""
    if not args.run_gaslit_analysis:
        log.info("Skipping GASLIT-AF analysis as requested")
        return True
    
    log.info(f"Running GASLIT-AF analysis")
    
    # Define input VCF file
    vcf_file = args.output_dir / f"{args.sample_name}.vcf.gz"
    
    # Check if input VCF file exists
    if not vcf_file.exists():
        log.error(f"Input VCF file not found: {vcf_file}")
        return False
    
    try:
        # Run direct system analysis
        cmd = [
            "python", "direct_system_analysis.py",
            str(vcf_file),
            "--output-dir", str(args.output_dir / f"gaslit_analysis_{args.sample_name}")
        ]
        log.info(f"Running GASLIT-AF analysis: {' '.join(cmd)}")
        
        subprocess.run(cmd, check=True)
        
        log.info(f"GASLIT-AF analysis completed successfully")
        return True
    except subprocess.SubprocessError as e:
        log.error(f"Error running GASLIT-AF analysis: {e}")
        return False

def main():
    """Main function to run the FASTQ to VCF pipeline."""
    # Parse command-line arguments
    args = parse_args()
    
    log.info(f"{'='*80}")
    log.info(f"FASTQ to VCF Pipeline for GASLIT-AF Variant Analysis")
    log.info(f"{'='*80}")
    
    # Check dependencies
    if not check_dependencies(args):
        log.error("Dependency check failed, aborting pipeline")
        sys.exit(1)
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Index reference genome if needed
    if not args.skip_alignment and not index_reference_genome(args):
        log.error("Reference genome indexing failed, aborting pipeline")
        sys.exit(1)
    
    # Step 2: Align FASTQ to BAM
    if not align_fastq_to_bam(args):
        log.error("Alignment failed, aborting pipeline")
        sys.exit(1)
    
    # Step 3: Call variants
    if not call_variants(args):
        log.error("Variant calling failed, aborting pipeline")
        sys.exit(1)
    
    # Step 4: Run GASLIT-AF analysis
    if not run_gaslit_analysis(args):
        log.error("GASLIT-AF analysis failed")
        # Don't exit here, as the VCF file was still generated successfully
    
    log.info(f"{'='*80}")
    log.info(f"Pipeline completed successfully!")
    log.info(f"Output directory: {args.output_dir}")
    log.info(f"{'='*80}")

if __name__ == "__main__":
    main()
