#!/usr/bin/env python3
"""
Download and prepare reference genome for FASTQ to VCF pipeline.

This script downloads the hg38 or T2T-CHM13 reference genome and prepares it for use
with the FASTQ to VCF pipeline. It also installs the necessary tools if requested.
"""

import os
import sys
import subprocess
import argparse
import logging
from pathlib import Path
import urllib.request
import tarfile
import gzip
import shutil
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
log = logging.getLogger("gaslit-af-reference-downloader")

# Reference genome URLs
REFERENCE_URLS = {
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    "t2t-chm13": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz"
}

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Download and prepare reference genome for FASTQ to VCF pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Reference genome arguments
    parser.add_argument(
        "--reference",
        choices=["hg38", "t2t-chm13"],
        default="hg38",
        help="Reference genome to download"
    )
    parser.add_argument(
        "--output-dir",
        default="reference",
        help="Output directory for reference genome"
    )
    
    # Tool installation arguments
    parser.add_argument(
        "--install-tools",
        action="store_true",
        help="Install necessary tools (BWA-MEM2, Samtools, Docker for DeepVariant)"
    )
    parser.add_argument(
        "--conda-env",
        default="bioinformatics",
        help="Conda environment name for tool installation"
    )
    
    # Parse arguments
    args = parser.parse_args()
    
    # Convert string paths to Path objects
    args.output_dir = Path(args.output_dir)
    
    return args

def download_file(url, output_path, desc=None):
    """Download a file with progress reporting."""
    if not desc:
        desc = f"Downloading {url.split('/')[-1]}"
    
    # Create parent directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Check if file already exists
    if output_path.exists():
        log.info(f"File already exists: {output_path}")
        return True
    
    try:
        with Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("{task.completed}/{task.total}"),
            TextColumn("({task.percentage:.0f}%)"),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            # Get file size
            with urllib.request.urlopen(url) as response:
                file_size = int(response.info().get("Content-Length", 0))
                
                # Create download task
                task = progress.add_task(desc, total=file_size)
                
                # Download file
                with urllib.request.urlopen(url) as response, open(output_path, "wb") as out_file:
                    chunk_size = 1024 * 1024  # 1 MB
                    downloaded = 0
                    
                    while True:
                        chunk = response.read(chunk_size)
                        if not chunk:
                            break
                        
                        out_file.write(chunk)
                        downloaded += len(chunk)
                        progress.update(task, completed=downloaded)
        
        log.info(f"Downloaded {url} to {output_path}")
        return True
    except Exception as e:
        log.error(f"Error downloading {url}: {e}")
        return False

def extract_gzip(input_path, output_path):
    """Extract a gzipped file."""
    log.info(f"Extracting {input_path} to {output_path}")
    
    # Check if output file already exists
    if output_path.exists():
        log.info(f"Output file already exists: {output_path}")
        return True
    
    try:
        with Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TimeElapsedColumn(),
            console=console
        ) as progress:
            # Create extraction task
            task = progress.add_task(f"Extracting {input_path.name}", total=None)
            
            # Extract file
            with gzip.open(input_path, "rb") as f_in, open(output_path, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            
            progress.update(task, completed=1, total=1)
        
        log.info(f"Extracted {input_path} to {output_path}")
        return True
    except Exception as e:
        log.error(f"Error extracting {input_path}: {e}")
        return False

def install_tools_conda(args):
    """Install necessary tools using Conda."""
    log.info(f"Installing tools using Conda in environment: {args.conda_env}")
    
    try:
        # Check if conda is installed
        subprocess.run(["conda", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        
        # Check if environment exists
        result = subprocess.run(["conda", "env", "list"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        env_exists = args.conda_env in result.stdout.decode()
        
        if env_exists:
            log.info(f"Conda environment {args.conda_env} already exists")
        else:
            # Create environment
            log.info(f"Creating Conda environment: {args.conda_env}")
            subprocess.run(["conda", "create", "-y", "-n", args.conda_env, "python=3.9"], check=True)
        
        # Install tools
        log.info(f"Installing bioinformatics tools in {args.conda_env}")
        subprocess.run([
            "conda", "install", "-y", "-n", args.conda_env,
            "-c", "bioconda", "-c", "conda-forge",
            "bwa-mem2", "samtools", "picard", "gatk4"
        ], check=True)
        
        # Print activation instructions
        log.info(f"Tools installed successfully in Conda environment: {args.conda_env}")
        log.info(f"To activate the environment, run: conda activate {args.conda_env}")
        
        return True
    except subprocess.SubprocessError as e:
        log.error(f"Error installing tools with Conda: {e}")
        return False
    except FileNotFoundError:
        log.error("Conda not found. Please install Conda first.")
        return False

def install_tools_apt(args):
    """Install necessary tools using apt-get (for Ubuntu/Debian)."""
    log.info("Installing tools using apt-get")
    
    try:
        # Check if we're on Ubuntu/Debian
        result = subprocess.run(["lsb_release", "-a"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if "Ubuntu" not in result.stdout.decode() and "Debian" not in result.stdout.decode():
            log.warning("Not running Ubuntu or Debian, apt-get installation may fail")
        
        # Update package list
        log.info("Updating package list")
        subprocess.run(["sudo", "apt-get", "update"], check=True)
        
        # Install Docker
        log.info("Installing Docker")
        subprocess.run([
            "sudo", "apt-get", "install", "-y",
            "docker.io", "docker-compose"
        ], check=True)
        
        # Start Docker service
        log.info("Starting Docker service")
        subprocess.run(["sudo", "systemctl", "start", "docker"], check=True)
        subprocess.run(["sudo", "systemctl", "enable", "docker"], check=True)
        
        # Add current user to Docker group
        log.info("Adding current user to Docker group")
        subprocess.run(["sudo", "usermod", "-aG", "docker", os.environ["USER"]], check=True)
        
        # Install BWA-MEM2 and Samtools
        log.info("Installing BWA-MEM2 and Samtools")
        subprocess.run([
            "sudo", "apt-get", "install", "-y",
            "bwa", "samtools"
        ], check=True)
        
        log.info("Tools installed successfully using apt-get")
        log.info("You may need to log out and log back in for Docker group changes to take effect")
        
        return True
    except subprocess.SubprocessError as e:
        log.error(f"Error installing tools with apt-get: {e}")
        return False

def download_reference_genome(args):
    """Download and prepare reference genome."""
    log.info(f"Downloading {args.reference} reference genome")
    
    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get reference URL
    reference_url = REFERENCE_URLS.get(args.reference)
    if not reference_url:
        log.error(f"Unknown reference genome: {args.reference}")
        return False
    
    # Download reference genome
    reference_gz = args.output_dir / f"{args.reference}.fa.gz"
    if not download_file(reference_url, reference_gz, f"Downloading {args.reference} reference genome"):
        return False
    
    # Extract reference genome
    reference_fa = args.output_dir / f"{args.reference}.fa"
    if not extract_gzip(reference_gz, reference_fa):
        return False
    
    log.info(f"Reference genome downloaded and extracted successfully: {reference_fa}")
    
    # Print next steps
    log.info("\nNext steps:")
    log.info(f"1. Index the reference genome using BWA-MEM2:")
    log.info(f"   bwa-mem2 index {reference_fa}")
    log.info(f"2. Run the FASTQ to VCF pipeline:")
    log.info(f"   python fastq_to_vcf_pipeline.py --fastq1 data/KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.1.fq.gz --fastq2 data/KetanRaturi-SQ4T88M8-30x-WGS-Sequencing_com-08-24-24.2.fq.gz --reference {reference_fa} --output-dir pipeline_output --run-gaslit-analysis")
    
    return True

def main():
    """Main function to download reference genome and install tools."""
    # Parse command-line arguments
    args = parse_args()
    
    log.info(f"{'='*80}")
    log.info(f"Reference Genome Downloader for FASTQ to VCF Pipeline")
    log.info(f"{'='*80}")
    
    # Install tools if requested
    if args.install_tools:
        # Try Conda installation first
        if not install_tools_conda(args):
            # Fall back to apt-get installation
            log.warning("Conda installation failed, trying apt-get installation")
            if not install_tools_apt(args):
                log.error("Tool installation failed")
                sys.exit(1)
    
    # Download reference genome
    if not download_reference_genome(args):
        log.error("Reference genome download failed")
        sys.exit(1)
    
    log.info(f"{'='*80}")
    log.info(f"Reference genome preparation completed successfully!")
    log.info(f"{'='*80}")

if __name__ == "__main__":
    main()
