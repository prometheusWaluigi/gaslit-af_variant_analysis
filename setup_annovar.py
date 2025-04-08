#!/usr/bin/env python3
"""
ANNOVAR Setup Script for GASLIT-AF Variant Analysis

This script downloads and sets up ANNOVAR for use with the GASLIT-AF framework,
creating a quantum coherence bridge between genomic architecture and the
theoretical model parameters through advanced functional annotations.
"""

import os
import sys
import argparse
import subprocess
import logging
import shutil
import requests
import tarfile
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
log = logging.getLogger("gaslit-af")

# ANNOVAR constants
ANNOVAR_URL = "https://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz"
DEFAULT_INSTALL_DIR = "./tools/annovar"
DEFAULT_HUMANDB_DIR = "./tools/annovar/humandb"
DEFAULT_BUILD = "hg38"

# Essential databases for GASLIT-AF analysis
ESSENTIAL_DATABASES = [
    "refGene",
    "clinvar_20220320",
    "exac03",
    "gnomad211_exome",
    "dbnsfp42a"
]

def parse_args():
    """Parse command-line arguments for ANNOVAR setup."""
    parser = argparse.ArgumentParser(description="ANNOVAR Setup for GASLIT-AF Variant Analysis")
    
    parser.add_argument("--install-dir", type=str, default=DEFAULT_INSTALL_DIR,
                      help=f"ANNOVAR installation directory (default: {DEFAULT_INSTALL_DIR})")
    parser.add_argument("--humandb-dir", type=str, default=DEFAULT_HUMANDB_DIR,
                      help=f"ANNOVAR humandb directory (default: {DEFAULT_HUMANDB_DIR})")
    parser.add_argument("--build", type=str, default=DEFAULT_BUILD,
                      help=f"Genome build version (default: {DEFAULT_BUILD})")
    parser.add_argument("--download-dbs", action="store_true",
                      help="Download essential ANNOVAR databases")
    parser.add_argument("--force", action="store_true",
                      help="Force reinstallation even if ANNOVAR already exists")
    parser.add_argument("--email", type=str,
                      help="Email address for ANNOVAR registration (required for download)")
    
    return parser.parse_args()

def download_annovar(install_dir, email=None):
    """
    Download ANNOVAR from the official website.
    
    Args:
        install_dir: Installation directory
        email: Email address for registration
        
    Returns:
        Path to downloaded file or None if download failed
    """
    if not email:
        log.error("Email address is required for ANNOVAR download")
        log.info("Please register at https://www.openbioinformatics.org/annovar/annovar_download_form.php")
        log.info("Then run this script again with your registered email")
        return None
    
    # Create download directory
    os.makedirs(os.path.dirname(install_dir), exist_ok=True)
    
    # Download URL with email parameter
    download_url = f"{ANNOVAR_URL}/{email}"
    
    # Download file
    log.info(f"Downloading ANNOVAR from {ANNOVAR_URL}")
    log.info("This may take a few minutes...")
    
    try:
        with Progress(
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
            TimeElapsedColumn()
        ) as progress:
            task = progress.add_task("[cyan]Downloading ANNOVAR...", total=None)
            
            response = requests.get(download_url, stream=True)
            if response.status_code != 200:
                progress.update(task, completed=1, total=1)
                log.error(f"Download failed with status code {response.status_code}")
                log.error("Please check your email address and try again")
                return None
            
            # Get file size if available
            total_size = int(response.headers.get('content-length', 0))
            if total_size:
                progress.update(task, total=total_size)
            
            # Download path
            download_path = os.path.join(os.path.dirname(install_dir), "annovar.tar.gz")
            
            # Download file
            with open(download_path, 'wb') as f:
                for chunk in response.iter_content(chunk_size=8192):
                    if chunk:
                        f.write(chunk)
                        if total_size:
                            progress.update(task, advance=len(chunk))
            
            progress.update(task, completed=1, total=1)
        
        log.info(f"ANNOVAR downloaded successfully to {download_path}")
        return download_path
    
    except Exception as e:
        log.error(f"Error downloading ANNOVAR: {e}")
        return None

def extract_annovar(tar_path, install_dir):
    """
    Extract ANNOVAR tarball to installation directory.
    
    Args:
        tar_path: Path to ANNOVAR tarball
        install_dir: Installation directory
        
    Returns:
        True if extraction successful, False otherwise
    """
    try:
        log.info(f"Extracting ANNOVAR to {install_dir}")
        
        # Create installation directory
        os.makedirs(install_dir, exist_ok=True)
        
        # Extract tarball
        with tarfile.open(tar_path, 'r:gz') as tar:
            tar.extractall(path=install_dir)
        
        log.info("ANNOVAR extracted successfully")
        return True
    
    except Exception as e:
        log.error(f"Error extracting ANNOVAR: {e}")
        return False

def download_databases(annovar_dir, humandb_dir, build, databases):
    """
    Download ANNOVAR databases.
    
    Args:
        annovar_dir: ANNOVAR installation directory
        humandb_dir: ANNOVAR humandb directory
        build: Genome build version
        databases: List of databases to download
        
    Returns:
        True if all downloads successful, False otherwise
    """
    log.info(f"Downloading ANNOVAR databases for {build}")
    
    # Create humandb directory
    os.makedirs(humandb_dir, exist_ok=True)
    
    # Download each database
    success = True
    for db in databases:
        log.info(f"Downloading database: {db}")
        
        cmd = [
            "perl", os.path.join(annovar_dir, "annotate_variation.pl"),
            "-buildver", build,
            "-downdb", db,
            humandb_dir
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                log.error(f"Error downloading database {db}: {result.stderr}")
                success = False
            else:
                log.info(f"Successfully downloaded database: {db}")
        except Exception as e:
            log.error(f"Exception downloading database {db}: {e}")
            success = False
    
    return success

def check_annovar_installation(install_dir, humandb_dir):
    """
    Check if ANNOVAR is properly installed.
    
    Args:
        install_dir: ANNOVAR installation directory
        humandb_dir: ANNOVAR humandb directory
        
    Returns:
        Tuple of (is_installed, has_databases)
    """
    # Check for essential ANNOVAR scripts
    required_scripts = ["table_annovar.pl", "convert2annovar.pl", "annotate_variation.pl"]
    is_installed = True
    
    for script in required_scripts:
        script_path = os.path.join(install_dir, script)
        if not os.path.exists(script_path):
            is_installed = False
            break
    
    # Check for humandb directory
    has_databases = os.path.exists(humandb_dir) and len(os.listdir(humandb_dir)) > 0
    
    return is_installed, has_databases

def main():
    """Main entry point for ANNOVAR setup."""
    args = parse_args()
    
    # Convert paths to absolute paths
    install_dir = os.path.abspath(args.install_dir)
    humandb_dir = os.path.abspath(args.humandb_dir)
    
    # Check if ANNOVAR is already installed
    is_installed, has_databases = check_annovar_installation(install_dir, humandb_dir)
    
    if is_installed and not args.force:
        log.info(f"ANNOVAR is already installed at {install_dir}")
        
        # Update paths in configuration
        log.info("Updating ANNOVAR paths in GASLIT-AF configuration")
        log.info(f"ANNOVAR Path: {install_dir}")
        log.info(f"HumanDB Path: {humandb_dir}")
        
        # Check for databases
        if not has_databases and args.download_dbs:
            log.info("No databases found. Downloading essential databases...")
            download_databases(install_dir, humandb_dir, args.build, ESSENTIAL_DATABASES)
        elif not has_databases:
            log.warning("No databases found. Use --download-dbs to download essential databases")
        else:
            log.info("ANNOVAR databases found")
        
        return
    
    # Download and install ANNOVAR
    if args.force and os.path.exists(install_dir):
        log.info(f"Removing existing ANNOVAR installation at {install_dir}")
        shutil.rmtree(install_dir)
    
    # Download ANNOVAR
    tar_path = download_annovar(install_dir, args.email)
    if not tar_path:
        log.error("ANNOVAR download failed")
        return
    
    # Extract ANNOVAR
    if not extract_annovar(tar_path, install_dir):
        log.error("ANNOVAR extraction failed")
        return
    
    # Download databases if requested
    if args.download_dbs:
        log.info("Downloading essential ANNOVAR databases")
        download_databases(install_dir, humandb_dir, args.build, ESSENTIAL_DATABASES)
    else:
        log.info("Skipping database download")
        log.info("Run with --download-dbs to download essential databases")
    
    # Final instructions
    log.info("\n[bold green]ANNOVAR setup completed successfully![/]")
    log.info(f"ANNOVAR installed at: {install_dir}")
    log.info(f"HumanDB directory: {humandb_dir}")
    log.info("\nTo use ANNOVAR with GASLIT-AF, run:")
    log.info(f"python analyze_modular.py your_vcf_file.vcf.gz --use-annovar --annovar-path {install_dir} --humandb-path {humandb_dir}")

if __name__ == "__main__":
    main()
