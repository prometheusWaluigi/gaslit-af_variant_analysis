"""
ClinVar Downloader Module for GASLIT-AF Variant Analysis.

This module is responsible for downloading ClinVar data files from NCBI's FTP server.
"""

import logging
import requests
from datetime import datetime
from pathlib import Path
from typing import Optional, Dict

# Configure logging
log = logging.getLogger("gaslit-af")

class ClinVarDownloader:
    """Handles downloading ClinVar data files."""
    
    def __init__(self, cache_dir: Optional[Path] = None):
        """
        Initialize the ClinVar downloader.
        
        Args:
            cache_dir: Directory to cache downloaded files
        """
        self.cache_dir = Path(cache_dir) if cache_dir else Path("./cache/clinvar")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        
        # Base URLs
        self.ftp_base_url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar"
        
        # Available file types
        self.file_types = {
            "vcv_xml": f"{self.ftp_base_url}/xml/ClinVarVCVRelease_00-latest.xml.gz",
            "rcv_xml": f"{self.ftp_base_url}/xml/RCV_release/ClinVarRCVRelease_00-latest.xml.gz",
            "variant_summary": f"{self.ftp_base_url}/tab_delimited/variant_summary.txt.gz",
            "vcf_grch37": f"{self.ftp_base_url}/vcf_GRCh37/clinvar.vcf.gz",
            "vcf_grch38": f"{self.ftp_base_url}/vcf_GRCh38/clinvar.vcf.gz",
            "var_citations": f"{self.ftp_base_url}/tab_delimited/var_citations.txt",
            "cross_references": f"{self.ftp_base_url}/tab_delimited/cross_references.txt"
        }
    
    def download_file(self, file_type: str, force_download: bool = False) -> Path:
        """
        Download a ClinVar file.
        
        Args:
            file_type: Type of file to download (from self.file_types)
            force_download: Whether to force download even if file exists
            
        Returns:
            Path to the downloaded file
        """
        if file_type not in self.file_types:
            raise ValueError(f"Unknown file type: {file_type}. Available types: {list(self.file_types.keys())}")
        
        url = self.file_types[file_type]
        filename = url.split("/")[-1]
        output_path = self.cache_dir / filename
        
        # Check if file already exists
        if output_path.exists() and not force_download:
            file_age_days = (datetime.now() - datetime.fromtimestamp(output_path.stat().st_mtime)).days
            if file_age_days < 30:
                log.info(f"Using cached {file_type} file (age: {file_age_days} days)")
                return output_path
        
        # Download the file
        log.info(f"Downloading {file_type} from {url}")
        response = requests.get(url, stream=True)
        response.raise_for_status()
        
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        
        log.info(f"Downloaded {file_type} to {output_path}")
        return output_path
    
    def get_latest_release_date(self) -> str:
        """
        Get the release date of the latest ClinVar data.
        
        Returns:
            Date string in the format YYYY-MM-DD
        """
        try:
            url = self.file_types["variant_summary"]
            response = requests.head(url)
            response.raise_for_status()
            
            # Extract date from Last-Modified header
            last_modified = response.headers.get('Last-Modified')
            if last_modified:
                release_date = datetime.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
                return release_date.strftime('%Y-%m-%d')
            
            return "Unknown"
        except Exception as e:
            log.error(f"Error getting latest release date: {e}")
            return "Unknown"
