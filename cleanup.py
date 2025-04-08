#!/usr/bin/env python3
"""
GASLIT-AF Variant Analysis Cleanup Script

This script recursively cleans up the GASLIT-AF project directory,
removing old analysis results, cached files, and other temporal artifacts
while preserving the essential quantum coherence bridges.

The script uses a fractal pattern recognition approach to identify
redundant temporal patterns and collapse them into higher-order coherence.
"""

import os
import re
import sys
import json
import shutil
import logging
import argparse
from pathlib import Path
from datetime import datetime, timedelta
from typing import List, Dict, Set, Tuple, Any, Optional
from collections import defaultdict
import hashlib

# Configure logging with rich for enhanced visual coherence
try:
    from rich.console import Console
    from rich.logging import RichHandler
    from rich.progress import Progress, TextColumn, BarColumn, TimeElapsedColumn
    from rich.table import Table
    
    console = Console()
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        datefmt="[%X]",
        handlers=[RichHandler(rich_tracebacks=True, console=console)]
    )
    RICH_AVAILABLE = True
except ImportError:
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s: %(message)s",
        datefmt="%H:%M:%S"
    )
    RICH_AVAILABLE = False

log = logging.getLogger("gaslit-af-cleanup")

# Define fractal cleanup patterns
TIMESTAMP_PATTERN = re.compile(r'\d{8}-\d{6}')
DATE_PATTERN = re.compile(r'\d{8}')
HASH_PATTERN = re.compile(r'[a-f0-9]{32}')

# Define directories that should be preserved completely
PRESERVE_DIRS = {
    'src',
    'tests',
    'data',  # Personal genomic data should be preserved
}

# Define directories that should be cleaned thoroughly
CLEAN_DIRS = {
    'analysis_results',
    'results',
    'output',
    'cache',
    '__pycache__',
    '.pytest_cache',
    'benchmark_results',
}

# Define file extensions to clean up
CLEAN_EXTENSIONS = {
    '.pyc',
    '.pyo',
    '.pyd',
    '.so',
    '.o',
    '.a',
    '.log',
    '.tmp',
    '.temp',
}

# Define patterns for redundant files (files with timestamps)
REDUNDANT_PATTERNS = [
    r'.*_\d{8}-\d{6}\.(csv|json|html|md|txt|png|svg)$',
    r'.*_\d{8}\.(csv|json|html|md|txt|png|svg)$',
    r'.*\.\d{8}-\d{6}\.(csv|json|html|md|txt|png|svg)$',
    r'.*\.\d{8}\.(csv|json|html|md|txt|png|svg)$',
]

class FractalCleanup:
    """
    Fractal cleanup system for GASLIT-AF Variant Analysis project.
    
    Uses recursive pattern recognition to identify and remove redundant
    temporal artifacts while preserving the essential quantum coherence bridges.
    """
    
    def __init__(self, root_dir: str, dry_run: bool = True, 
                days_to_keep: int = 7, keep_latest: int = 3,
                interactive: bool = False):
        """
        Initialize the fractal cleanup system.
        
        Args:
            root_dir: Root directory of the project
            dry_run: If True, only show what would be deleted without actually deleting
            days_to_keep: Number of days of data to keep
            keep_latest: Number of latest files to keep for each pattern
            interactive: If True, ask for confirmation before deleting each file
        """
        self.root_dir = Path(root_dir)
        self.dry_run = dry_run
        self.days_to_keep = days_to_keep
        self.keep_latest = keep_latest
        self.interactive = interactive
        
        self.cutoff_date = datetime.now() - timedelta(days=days_to_keep)
        
        # Statistics
        self.stats = {
            'dirs_scanned': 0,
            'files_scanned': 0,
            'dirs_removed': 0,
            'files_removed': 0,
            'bytes_freed': 0,
            'redundant_sets': 0,
            'preserved_files': 0,
        }
        
        # Track deleted files
        self.deleted_files = []
        
        # Track redundant file sets
        self.redundant_sets = defaultdict(list)
        
        # Track file hashes for duplicate detection
        self.file_hashes = {}
        
        log.info(f"Initializing fractal cleanup for {self.root_dir}")
        log.info(f"Mode: {'Dry run' if self.dry_run else 'Actual deletion'}")
        log.info(f"Keeping files newer than {self.cutoff_date.strftime('%Y-%m-%d')}")
        log.info(f"Keeping {self.keep_latest} latest files for each pattern")
    
    def should_preserve_path(self, path: Path) -> bool:
        """
        Check if a path should be preserved based on predefined patterns.
        
        Args:
            path: Path to check
            
        Returns:
            True if the path should be preserved
        """
        # Check if path is in preserved directories
        for preserve_dir in PRESERVE_DIRS:
            if preserve_dir in path.parts:
                # Special case for __pycache__ in preserved directories
                if '__pycache__' in path.parts:
                    return False
                return True
        
        # Always preserve .git directory
        if '.git' in path.parts:
            return True
        
        # Always preserve poetry.lock and pyproject.toml
        if path.name in ('poetry.lock', 'pyproject.toml', '.gitignore', 'README.md', 'PRD.md'):
            return True
        
        # Always preserve Python source files in the root directory
        if path.parent == self.root_dir and path.suffix == '.py':
            return True
        
        return False
    
    def should_clean_path(self, path: Path) -> bool:
        """
        Check if a path should be cleaned based on predefined patterns.
        
        Args:
            path: Path to check
            
        Returns:
            True if the path should be cleaned
        """
        # Check if path is in directories to clean
        for clean_dir in CLEAN_DIRS:
            if clean_dir in path.parts:
                return True
        
        # Check file extensions
        if path.suffix in CLEAN_EXTENSIONS:
            return True
        
        # Check for redundant file patterns
        for pattern in REDUNDANT_PATTERNS:
            if re.match(pattern, path.name):
                return True
        
        return False
    
    def extract_timestamp(self, filename: str) -> Optional[datetime]:
        """
        Extract timestamp from filename.
        
        Args:
            filename: Filename to extract timestamp from
            
        Returns:
            Datetime object if timestamp found, None otherwise
        """
        # Try to extract timestamp in format YYYYMMDD-HHMMSS
        match = TIMESTAMP_PATTERN.search(filename)
        if match:
            timestamp_str = match.group(0)
            try:
                return datetime.strptime(timestamp_str, '%Y%m%d-%H%M%S')
            except ValueError:
                pass
        
        # Try to extract date in format YYYYMMDD
        match = DATE_PATTERN.search(filename)
        if match:
            date_str = match.group(0)
            try:
                return datetime.strptime(date_str, '%Y%m%d')
            except ValueError:
                pass
        
        return None
    
    def compute_file_hash(self, file_path: Path) -> str:
        """
        Compute MD5 hash of a file.
        
        Args:
            file_path: Path to the file
            
        Returns:
            MD5 hash of the file
        """
        try:
            hasher = hashlib.md5()
            with open(file_path, 'rb') as f:
                buf = f.read(65536)
                while len(buf) > 0:
                    hasher.update(buf)
                    buf = f.read(65536)
            return hasher.hexdigest()
        except Exception as e:
            log.warning(f"Error computing hash for {file_path}: {e}")
            return "error"
    
    def group_redundant_files(self):
        """
        Group redundant files based on name patterns and timestamps.
        """
        log.info("Grouping redundant files...")
        
        # Group files by base pattern (removing timestamps)
        pattern_groups = defaultdict(list)
        
        for root, dirs, files in os.walk(self.root_dir):
            root_path = Path(root)
            
            for filename in files:
                file_path = root_path / filename
                
                # Skip if file should be preserved
                if self.should_preserve_path(file_path):
                    continue
                
                # Check for redundant patterns
                for pattern in REDUNDANT_PATTERNS:
                    if re.match(pattern, filename):
                        # Extract base name by removing timestamp
                        base_name = re.sub(r'_\d{8}-\d{6}', '', filename)
                        base_name = re.sub(r'_\d{8}', '', base_name)
                        base_name = re.sub(r'\.\d{8}-\d{6}', '', base_name)
                        base_name = re.sub(r'\.\d{8}', '', base_name)
                        
                        # Group by directory and base name
                        group_key = f"{root}:{base_name}"
                        pattern_groups[group_key].append(file_path)
                        break
        
        # Process each group to keep only the latest N files
        for group_key, file_paths in pattern_groups.items():
            if len(file_paths) <= self.keep_latest:
                continue
            
            # Sort files by modification time (newest first)
            sorted_files = sorted(file_paths, 
                                key=lambda p: p.stat().st_mtime if p.exists() else 0, 
                                reverse=True)
            
            # Keep the latest N files
            keep_files = sorted_files[:self.keep_latest]
            remove_files = sorted_files[self.keep_latest:]
            
            # Add to redundant sets
            if remove_files:
                self.redundant_sets[group_key] = {
                    'keep': [str(p) for p in keep_files],
                    'remove': [str(p) for p in remove_files]
                }
                self.stats['redundant_sets'] += 1
        
        log.info(f"Found {self.stats['redundant_sets']} redundant file sets")
    
    def find_duplicate_files(self):
        """
        Find duplicate files based on content hash.
        """
        log.info("Finding duplicate files...")
        
        # Reset hash dictionary
        self.file_hashes = {}
        
        # Track duplicates
        duplicates = defaultdict(list)
        
        for root, dirs, files in os.walk(self.root_dir):
            root_path = Path(root)
            
            for filename in files:
                file_path = root_path / filename
                
                # Skip if file should be preserved
                if self.should_preserve_path(file_path):
                    continue
                
                # Skip very large files
                try:
                    file_size = file_path.stat().st_size
                    if file_size > 10 * 1024 * 1024:  # Skip files larger than 10MB
                        continue
                except Exception:
                    continue
                
                # Compute hash
                file_hash = self.compute_file_hash(file_path)
                
                if file_hash in self.file_hashes:
                    # This is a duplicate
                    duplicates[file_hash].append(str(file_path))
                else:
                    self.file_hashes[file_hash] = str(file_path)
        
        # Add duplicates to redundant sets
        for file_hash, paths in duplicates.items():
            if len(paths) > 1:
                original = self.file_hashes[file_hash]
                self.redundant_sets[f"duplicate:{file_hash}"] = {
                    'keep': [original],
                    'remove': paths
                }
                self.stats['redundant_sets'] += 1
        
        log.info(f"Found {len(duplicates)} duplicate file sets")
    
    def clean_directory(self, directory: Path):
        """
        Clean a directory recursively.
        
        Args:
            directory: Directory to clean
        """
        if not directory.exists() or not directory.is_dir():
            return
        
        self.stats['dirs_scanned'] += 1
        
        # First, process all files in this directory
        for item in directory.iterdir():
            if item.is_file():
                self.stats['files_scanned'] += 1
                
                # Skip if file should be preserved
                if self.should_preserve_path(item):
                    self.stats['preserved_files'] += 1
                    continue
                
                # Check if file should be cleaned
                should_clean = False
                
                # Check if file is older than cutoff date
                try:
                    mtime = datetime.fromtimestamp(item.stat().st_mtime)
                    if mtime < self.cutoff_date:
                        should_clean = True
                except Exception:
                    pass
                
                # Check if file matches clean patterns
                if self.should_clean_path(item):
                    should_clean = True
                
                # Clean file if needed
                if should_clean:
                    self.remove_file(item)
        
        # Then, recursively process subdirectories
        for item in list(directory.iterdir()):
            if item.is_dir():
                # Skip if directory should be preserved
                if self.should_preserve_path(item):
                    continue
                
                # Recursively clean subdirectory
                self.clean_directory(item)
                
                # Remove directory if empty
                try:
                    if item.exists() and not any(item.iterdir()):
                        self.remove_directory(item)
                except Exception as e:
                    log.warning(f"Error checking if directory is empty: {e}")
    
    def remove_file(self, file_path: Path):
        """
        Remove a file and update statistics.
        
        Args:
            file_path: Path to the file to remove
        """
        try:
            # Get file size for statistics
            file_size = file_path.stat().st_size
            
            # Ask for confirmation if interactive mode is enabled
            if self.interactive:
                if RICH_AVAILABLE:
                    response = console.input(f"Delete file [cyan]{file_path}[/]? (y/n): ")
                else:
                    response = input(f"Delete file {file_path}? (y/n): ")
                
                if response.lower() != 'y':
                    log.info(f"Skipping file: {file_path}")
                    return
            
            # Remove file
            if not self.dry_run:
                file_path.unlink()
                log.info(f"Removed file: {file_path}")
            else:
                log.info(f"Would remove file: {file_path}")
            
            # Update statistics
            self.stats['files_removed'] += 1
            self.stats['bytes_freed'] += file_size
            self.deleted_files.append(str(file_path))
            
        except Exception as e:
            log.warning(f"Error removing file {file_path}: {e}")
    
    def remove_directory(self, dir_path: Path):
        """
        Remove a directory and update statistics.
        
        Args:
            dir_path: Path to the directory to remove
        """
        try:
            # Ask for confirmation if interactive mode is enabled
            if self.interactive:
                if RICH_AVAILABLE:
                    response = console.input(f"Delete directory [cyan]{dir_path}[/]? (y/n): ")
                else:
                    response = input(f"Delete directory {dir_path}? (y/n): ")
                
                if response.lower() != 'y':
                    log.info(f"Skipping directory: {dir_path}")
                    return
            
            # Remove directory
            if not self.dry_run:
                shutil.rmtree(dir_path)
                log.info(f"Removed directory: {dir_path}")
            else:
                log.info(f"Would remove directory: {dir_path}")
            
            # Update statistics
            self.stats['dirs_removed'] += 1
            
        except Exception as e:
            log.warning(f"Error removing directory {dir_path}: {e}")
    
    def clean_redundant_files(self):
        """
        Clean redundant files based on grouped patterns.
        """
        log.info("Cleaning redundant files...")
        
        for group_key, file_set in self.redundant_sets.items():
            log.info(f"Processing redundant set: {group_key}")
            
            # Log files to keep
            for keep_file in file_set['keep']:
                log.info(f"Keeping file: {keep_file}")
            
            # Remove redundant files
            for remove_file in file_set['remove']:
                self.remove_file(Path(remove_file))
    
    def run(self):
        """
        Run the fractal cleanup process.
        """
        log.info("Starting fractal cleanup process...")
        
        # Group redundant files
        self.group_redundant_files()
        
        # Find duplicate files
        self.find_duplicate_files()
        
        # Clean redundant files
        self.clean_redundant_files()
        
        # Clean directories recursively
        self.clean_directory(self.root_dir)
        
        # Print summary
        self.print_summary()
        
        # Save cleanup report
        self.save_report()
        
        log.info("Fractal cleanup process completed.")
    
    def print_summary(self):
        """
        Print summary of cleanup operations.
        """
        if RICH_AVAILABLE:
            table = Table(title="Fractal Cleanup Summary")
            table.add_column("Metric", style="cyan")
            table.add_column("Value", style="green")
            
            table.add_row("Directories Scanned", str(self.stats['dirs_scanned']))
            table.add_row("Files Scanned", str(self.stats['files_scanned']))
            table.add_row("Directories Removed", str(self.stats['dirs_removed']))
            table.add_row("Files Removed", str(self.stats['files_removed']))
            table.add_row("Bytes Freed", f"{self.stats['bytes_freed'] / (1024*1024):.2f} MB")
            table.add_row("Redundant Sets Found", str(self.stats['redundant_sets']))
            table.add_row("Files Preserved", str(self.stats['preserved_files']))
            
            console.print(table)
        else:
            log.info("=== Fractal Cleanup Summary ===")
            log.info(f"Directories Scanned: {self.stats['dirs_scanned']}")
            log.info(f"Files Scanned: {self.stats['files_scanned']}")
            log.info(f"Directories Removed: {self.stats['dirs_removed']}")
            log.info(f"Files Removed: {self.stats['files_removed']}")
            log.info(f"Bytes Freed: {self.stats['bytes_freed'] / (1024*1024):.2f} MB")
            log.info(f"Redundant Sets Found: {self.stats['redundant_sets']}")
            log.info(f"Files Preserved: {self.stats['preserved_files']}")
    
    def save_report(self):
        """
        Save cleanup report to a file.
        """
        timestamp = datetime.now().strftime('%Y%m%d-%H%M%S')
        report_file = self.root_dir / f"cleanup_report_{timestamp}.json"
        
        report = {
            'timestamp': timestamp,
            'mode': 'dry_run' if self.dry_run else 'actual',
            'stats': self.stats,
            'deleted_files': self.deleted_files,
            'redundant_sets': dict(self.redundant_sets),
        }
        
        try:
            with open(report_file, 'w') as f:
                json.dump(report, f, indent=2)
            
            log.info(f"Saved cleanup report to {report_file}")
        except Exception as e:
            log.warning(f"Error saving cleanup report: {e}")

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="GASLIT-AF Variant Analysis Cleanup Script",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        "--root-dir", 
        type=str, 
        default=".",
        help="Root directory of the GASLIT-AF project"
    )
    
    parser.add_argument(
        "--dry-run", 
        action="store_true",
        help="Perform a dry run without actually deleting files"
    )
    
    parser.add_argument(
        "--days-to-keep", 
        type=int, 
        default=7,
        help="Number of days of data to keep"
    )
    
    parser.add_argument(
        "--keep-latest", 
        type=int, 
        default=3,
        help="Number of latest files to keep for each pattern"
    )
    
    parser.add_argument(
        "--interactive", 
        action="store_true",
        help="Ask for confirmation before deleting each file"
    )
    
    return parser.parse_args()

def main():
    """
    Main entry point for the fractal cleanup script.
    """
    args = parse_args()
    
    # Create and run fractal cleanup
    cleanup = FractalCleanup(
        root_dir=args.root_dir,
        dry_run=args.dry_run,
        days_to_keep=args.days_to_keep,
        keep_latest=args.keep_latest,
        interactive=args.interactive
    )
    
    cleanup.run()

if __name__ == "__main__":
    main()
