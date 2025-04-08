"""
Device initialization module for GASLIT-AF Variant Analysis.
Handles GPU/CPU device selection and initialization using Intel oneAPI.
"""

import os
import sys
import logging
import dpctl

# Configure logging
log = logging.getLogger("gaslit-af")

def initialize_device():
    """
    Initialize Intel oneAPI device queue with GPU preference and CPU fallback.
    
    Returns:
        dpctl.SyclQueue: Initialized device queue
    """
    try:
        # Request maximum GPU performance
        os.environ["SYCL_CACHE_PERSISTENT"] = "1"
        
        # Create a queue with GPU selection
        queue = dpctl.SyclQueue("gpu")
        log.info(f"🚀 Using GPU device: {queue.sycl_device}")
        
        # Check if we're actually using the Intel Arc GPU
        if "Arc" in str(queue.sycl_device):
            log.info("✅ Successfully connected to Intel Arc GPU")
        
        return queue
    
    except Exception as e:
        log.warning(f"⚠️ Could not initialize GPU device: {e}")
        log.info("⚙️ Falling back to CPU")
        
        try:
            queue = dpctl.SyclQueue()
            log.info(f"🖥️ Using default device: {queue.sycl_device}")
            return queue
        
        except Exception as e2:
            log.error(f"Could not initialize any SYCL device: {e2}")
            sys.exit(1)

def get_memory_usage():
    """
    Get current memory usage in GB.
    
    Returns:
        float: Current memory usage in GB
    """
    try:
        import psutil
        process = psutil.Process(os.getpid())
        memory_info = process.memory_info()
        return memory_info.rss / (1024 ** 3)  # Convert bytes to GB
    except ImportError:
        log.warning("psutil not available, cannot monitor memory usage")
        return 0.0

def check_memory_limits(current_usage, max_ram_usage=64, ram_buffer=16):
    """
    Check if memory usage is approaching limits.
    
    Args:
        current_usage: Current memory usage in GB
        max_ram_usage: Maximum allowed RAM usage in GB
        ram_buffer: RAM buffer in GB
        
    Returns:
        bool: True if memory usage is within limits, False otherwise
    """
    if current_usage > (max_ram_usage - ram_buffer):
        log.warning(f"⚠️ Memory usage ({current_usage:.2f} GB) is approaching limit ({max_ram_usage} GB)")
        return False
    return True
