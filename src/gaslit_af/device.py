"""
Device initialization module for GASLIT-AF Variant Analysis.
Handles GPU/CPU device selection and initialization using Intel oneAPI.
Provides enhanced memory management with adaptive batch sizing and throttling.
"""

import os
import sys
import time
import logging
import dpctl
from typing import Tuple, Optional

from src.gaslit_af.exceptions import DeviceError, MemoryError

# Configure logging
log = logging.getLogger("gaslit-af")

# Memory usage thresholds
MEMORY_CRITICAL = 0.9  # 90% of max memory
MEMORY_HIGH = 0.8      # 80% of max memory
MEMORY_MODERATE = 0.7  # 70% of max memory
MEMORY_LOW = 0.5       # 50% of max memory

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
            error_msg = f"Could not initialize any SYCL device: {e2}"
            log.error(error_msg)
            raise DeviceError(error_msg)

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
        error_msg = "psutil not available, cannot monitor memory usage"
        log.warning(error_msg)
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
        warning_msg = f"⚠️ Memory usage ({current_usage:.2f} GB) is approaching limit ({max_ram_usage} GB)"
        log.warning(warning_msg)
        return False
    return True

def get_memory_pressure(current_usage, max_ram_usage):
    """
    Calculate memory pressure as a ratio of current usage to max allowed.
    
    Args:
        current_usage: Current memory usage in GB
        max_ram_usage: Maximum allowed RAM usage in GB
        
    Returns:
        float: Memory pressure ratio (0.0 to 1.0)
    """
    return current_usage / max_ram_usage

def get_optimal_batch_size(current_batch_size, max_ram_usage, current_usage, min_batch_size=1000):
    """
    Dynamically adjust batch size based on memory pressure.
    
    Args:
        current_batch_size: Current batch size
        max_ram_usage: Maximum allowed RAM usage in GB
        current_usage: Current memory usage in GB
        min_batch_size: Minimum allowed batch size
        
    Returns:
        int: Optimal batch size
    """
    memory_pressure = get_memory_pressure(current_usage, max_ram_usage)
    
    if memory_pressure > MEMORY_CRITICAL:
        # Critical memory pressure - reduce batch size significantly
        new_batch_size = int(current_batch_size * 0.5)
        log.warning(f"⚠️ Critical memory pressure ({memory_pressure:.2f}), reducing batch size by 50%")
    elif memory_pressure > MEMORY_HIGH:
        # High memory pressure - reduce batch size moderately
        new_batch_size = int(current_batch_size * 0.7)
        log.warning(f"⚠️ High memory pressure ({memory_pressure:.2f}), reducing batch size by 30%")
    elif memory_pressure > MEMORY_MODERATE:
        # Moderate memory pressure - reduce batch size slightly
        new_batch_size = int(current_batch_size * 0.9)
        log.info(f"ℹ️ Moderate memory pressure ({memory_pressure:.2f}), reducing batch size by 10%")
    elif memory_pressure < MEMORY_LOW:
        # Low memory pressure - increase batch size
        new_batch_size = int(current_batch_size * 1.2)
        log.info(f"ℹ️ Low memory pressure ({memory_pressure:.2f}), increasing batch size by 20%")
    else:
        # Normal memory pressure - keep current batch size
        return current_batch_size
    
    # Ensure batch size doesn't go below minimum
    return max(new_batch_size, min_batch_size)

def throttle_processing(memory_pressure):
    """
    Throttle processing based on memory pressure by introducing delays.
    
    Args:
        memory_pressure: Memory pressure ratio (0.0 to 1.0)
    """
    if memory_pressure > MEMORY_CRITICAL:
        # Critical memory pressure - pause processing
        log.warning(f"⚠️ Critical memory pressure ({memory_pressure:.2f}), pausing for 2 seconds")
        time.sleep(2.0)
    elif memory_pressure > MEMORY_HIGH:
        # High memory pressure - slow down processing
        log.warning(f"⚠️ High memory pressure ({memory_pressure:.2f}), pausing for 1 second")
        time.sleep(1.0)
    elif memory_pressure > MEMORY_MODERATE:
        # Moderate memory pressure - slight delay
        log.info(f"ℹ️ Moderate memory pressure ({memory_pressure:.2f}), pausing for 0.5 seconds")
        time.sleep(0.5)

def manage_memory_proactively(current_usage, max_ram_usage, current_batch_size, min_batch_size=1000) -> Tuple[bool, int]:
    """
    Proactively manage memory by adjusting batch size and throttling.
    
    Args:
        current_usage: Current memory usage in GB
        max_ram_usage: Maximum allowed RAM usage in GB
        current_batch_size: Current batch size
        min_batch_size: Minimum allowed batch size
        
    Returns:
        Tuple[bool, int]: (continue_processing, new_batch_size)
    """
    memory_pressure = get_memory_pressure(current_usage, max_ram_usage)
    
    # Adjust batch size
    new_batch_size = get_optimal_batch_size(
        current_batch_size, 
        max_ram_usage, 
        current_usage, 
        min_batch_size
    )
    
    # Apply throttling if needed
    throttle_processing(memory_pressure)
    
    # Determine if processing should continue
    continue_processing = memory_pressure < MEMORY_CRITICAL + 0.05  # Allow a small buffer
    
    if not continue_processing:
        error_msg = f"🛑 Memory pressure too high ({memory_pressure:.2f}), stopping processing"
        log.error(error_msg)
        # Don't raise an exception here, just return False to allow graceful shutdown
    
    return continue_processing, new_batch_size
