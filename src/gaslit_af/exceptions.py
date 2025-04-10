"""
Exceptions module for GASLIT-AF Variant Analysis.
Defines custom exception classes for better error handling.
"""

class GaslitAFError(Exception):
    """Base exception class for all GASLIT-AF errors."""
    
    def __init__(self, message="An error occurred in GASLIT-AF", details=None):
        self.message = message
        self.details = details
        super().__init__(self.message)
    
    def __str__(self):
        if self.details:
            return f"{self.message}: {self.details}"
        return self.message


class DataProcessingError(GaslitAFError):
    """Exception raised for errors during data processing."""
    
    def __init__(self, message="Error during data processing", details=None):
        super().__init__(message, details)


class AnnotationError(GaslitAFError):
    """Exception raised for errors during variant annotation."""
    
    def __init__(self, message="Error during variant annotation", details=None):
        super().__init__(message, details)


class DeviceError(GaslitAFError):
    """Exception raised for errors related to device initialization or operation."""
    
    def __init__(self, message="Error with device initialization or operation", details=None):
        super().__init__(message, details)


class MemoryError(GaslitAFError):
    """Exception raised for memory-related errors."""
    
    def __init__(self, message="Memory error", details=None):
        super().__init__(message, details)


class CacheError(GaslitAFError):
    """Exception raised for errors related to caching."""
    
    def __init__(self, message="Error with cache operations", details=None):
        super().__init__(message, details)


class APIError(GaslitAFError):
    """Exception raised for errors related to API operations."""
    
    def __init__(self, message="Error with API operations", details=None):
        super().__init__(message, details)


class FileError(GaslitAFError):
    """Exception raised for errors related to file operations."""
    
    def __init__(self, message="Error with file operations", details=None):
        super().__init__(message, details)


class ConfigurationError(GaslitAFError):
    """Exception raised for errors related to configuration."""
    
    def __init__(self, message="Error with configuration", details=None):
        super().__init__(message, details)


class VisualizationError(GaslitAFError):
    """Exception raised for errors related to visualization."""
    
    def __init__(self, message="Error generating visualization", details=None):
        super().__init__(message, details)


class ReportingError(GaslitAFError):
    """Exception raised for errors related to report generation."""
    
    def __init__(self, message="Error generating report", details=None):
        super().__init__(message, details)


# Utility function for retrying operations
def retry_operation(max_attempts=3, retry_delay=1, retry_exceptions=(APIError, ConnectionError)):
    """
    Decorator for retrying operations that might fail transiently.
    
    Args:
        max_attempts: Maximum number of attempts
        retry_delay: Delay between attempts in seconds
        retry_exceptions: Tuple of exceptions to catch and retry
        
    Returns:
        Decorator function
    """
    import time
    import functools
    import logging
    
    log = logging.getLogger("gaslit-af")
    
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            attempts = 0
            while attempts < max_attempts:
                try:
                    return func(*args, **kwargs)
                except retry_exceptions as e:
                    attempts += 1
                    if attempts == max_attempts:
                        log.error(f"Operation failed after {max_attempts} attempts: {e}")
                        raise
                    log.warning(f"Operation failed, retrying ({attempts}/{max_attempts}): {e}")
                    time.sleep(retry_delay)
        return wrapper
    return decorator
