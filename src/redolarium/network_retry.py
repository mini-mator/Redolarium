import time
import logging
import functools

logger = logging.getLogger("redolarium.network_retry")

def with_retry(max_retries=5, base_delay=2.0, max_delay=32.0, exceptions=(Exception,)):
    """
    Decorator for wrapping external network API calls (Entrez, BLAST, KEGG, MIBiG) with an 
    exponential backoff retry queue to prevent pipeline failure due to transient network drops.
    """
    def decorator(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            delay = base_delay
            for attempt in range(max_retries):
                try:
                    return func(*args, **kwargs)
                except exceptions as e:
                    if attempt == max_retries - 1:
                        logger.error(f"CRITICAL: Function '{func.__name__}' failed after {max_retries} attempts. Final error: {e}")
                        raise
                    logger.warning(f"Network error in '{func.__name__}': {e}. Retrying in {delay} seconds... (Attempt {attempt+1}/{max_retries})")
                    time.sleep(delay)
                    delay = min(delay * 2, max_delay)
        return wrapper
    return decorator
