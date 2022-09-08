import functools
import time

from threading import Lock

from dxcint.Logger import Logger

logger = Logger.make(name=__name__, verbosity="info")

DEFAULT_INSTANCE_TYPE = "mem1_ssd1_v2_x4"


def rm_suffix(original_string: str, suffix: str) -> str:
    if not original_string:
        return ""
    if not suffix:
        return original_string
    if original_string[-1] == suffix[-1]:
        return rm_suffix(original_string[:-1], suffix[:-1])
    else:
        raise ValueError("rm_suffix(): suffix is not present in the original string")


def rm_prefix(original_string: str, prefix: str) -> str:
    if not original_string:
        return ""
    if not prefix:
        return original_string
    if original_string[0] == prefix[0]:
        return rm_prefix(original_string[1:], prefix[1:])
    else:
        raise ValueError("rm_suffix(): prefix is not present in the original string")


def async_retry(max_retries: int = 5, delay: int = 5):
    """
    A decorator function to perform async retry of a decorated method.
    Args:
        max_retries: int. Number of maximum retries.
        delay: int. Amount of time to sleep in seconds between retries.
    Returns: A result of a decorated callable.
    """

    def async_retry_inner(func):
        @functools.wraps(func)
        def async_retry_wrapper(*args, **kwargs):
            try:
                lock = args[
                    0
                ].context.lock  # For methods of a class with Context property
            except AttributeError:
                lock = Lock()  # For other functions.
            for i in range(0, max_retries):
                try:
                    with lock:
                        logger.info(
                            f"Retry {i} for function `{func.__name__}({args}, {kwargs})`"
                        )
                    ret_value = func(*args, **kwargs)
                    return ret_value
                except Exception as e:
                    with lock:
                        logger.info(
                            f"Error when running an async retry for function `{func.__name__}`\n"
                            f"Error CONTENT: {e}\n"
                            f"With ARGS: {args}\n"
                            f"With KWARGS: {kwargs}\n"
                            f"Retry in {delay} sec"
                        )
                    time.sleep(delay)
            else:
                raise Exception(
                    f"Failed after {max_retries} retries for function `{func.__name__}`\n"
                    f"With ARGS: {args}\n"
                    f"With KWARGS: {kwargs}"
                )

        return async_retry_wrapper

    return async_retry_inner
