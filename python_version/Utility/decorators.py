import functools
import time

def retry_and_log(retry_attempts):
    # type: (int, Callable[..., Any]) -> Callable[..., Any]
    """
    A function wrapper for catching all exceptions and logging them, with retry attempts
    """
    def wrapper(f):
        # type: (Callable[..., Any]) -> Callable[..., Any]
        @functools.wraps(f)
        def inner(*args, **kwargs):
            # type: (*Any, **Any) -> Any
            for i in range(0, retry_attempts):
                try:
                    return f(*args, **kwargs)
                except Exception as ex:
                    print(ex)
                    time.sleep(i+1)
        return inner

    return wrapper

