from functools import wraps

# FROM: https://stackoverflow.com/a/55896748
# SEE ALSO: https://stackoverflow.com/questions/29402606/possible-to-create-a-synchronized-decorator-thats-aware-of-a-methods-object
def synchronized(lock):
    def wrapper(f):
        @wraps(f)
        def inner_wrapper(*args, **kwargs):
            with lock:
                return f(*args, **kwargs)

        return inner_wrapper

    return wrapper
