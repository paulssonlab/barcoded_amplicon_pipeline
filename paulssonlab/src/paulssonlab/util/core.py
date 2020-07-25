from itertools import zip_longest


def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


# FROM: https://docs.python.org/3/library/itertools.html#recipes
def grouper(iterable, n, fillvalue=None):
    """Collect data into fixed-length chunks or blocks."""
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)
