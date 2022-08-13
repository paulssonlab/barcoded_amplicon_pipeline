import numpy as np
import warnings


def silent_nanquantile(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        return np.nanquantile(*args, **kwargs)
