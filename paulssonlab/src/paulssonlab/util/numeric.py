import warnings

import numpy as np


def silent_nanquantile(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        return np.nanquantile(*args, **kwargs)
