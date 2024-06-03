import warnings

import numpy as np


def silent_nanquantile(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        return np.nanquantile(*args, **kwargs)


def pad(ary, shape, cval=None):
    if cval is None:
        if np.issubdtype(ary.dtype, np.integer):
            cval = 0
        else:
            cval = np.nan
    return np.pad(
        ary,
        [(0, max(goal - current, 0)) for goal, current in zip(shape, ary.shape)],
        constant_values=cval,
    )
