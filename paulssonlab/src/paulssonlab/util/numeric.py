import warnings

import numpy as np


def scale_image(img, scale=False, downsample=None):
    # needed for downsample to work if we have extraneous singleton dimensions up front
    img = img.squeeze()
    if downsample is not None:
        img = img[::downsample, ::downsample, ...]
    if scale is True:
        img_min = np.nanmin(img)
        img = (img - img_min) / (np.nanmax(img) - img_min)
    elif scale:
        img = img - np.nanmin(img)
        img = img / np.nanpercentile(img, scale * 100)
    img = np.clip(img, 0, 1)
    return img


def silent_nanquantile(*args, **kwargs):
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="All-NaN slice encountered")
        return np.nanquantile(*args, **kwargs)
