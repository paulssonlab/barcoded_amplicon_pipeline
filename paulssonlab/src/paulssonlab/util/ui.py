import imageio.v3 as iio
import numpy as np
from IPython.display import Image


def display_image(img, scale=False, downsample=1, format="jpg"):
    if downsample != 1:
        img = img[::downsample, ::downsample, ...]
    if scale is True:
        img_min = img.min()
        img = (img - img_min) / (img.max() - img_min)
    elif scale:
        img = img - img.min()
        img = img / np.percentile(img, scale * 100)
    img = (np.clip(img, 0, 1) * 255).astype(np.uint8)
    if format.lower() in ("jpg", "jpeg"):
        bytes_ = iio.imwrite("<bytes>", img, extension=".jpeg", quality=95)
    elif format.lower() == "png":
        bytes_ = iio.imwrite("<bytes>", img, extension=".png")
    else:
        raise ValueError(f"unknown format: {format}")
    return Image(bytes_)
