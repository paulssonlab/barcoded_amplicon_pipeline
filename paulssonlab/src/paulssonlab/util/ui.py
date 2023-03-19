import numpy as np
import imageio.v3 as iio
from IPython.display import Image


def display_image(img, format="jpg"):
    img = (img * 255).astype(np.uint8)
    if format.lower() in ("jpg", "jpeg"):
        bytes_ = iio.imwrite("<bytes>", img, extension=".jpeg", quality=95)
    elif format.lower() == "png":
        bytes_ = iio.imwrite("<bytes>", img, extension=".png")
    else:
        raise ValueError(f"unknown format: {format}")
    return Image(bytes_)
