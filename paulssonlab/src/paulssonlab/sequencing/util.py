from glob import glob as glob_


def detect_format(format, filename, formats, glob=False):
    if format is None:
        if glob:
            # get first globbed filename
            filenames = list(glob_(filename))
            if not filenames:
                raise ValueError(
                    f"did not find any matches for glob pattern: {filename}"
                )
            filename = filenames[0]
        for f in formats:
            if filename.endswith(f".{f}"):
                format = f
                break
        else:
            raise ValueError(f"unknown file extension: {filename}")
    return format
