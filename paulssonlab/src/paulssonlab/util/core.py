from cytoolz import dissoc


def extract_keys(d, keys):
    return {k: d[k] for k in d.keys() & set(keys)}


def pop_keys(d, keys):
    return extract_keys(d, keys), dissoc(d, *keys)
