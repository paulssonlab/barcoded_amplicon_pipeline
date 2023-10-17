_RC_BASES = ("ACRMBH", "TGYKVD")
RC_MAP = {
    k: v
    for lower in (False, True)
    for a, b in zip(*[bases.lower() if lower else bases for bases in _RC_BASES])
    for k, v in [(a, b), (b, a)]
}


def reverse_complement(seq):
    return "".join(RC_MAP.get(base, base) for base in reversed(seq))
