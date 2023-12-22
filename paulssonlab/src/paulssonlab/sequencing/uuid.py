import numba
import numpy as np

IDXS_WITHOUT_HYPHENS = list(set(range(36)) - set([8, 13, 18, 23]))

hex_to_nibble = np.zeros(256, dtype=np.uint8)
for char in range(ord("0"), ord("9") + 1):
    hex_to_nibble[char] = char - ord("0")
for char in range(ord("a"), ord("f") + 1):
    hex_to_nibble[char] = char - ord("a") + 10
for char in range(ord("A"), ord("F") + 1):
    hex_to_nibble[char] = char - ord("A") + 10


def remove_hyphens(ary):
    return ary[:, IDXS_WITHOUT_HYPHENS]


@numba.njit(nogil=True)
def parse_uuids(ary):
    num_bytes, remainder = divmod(len(ary[0]), 2)
    if remainder != 0:
        return None
    out = np.empty((len(ary), num_bytes), dtype=np.uint8)
    for idx in range(len(ary)):
        for num_byte in range(num_bytes):
            out[idx, num_byte] = (
                hex_to_nibble[ary[idx][2 * num_byte]] << 4
            ) + hex_to_nibble[ary[idx][2 * num_byte + 1]]
    return out


def parse_uuid(s):
    if isinstance(s, str):
        s = s.encode()
    ary = np.array([list(s)], dtype=np.uint8)
    return parse_uuids(ary)[0]
