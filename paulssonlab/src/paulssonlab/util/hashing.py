import collections.abc
import functools
import hashlib
import json

from frozendict import frozendict


# FROM: https://stackoverflow.com/a/53035002
@functools.singledispatch
def freeze(x):
    raise TypeError("Don't know how to freez {} objects".format(type(x).__name__))


@freeze.register(collections.abc.Hashable)
@freeze.register(str)  # otherwise Sequence takes precedence over Hashable
def _(x):
    return x


@freeze.register
def _(x: collections.abc.Sequence):
    return tuple(freeze(xi) for xi in x)


@freeze.register
def _(x: collections.abc.Set):
    return frozenset(freeze(xi) for xi in x)


@freeze.register
def _(x: collections.abc.Mapping):
    return frozendict({freeze(k): freeze(v) for k, v in x.items()})


# TODO: put snakenbake's ndarray hashing code here?


# SEE: https://stackoverflow.com/a/64572412
# SEE: https://stackoverflow.com/a/22003440
def hash_json(obj, digest_size=4, hash_func=hashlib.blake2b):
    s = json.dumps(obj, sort_keys=True, ensure_ascii=True, default=str)
    return hash_str(s, digest_size=digest_size, hash_func=hash_func)


def hash_str(s, digest_size=4, hash_func=hashlib.blake2b):
    h = hash_func(digest_size=digest_size)
    h.update(s.encode())
    return h.hexdigest()
