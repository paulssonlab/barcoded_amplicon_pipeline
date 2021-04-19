import collections.abc
from frozendict import frozendict

# FROM: https://stackoverflow.com/a/53035002
@functools.singledispatch
def hashify(x):
    raise TypeError(
        "Don't know how to make {} objects hashable".format(type(x).__name__)
    )


@hashify.register
def _(x: collections.abc.Hashable):
    return x


@hashify.register
def _(x: collections.abc.Sequence):
    return tuple(hashify(xi) for xi in x)


@hashify.register
def _(x: collections.abc.Set):
    return frozenset(hashify(xi) for xi in x)


@hashify.register
def _(x: collections.abc.Mapping):
    return frozendict({k: hashify(v) for k, v in x.items()})
