from itertools import zip_longest
from collections import Sequence, Mapping
from numbers import Number


def any_none(*args):
    return any(a is None for a in args)


def any_not_none(*args):
    return any(a is not None for a in args)


def sign(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


def format_sign(sgn):
    if sgn > 0:
        return "+"
    elif sgn < 0:
        return "-"
    else:
        return "|"  # TODO: is this the most comprehensible symbol?


def format_number(fmt, x):
    if isinstance(x, Number):
        return fmt.format(x)
    else:
        return "{}".format(x)


UNEVEN_GROUPS = object()

# TODO: what does this fancy version even do?? document!!!
# FROM: https://docs.python.org/3/library/itertools.html#recipes
def grouper(iterable, n, fillvalue=UNEVEN_GROUPS):
    """Collect data into fixed-length chunks or blocks."""
    args = [iter(iterable)] * n
    filter_filler = False
    if fillvalue == UNEVEN_GROUPS:
        fillvalue = object()
        filter_filler = True
    groups = zip_longest(*args, fillvalue=fillvalue)
    if filter_filler:
        groups = map(compose(tuple, partial(filter, lambda x: x != fillvalue)), groups)
    return groups


# TODO: simpler version of matriarch.util.get_one
def first(obj):
    return next(iter(obj))


def only(obj, msg="expecting a length-one iterable"):
    iterator = iter(obj)
    first_ = next(iterator)
    try:
        second_ = next(iterator)
    except:
        return first_
    else:
        raise ValueError(msg)


# FROM: https://stackoverflow.com/questions/42095393/python-map-a-function-over-recursive-iterables/42095505
# TODO: document!!!
def recursive_map(
    func,
    data,
    shortcircuit=(),
    ignore=(),
    keys=False,
    max_level=None,
    predicate=None,
    key_predicate=None,
):
    if max_level is not None and max_level is not False:
        if max_level == 0:
            return func(data)
        max_level -= 1
    apply = lambda x: recursive_map(
        func,
        x,
        shortcircuit=shortcircuit,
        ignore=ignore,
        keys=keys,
        max_level=max_level,
        predicate=predicate,
        key_predicate=key_predicate,
    )
    if isinstance(data, shortcircuit):
        return func(data)
    # to avoid an infinite recursion error
    # we need to hardcode that strs are ignored, because str[0] is a str, and hence a Sequence
    elif isinstance(data, str) or ignore is not True and isinstance(data, ignore):
        return data
    elif isinstance(data, Mapping):
        if keys:
            values = {apply(k): apply(v) for k, v in data.items()}
        else:
            values = {k: apply(v) for k, v in data.items()}
        if key_predicate is not None:
            values = {k: v for k, v in values.items() if key_predicate(k)}
        if predicate is not None:
            values = {k: v for k, v in values.items() if predicate(v)}
        return type(data)(values)
    elif isinstance(data, Sequence):
        values = [apply(v) for v in data]
        if predicate is not None:
            values = [v for v in values if predicate(v)]
        return type(data)(values)
    elif ignore is True or True in ignore:
        return data
    else:
        return func(data)
