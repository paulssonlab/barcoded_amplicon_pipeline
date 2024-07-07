import string
from collections.abc import Iterable, KeysView, Mapping, Sequence, ValuesView
from functools import partial
from itertools import chain, zip_longest
from numbers import Number

from cytoolz import compose, dissoc


def extract_keys(d, keys):
    return {k: d[k] for k in d.keys() & set(keys)}


def pop_keys(d, keys):
    return extract_keys(d, keys), dissoc(d, *keys)


def any_none(*args):
    return any(a is None for a in args)


def any_not_none(*args):
    return any(a is not None for a in args)


def parse_number(s):
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            pass
    return s


# FROM: https://stackoverflow.com/a/53718454
def str_placeholders(format_string):
    return [
        parse_number(x[1]) if x[1] else idx
        for idx, x in enumerate(string.Formatter().parse(format_string))
    ]


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
# TODO: steal better version from recipes, or find in a library
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
    try:
        first_ = next(iterator)
    except:
        raise ValueError(msg)
    try:
        second_ = next(iterator)
    except:
        return first_
    else:
        raise ValueError(msg)


# FROM: https://stackoverflow.com/questions/42095393/python-map-a-function-over-recursive-iterables/42095505
# TODO: document!!!
def map_recursive(
    func,
    data,
    shortcircuit=(),
    ignore=(str,),
    keys=False,
    max_level=None,
    predicate=None,
    key_predicate=None,
):
    if max_level is not None and max_level is not False:
        if max_level == 0:
            return func(data)
        max_level -= 1

    def apply(x):
        return map_recursive(
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
    # we need ignore str by default, because str[0] is a str, and hence a Sequence
    elif isinstance(data, ignore):
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
    else:
        return func(data)


def iter_recursive(
    obj,
    ignore=(str,),
    recurse=(tuple, list, dict, KeysView, ValuesView),
):
    # to avoid an infinite recursion error
    # we need ignore str by default, because str[0] is a str, and hence a Sequence
    should_recurse = not recurse or isinstance(obj, recurse)
    if isinstance(obj, ignore):
        return iter((obj,))
    elif should_recurse and isinstance(obj, Mapping):
        return chain(iter_recursive(obj.keys()), iter_recursive(obj.values()))
    elif should_recurse and isinstance(obj, Iterable):
        return chain.from_iterable(iter_recursive(x) for x in obj)
    else:
        return iter((obj,))


# FROM: https://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
def flatten_dict(
    d, parent_key=None, sep=None, predicate=None, lookahead=None, concatenate_keys=False
):
    if parent_key is None:
        if sep is not None:
            parent_key = ""
        else:
            parent_key = ()
    items = []
    for k, v in d.items():
        if sep is not None:
            new_key = parent_key + sep + str(k) if parent_key else str(k)
        else:
            if concatenate_keys and isinstance(k, tuple):
                new_key = parent_key + k
            else:
                new_key = parent_key + (k,)
        if isinstance(v, Mapping) and (lookahead is None or lookahead(v)):
            items.extend(
                flatten_dict(
                    v,
                    parent_key=new_key,
                    sep=sep,
                    predicate=predicate,
                    lookahead=lookahead,
                    concatenate_keys=concatenate_keys,
                ).items()
            )
        else:
            if predicate is None or predicate(k, v):
                items.append((new_key, v))
    return d.__class__(items)


class ItemProxy(object):
    def __init__(self, obj, name):
        self._obj = obj
        self._name = name

    def __contains__(self, key):
        return getattr(self._obj, f"_{self._name}_contains")(key)

    def __delitem__(self, key):
        return getattr(self._obj, f"_{self._name}_delitem")(key)

    def __iter__(self):
        return getattr(self._obj, f"_{self._name}_iter")()

    def items(self):
        return getattr(self._obj, f"_{self._name}_items")()

    def __getitem__(self, key):
        return getattr(self._obj, f"_{self._name}_getitem")(key)

    def __setitem__(self, key, value):
        return getattr(self._obj, f"_{self._name}_setitem")(key, value)

    def update(self, other):
        setitem = getattr(self._obj, f"_{self._name}_setitem")
        for key, value in other.items():
            setitem(key, value)
