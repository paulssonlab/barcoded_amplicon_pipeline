def split_delimited_list(ctx, param, value, delimiter=","):
    if value:
        return [c.strip() for elem in value for c in elem.split(delimiter)]
    else:
        return []


def parse_kv(ctx, param, value):
    d = {}
    for k, v in value:
        try:
            v = int(v)
        except:
            try:
                v = float(v)
            except:
                pass
        d[k] = v
    return d
