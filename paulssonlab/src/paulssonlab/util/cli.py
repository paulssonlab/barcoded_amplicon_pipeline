def split_delimited_list(ctx, param, value, delimiter=","):
    if value:
        return [c.strip() for elem in value for c in elem.split(delimiter)]
    else:
        return []
