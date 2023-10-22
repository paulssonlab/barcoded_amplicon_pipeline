def split_delimited_list(ctx, param, value, delimiter=","):
    return [c.strip() for elem in value for c in elem.split(delimiter)]
