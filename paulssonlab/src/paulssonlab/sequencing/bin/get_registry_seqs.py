#!/usr/bin/env python
import time
import sys
import os
import re
import json
import hashlib
from frozendict import frozendict
import toml
from pathlib import Path
import pygsheets

sys.path.append(str(Path(os.environ["src"]).parent))
import paulssonlab.cloning.registry as registry
from paulssonlab.cloning.commands.parser import expr_list_parser, unparse_expr
from paulssonlab.cloning.commands.semantics import eval_expr
from paulssonlab.util.hashing import hash_json, freeze

RETRIES = 3
RETRY_DELAY = 0.3


def get_registry(config_dir):
    config_dir = Path(config_dir)
    config = toml.load(config_dir / "config.toml")
    gc = pygsheets.authorize(service_account_file=config_dir / "credentials.json")
    reg = registry.Registry(gc, config["registry"]["folder"])
    return reg


def get_registry_seq(reg, name, filename):
    res = eval_expr(name, reg.get)
    if "_seq" not in res:
        raise ValueError(f"did not find sequence for {name}")
    seq = res["_seq"]
    with open(filename, "w") as f:
        f.write(seq.format("gb"))


def parse_references(refs):
    if not isinstance(refs, list):
        refs = [refs]
    exprs = sum([list(expr_list_parser.parse(s)) for s in refs], [])
    exprs = [freeze(expr) for expr in exprs]
    return exprs


def filename_for_expr(expr):
    if isinstance(expr, (list, tuple)):
        raise ValueError("expecting a single expression, not a list")
    elif not isinstance(expr, (dict, frozendict)) or "_type" not in expr:
        return hash_json(expr)
    type_ = expr["_type"]
    if type_ == "name":
        return f"{expr['name']}"
    elif type_ == "pcr":
        return f"pcr_{expr['template']['name']}_with_{expr['primer1']['name']}_and_{expr['primer2']['name']}"
    elif type_ == "digest":
        return f"digest_{expr['input']['name']}_with_{expr['enzyme']['name']}"
    elif type_ == "anneal":
        return f"anneal_{expr['strand1']['name']}_with_{expr['strand2']['name']}"
    else:
        raise ValueError(f"unknown expression type: {type_}")


config_dir = Path(sys.argv[1])
if not config_dir.is_dir():
    raise ValueError(f"expecting config directory: {config_dir}")

output_dir = Path(sys.argv[2])
# ensure output directory exists
output_dir.mkdir(parents=True, exist_ok=True)

samples = json.load(sys.stdin)

keys_to_exprs = {}
exprs_to_get = set()

for key, references in samples.items():
    exprs = parse_references(references)
    keys_to_exprs[key] = exprs
    exprs_to_get.update(exprs)

names = {e: unparse_expr(e) for e in exprs_to_get}
filenames = {e: f"{filename_for_expr(e)}.gb" for e in exprs_to_get}

exprs_to_get = list(exprs_to_get)

reg = None
retry = 1
while exprs_to_get:
    try:
        expr = exprs_to_get[-1]
        name = names[expr]
        filename = filenames[expr]
        output_path = output_dir / filename
        if output_path.exists():
            exprs_to_get.pop()
            continue
        print(f"Fetching registry sequence '{name}'... ", end=None, file=sys.stderr)
        if reg is None:
            reg = get_registry(config_dir)
        get_registry_seq(reg, name, output_path)
        exprs_to_get.pop()
        print("done.", file=sys.stderr)
    except Exception as e:
        print(f"Got exception: {e}", file=sys.stderr)
        if retry == RETRIES:
            print("Maximum retries reached, aborting.", file=sys.stderr)
            raise
            # sys.exit(1)
        retry += 1
        print(f"Retrying (attempt {retry})...", file=sys.stderr)
        time.sleep(RETRY_DELAY * retry)

output = {}
for key, exprs in keys_to_exprs.items():
    output[key] = {
        "reference_names": [names[e] for e in exprs],
        "references": [filenames[e] for e in exprs],
    }

json.dump(output, sys.stdout)
print()  # add newline
