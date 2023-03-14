#!/usr/bin/env python
import time
import sys
import os
import shutil
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
from paulssonlab.util.hashing import hash_json, hash_str, freeze

RETRIES = 3
RETRY_DELAY = 0.3


def get_registry(config_dir):
    config_dir = Path(config_dir)
    config = toml.load(config_dir / "config.toml")
    gc = pygsheets.authorize(service_account_file=config_dir / "credentials.json")
    reg = registry.Registry(gc, config["registry"]["folder"])
    return reg


def get_registry_seq(reg, expr):
    return eval_expr(expr, reg.get)


def download_registry_seq(reg, expr, filename_base):
    res = get_registry_seq(reg, expr)
    if "_seq" not in res:
        raise ValueError(f"did not find sequence for {expr}")
    seq = res["_seq"]
    if isinstance(seq, (Seq, SeqRecord, DsSeqRecord)):
        filename = f"{filename_base}.gb"
        output = seq.format("gb")
    else:
        raise ValueError(f"unexpected sequence object {seq!r}")
    with open(filename, "w") as f:
        f.write(output)
    return filename


def get_expr_seq():
    pass


def parse_references(refs):
    if not isinstance(refs, list):
        refs = [refs]
    exprs = sum([list(expr_list_parser.parse(s)) for s in refs], [])
    exprs = [freeze(expr) for expr in exprs]
    return exprs


def name_for_expr(expr):
    """If expression has only name-type expressions as inputs (so no recursive
    expression tree), then use a human-readable name for the filename.

    Otherwise uses the hash of the expression tree.
    """
    if isinstance(expr, (list, tuple)):
        raise ValueError("expecting a single expression, not a list")
    elif isinstance(expr, (dict, frozendict)) and "_type" in expr:
        type_ = expr["_type"]
        if type_ == "name":
            return f"{expr['name']}"
        else:
            input_keys = sorted([k for k in expr.keys() if not k.startswith("_")])
            if all(
                expr[k].get("_type") == "name" and expr[k]["name"] for k in input_keys
            ):
                if type_ == "pcr":
                    return f"pcr_{expr['template']['name']}_with_{expr['primer1']['name']}_and_{expr['primer2']['name']}"
                elif type_ == "digest":
                    return (
                        f"digest_{expr['input']['name']}_with_{expr['enzyme']['name']}"
                    )
                elif type_ == "anneal":
                    return f"anneal_{expr['strand1']['name']}_with_{expr['strand2']['name']}"
                else:
                    raise ValueError(f"unknown expression type: {type_}")
    return hash_json(expr)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise ValueError("expecting two arguments: config_dir and output_dir")
    config_dir = Path(sys.argv[1])
    if not config_dir.is_dir():
        raise ValueError(f"expecting config directory: {config_dir}")
    output_dir = Path(sys.argv[2])
    # ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    samples = json.load(sys.stdin)
    output = {}
    reg = None
    retry = 1
    for key, references in samples.items():
        output[key] = {"references": [], "reference_exprs": []}
        if re.match(r"^\s*(?:/|\.|\.\.)", references):
            # SEE: https://stackoverflow.com/a/21107911
            for source_path in re.split(r"(?<!\\),", references):
                source_path = Path(source_path.strip()).resolve()
                dest_filename = f"{source_path.stem}_{hash_str(str(source_path))}{source_path.suffix}"
                dest_path = output_dir / dest_filename
                # print("!", dest_path, "|", output_dir, ">", sys.argv[2])
                if not dest_path.exists():
                    shutil.copyfile(source_path, dest_path)
                output[key]["references"].append(str(dest_path))
                output[key]["reference_exprs"].append(str(source_path))
        else:
            continue
            exprs = parse_references(references)
            for expr in exprs:
                while True:
                    try:
                        expr = exprs_to_get[-1]
                        expr_str = expr_strs[expr]
                        filename = filenames[expr]
                        output_path = output_dir / filename
                        if output_path.exists():
                            exprs_to_get.pop()
                            continue
                        print(
                            f"Fetching registry sequence '{expr_str}'... ",
                            end=None,
                            file=sys.stderr,
                        )
                        if reg is None:
                            reg = get_registry(config_dir)
                        get_registry_seq(reg, name, output_path)
                        exprs_to_get.pop()
                        dest_filename = 0
                        output[key].append(dest_filename)
                        print("done.", file=sys.stderr)
                        break
                    except Exception as e:
                        print(f"Got exception: {e}", file=sys.stderr)
                        if retry == RETRIES:
                            print("Maximum retries reached, aborting.", file=sys.stderr)
                            raise
                            # sys.exit(1)
                        retry += 1
                        print(f"Retrying (attempt {retry})...", file=sys.stderr)
                        time.sleep(RETRY_DELAY * retry)

    ########
    # keys_to_exprs = {}
    # exprs_to_get = set()
    # for key, references in samples.items():
    #     exprs = parse_references(references)
    #     keys_to_exprs[key] = exprs
    #     exprs_to_get.update(exprs)
    # expr_strs = {e: unparse_expr(e) for e in exprs_to_get}
    # names = {e: f"{name_for_expr(e)}.gb" for e in exprs_to_get}
    # exprs_to_get = list(exprs_to_get)
    #####
    json.dump(output, sys.stdout)
    print()  # add newline
