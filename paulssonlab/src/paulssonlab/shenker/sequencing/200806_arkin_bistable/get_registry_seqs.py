import time
import sys
import os
import re
import toml
from pathlib import Path
import pygsheets

sys.path.append(str(Path(os.environ["src"]).parent))
import paulssonlab.cloning.registry as registry

RETRIES = 3


def get_registry(config_dir):
    config_dir = Path(config_dir)
    config = toml.load(config_dir / "config.toml")
    gc = pygsheets.authorize(service_account_file=config_dir / "credentials.json")
    reg = registry.Registry(gc, config["registry"]["folder"])
    return reg


def get_registry_seq(reg, name, filename):
    entry = reg.get(name)
    if "_seq" not in entry:
        raise ValueError(f"did not find sequence for {name}")
    seq = entry["_seq"]
    with open(filename, "w") as f:
        f.write(seq.format("gb"))


ids = re.split(r"\s*,\s*", sys.stdin.read().rstrip())

print(f"Fetching registry sequences: {ids}")
print()

reg = None
for retry in range(1, RETRIES + 1):
    while ids:
        try:
            if reg is None:
                reg = get_registry(sys.argv[1])
            id = ids[-1]
            print(f"Fetching registry sequence '{id}'... ", end=None)
            get_registry_seq(reg, id, f"{id}.gb")
            ids.pop()
            print("done.")
        except Exception as e:
            print()
            print(f"Got exception: {e}")
            if retry == RETRIES:
                print("Maximum retries reached, aborting.")
                sys.exit(1)
            print(f"Retrying (attempt {retry})...")
            time.sleep(0.3 * retry)
        else:
            break
