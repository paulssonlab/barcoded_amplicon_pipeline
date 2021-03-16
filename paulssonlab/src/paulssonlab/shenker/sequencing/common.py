import os
from pathlib import Path
from functools import lru_cache

REMOTE_HOST = "transfer.rc.hms.harvard.edu"
REMOTE_BASE_DIR = "/n/files/SysBio/PAULSSON LAB/Personal Folders"
SCRATCH_SUBDIR = "sequencing"
CLONING_CONFIG_DIR = Path(__file__).parent.parent / "cloning"


def get_workdir():
    if "PAULSSONLAB_SCRATCH" in os.environ:
        project_name = Path.cwd().name
        return Path(os.environ["PAULSSONLAB_SCRATCH"]) / SCRATCH_SUBDIR / project_name
    else:
        return "."


@lru_cache
def get_registry(config_dir=CLONING_CONFIG_DIR):
    import toml
    import pygsheets
    import paulssonlab.cloning.registry as registry

    config_dir = Path(config_dir)
    config = toml.load(config_dir / "config.toml")
    gc = pygsheets.authorize(service_account_file=config_dir / "credentials.json")
    reg = registry.Registry(gc, config["registry"]["folder"])
    return reg


def get_registry_seq(name, filename, config_dir=CLONING_CONFIG_DIR):
    reg = get_registry(config_dir)
    entry = reg.get(name)
    if "_seq" not in entry:
        raise ValueError(f"did not find sequence for {name}")
    seq = entry["_seq"]
    # filename = output_dir / f"{name}.gb"
    with open(filename, "w") as f:
        f.write(seq.format("gb"))
