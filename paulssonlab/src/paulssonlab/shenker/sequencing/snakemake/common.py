import os
from functools import lru_cache
from pathlib import Path

DIRS_TO_LINK = ["data", "references", "output", "logs"]
REMOTE_HOST = "transfer.rc.hms.harvard.edu"
REMOTE_BASE_DIR = "/n/files/SysBio/PAULSSON LAB/Personal Folders"
SCRATCH_SUBDIR = "sequencing"
CLONING_CONFIG_DIR = Path(__file__).parent.parent / "cloning"


def get_workdir(make_links=False, dirs_to_link=DIRS_TO_LINK):
    workdir = _get_workdir()
    if workdir != "." and make_links:
        cwd = Path.cwd()
        for name in dirs_to_link:
            src = cwd / name
            if not src.is_symlink():
                src.symlink_to(workdir / name, target_is_directory=True)
    return workdir


def _get_workdir():
    if "PAULSSONLAB_SCRATCH" in os.environ:
        project_name = Path.cwd().name
        return Path(os.environ["PAULSSONLAB_SCRATCH"]) / SCRATCH_SUBDIR / project_name
    else:
        return "."


# TODO: Instead of reloading the full registry, this should be replaced
#       with a lightweight query to google drive for any files named "pLIBxx.gb"
#       and no calls to gsheets. This would also allow using it for non-registry-compliant
#       directories of plasmid maps
# TODO: this caching doesn't work because snakemake forks for each job
@lru_cache
def get_registry(config_dir=CLONING_CONFIG_DIR):
    import pygsheets
    import toml

    import paulssonlab.cloning.registry as registry

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
