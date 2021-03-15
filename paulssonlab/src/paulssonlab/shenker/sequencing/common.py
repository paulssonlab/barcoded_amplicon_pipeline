import os
from pathlib import Path

REMOTE_HOST = "transfer.rc.hms.harvard.edu"

REMOTE_BASE_DIR = "/n/files/SysBio/PAULSSON LAB/Personal Folders"


def get_workdir():
    if "PAULSSONLAB_SCRATCH" in os.environ:
        subdir = Path.cwd().parent.name
        return Path(os.environ["PAULSSONLAB_SCRATCH"]) / subdir
    else:
        return "."
