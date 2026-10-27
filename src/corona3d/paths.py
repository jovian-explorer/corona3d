"""The working directory: --home, the environment variable CORONA3D_HOME, or the current folder
when it holds a config.json written by `corona3d init`."""
import json
import os
import sys
from importlib.resources import files
from pathlib import Path


def default_config():
    return json.loads(files("corona3d").joinpath("default_config.json").read_text(encoding="utf-8"))


def home(arg=None, must_exist=True):
    h = Path(arg or os.environ.get("CORONA3D_HOME") or Path.cwd()).resolve()
    if must_exist and not (h / "config.json").exists():
        sys.exit(f"{h} is not a corona3d working directory (no config.json); "
                 "run `corona3d init` there first, or pass --home")
    return h


def merge(base, new):
    """Copy of base with the keys of new; nested dictionaries are merged key by key."""
    out = dict(base)
    for k, v in (new or {}).items():
        out[k] = merge(out[k], v) if isinstance(v, dict) and isinstance(out.get(k), dict) else v
    return out


def config(h):
    """The default event, updated key by key with the working directory's config.json (if any)."""
    f = Path(h) / "config.json"
    user = json.loads(f.read_text(encoding="utf-8")) if f.exists() else {}
    return merge(default_config(), user)


def model_dir(h, cfg):
    return Path(h) / "data" / cfg["run"].replace("/", "_")


def wl_dir(h):
    return Path(h) / "data" / "wl"
