"""The working directory: --home, the environment variable CORONA3D_HOME, or the current folder
when it holds a config.json written by `corona3d init`."""
import json
import os
import sys
from importlib.resources import files
from pathlib import Path


def default_config():
    return json.loads(files("corona3d").joinpath("default_config.json").read_text(encoding="utf-8"))


def events():
    """Names of the event presets: the default event and the files in corona3d/events/."""
    return ["2013-10-25"] + sorted(f.name[:-5] for f in files("corona3d").joinpath("events").iterdir()
                                   if f.name.endswith(".json"))


def event(name):
    return json.loads(files("corona3d").joinpath("events", f"{name}.json").read_text(encoding="utf-8"))


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


def legacy(user):
    """Keys of a config.json written before the STEREO-A support and the OVRO-LWA layer (version
    0.1.0 and the commits after it) in the current names: `stereo_b` becomes `stereo` with
    `stereo_sc` "B", the views from `stereo_a` and `stereo_b` become views from `stereo`, and
    `radio_nrh` becomes an entry of `radio_images`. Returns a new dictionary."""
    user = dict(user)
    if "stereo_b" in user:
        sb = user.pop("stereo_b")
        if "stereo" not in user:
            user["stereo"] = sb
            user.setdefault("stereo_sc", "B")
    views = []
    for v in user.get("views", []):
        if v.get("from") in ("stereo_a", "stereo_b"):
            user.setdefault("stereo_sc", v["from"][-1].upper())
            v = dict(v, **{"from": "stereo"})
        views.append(v)
    if "views" in user:
        user["views"] = views
    if "radio_nrh" in user:
        nrh = dict(user.pop("radio_nrh"))
        if "colormaps" in nrh:
            user.setdefault("radio_colormaps", nrh.pop("colormaps"))
        if "radio_images" not in user:
            user["radio_images"] = [dict(nrh, instrument="nrh")]
    return user


def config(h):
    """The default event, updated key by key with the working directory's config.json (if any);
    keys of older versions are converted (see `legacy`)."""
    f = Path(h) / "config.json"
    user = legacy(json.loads(f.read_text(encoding="utf-8"))) if f.exists() else {}
    return merge(default_config(), user)


def model_dir(h, cfg):
    return Path(h) / "data" / cfg["run"].replace("/", "_")


def wl_dir(h):
    return Path(h) / "data" / "wl"
