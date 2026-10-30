"""The working directory (home): --home, the environment variable CORONA3D_HOME, the current folder
when it holds a config.json written by `corona3d init`, or else the platform cache directory
(platformdirs), e.g. ~/.cache/corona3d on Linux and %LOCALAPPDATA%/corona3d/Cache on Windows."""
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


def default_home():
    """The platform cache directory for corona3d (used when no other home is chosen)."""
    try:
        from platformdirs import user_cache_dir
        return Path(user_cache_dir("corona3d", appauthor=False))
    except ImportError:
        return Path.home() / ".corona3d"


def home_source(arg=None):
    """(path, where it came from) of the working directory chosen by home()."""
    if arg:
        return Path(arg).resolve(), "--home"
    if os.environ.get("CORONA3D_HOME"):
        return Path(os.environ["CORONA3D_HOME"]).resolve(), "CORONA3D_HOME"
    if (Path.cwd() / "config.json").exists():
        return Path.cwd().resolve(), "current folder (config.json)"
    return default_home().resolve(), "platform cache directory"


def home(arg=None, must_exist=True):
    h = home_source(arg)[0]
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
    return colormap_rule(user)


def colormap_rule(user):
    """Colour maps chosen by the user (not those that `init` and `new` write) replace the default radio colours:
    `radio_colors` becomes None unless it is given. Returns a new dictionary."""
    user = dict(user)
    old = ["spring", "winter", "autumn", "summer"]  # the colour maps that configs written by `init` and `new` carry
    cm_ = user.get("radio_colormaps")
    if cm_ and list(cm_) != old[:len(cm_)]:
        user.setdefault("radio_colors", None)
    return user


def config(h):
    """The default event, updated key by key with the working directory's config.json (if any);
    keys of older versions are converted (see `legacy`), and `density_scale` "pb" is replaced by the scale of
    `pb_fit` (corona3d.pb.resolve)."""
    f = Path(h) / "config.json"
    user = legacy(json.loads(f.read_text(encoding="utf-8"))) if f.exists() else {}
    from .pb import resolve
    cfg = merge(default_config(), user)
    try:
        return resolve(cfg)                               # density_scale "pb" becomes the scale fitted by corona3d pb
    except ValueError:                                    # no pb_fit yet: every command must still start (corona3d pb makes it)
        import warnings
        warnings.warn('config.json has density_scale "pb" but no pb_fit; the default density scale is used until '
                      '`corona3d pb` has been run (reconstruct, shock, sources and spectra stop until then)', stacklevel=2)
        return dict(cfg, density_scale=default_config()["density_scale"], density_scale_pending_pb=True)


def model_dir(h, cfg):
    return Path(h) / "data" / cfg["run"].replace("/", "_")


def wl_dir(h):
    return Path(h) / "data" / "wl"
