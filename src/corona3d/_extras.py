"""Optional dependencies: each heavy package is imported where it is first needed, and a missing one
stops with the name of the extra that installs it (`pip install "corona3d[render]"`)."""
from importlib import import_module

EXTRAS = {                                   # extra: the packages (import name, pip name) it installs
    "render": [("pyvista", "pyvista"), ("PyThea", "PyThea")],
    "gcs": [("PyThea", "PyThea"), ("skimage", "scikit-image")],
    "radio": [("h5py", "h5py")],
}


class MissingExtra(ImportError):
    """An optional package is not installed; the message names the extra."""


def need(extra, module, what):
    """Import `module` for `what`; on failure raise MissingExtra naming the extra to install."""
    try:
        return import_module(module)
    except ImportError as e:
        raise MissingExtra(f"{what} needs the package {module}, which is not installed. "
                           f"Install it with: pip install \"corona3d[{extra}]\"  (or corona3d[all])") from e


def status():
    """{extra: {package: version or None}} for the optional extras."""
    from importlib import metadata
    out = {}
    for extra, pk in EXTRAS.items():
        out[extra] = {}
        for mod, pip in pk:
            try:
                import_module(mod)
                out[extra][pip] = metadata.version(pip)
            except Exception:
                out[extra][pip] = None
    return out
