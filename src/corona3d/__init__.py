"""corona3d: three-dimensional reconstruction of a coronal mass ejection and the corona from public
MHD, coronagraph, and radio data."""
__version__ = "0.1.0"

from .layers import LAYERS, ON  # noqa: E402
from .paths import config as load_config  # noqa: E402

__all__ = ["LAYERS", "ON", "make_figure", "reconstruct", "newkirk_r", "load_config", "build", "export",
           "Reconstruction", "__version__"]

_LAZY = {"make_figure": "reconstruction", "reconstruct": "reconstruction", "newkirk_r": "reconstruction",
         "build": "api", "Reconstruction": "api", "export": "exporting"}


def __getattr__(name):
    """The functions that need PyVista, SunPy, or VTK are imported only when first used."""
    if name in _LAZY:
        from importlib import import_module
        mod = import_module(f".{_LAZY[name]}", __name__)
        return getattr(mod, name)
    raise AttributeError(name)
