"""corona3d: physics-based, reproducible three-dimensional reconstruction of the solar corona,
coronal mass ejections, and radio sources for any event from public MHD, coronagraph, and radio
data, following the method described in Morosan et al. (2022)."""
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
