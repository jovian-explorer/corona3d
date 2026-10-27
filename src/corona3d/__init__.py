"""corona3d: three-dimensional reconstruction of a coronal mass ejection and the corona from public
MHD, coronagraph, and radio data."""
__version__ = "0.1.0"

from .layers import LAYERS, ON  # noqa: E402

__all__ = ["LAYERS", "ON", "make_figure", "reconstruct", "newkirk_r", "__version__"]


def __getattr__(name):
    """make_figure, reconstruct, and newkirk_r load PyVista and VTK only when first used."""
    if name in ("make_figure", "reconstruct", "newkirk_r"):
        from . import reconstruction
        return getattr(reconstruction, name)
    raise AttributeError(name)
