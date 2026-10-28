"""Opacity of every layer: a number, a list per contour level or per blob, or a transfer function
of the pixel value for the surface and the coronagraph images."""
import json

import numpy as np
import pytest

from corona3d.cli import main, parse_opacity
from corona3d.reconstruction import check_opacity, opacities


def test_opacities_are_spread_over_the_items():
    assert opacities(0.4, 3).tolist() == [0.4, 0.4, 0.4]
    assert opacities([0.2, 0.5, 0.8], 3).tolist() == [0.2, 0.5, 0.8]
    assert opacities([0.2, 1.0], 5) == pytest.approx([0.2, 0.4, 0.6, 0.8, 1.0])


def test_wrong_opacities_are_refused():
    check_opacity({"radio_contours": [0.2, 1.0], "surface": "sigmoid", "wl": [0.0, 0.8], "gcs": 0.5})
    for bad in ({"gcs": 1.5}, {"gcs": [0.2, 0.4]}, {"radio": "linear"}, {"wl": []}, {"open": True}):
        with pytest.raises(ValueError, match="opacity of"):
            check_opacity(bad)


def test_opacity_option_is_read_as_number_list_or_name():
    assert parse_opacity("gcs=0.5") == ("gcs", 0.5)
    assert parse_opacity("radio_contours=0.2,0.6,1") == ("radio_contours", [0.2, 0.6, 1.0])
    assert parse_opacity("wl=linear") == ("wl", "linear")
    with pytest.raises(SystemExit):
        parse_opacity("flares=0.5")


def test_render_with_semitransparent_layers(tmp_path):
    """Sphere and shell half transparent, and radio blobs fading in along the chain, give a lighter
    panel than the opaque layers."""
    from PIL import Image
    main(["--home", str(tmp_path), "init"])
    off = [f"--no-{k}" for k in ("surface", "open", "closed", "streamer")]
    src = json.dumps([{"freq": 173, "x": [-1.19, -1.09, -0.99], "y": [-0.15, -0.08, 0.0]}])
    common = ["--home", str(tmp_path), "render", *off, "--radio", "--radio-sources", src]
    main(common + ["--out", "figures/opaque.png"])
    opaque = np.asarray(Image.open(tmp_path / "figures" / "view_a.png").convert("L"), float).mean()
    main(common + ["--opacity", "sphere=0.4", "--opacity", "gcs=0.5", "--opacity", "radio=0.2,1",
                   "--out", "figures/semi.png"])
    semi = np.asarray(Image.open(tmp_path / "figures" / "view_a.png").convert("L"), float).mean()
    assert semi > opaque + 1
    with pytest.raises(ValueError, match="opacity of gcs"):
        main(common + ["--opacity", "gcs=2", "--out", "figures/bad.png"])
