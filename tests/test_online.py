"""End to end from public data (CORONA3D_ONLINE=1; about 150 MB and 5-10 minutes): download the
MAS/MAST run and the white-light images of 2013-10-25 and draw both figures of the README."""
import os

import numpy as np
import pytest

from corona3d.cli import main

pytestmark = pytest.mark.skipif(os.environ.get("CORONA3D_ONLINE") != "1",
                                reason="set CORONA3D_ONLINE=1 to download public data")


def test_fetch_and_render_from_public_data(tmp_path):
    main(["--home", str(tmp_path), "init"])
    main(["--home", str(tmp_path), "fetch"])
    names = {p.name for p in (tmp_path / "data" / "wl").iterdir()}
    assert any(n.endswith("d4c2B.fts") for n in names) and any(n.endswith("s4c1B.fts") for n in names)
    main(["--home", str(tmp_path), "render", "--out", "figures/corona3d.png"])
    main(["--home", str(tmp_path), "render", "--c2", "--cor1", "--cor2", "--out", "figures/wl.png"])
    from PIL import Image
    a = np.asarray(Image.open(tmp_path / "figures" / "view_earth.png").convert("L"), float)
    assert a.std() > 20                                    # a drawn panel, not a blank one
