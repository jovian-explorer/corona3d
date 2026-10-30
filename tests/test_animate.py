"""Panel selection and layouts of corona3d render, and corona3d animate: height laws, the shared GIF palette, the MP4,
and the frames of a synthetic reconstruction (no download)."""
import shutil
from datetime import datetime

import numpy as np
import pytest
from PIL import Image

from corona3d import animate as an
from corona3d.cli import main
from corona3d.whitelight import utc

OFF = [f"--no-{k}" for k in ("surface", "open", "closed", "streamer")]
P = {"time": "2013-10-25T15:00", "gcs": {"height": 2.0}, "gcs_wl": {"time": "2013-10-25T15:24", "height": 5.0}}


def test_height_law_line_speed_and_measured_points():
    t0 = utc("2013-10-25T15:00")
    f = an.height_law(P)                                           # the line through the two shells of the configuration
    assert f(utc("2013-10-25T15:12")) == pytest.approx(3.5)
    assert f(utc("2013-10-25T15:36")) == pytest.approx(6.5)        # a straight line also beyond the two fits
    f = an.height_law(P, fits=an.parse_height_points("2013-10-25T15:00=1.0,2013-10-25T16:00=4.0"))
    assert f(utc("2013-10-25T15:20")) == pytest.approx(2.0)
    f = an.height_law(P, speed_km_s=695.7)                         # 0.001 R_sun per second
    assert f(utc("2013-10-25T15:01")) == pytest.approx(2.06)
    f = an.height_law(P, speed_km_s=695.7, anchor=(t0, 3.0))
    assert f(utc("2013-10-25T15:00:10")) == pytest.approx(3.01)
    f = an.height_law(P, table=[(utc("2013-10-25T15:02"), 3.0), (t0, 2.0)])
    assert f(utc("2013-10-25T15:01")) == pytest.approx(2.5)
    with pytest.raises(ValueError, match="outside the measured"):
        f(utc("2013-10-25T15:03"))
    with pytest.raises(ValueError, match="one shell only"):
        an.height_law({"time": "2013-10-25T15:00", "gcs": {"height": 2.0}})


def test_heights_file_and_frame_times(tmp_path):
    (tmp_path / "h.csv").write_text("time,height\n# EUVI-B\n2013-10-25T15:00:00,1.5\n2013-10-25T15:01:00,2.0\n")
    assert an.read_heights(tmp_path / "h.csv")[1] == (utc("2013-10-25T15:01:00"), 2.0)
    ts = an.frame_times(utc("2013-10-25T14:59:00"), utc("2013-10-25T15:03:00"), 40)
    assert [f"{t:%H:%M:%S}" for t in ts] == ["14:59:00", "14:59:40", "15:00:20", "15:01:00", "15:01:40", "15:02:20", "15:03:00"]
    with pytest.raises(ValueError):
        an.frame_times(utc("2013-10-25T15:00"), utc("2013-10-25T14:00"), 40)


def _frames(reserved, n=3, size=(210, 120)):
    """Frames of coarse noise with every reserved colour in the first column."""
    rng = np.random.default_rng(1)
    out = []
    for i in range(n):
        f = rng.integers(0, 256, (*size, 3), dtype=np.uint8) // 8 * 8 + 100
        f[:, 1:40] = 200
        f[: len(reserved), 0] = reserved
        out.append(f)
    return out


def test_reserved_colours_survive_the_gif_palette_exactly(tmp_path):
    p = {"mesh_rings": {"color": "#c23ad8"}, "radio_colors": ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#8c564b", "#17becf"]}
    res = an.reserved_colours(p)
    assert res.shape[1] == 3 and len(res) <= 200
    assert {(0, 0, 0), (255, 255, 255), (0xc2, 0x3a, 0xd8), (0x9a, 0x6f, 0xc0)} <= {tuple(c) for c in res}
    frames = _frames(res)
    an.save_gif(frames, tmp_path / "a.gif", 2, res, scale=1)
    g = Image.open(tmp_path / "a.gif")
    assert g.n_frames == 3
    for i in range(g.n_frames):
        g.seek(i)
        px = np.asarray(g.convert("RGB")).reshape(-1, 3)
        found = {tuple(c) for c in px}
        assert {tuple(c) for c in res} <= found                       # every reserved colour is there, exactly
    assert len({tuple(c) for c in px}) <= 256


def test_gif_palette_is_the_same_in_every_frame_and_has_no_dithering(tmp_path):
    res = an.reserved_colours({"mesh_rings": {"color": "#c23ad8"}, "radio_colors": ["#1f77b4"]})
    ramp = np.linspace(0, 255, 100).astype(np.uint8)
    f = np.repeat(np.repeat(ramp[None, :, None], 3, 2), 50, 0)          # a grey gradient, 100 columns
    an.save_gif([f, f[:, ::-1]], tmp_path / "g.gif", 1, res, scale=1)
    g = Image.open(tmp_path / "g.gif")
    a = np.asarray(g.convert("RGB"))
    assert (a[:, :, 0] == a[:, :, 1]).all()                              # no dither noise in the colour channels
    assert (a == a[:1]).all()                                            # the gradient is constant down the columns
    g.seek(1)
    assert (np.asarray(g.convert("RGB")) == a[:, ::-1]).all()            # the same palette: a colour is a colour in every frame


def test_mp4_has_even_size_and_h264(tmp_path):
    try:
        an.ffmpeg_exe()
    except RuntimeError:
        pytest.skip("no ffmpeg")
    f = np.zeros((61, 81, 3), np.uint8)
    f[:, :, 0] = 200
    out = an.save_mp4([f, f, f], tmp_path / "a.mp4", 2)
    data = out.read_bytes()
    assert data[4:8] == b"ftyp" and b"avc1" in data[:200] and out.stat().st_size > 500


def test_panels_and_layouts_of_render(tmp_path):
    main(["--home", str(tmp_path), "init"])
    r = lambda *x: main(["--home", str(tmp_path), "render", *OFF, *x])
    r("--out", "f/a.png")
    w, h = Image.open(tmp_path / "f" / "a.png").size
    assert (w, h) == (2400, 735)                                         # the layout of Fig. 7
    r("--panel", "north", "--panel", "earth", "--panel", "300,20", "--out", "f/b.png")
    w2, h2 = Image.open(tmp_path / "f" / "b.png").size
    assert (w2, h2) == (2400, 735)                                       # one or two panels keep the figure width
    r("--panel", "earth,north,stereo", "--panel", "120,-10,45", "--out", "f/c.png")
    assert Image.open(tmp_path / "f" / "c.png").size == (3200, 735)      # four panels: wider
    r("--panel", "earth,north,stereo", "--layout", "portrait", "--out", "f/d.png")
    w3, h3 = Image.open(tmp_path / "f" / "d.png").size
    assert h3 > 3 * w3 / 2 and w3 == 1050
    r("--panel", "earth", "--panel", "40,5", "--layout", "separate", "--out", "f/e.png")
    for k in ("earth", "lon40_lat5"):
        assert Image.open(tmp_path / "f" / f"e_{k}.png").size[0] == 1050
    assert not (tmp_path / "f" / "e.png").exists()
    with pytest.raises(ValueError, match="panel"):
        r("--panel", "mars")
    with pytest.raises(ValueError, match="lon,lat"):
        r("--panel", "1,2,3,4")


def test_selected_view_keeps_the_framing_of_the_configuration(tmp_path):
    from corona3d.reconstruction import select_views
    main(["--home", str(tmp_path), "init"])
    from corona3d.paths import config
    p = config(tmp_path)
    (k, v), = select_views(p, ["stereo"])
    assert k == "stereo" and v["frame"] == p["views"][2]["frame"]
    assert [k for k, _ in select_views(p)] == ["earth", "north", "stereo"]


def test_animation_computes_the_model_once_and_writes_gif_and_mp4(synthetic_home, monkeypatch):
    import json
    from corona3d import geometry, mast_model
    cfg = json.loads((synthetic_home / "config.json").read_text())
    (synthetic_home / "config.json").write_text(json.dumps(dict(cfg, n_seeds=40)))
    n, lines = [], []
    real, real_lines = mast_model.load, geometry.FieldTracer.field_lines
    monkeypatch.setattr(mast_model, "load", lambda *a, **k: n.append(1) or real(*a, **k))
    monkeypatch.setattr(geometry.FieldTracer, "field_lines", lambda self, x: lines.append(1) or real_lines(self, x))
    args = ["--home", str(synthetic_home), "animate", "--no-surface", "--panel", "earth,north",
            "--start", "2013-10-25T15:00:00", "--end", "2013-10-25T15:01:00", "--step", "30", "--fps", "2",
            "--height-fits", "2013-10-25T15:00=1.5,2013-10-25T15:02=2.5", "--gif", "--out", "movie/m"]
    if shutil.which("ffmpeg") or _has_imageio():
        args += ["--mp4"]
    main(args)
    assert len(n) == 1 and len(lines) == 2                                # one model and two traces (the field, the streamer) for three frames
    g = Image.open(synthetic_home / "movie" / "m.gif")
    assert g.n_frames == 3 and g.size == (1200, 368)
    if "--mp4" in args:
        assert (synthetic_home / "movie" / "m.mp4").stat().st_size > 1000


def _has_imageio():
    try:
        import imageio_ffmpeg  # noqa: F401
        return True
    except ImportError:
        return False


def test_animation_frame_time_and_height_are_on_the_figure(synthetic_home):
    from corona3d.animate import animate
    files = animate("m", synthetic_home, "2013-10-25T15:00:00", "2013-10-25T15:00:40", 40, 1, True, False, 3.0,
                    None, ["earth"], "portrait", an.height_law(P), 1.0, surface=False, open=False, closed=False, streamer=False)
    assert [f.name for f in files] == ["m.gif"]
    g = Image.open(files[0])
    assert g.n_frames == 2
    a = np.asarray(g.convert("RGB")).astype(int)
    g.seek(1)
    b = np.asarray(g.convert("RGB")).astype(int)
    assert (a != b).any()                                                # the clock and the shell change


def test_animation_stops_below_the_solar_surface(synthetic_home):
    from corona3d.animate import animate
    with pytest.raises(ValueError, match="below the solar surface"):
        animate("m", synthetic_home, "2013-10-25T15:00:00", "2013-10-25T15:00:40", 40, 1, True, False, 3.0, None, None,
                "landscape", lambda t: 0.9, surface=False, open=False, closed=False, streamer=False)
