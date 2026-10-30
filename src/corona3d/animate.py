"""Animations of the reconstruction (`corona3d animate`): one frame per time, as a GIF and an MP4.

Every frame is the figure of `corona3d render` at the frame time: the shell has the apex height that a height
law gives for that time, and the radio sources accumulate up to it. The model surface and the field lines are
computed once, the cameras, the colour limits, the panel scale, and the radio time range are fixed for the whole
animation, and the frame time is written on the figure. A GIF holds 256 colours, so one palette is built for all
frames (`shared_palette`).
"""
import csv
import shutil
import subprocess
from datetime import timedelta
from pathlib import Path

import numpy as np

R_SUN_KM = 695700.0


def parse_height_points(text):
    """'2013-10-25T15:00=1.931,2013-10-25T15:24=5.538' -> [(datetime, height [R_sun])]."""
    from .whitelight import utc
    try:
        return [(utc(t), float(h)) for t, h in (x.rsplit("=", 1) for x in text.split(","))]
    except ValueError:
        raise ValueError(f"{text!r}: give TIME=HEIGHT pairs, e.g. 2013-10-25T15:00=1.931,2013-10-25T15:24=5.538")


def read_heights(path):
    """Measured heights from a CSV file of two columns, time (ISO) and apex height [R_sun]; a header line and
    lines starting with # are skipped."""
    from .whitelight import utc
    pts = []
    for row in csv.reader(l for l in Path(path).read_text(encoding="utf-8").splitlines() if l.strip() and l[0] != "#"):
        try:
            pts.append((utc(row[0].strip()), float(row[1])))
        except (ValueError, IndexError):
            if pts:
                raise ValueError(f"{path}: cannot read the line {row}")        # only the header may be unreadable
    if len(pts) < 2:
        raise ValueError(f"{path}: at least two (time, height) lines are needed")
    return sorted(pts)


def height_law(p, fits=None, speed_km_s=None, anchor=None, table=None):
    """Function from a time (datetime) to the apex height of the shell [R_sun].

    table: measured (time, height) points (EUV or Solar Orbiter), interpolated linearly, valid only between the
    first and the last of them. speed_km_s: constant speed from anchor, one (time, height) point (default: the
    shell of the configuration at its time). fits: two (time, height) points, for example the two GCS fits, with
    a straight line through them (also beyond them). Default: the line through the shell `gcs` at the time of
    the reconstruction and the shell `gcs_wl` at its own time, the two fits of the configuration."""
    from .whitelight import utc
    if table:
        t, h = zip(*sorted(table))
        s = np.array([(x - t[0]).total_seconds() for x in t])
        def law(time):
            x = (time - t[0]).total_seconds()
            if not s[0] <= x <= s[-1]:
                raise ValueError(f"{time:%Y-%m-%d %H:%M:%S} UT is outside the measured heights ({t[0]:%H:%M:%S} to "
                                 f"{t[-1]:%H:%M:%S} UT); change --start and --end, or use a line or a speed")
            return float(np.interp(x, s, h))
        return law
    if speed_km_s is not None:
        t0, h0 = anchor or (utc((p.get("gcs_fit") or {}).get("t_out") or p["time"]), p["gcs"]["height"])
        return lambda time: h0 + speed_km_s * (time - t0).total_seconds() / R_SUN_KM
    if fits is None:
        g, w = p["gcs"], p.get("gcs_wl")
        if not w or "time" not in w:
            raise ValueError("the configuration has one shell only: give --height-fits, --speed, or --heights")
        fits = [(utc((p.get("gcs_fit") or {}).get("t_out") or p["time"]), g["height"]), (utc(w["time"]), w["height"])]
    (t0, h0), (t1, h1) = fits
    if t0 == t1:
        raise ValueError("the two fits are at the same time")
    return lambda time: h0 + (h1 - h0) * (time - t0).total_seconds() / (t1 - t0).total_seconds()


def frame_times(start, end, step_s):
    """Times from start to end (inclusive) in steps of step_s seconds."""
    if step_s <= 0 or end < start:
        raise ValueError("give --step above zero and --end not before --start")
    n = int((end - start).total_seconds() // step_s)
    return [start + timedelta(seconds=step_s * i) for i in range(n + 1)]


def reserved_colours(p, n_max=200):
    """The colours that every frame needs exactly, as uint8 (n x 3): the v_A colour map (Spectral_r, 40 steps)
    and two copies in the shade of the surface (0.75 and 0.5 of the brightness), the light-to-dark ramp of every
    radio colour (12 steps), the shell colour, the purple of the streamer lines, black, and white."""
    from matplotlib import colormaps, colors
    from .reconstruction import chain_colours
    sp = colormaps["Spectral_r"](np.linspace(0, 1, 40))[:, :3]
    single = np.array([colors.to_rgb(c) for c in (p["mesh_rings"]["color"], "#9a6fc0", "#000000", "#ffffff")])
    base = p.get("radio_colors") or []
    n_ramp = min(12, max(2, (n_max - 120 - len(single)) // max(len(base), 1)))
    ramps = [chain_colours(c, np.linspace(0, 1, n_ramp)) for c in base]
    return (np.vstack([sp, sp * 0.75, sp * 0.5, *ramps, single]) * 255).round().astype(np.uint8)


def shared_palette(frames, reserved, n_total=256):
    """One 256-colour palette image for all frames: the reserved colours, then a median cut of the frames for
    the rest. Mostly grey frames would make a palette of each frame spend its colours on greys and shift the
    v_A and radio colours from frame to frame."""
    from PIL import Image
    n_rest = n_total - len(reserved)
    sample = np.vstack([f[::2, ::2] if f.shape[0] * f.shape[1] > 10**6 else f for f in frames])
    both = Image.fromarray(np.ascontiguousarray(sample).reshape(-1, sample.shape[1], 3))
    rest = both.quantize(colors=n_rest, method=Image.Quantize.MEDIANCUT, dither=Image.Dither.NONE).convert("RGB")
    rest = np.unique(np.asarray(rest).reshape(-1, 3), axis=0)[:n_rest]
    pal = np.vstack([reserved, rest])
    pimg = Image.new("P", (1, 1))
    pimg.putpalette(list(pal.ravel()) + [0] * (768 - pal.size))
    return pimg


def save_gif(frames, path, fps, reserved, scale=0.5):
    """Write the RGB frames (uint8 arrays) as a GIF with one shared palette and no dithering."""
    from PIL import Image
    imgs = [Image.fromarray(f) for f in frames]
    if scale != 1:
        imgs = [i.resize((max(1, round(i.width * scale)), max(1, round(i.height * scale))), Image.Resampling.LANCZOS)
                for i in imgs]
    pimg = shared_palette([np.asarray(i) for i in imgs], reserved)
    q = [i.quantize(palette=pimg, dither=Image.Dither.NONE) for i in imgs]
    q[0].save(path, save_all=True, append_images=q[1:], duration=round(1000 / fps), loop=0, optimize=False)
    return Path(path)


def ffmpeg_exe():
    """The ffmpeg of the package imageio-ffmpeg, or else the one on the path."""
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except Exception:
        exe = shutil.which("ffmpeg")
        if exe:
            return exe
    raise RuntimeError("an MP4 needs ffmpeg: install it with pip install \"corona3d[animate]\" (the package "
                       "imageio-ffmpeg) or put ffmpeg on the path")


def save_mp4(frames, path, fps):
    """Write the RGB frames as an H.264 MP4 (libx264, yuv420p; the size is cut to even numbers)."""
    h, w = frames[0].shape[0] // 2 * 2, frames[0].shape[1] // 2 * 2
    cmd = [ffmpeg_exe(), "-y", "-loglevel", "error", "-f", "rawvideo", "-pix_fmt", "rgb24", "-s", f"{w}x{h}",
           "-r", str(fps), "-i", "-", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
           "-movflags", "+faststart", str(path)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        for f in frames:
            proc.stdin.write(np.ascontiguousarray(f[:h, :w]).tobytes())
        proc.stdin.close()
        err = proc.stderr.read().decode(errors="replace")
    finally:
        proc.wait()
    if proc.returncode:
        raise RuntimeError(f"ffmpeg failed: {err[-500:]}")
    return Path(path)


def animate(out="figures/corona3d_animation", data_root=".", start=None, end=None, step_s=None, fps=2.0,
            gif=True, mp4=False, half=None, params=None, panels=None, layout="landscape", law=None,
            gif_scale=0.5, **layers):
    """Render one frame per time from start to end (ISO texts; default: the radio images of the configuration)
    and write out.gif and/or out.mp4 (for layout "separate", one animation per view, out_VIEW.gif).
    law: function from a time to the apex height (height_law; default: the line through the two shells of
    the configuration). The other arguments are those of make_figure. Returns the list of files."""
    from . import reconstruction as rc
    from .layers import ON
    from .paths import config, merge
    from .whitelight import utc
    bad = set(layers) - set(rc.LAYERS)
    if bad:
        raise ValueError(f"unknown layers {sorted(bad)}; choose from {', '.join(rc.LAYERS)}")
    if half is not None and half <= 0:
        raise ValueError("half must be positive")
    on = dict(ON, **layers)
    p0 = merge(config(data_root), params)
    rc.check_opacity(dict(p0["opacity"], **({"wl": p0["wl_opacity"]} if "wl_opacity" in p0 else {})))
    ri = p0["radio_images"] if on["radio"] else []
    if start is None or end is None:
        if not ri:
            raise ValueError("give --start and --end")
        start = start or min((rs["start"] for rs in ri), key=utc)
        end = end or max((rs["end"] for rs in ri), key=utc)
    step_s = step_s or (ri[0]["step_s"] if ri else 60.0)
    times = frame_times(utc(start), utc(end), step_s)
    law = law or height_law(p0)
    rc.select_views(p0, panels)
    out = Path(out) if Path(out).is_absolute() else Path(data_root) / out
    out.parent.mkdir(parents=True, exist_ok=True)
    # one colour scale of the radio time for every frame, from the first radio image to the last frame
    rng = p0.get("radio_time_range") or ([min((rs["start"] for rs in ri), key=utc), times[-1].isoformat()] if ri else None)
    cache, frames = {}, []
    for t in times:                                        # pass 1: the layers of every frame
        iso = t.isoformat()
        h = law(t)
        if h <= 1.0:
            raise ValueError(f"{iso}: the height law gives {h:.2f} R_sun, below the solar surface")
        pf = dict(p0, time=iso, anim=dict(height=h), radio_time_range=rng, _png_dir=None,
                  radio_images=[dict(rs, end=iso) for rs in p0["radio_images"]])
        sc = rc.reconstruct(pf, on, data_root, cache)
        sc["times"]["frame"] = iso
        frames.append((pf, sc))
        print(f"{iso[11:19]} UT: shell height {h:.3f} R_sun")
    # a coronagraph image keeps one grey scale: the largest 98th percentile of all frames
    clim = {i: max(np.nanpercentile(np.abs(sc[i]["wl"]), 98) for _, sc in frames) for i in rc.WL if on[i]}
    pix = {}
    for pf, sc in frames:                                  # pass 2: the panels, with fixed cameras and limits
        pf["wl_clim"] = clim
        for key, fig in rc.compose(sc, on, pf, half, panels, layout):
            pix.setdefault(key, []).append(rc.fig_array(fig))
    res = reserved_colours(p0)
    files = []
    for key, arrs in pix.items():
        base = out if key is None else out.with_name(f"{out.name}_{key}")
        if gif or not mp4:
            files.append(save_gif(arrs, base.with_name(base.name + ".gif"), fps, res, gif_scale))
        if mp4:
            files.append(save_mp4(arrs, base.with_name(base.name + ".mp4"), fps))
    for f in files:
        print("saved", f)
    return files
