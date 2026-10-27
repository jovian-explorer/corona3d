"""The burst interval of the radio images, chosen from the images (`corona3d radio-interval`).

For every frequency of a `radio_images` entry, the brightest pixel of each image in the search
window is compared with `min_tb`. The frames at or above it form runs, separated where two such
frames are more than three times the cadence of the images apart (at least 60 s); the run that
holds the brightest frame of the frequency is its burst interval. The entry takes the earliest
start and the latest end of its frequencies. Only entries marked `"interval_source": "default"`
(written by `corona3d new`) are changed; an entry without that key, whose interval was given by
hand, and an entry already chosen from the images (`"data"`) are kept unless `force` is set.
"""
import json
import math
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

ISO = "%Y-%m-%dT%H:%M:%S"


def burst_run(series, min_tb):
    """series: [(time, peak [K])] in time order. Returns dict(start, end, n, peak, t_peak) of the run
    of frames with peak >= min_tb that holds the brightest frame, or None when there is none."""
    series = sorted(series)
    if len(series) < 1 or max(p for _, p in series) < min_tb:
        return None
    dts = [(b[0] - a[0]).total_seconds() for a, b in zip(series, series[1:])]
    gap = max(60.0, 3 * float(np.median(dts))) if dts else 60.0
    runs, cur = [], []
    for t, p in series:
        if p < min_tb:
            continue
        if cur and (t - cur[-1][0]).total_seconds() > gap:
            runs.append(cur)
            cur = []
        cur.append((t, p))
    runs.append(cur)
    best = max(runs, key=lambda r: max(p for _, p in r))
    t_peak, peak = max(best, key=lambda x: x[1])
    return dict(start=best[0][0], end=best[-1][0], n=len(best), peak=float(peak), t_peak=t_peak)


def peaks(instrument, folder, freqs, t0, t1):
    """{freq: [(time, brightest pixel [K])]} of the images between t0 and t1 in the folder."""
    from . import nrh, ovro_lwa
    out = {f: [] for f in freqs}
    if instrument == "nrh":
        for f in freqs:
            for p in nrh.files(folder).get(round(f, 1), []):
                r = nrh.read(p)
                out[f] += [(t, float(np.nanmax(im))) for t, im in zip(r["times"], r["images"]) if t0 <= t <= t1]
    else:
        for tf, p in ovro_lwa.files(folder):
            if not t0 <= tf <= t1:
                continue
            _, bands, imgs = ovro_lwa.read(p)
            for f in freqs:
                k = int(np.argmin(abs(bands - f)))
                if abs(bands[k] - f) <= 3.0:
                    out[f].append((tf, float(np.nanmax(imgs[k]))))
    return out


def entry_interval(rs, folder, t0, t1, reader=peaks):
    """Burst interval of one radio_images entry: ({freq: run or None}, (start, end) or None)."""
    series = reader(rs["instrument"], folder, rs["freqs"], t0, t1)
    runs = {f: burst_run(series[f], rs["min_tb"]) for f in rs["freqs"]}
    found = [r for r in runs.values() if r]
    return runs, ((min(r["start"] for r in found), max(r["end"] for r in found)) if found else None)


def select(home, cfg, minutes=10.0, force=False, reader=peaks, write=True):
    """Set start, end, and step_s of every radio_images entry from the images in its folder,
    within `minutes` of the event time; returns [(entry index, message)] for the report.
    Writes config.json (unless write is false) and marks the entry `"interval_source": "data"`."""
    from .whitelight import utc
    t = utc(cfg["time"])
    t0, t1 = t - timedelta(minutes=minutes), t + timedelta(minutes=minutes)
    entries = [dict(e) for e in cfg["radio_images"]]
    notes, changed = [], False
    for i, rs in enumerate(entries):
        label = f"{rs['instrument']} entry {i}"
        folder = Path(home) / rs["folder"]
        src = rs.get("interval_source", "manual")
        if src != "default" and not force:
            notes.append((i, f"{label}: interval {rs['start']} to {rs['end']} " +
                          ("was already chosen from the images" if src == "data" else "was given by hand") + "; kept"))
            continue
        if not folder.exists() or not any(folder.iterdir()):
            notes.append((i, f"{label}: no images in {folder}; interval left at {rs['start']} to {rs['end']}"))
            continue
        runs, span = entry_interval(rs, folder, t0, t1, reader)
        for f, r in runs.items():
            notes.append((i, f"{label} {f:g} MHz: " + (
                f"{r['start']:%H:%M:%S} to {r['end']:%H:%M:%S} UT, {r['n']} images, peak {r['peak']:.2e} K at "
                f"{r['t_peak']:%H:%M:%S}" if r else f"no image at or above {rs['min_tb']:.1e} K")))
        if span is None:
            notes.append((i, f"{label}: no frequency reaches min_tb within {minutes:g} min of {cfg['time']}; "
                             f"interval left at {rs['start']} to {rs['end']}"))
            continue
        rs["start"], rs["end"] = span[0].strftime(ISO), span[1].strftime(ISO)
        length = (span[1] - span[0]).total_seconds()
        rs["step_s"] = max(int(rs.get("step_s", 60)), int(math.ceil(length / 12 / 10)) * 10)   # at most about 13 blobs
        rs["interval_source"] = "data"
        changed = True
        notes.append((i, f"{label}: interval {rs['start']} to {rs['end']}, one image every {rs['step_s']} s"))
    if changed and write:
        f = Path(home) / "config.json"
        user = json.loads(f.read_text(encoding="utf-8"))
        user["radio_images"] = entries
        f.write_text(json.dumps(user, indent=2) + "\n", encoding="utf-8")
    return notes


def run(home, cfg, a):
    """The `corona3d radio-interval` command."""
    for _, msg in select(home, cfg, a.search, a.force):
        print(msg)
