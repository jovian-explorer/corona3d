"""Positions of the radio sources under every placement (`corona3d sources`): the sphere of the four-fold Newkirk
height and the plasma-frequency surface of the scaled MAS density, for fundamental (F) and harmonic (H) emission,
with and without the correction for scattering, together with the polarisation that decides between F and H, the
observed and the scattering size of each source, and the OVRO-LWA position check on the quiet Sun.

The table `analysis/radio_sources.csv` has one row per centroid of the radio layer (as `corona3d export`), with the
observed sky-plane position (x_obs, y_obs, R_sun, seen from Earth), the offset subtracted from it (OVRO-LWA
`position_offset`), the degree of circular polarisation (NRH Stokes V), the emission chosen for the row (`emission`
of the entry, or from the polarisation with "auto") and its reason, and for each placement P in sphere_H, sphere_F,
model_H, model_F the distance r_P and the position x_P, y_P, z_P (Carrington Cartesian). The columns *_scat repeat
the chosen emission with the scattering correction: the radial shift dr and its range dr_lo to dr_hi, the corrected
sky-plane position, and the new three-dimensional position on each surface. fwhm_obs_deg is the diameter of the
circle with the area of the pixels above `level` times the peak, and fwhm_scat_deg the size expected from scattering
alone (Kontar et al. 2019).
"""
import csv
import json
from pathlib import Path

import numpy as np


def observed(home, p):
    """Radio table of the centroids placed by the default rule (no offset, H, sphere), from corona3d.api.build."""
    from .api import build
    params = dict(radio_deproject="sphere", radio_emission="H", scattering=dict(apply=False),
                  radio_images=[{k: v for k, v in e.items() if k not in ("emission", "position_offset")}
                                for e in p["radio_images"]])
    rec = build(home, params, surface=False, open=False, closed=False, streamer=False, gcs=False, radio=True)
    return rec.radio


def observed_size(home, p, inst, f, t):
    """Equivalent diameter [deg] of the pixels above `level` times the peak in the image of the source."""
    from . import nrh, ovro_lwa
    from .whitelight import utc
    e = next((e for e in p["radio_images"] if e["instrument"] == inst and any(abs(f - x) < 0.05 for x in e["freqs"])), None)
    if e is None or not t:
        return np.nan
    mod = nrh if inst == "nrh" else ovro_lwa
    img, _, x, y = mod.image_at(Path(home) / e["folder"], f, utc(t))
    area = (img >= e["level"] * np.nanmax(img)).sum() * abs((x[1] - x[0]) * (y[1] - y[0]))
    from astropy import units as u
    from sunpy.coordinates import sun
    rs_deg = sun.angular_radius(p["time"]).to_value(u.deg)
    return float(2 * np.sqrt(area / np.pi) * rs_deg)


def analyse(home, p=None, ovro_offset=False, log=print):
    """Rows of radio_sources.csv and the summary per instrument and frequency."""
    from . import heights, mast_model, nrh, ovro_lwa, scattering
    from .paths import config, merge, model_dir
    from .pb import resolve
    p = resolve(merge(config(home), p))
    rt = observed(home, p)
    summ = dict(time=p["time"], radio_deproject=p.get("radio_deproject", "sphere"), density_scale=p["density_scale"],
                eps2_h_km=list((p.get("scattering") or {}).get("eps2_h_km") or scattering.EPS2_H), groups=[])
    if ovro_offset:
        try:
            summ["ovro_lwa_offset"] = ovro_check(home, p, log)
        except Exception as e:                             # a failed download must not stop the source positions
            log(f"sources: OVRO-LWA position check failed ({e})")
            summ["ovro_lwa_offset"] = dict(failed=str(e))
    if not rt or not len(rt.get("x_rsun", [])):
        log("sources: no radio source in the radio layer")
        return [], summ
    m = mast_model.load(model_dir(home, p), density_scale=p["density_scale"], va_density=p["va_density"])
    ne_at = heights.ne_sampler(m, p)
    r_max = float(m["r"][-1])
    pol = nrh.polarisation_rows(home, p, rt)
    sc = p.get("scattering") or {}
    eps = tuple(sc.get("eps2_h_km") or scattering.EPS2_H)
    rows = []
    for j, (inst, f, t, x, y, tb) in enumerate(zip(rt["instrument"], rt["freq_mhz"], rt["time_utc"], rt["x_rsun"],
                                                   rt["y_rsun"], rt["tb_k"])):
        f = float(f)
        e = next((e for e in p["radio_images"] if e["instrument"] == inst and any(abs(f - q) < 0.05 for q in e["freqs"])), {})
        dx, dy = ovro_lwa.offset_for(home, e, f) if inst == "ovro_lwa" and e else (0.0, 0.0)
        xo, yo = float(x) - dx, float(y) - dy
        dcp, valid = float(pol["dcp_source"][j]), bool(pol["dcp_valid"][j])
        em = heights.emission_of(p, e, f) if inst != "user" else heights.emission_of(p)
        if em == "auto":
            em, why = heights.emission_from_polarisation(dcp, valid, p.get("dcp_fundamental", 0.3))
        else:
            why = "given in the configuration" + (f"; V/I = {dcp:+.2f}" if np.isfinite(dcp) and valid else "")
        row = dict(instrument=str(inst), freq_mhz=f, time_utc=str(t), x_obs=float(x), y_obs=float(y), tb_k=float(tb),
                   offset_x=dx, offset_y=dy, dcp_source=dcp, dcp_valid=valid, emission=em, emission_reason=why)
        for mode in ("sphere", "model"):
            for e_ in ("H", "F"):
                q = heights.locate(xo, yo, f, p, e_ == "H", mode, ne_at, r_max)
                k = f"{mode}_{e_}"
                row.update({f"r_{k}": q["r"], f"x_{k}": q["xyz"][0], f"y_{k}": q["xyz"][1], f"z_{k}": q["xyz"][2],
                            f"surface_{k}": q["surface"]})
            if em == "H":                                  # the shift that fundamental emission would have
                qf = heights.locate(xo, yo, f, p, False, mode, ne_at, r_max, True, eps)
                row.update({f"dr_{mode}_F": qf["dr"], f"dr_lo_{mode}_F": qf["dr_lo"], f"dr_hi_{mode}_F": qf["dr_hi"]})
            q = heights.locate(xo, yo, f, p, em == "H", mode, ne_at, r_max, True, eps)
            if em == "F":
                row.update({f"dr_{mode}_F": q["dr"], f"dr_lo_{mode}_F": q["dr_lo"], f"dr_hi_{mode}_F": q["dr_hi"]})
            k = f"{mode}_scat"
            row.update({f"dr_{k}": q["dr"], f"dr_lo_{k}": q["dr_lo"], f"dr_hi_{k}": q["dr_hi"], f"xsky_{k}": q["x"],
                        f"ysky_{k}": q["y"], f"r_{k}": q["r"], f"x_{k}": q["xyz"][0], f"y_{k}": q["xyz"][1],
                        f"z_{k}": q["xyz"][2]})
        row["fwhm_obs_deg"] = observed_size(home, p, inst, f, t) if inst != "user" else np.nan
        row["fwhm_scat_deg"] = float(scattering.fwhm_deg(f))
        rows.append(row)
    keys = sorted({(r["instrument"], r["freq_mhz"]) for r in rows})
    for inst, f in keys:
        R = [r for r in rows if r["instrument"] == inst and r["freq_mhz"] == f]
        med = lambda k: float(np.nanmedian([r[k] for r in R])) if np.isfinite([r[k] for r in R]).any() else np.nan
        on = lambda k: float(np.nanmedian([r[f"r_{k}"] for r in R if r[f"surface_{k}"] == "model"]))             if any(r[f"surface_{k}"] == "model" for r in R) else np.nan
        ems = sorted({r["emission"] for r in R})
        g = dict(instrument=inst, freq_mhz=f, n=len(R), emission="/".join(ems),
                 dcp_median=float(np.nanmedian([r["dcp_source"] if r["dcp_valid"] else np.nan for r in R])) if any(r["dcp_valid"] for r in R) else None,
                 **{f"r_{k}": med(f"r_{k}") for k in ("sphere_H", "sphere_F", "sphere_scat", "model_scat")},
                 r_model_H=on("model_H"), r_model_F=on("model_F"),
                 dr_sphere=med("dr_sphere_scat"), dr_sphere_lo=med("dr_lo_sphere_scat"), dr_sphere_hi=med("dr_hi_sphere_scat"),
                 dr_model=med("dr_model_scat"), dr_model_lo=med("dr_lo_model_scat"), dr_model_hi=med("dr_hi_model_scat"),
                 dr_sphere_F=med("dr_sphere_F"), dr_model_F=med("dr_model_F"),
                 model_fallback=int(sum(r["surface_model_H"] != "model" for r in R)),
                 model_fallback_F=int(sum(r["surface_model_F"] != "model" for r in R)),
                 fwhm_obs_deg=med("fwhm_obs_deg"), fwhm_scat_deg=med("fwhm_scat_deg"),
                 offset=[R[0]["offset_x"], R[0]["offset_y"]])
        summ["groups"].append(g)
        log(f"sources: {inst} {f:g} MHz ({len(R)}): r sphere H {g['r_sphere_H']:.2f}, F {g['r_sphere_F']:.2f}; model H "
            f"{g['r_model_H']:.2f}, F {g['r_model_F']:.2f}; scattering dr {g['dr_sphere']:.2f} (sphere), {g['dr_model']:.2f} (model)")
    return rows, summ


def ovro_check(home, p, log=print, step_s=600, hours=(-2.0, 1.0)):
    """Position check of OVRO-LWA on the quiet Sun: image files one every step_s seconds from hours[0] to hours[1]
    around the event are downloaded (corona3d.ovro_lwa.fetch_quiet), the offset of the radio disc from Sun centre is
    measured in the quiet images (ovro_lwa.quiet_offsets), and the result is written to analysis/ovro_lwa_offset.json,
    which `position_offset: "quiet_sun"` uses."""
    from datetime import timedelta
    from . import ovro_lwa
    from .whitelight import utc
    entries = [e for e in p["radio_images"] if e["instrument"] == "ovro_lwa"]
    if not entries:
        return None
    e = entries[0]
    t = utc(p["time"])
    t0, t1 = t + timedelta(hours=hours[0]), t + timedelta(hours=hours[1])
    n = ovro_lwa.fetch_quiet(home, e, t0, t1, step_s)
    res = ovro_lwa.quiet_offsets(ovro_lwa.quiet_dir(home, e), e["freqs"], t0, t1)
    fit = ovro_lwa.refraction_fit(res)
    rec = dict(window=[t0.isoformat(), t1.isoformat()], step_s=step_s, failed_downloads=n, criteria=ovro_lwa.QUIET,
               bands={f"{f:g}": {k: v for k, v in r.items() if k != "frames"} | dict(frames=r["frames"]) for f, r in res.items()},
               refraction_fit=fit)
    f = Path(home) / ovro_lwa.OFFSET_FILE
    f.parent.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(rec, indent=1, default=str) + "\n", encoding="utf-8")
    for fr, r in res.items():
        log(f"sources: OVRO-LWA {fr:g} MHz: {r['n']} quiet of {r['n_frames']} images" +
            (f", offset ({r['dx']:+.3f}, {r['dy']:+.3f}) R_sun" if r["n"] else ""))
    return {k: v for k, v in rec.items() if k != "bands"} | dict(bands={k: {q: w for q, w in v.items() if q != "frames"}
                                                                         for k, v in rec["bands"].items()})


def write(out, rows, summ):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    if rows:
        with open(out / "radio_sources.csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, list(rows[0]))
            w.writeheader()
            w.writerows(rows)
    (out / "radio_sources.json").write_text(json.dumps(summ, indent=1, default=str) + "\n", encoding="utf-8")
    return out


def _f(x, fmt=".2f"):
    return "n/a" if x is None or not np.isfinite(x) else format(x, fmt)


def report_lines(summ):
    L = [f"- Radio sources placed on the {summ['radio_deproject']} surface in the figures; density_scale {summ['density_scale']}; "
         f"scattering with eps^2/h from {summ['eps2_h_km'][0]:g} to {summ['eps2_h_km'][1]:g} per km (Chrysaphi et al. 2018)."]
    for g in summ["groups"]:
        pol = f", median V/I {g['dcp_median']:+.2f}" if g.get("dcp_median") is not None else ", no Stokes V"
        L.append(f"- {g['instrument']} {g['freq_mhz']:g} MHz ({g['n']} positions, emission {g['emission']}{pol}): distance "
                 f"{_f(g['r_sphere_H'])} (H) and {_f(g['r_sphere_F'])} R_sun (F) on the Newkirk sphere, {_f(g['r_model_H'])} (H) and "
                 f"{_f(g['r_model_F'])} R_sun (F) on the MAS plasma surface" +
                 (f" (the lines of sight of {g['model_fallback']} (H) and {g['model_fallback_F']} (F) positions miss the surface; they "
                  "are left out of these medians and kept on the sphere)" if g["model_fallback"] or g["model_fallback_F"] else "")
                 + f"; the scattering shift is {_f(g['dr_sphere'])} ({_f(g['dr_sphere_lo'])} to {_f(g['dr_sphere_hi'])}) R_sun in the "
                 f"Newkirk density and {_f(g['dr_model'])} ({_f(g['dr_model_lo'])} to {_f(g['dr_model_hi'])}) R_sun in the MAS density, which "
                 f"moves the sources to {_f(g['r_sphere_scat'])} and {_f(g['r_model_scat'])} R_sun (for fundamental emission the shift "
                 f"would be {_f(g['dr_sphere_F'])} and {_f(g['dr_model_F'])} R_sun); observed size {_f(g['fwhm_obs_deg'])} deg "
                 f"against {_f(g['fwhm_scat_deg'])} deg from scattering alone.")
    o = summ.get("ovro_lwa_offset")
    if o and "failed" in o:
        L.append(f"- OVRO-LWA position check failed ({o['failed']}); no offset was applied.")
    elif o:
        good = {k: v for k, v in o["bands"].items() if v["n"]}
        if good:
            L.append("- OVRO-LWA position check on the quiet Sun: " + "; ".join(
                f"{k} MHz ({v['n']} images) offset ({v['dx']:+.3f}, {v['dy']:+.3f}) R_sun" for k, v in good.items()) + ".")
        else:
            n = sum(v["n_frames"] for v in o["bands"].values())
            L.append(f"- OVRO-LWA position check: none of the {n} band images between {o['window'][0][11:16]} and "
                     f"{o['window'][1][11:16]} UT shows a quiet disc (equivalent radius {o['criteria']['radius'][0]:g} to "
                     f"{o['criteria']['radius'][1]:g} R_sun, peak below {o['criteria']['contrast']:g} times and disc brightness below "
                     f"{o['criteria']['max_tb']:.0e} K), so that the offset could not be measured and no correction is applied.")
    return L


def run(home, cfg, a):
    """The `corona3d sources` command."""
    from .report_sections import update
    rows, summ = analyse(home, ovro_offset=a.ovro_offset)
    out = write(Path(home) / a.out, rows, summ)
    lines = report_lines(summ)
    print("\n".join(lines))
    update(home, "Radio source positions", lines)
    print(f"wrote {out}")
