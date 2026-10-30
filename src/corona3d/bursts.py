"""Type II and type III radio bursts of a day from the edited event list of NOAA SWPC.

When `corona3d event` is given only a date, the burst times of that day are listed and the user
chooses one with `--burst N`. The source is the daily file `YYYYMMDDevents.txt` of the SWPC edited
events (patrol reports of ground radio observatories, code RSP), kept by NOAA NCEI since 1996 and reachable without a login at
https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/solar_event_reports/ ;
the recent files are also on ftp://ftp.swpc.noaa.gov/pub/indices/events/ (used as a fallback).
A line of type RSP gives the begin and end time [UT], the observatory, the frequency range [MHz],
the burst type (II, III, IV, V, VI, CTM with an intensity class) and, for type II, the shock speed [km/s].
"""
import re
import urllib.request
from datetime import datetime, timedelta
from pathlib import Path

NCEI = ("https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/solar_event_reports/"
        "{y}/{m}/{d}events.txt")
FTP = "ftp://ftp.swpc.noaa.gov/pub/indices/events/{d}events.txt"
TYPES = ("II", "III")
DATE = re.compile(r"\d{4}-\d{2}-\d{2}")


def is_date(text):
    """True for a date without a time, e.g. 2024-05-09."""
    return bool(DATE.fullmatch(text.strip()))


def parse(text, day):
    """Burst rows [dict(type, begin, end, obs, freq_mhz, speed_km_s)] of the RSP lines of a SWPC
    event file for the date `day` (datetime.date or 'YYYY-MM-DD'), in time order. A begin or end
    time may carry a flag letter (A, B, D, U) and an end of '////' is missing; the end then is None."""
    day = datetime.fromisoformat(str(day)).date()
    rows = []
    for line in text.splitlines():
        p = line.split()
        if len(p) < 3 or line.startswith(("#", ":")) or "RSP" not in p:
            continue
        i = p.index("RSP")
        if i < 5 or i + 2 >= len(p):                                # ..., begin, max, end, observatory, quality, RSP, range, type
            continue
        begin, end = (re.sub(r"^[A-Z]", "", p[i - k]) for k in (5, 3))
        if not re.fullmatch(r"\d{4}", begin):
            continue
        midnight = datetime.combine(day, datetime.min.time())
        clock = lambda t: midnight + timedelta(hours=int(t[:2]), minutes=int(t[2:]))
        kind = p[i + 2].split("/")[0]
        speed = float(p[i + 3]) if kind == "II" and len(p) > i + 3 and re.fullmatch(r"\d+(\.\d+)?", p[i + 3]) else None
        t_end = clock(end) if re.fullmatch(r"\d{4}", end) else None
        if t_end is not None and t_end < clock(begin):               # the burst crosses midnight
            t_end += timedelta(days=1)
        rows.append(dict(type=kind, begin=clock(begin), end=t_end,
                         obs=p[i - 2], freq_mhz=p[i + 1], speed_km_s=speed))
    return sorted(rows, key=lambda r: r["begin"])


def select(rows, types=TYPES):
    """The rows whose burst type is in `types` (default II and III)."""
    return [r for r in rows if r["type"] in types]


def fetch(day, home=None, opener=urllib.request.urlopen):
    """Text of the SWPC event file of `day`, from the cache folder data/bursts of `home` when present;
    otherwise from NCEI and, if that fails, from the SWPC FTP server. Raises RuntimeError with both reasons."""
    day = datetime.fromisoformat(str(day)).date()
    d = f"{day:%Y%m%d}"
    cache = Path(home) / "data" / "bursts" / f"{d}events.txt" if home else None
    if cache and cache.exists():
        return cache.read_text(encoding="utf-8", errors="replace")
    why = []
    for url in (NCEI.format(y=f"{day:%Y}", m=f"{day:%m}", d=d), FTP.format(d=d)):
        try:
            with opener(url, timeout=60) as r:
                text = r.read().decode("utf-8", errors="replace")
            if "Edited Events" not in text:
                raise ValueError("not an event file")
            if cache:
                cache.parent.mkdir(parents=True, exist_ok=True)
                cache.write_text(text, encoding="utf-8")
            return text
        except Exception as e:
            why.append(f"{url}: {e}")
    raise RuntimeError(f"the SWPC event list of {day} could not be fetched ({'; '.join(why)}); "
                       "give --time with the time of the burst instead")


def listing(rows):
    """Numbered text lines of burst rows, for the user to choose from with --burst."""
    out = []
    for n, r in enumerate(rows, 1):
        end = f"{r['end']:%H:%M}" if r["end"] else "----"
        sp = f", shock speed {r['speed_km_s']:g} km/s" if r["speed_km_s"] else ""
        out.append(f"  {n:2d}  type {r['type']:<3s} {r['begin']:%H:%M} to {end} UT  {r['freq_mhz']} MHz  {r['obs']}{sp}")
    return out


def choose(date, burst, home=None, types=TYPES, opener=urllib.request.urlopen, log=print):
    """Time (ISO text, minutes) of burst number `burst` (1-based) of `date`; with burst None the
    candidates are listed and None is returned."""
    rows = select(parse(fetch(date, home, opener), date), types)
    if not rows:
        raise RuntimeError(f"no type {'/'.join(types)} burst in the SWPC list of {date}; give --time")
    if burst is None:
        log(f"Type {'/'.join(types)} bursts of {date} (SWPC edited events; times are UT):")
        for t in listing(rows):
            log(t)
        log("Choose one with --burst N; its begin time becomes --time (use --time to give another).")
        return None
    if not 1 <= burst <= len(rows):
        raise RuntimeError(f"--burst {burst}: the list of {date} has {len(rows)} bursts")
    return f"{rows[burst - 1]['begin']:%Y-%m-%dT%H:%M}"
