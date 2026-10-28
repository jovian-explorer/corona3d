"""Type II and III bursts of a day from the SWPC edited event list: the parser on a stored file of
2013-10-25 (public NOAA product), the choice with --burst, the cache, and the fallback and error paths.
The network is replaced by an opener that returns the stored file."""
import io
from datetime import datetime
from pathlib import Path

import pytest

from corona3d import bursts
from corona3d.cli import main

SAMPLE = Path(__file__).resolve().parent / "data" / "20131025events.txt"
needs_sample = pytest.mark.skipif(not SAMPLE.exists(), reason="source tree not present")


def opener_of(text, urls=None, fail=()):
    def op(url, timeout=0):
        if urls is not None:
            urls.append(url)
        if any(url.startswith(f) for f in fail):
            raise OSError("unreachable")
        return io.BytesIO(text.encode())
    return op


@needs_sample
def test_the_parser_finds_every_type_ii_and_iii_burst_of_the_day():
    rows = bursts.select(bursts.parse(SAMPLE.read_text(), "2013-10-25"))
    assert len(rows) == 25 and [r["type"] for r in rows].count("II") == 5
    ii = [r for r in rows if r["type"] == "II"]
    assert [r["begin"] for r in ii][-1] == datetime(2013, 10, 25, 14, 58)           # the event of the example
    assert ii[-1]["end"] == datetime(2013, 10, 25, 15, 14) and ii[-1]["speed_km_s"] == 2078.0
    assert ii[-1]["obs"] == "SAG" and ii[-1]["freq_mhz"] == "025-180"
    assert rows == sorted(rows, key=lambda r: r["begin"])
    assert all(r["type"] in ("II", "III") for r in rows)


def test_flags_in_front_of_times_and_a_missing_end_are_handled():
    text = (" 10      B0850  U0851      ////  SVI  C   RSP  045-145   II/1    686       1882\n"
            " 11 +     0903   ////      0905  SVI  C   RSP  025-180   III/1                 \n"
            " 12       0905   0905      0905  SVI  G   RBR  410       160                   \n")
    rows = bursts.parse(text, "2013-10-25")
    assert [r["type"] for r in rows] == ["II", "III"]
    assert rows[0]["begin"] == datetime(2013, 10, 25, 8, 50) and rows[0]["end"] is None and rows[0]["speed_km_s"] == 686.0


@needs_sample
def test_choose_lists_the_candidates_then_returns_the_begin_time_of_the_chosen_one(tmp_path):
    op = opener_of(SAMPLE.read_text())
    lines = []
    assert bursts.choose("2013-10-25", None, tmp_path, opener=op, log=lines.append) is None
    assert any("type II  14:58 to 15:14 UT" in x and "2078" in x for x in lines) and "--burst" in lines[-1]
    n = next(i for i, x in enumerate(lines[1:], 1) if "14:58" in x)
    assert bursts.choose("2013-10-25", n, tmp_path, opener=op) == "2013-10-25T14:58"
    with pytest.raises(RuntimeError, match="has 25 bursts"):
        bursts.choose("2013-10-25", 26, tmp_path, opener=op)


@needs_sample
def test_the_file_is_cached_and_the_ftp_server_is_the_fallback(tmp_path):
    urls = []
    text = bursts.fetch("2013-10-25", tmp_path, opener_of(SAMPLE.read_text(), urls, fail=("https",)))
    assert "Edited Events for 2013 Oct 25" in text and urls[0].startswith("https://www.ngdc.noaa.gov/") and urls[-1].startswith("ftp://")
    assert (tmp_path / "data" / "bursts" / "20131025events.txt").exists()
    urls.clear()
    bursts.fetch("2013-10-25", tmp_path, opener_of("", urls))                      # from the cache: no request
    assert urls == []


def test_an_unreachable_source_gives_a_clear_error(tmp_path):
    with pytest.raises(RuntimeError, match="could not be fetched.*give --time"):
        bursts.fetch("2013-10-25", tmp_path, opener_of("", fail=("https", "ftp")))
    with pytest.raises(RuntimeError, match="could not be fetched"):                  # a page that is not an event file
        bursts.fetch("2013-10-25", tmp_path, opener_of("<html>404</html>"))


def test_a_day_without_bursts_asks_for_a_time(tmp_path):
    text = ":Product: x\n#  Edited Events for 2013 Oct 25\n 10       0613   0613     A0644  LEA  2   FLA  N06W24    SF\n"
    with pytest.raises(RuntimeError, match="no type II/III burst"):
        bursts.choose("2013-10-25", None, tmp_path, opener=opener_of(text))


def test_date_without_time_is_recognised():
    assert bursts.is_date("2024-05-09") and not bursts.is_date("2024-05-09T17:32") and not bursts.is_date("17:32")


@needs_sample
def test_event_with_only_a_date_lists_the_bursts_and_burst_n_sets_the_time(tmp_path, monkeypatch, capsys):
    from corona3d import pipeline
    monkeypatch.setattr(bursts.urllib.request, "urlopen", opener_of(SAMPLE.read_text()))
    main(["--home", str(tmp_path), "event", "--time", "2013-10-25", "--flare", "S06E69"])        # list only
    out = capsys.readouterr().out
    assert "type II  14:58 to 15:14 UT" in out and not (tmp_path / "config.json").exists()
    seen = {}
    monkeypatch.setattr(pipeline, "run", lambda h, a: seen.update(time=a.time) or 0)
    with pytest.raises(SystemExit) as e:
        main(["--home", str(tmp_path), "event", "--time", "2013-10-25", "--burst", "18", "--flare", "S06E69"])
    assert e.value.code == 0 and seen["time"] == "2013-10-25T14:58"
    with pytest.raises(SystemExit, match="--burst needs a date"):
        main(["--home", str(tmp_path), "event", "--time", "2013-10-25T15:00", "--burst", "1", "--flare", "S06E69"])
