"""corona3d new: flare positions, the PSI run of the Carrington rotation, the STEREO spacecraft
of the date, and the OVRO-LWA entry from 2024 April onward, without network access."""
import pytest

from corona3d.events import new_config, parse_flare, psi_run


def test_flare_positions():
    assert parse_flare("S14W28") == (28.0, -14.0)
    assert parse_flare("n06e69") == (-69.0, 6.0)
    with pytest.raises(ValueError):
        parse_flare("14S28W")


def test_psi_run_is_that_of_the_carrington_rotation():
    listing = "cr2283-high/ cr2284-high/ cr2143-medium/ cr2143-high/"
    assert psi_run("2024-05-09T17:32", listing) == "cr2284-high/hmi_mast_mas_std_0201"   # CR 2284 from May 6
    assert psi_run("2013-10-25T15:00", listing) == "cr2143-medium/hmi_mast_mas_std_0201"
    with pytest.raises(ValueError, match="2290"):
        psi_run("2024-10-20T12:00", listing)


def test_new_event_configuration():
    obs = ([317.1, -3.2], [329.6, -1.7])
    c = new_config("2024-05-09T17:32", "S14W28", earth_stereo=obs, run="cr2284-high/x")
    assert c["stereo_sc"] == "A" and c["gcs"]["lon_stonyhurst"] == 28.0 and c["gcs"]["lat_deg"] == -14.0
    assert [r["instrument"] for r in c["radio_images"]] == ["nrh", "ovro_lwa"]
    assert c["wl_window"] == ["2024-05-09T16:32:00", "2024-05-09T20:02:00"]
    old = new_config("2013-10-25T15:00", "S06E69", earth_stereo=obs, run="r")
    assert old["stereo_sc"] == "B" and [r["instrument"] for r in old["radio_images"]] == ["nrh"]
