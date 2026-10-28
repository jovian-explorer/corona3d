"""Installation: the version and doctor commands, the default home, lazy optional imports, package data."""
import pytest

import corona3d
from corona3d import _extras, paths
from corona3d.cli import main


def test_version_command(capsys):
    with pytest.raises(SystemExit) as e:
        main(["--version"])
    assert e.value.code == 0 and f"corona3d {corona3d.__version__}" in capsys.readouterr().out


def test_doctor_offline_reports_extras_and_home(tmp_path, capsys):
    with pytest.raises(SystemExit) as e:
        main(["--home", str(tmp_path), "doctor", "--offline"])
    out = capsys.readouterr().out
    assert e.value.code == 0
    for word in ("Core packages", "[render]", "[gcs]", "[radio]", str(tmp_path), "--offline", "2024-05-09"):
        assert word in out


def test_home_order_arg_env_cwd_default(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "default_home", lambda: tmp_path / "cache")
    monkeypatch.delenv("CORONA3D_HOME", raising=False)
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    assert paths.home_source() == ((tmp_path / "cache").resolve(), "platform cache directory")
    (cwd / "config.json").write_text("{}")
    assert paths.home_source()[1].startswith("current folder")
    monkeypatch.setenv("CORONA3D_HOME", str(tmp_path / "env"))
    assert paths.home_source()[1] == "CORONA3D_HOME"
    assert paths.home_source(str(tmp_path / "arg"))[1] == "--home"


def test_default_home_is_a_platform_cache_directory_outside_the_working_folder():
    assert "corona3d" in str(paths.default_home()).lower()


def test_missing_extra_names_the_extra(monkeypatch):
    real = _extras.import_module

    def fake(name):
        if name.split(".")[0] in ("pyvista", "h5py"):
            raise ImportError(name)
        return real(name)
    monkeypatch.setattr(_extras, "import_module", fake)
    with pytest.raises(_extras.MissingExtra, match=r"corona3d\[radio\]"):
        _extras.need("radio", "h5py", "reading OVRO-LWA images")
    with pytest.raises(_extras.MissingExtra, match=r"corona3d\[render\]"):
        _extras.need("render", "pyvista", "the drawing")


def test_package_data_is_found_through_the_package():
    assert "2024-05-09" in paths.events() and paths.default_config()["run"]
    assert paths.event("2024-05-09")["time"].startswith("2024-05-09")
