<p align="center"><img src="docs/logo.png" width="180" alt="corona3d logo"></p>

# corona3d

> **corona3d is private** (pre-release 0.0.1). It is ready for use and for testing, and its results can still change between versions. Please report problems, wrong results, and suggestions on the [issue tracker](https://github.com/jovian-explorer/corona3d/issues/new/choose). To cite corona3d, use [CITATION.cff](CITATION.cff); a DOI will be added when the repository is made public.

Coronal mass ejections (CMEs) drive shocks through a corona whose density, magnetic field, and Alfvén speed vary strongly with height, latitude, and longitude, and the radio bursts that these shocks excite, such as the herringbones of type II bursts, are emitted where the shock meets plasma of a particular density and field geometry. Whether a shock becomes super-Alfvénic, and whether it accelerates electrons along closed streamer loops or along open field, therefore depends on the three-dimensional structure of the corona at the moment the CME passes, which no single image can show. Morosan et al. ([2022](https://doi.org/10.1051/0004-6361/202244432)) combined a magnetohydrodynamic (MHD) model of the corona, a graduated cylindrical shell (GCS) of the CME, and radio-source positions in one three-dimensional reconstruction seen from three viewpoints (their Fig. 7), and showed that the herringbones originate where the CME flank meets the closed loops of a streamer. corona3d builds this reconstruction from public data alone, with every layer switchable, so that the same analysis can be repeated, examined layer by layer, and extended with coronagraph images and with radio sources from other observations.

## What can I do with it?

* **Draw the corona at a plasma level**: the iso-surface of the electron density of a Predictive Science Inc. (PSI) MAS/MAST model at any chosen density, for example the plasma level of a radio frequency, coloured by the Alfvén speed.
* **Trace the coronal magnetic field**: open and closed field lines from a dense set of footpoints, and the closed loops of a streamer from a seed line along the flank of the CME.
* **Place the CME in three dimensions**: the GCS shell of the CME (Thernisien [2011](https://doi.org/10.1088/0067-0049/194/2/33)), built with PyThea (Kouloumvakos et al. [2022](https://doi.org/10.3389/fspas.2022.974137)) from its height, aspect ratio, half-angle, direction, and tilt.
* **Add radio sources**: sky-plane positions of radio sources, each placed at the height that a four-fold Newkirk ([1961](https://doi.org/10.1086/147104)) density model gives for its frequency.
* **Overlay coronagraph images**: base-difference images of the Large Angle and Spectrometric Coronagraph (LASCO) C2 and of the STEREO-B COR1 and COR2 coronagraphs, placed in the plane of the sky of the spacecraft that recorded them, so that the model and the observed CME can be compared directly.
* **Switch every layer on or off** from the command line or from Python, and change the event, the model run, or any parameter in one configuration file.

## Quick start

1. Install corona3d with pip (Python 3.10 or newer; details in [How to Install and Run](#-how-to-install-and-run)):
   ```
   pip install "git+https://github.com/jovian-explorer/corona3d.git"
   ```
2. Create an empty working directory, enter it, write the configuration of the event of 2013 October 25, and download its public data (about 130 MB):
   ```
   mkdir corona3d_work
   cd corona3d_work
   corona3d init
   corona3d fetch
   ```
3. Draw the figure, and the same figure with the coronagraph images:
   ```
   corona3d render
   corona3d render --c2 --cor1 --cor2 --out figures/corona3d_wl.png
   ```

The figures go to `figures/`, together with the image of each panel (`view_earth.png`, `view_north.png`, `view_stb.png`). The [user guide](docs/USER_GUIDE.md) describes every layer, option, and parameter, and [METHODS.md](METHODS.md) describes how each element of the reconstruction is computed.

### Worked example: the CME of 2013 October 25

The default configuration is the event studied by Morosan et al. (2022): a CME launched from the eastern hemisphere at about 15:00 UT and observed from Earth and from STEREO-B, 142° apart, during which the Nançay Radioheliograph recorded herringbones at 173 and 228 MHz. On an ordinary computer, `corona3d render` takes about 45 s and draws the electron-density surface at 9 × 10⁷ cm⁻³ (the plasma level of 173 MHz emitted at the harmonic) coloured by the Alfvén speed between 100 and 700 km s⁻¹, about 4,500 open and closed field lines, the streamer loops above the CME flank in purple, the 1.2 R⊙ reference sphere, and the GCS shell in magenta, seen from Earth, from above the solar north pole, and from STEREO-B.

<p align="center"><img src="docs/corona3d_2013-10-25.png" width="100%" alt="The reconstruction of 2013 October 25 seen from Earth, the solar north pole, and STEREO-B"></p>

With `--c2 --cor1 --cor2` the LASCO C2 image of 15:24 UT appears in the Earth view and the COR1 and COR2 images of STEREO-B in the STEREO-B view, each minus its first image of the day. The CME front lies on the east limb in C2 and on the west limb in COR2-B, on the same side as the GCS shell in each view.

<p align="center"><img src="docs/corona3d_2013-10-25_wl.png" width="100%" alt="The same reconstruction with LASCO C2 and STEREO-B COR1 and COR2 base-difference images"></p>

The public model differs from the run used by Morosan et al. (2022), which PSI no longer serves, so that the surface and the field lines agree with their figure in outline but not in every detail; the user guide lists the differences.

---

## 📡 Data Sources & Provenance

All data are public and are downloaded by `corona3d fetch` into the working directory.

| Data | Source | Used for |
|---|---|---|
| MAS/MAST thermodynamic MHD run `cr2143-medium/hmi_mast_mas_std_0201` (density, magnetic field, Alfvén speed) | [PSI model archive](https://www.predsci.com/data/runs/) (Lionello et al. [2009](https://doi.org/10.1088/0004-637X/690/1/902); Mikić et al. [2018](https://doi.org/10.1038/s41550-018-0562-5)) | density surface, Alfvén speed, field lines |
| SOHO/LASCO C2 total-brightness images, level 0.5 | Virtual Solar Observatory through SunPy (Brueckner et al. [1995](https://doi.org/10.1007/BF00733434)) | Earth-view overlay |
| STEREO-B/SECCHI COR1 polarisation triplets and COR2 total-brightness images, level 0.5 | Virtual Solar Observatory through SunPy (Howard et al. [2008](https://doi.org/10.1007/s11214-008-9341-4)) | STEREO-B-view overlay |
| GCS parameters of the CME | fitted to LASCO C2 and COR1-A and -B at 15:24 UT, and to SDO/AIA 193 Å and STEREO/EUVI 195 Å at 15:00 UT (METHODS.md) | CME shell |
| Radio-source positions | supplied by the user | radio layer |

The radio images of the Nançay Radioheliograph are not downloaded: the archive of the Paris Observatory (BASS2000) delivers them on request as visibilities that must be imaged with SolarSoft, so that corona3d takes the positions of the sources from the user (user guide, Section 5). The references are collected in [docs/refs.bib](docs/refs.bib).

---

## 🚀 How to Install and Run

### Install with pip (Python 3.10+)

corona3d is not on PyPI; pip installs it from GitHub (access to the private repository is needed):

```
pip install "git+https://github.com/jovian-explorer/corona3d.git"
```

To work on the code, or to run the tests, install it from a copy of the repository:

```
git clone https://github.com/jovian-explorer/corona3d.git
cd corona3d
pip install -e ".[test]"
```

The dependencies are NumPy, SciPy, Matplotlib, PyVista (Sullivan & Kaszynski [2019](https://doi.org/10.21105/joss.01450)), pyhdf, Astropy, SunPy ([2020](https://doi.org/10.3847/1538-4357/ab4f7a)), and PyThea; pip installs all of them.

### Command-line options

```
corona3d init                                   # write config.json (the event of 2013-10-25) in the working directory
corona3d fetch                                  # download the MAS/MAST run and the white-light images
corona3d fetch --model                          # only the model run (or --wl: only the images)
corona3d render                                 # the three-panel figure with the default layers
corona3d render --c2 --cor1 --cor2              # with the coronagraph images
corona3d render --no-open --no-closed           # without the open and closed field lines
corona3d render --radio --radio-sources '[{"freq": 173, "x": -1.1, "y": -0.1}]'
corona3d render --half 3 --wl-time 2013-10-25T15:12 --out figures/wide.png
corona3d --version
```

| Layer | Switch | Default |
|---|---|---|
| electron-density surface coloured by the Alfvén speed | `--surface` / `--no-surface` | on |
| 1.2 R⊙ reference sphere | `--sphere` / `--no-sphere` | on |
| open field lines | `--open` / `--no-open` | on |
| closed field lines | `--closed` / `--no-closed` | on |
| streamer loops above the CME flank | `--streamer` / `--no-streamer` | on |
| GCS shell of the CME | `--gcs` / `--no-gcs` | on |
| radio sources | `--radio` / `--no-radio` | off |
| LASCO C2 (Earth view) | `--c2` / `--no-c2` | off |
| STEREO-B COR1 (STEREO-B view) | `--cor1` / `--no-cor1` | off |
| STEREO-B COR2 (STEREO-B view) | `--cor2` / `--no-cor2` | off |

From Python, the same figure is drawn with

```python
from corona3d import make_figure
make_figure("figures/corona3d_wl.png", data_root="corona3d_work", c2=True, cor1=True, cor2=True, open=False)
```

### Where corona3d keeps its files

Every command works inside one working directory, which holds `config.json`, the downloaded data (`data/`, with one folder for the model run and `data/wl/` for the images), and the figures (`figures/`). It is chosen with `--home`, with the environment variable `CORONA3D_HOME`, or, by default, as the current directory when that folder was created with `corona3d init`.

---

## Licence

corona3d is free software under the MIT licence (`LICENSE`).

---

## 👨‍🔬 Creator & Lead Architect

* **Keshav Aggarwal**
  * **Current Affiliation:** Student Visitor, Space Physics Laboratory (SPL), Vikram Sarabhai Space Centre (VSSC), Indian Space Research Organisation (ISRO), Thiruvananthapuram.
  * **Doctoral Education:** Ph.D. in Astronomy and Space Engineering, Department of Astronomy, Astrophysics and Space Engineering (DAASE), Indian Institute of Technology (IIT) Indore (Thesis Submitted).
  * **Fellowship:** Prime Minister's Research Fellow (PMRF).
  * **Research Focus:** Planetary/exoplanetary atmospheres, ionospheres, planetary radio occultation, solar-planetary space weather interactions across MOM, Chandrayaan-2, Chandrayaan-3, Akatsuki, and Aditya-L1.

---

## 🧪 Automated Test Suite

```
pip install -e ".[test]"
python -m pytest                          # offline checks
CORONA3D_ONLINE=1 python -m pytest        # also downloads the public data of 2013-10-25 and draws both figures (about 130 MB, 5-10 minutes); on Windows: set CORONA3D_ONLINE=1 first
```

The offline checks require that the field-line tracer follow the analytic field lines of a dipole (apex height within 0.02 R⊙, the dipole constant r/sin²θ constant to 1.5 %), that the radio heights reproduce the 1.576 and 1.449 R⊙ of Morosan et al. (2022) at 173 and 228 MHz, that a feature placed in a synthetic coronagraph image at any roll angle appear at its true position in three dimensions, and that the layers which need no data render. The workflow `.github/workflows/ci.yml` runs the offline checks on Windows, macOS, and Linux for every push to `main` and every pull request.
