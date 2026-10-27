<p align="center"><img src="docs/logo.png" width="180" alt="corona3d logo"></p>

# corona3d

> **corona3d is private** (pre-release 0.1.0). It is ready for use and for testing, although its results can still change between versions, and problems, wrong results, and suggestions are welcome on the [issue tracker](https://github.com/jovian-explorer/corona3d/issues/new/choose). To cite corona3d, please use [CITATION.cff](CITATION.cff), to which a DOI will be added when the repository is made public.

Coronal mass ejections (CMEs) drive shocks through a corona whose density, magnetic field, and Alfvén speed vary strongly with height, latitude, and longitude. The radio bursts that these shocks excite, such as the herringbones of type II bursts, are emitted where the shock meets plasma of a particular density and field geometry. Whether a shock becomes super-Alfvénic, and whether it accelerates electrons along closed streamer loops or along open field, therefore depends on the three-dimensional structure of the corona at the moment the CME passes, which no single image can show. Morosan et al. ([2022](https://doi.org/10.1051/0004-6361/202244432)) combined a magnetohydrodynamic (MHD) model of the corona, a graduated cylindrical shell (GCS) of the CME, and radio-source positions in one three-dimensional reconstruction seen from three viewpoints (their Fig. 7), and from it they concluded that the herringbones originate where the CME flank meets the closed loops of a streamer. corona3d follows the method described in Morosan et al. (2022), and the layout of their Fig. 7, to build such a reconstruction from a public MHD run, public coronagraph images, and radio images of the Nançay Radioheliograph (NRH; Kerdraon & Delouis [1997](https://doi.org/10.1007/BFb0106458)), with every layer switchable and every viewpoint selectable. Since their model run is no longer served and their CME shell and radio positions were not published, the figure that corona3d draws follows their approach without reproducing their Fig. 7 exactly. It can therefore be repeated from the data, examined layer by layer, and applied to other events. The GCS parameters of the default event were fitted with the script in `scripts/`, which is distributed with the code but is not part of the installed package.

## What can I do with it?

* **Draw the corona at a plasma level**: the iso-surface of the electron density of a Predictive Science Inc. (PSI) MAS/MAST model at any chosen density, for example the plasma level of a radio frequency, coloured by the Alfvén speed.
* **Trace the coronal magnetic field**: open and closed field lines from a dense set of footpoints, traced through the whole model out to 30 R⊙, and the closed field lines that pass through the shell of the CME, which are the loops of the streamer that the CME meets.
* **Place the CME in three dimensions**: the GCS shell of the CME (Thernisien [2011](https://doi.org/10.1088/0067-0049/194/2/33)), built with PyThea (Kouloumvakos et al. [2022](https://doi.org/10.3389/fspas.2022.974137)) from its height, aspect ratio, half-angle, direction, and tilt, at the time of the model and at the time of the coronagraph images.
* **Add radio sources**: the sources of a radio burst measured in images of the NRH (150 to 450 MHz) or of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA, 32 to 87 MHz, from 2024 April onward), drawn as a chain of coloured blobs for each frequency, one blob per time step, at the height that a four-fold Newkirk ([1961](https://doi.org/10.1086/147104)) density model gives for that frequency, in the manner of the blobs in Fig. 7 of Morosan et al. (2022); sky-plane positions from other instruments can be added as well.
* **Overlay coronagraph images**: images of the Large Angle and Spectrometric Coronagraph (LASCO) C2 and of the STEREO-B COR1 and COR2 coronagraphs, divided by an image taken before the eruption, placed in the plane of the sky of the spacecraft that recorded them and drawn with an orthographic projection, so that the model, the shell, and the observed CME share one scale.
* **Choose the viewpoint of every panel**: Earth, STEREO-B, above the solar north pole, or any Carrington longitude and latitude with a roll angle.
* **Switch every layer on or off** from the command line or from Python, and change the event, the model run, or any parameter in one configuration file.

## Quick start

1. Install corona3d with pip (Python 3.12 or newer; details in [How to Install and Run](#-how-to-install-and-run)):
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

The figures are written to `figures/` in the working directory, together with the image of each panel (`view_a.png`, `view_b.png`, `view_c.png`). Moreover, the [user guide](docs/USER_GUIDE.md) describes every layer, option, and parameter, [METHODS.md](METHODS.md) describes how each element of the reconstruction is computed, and [docs/MODEL_AUDIT.md](docs/MODEL_AUDIT.md) lists every method and number with its source and its test.

### Worked example: the CME of 2013 October 25

The default configuration is the event studied by Morosan et al. (2022), in which a CME was launched from the eastern hemisphere at about 15:00 UT and observed from Earth and from STEREO-B, 142° apart, while the Nançay Radioheliograph recorded herringbones at 173 and 228 MHz. On an ordinary computer, `corona3d render` takes about 40 s and draws this event as seen from Earth, from above the solar north pole, and from STEREO-B. The electron-density surface lies at 9 × 10⁷ cm⁻³, the level that Morosan et al. used for the 173 MHz harmonic emission, and is coloured by the Alfvén speed between 100 and 700 km s⁻¹. The open and closed field lines are traced from 4,500 footpoints, of which 359 lines are open and 4,141 closed in this model. In particular, the closed lines that pass through the CME shell are drawn in purple, and they number 478 of the 943 lines through the shell. Together with these layers, the figure shows the 1.2 R⊙ reference sphere and the GCS shell at 15:00 UT in magenta, and each panel states the time of the model and of the shell.

<p align="center"><img src="docs/corona3d_2013-10-25.png" width="100%" alt="The reconstruction of 2013 October 25 seen from Earth, the solar north pole, and STEREO-B"></p>

With `--c2 --cor1 --cor2`, the LASCO C2 image of 15:24 UT appears in the Earth view and the COR1 and COR2 images of STEREO-B appear in the STEREO-B view, each divided by the first image of the download window (14:35 to 14:39 UT). At the same time, the shell is drawn at 15:24 UT, the time of the images, and the bright CME is observed inside the shell on the east limb in C2 and on the west limb in COR2-B.

<p align="center"><img src="docs/corona3d_2013-10-25_wl.png" width="100%" alt="The same reconstruction with LASCO C2 and STEREO-B COR1 and COR2 images"></p>

With `--radio` and the NRH images of the event copied into `data/nrh/`, the sources of the burst are added as coloured blobs (user guide, Section 5.7). The images available for the afternoon of 2013 October 25 are those at 150.9 and 432 MHz, and the burst is bright in them from 14:58:55 UT. Between 14:59 and 15:03 UT, one image per minute gives five blobs at 150.9 MHz, coloured from magenta to yellow, and three at 432 MHz, coloured from blue to green. They lie at 1.647 and 1.221 R⊙, respectively, and the 150.9 MHz sources lie on the eastern flank of the CME shell. The frequencies differ from the 173 and 228 MHz of Morosan et al. (2022), so that the blobs follow their figure in kind but not in position.

<p align="center"><img src="docs/corona3d_2013-10-25_radio.png" width="100%" alt="The reconstruction with the NRH radio sources at 150.9 and 432 MHz as coloured blobs"></p>

However, the public model differs from the run used by Morosan et al. (2022), which PSI no longer serves, so that the surface and the field lines agree with their figure in outline but not in every detail. The default figure therefore does not test their conclusion independently, and the user guide lists the differences.

---

## 📡 Data Sources & Provenance

The model run and the images are public, and `corona3d fetch` downloads them into the working directory.

| Data | Source | Used for |
|---|---|---|
| MAS/MAST thermodynamic MHD run `cr2143-medium/hmi_mast_mas_std_0201` (density and magnetic field; the mesh of the Alfvén-speed file) | [PSI model archive](https://www.predsci.com/data/runs/) (Lionello et al. [2009](https://doi.org/10.1088/0004-637X/690/1/902); Mikić et al. [2018](https://doi.org/10.1038/s41550-018-0562-5)) | density surface, Alfvén speed, field lines |
| SOHO/LASCO C2 images with the clear polariser (orange filter), level 0.5 | Virtual Solar Observatory through SunPy (Brueckner et al. [1995](https://doi.org/10.1007/BF00733434)) | Earth-view overlay, GCS fit |
| STEREO-B/SECCHI COR1 polarisation triplets and COR2 total-brightness images, level 0.5; STEREO-A COR1 with `--cor1a` | Virtual Solar Observatory through SunPy (Howard et al. [2008](https://doi.org/10.1007/s11214-008-9341-4)) | STEREO-B-view overlay, GCS fit |
| GCS parameters of the CME at 15:24 and 15:00 UT | fitted to LASCO C2 and COR1-A and -B with `corona3d gcs` (METHODS.md, Section 7) | CME shell |
| NRH2 images (Stokes I, about 10 s) at 150.9 and 432 MHz | BASS2000 archive of the Paris Observatory, delivered on request; copied by the user into `data/nrh/` | radio layer |
| OVRO-LWA spectral images (Stokes I, 10 s, 12 bands at 34 to 85 MHz), from 2024 April | [OVRO-LWA solar data](https://ovsa.njit.edu/docs/pages/ovro-lwa-solar-data-products/), public bucket on Amazon Web Services (CC BY 4.0); `corona3d fetch --ovro-lwa` | radio layer |
| Other radio-source positions (optional) | supplied by the user | radio layer |

The NRH images are not downloaded by `corona3d fetch`, since the archive of the Paris Observatory (BASS2000) delivers them only on request; the image files (`nrh2_*_i.fts`) are copied into `data/nrh/` of the working directory, where the radio layer reads them (user guide, Section 5.7). All references are collected in [docs/refs.bib](docs/refs.bib).

---

## 🚀 How to Install and Run

### Install with pip (Python 3.12+)

Since corona3d is not on PyPI, pip installs it from GitHub, which requires access to the private repository:

```
pip install "git+https://github.com/jovian-explorer/corona3d.git"
```

To work on the code, to run the tests, or to run the GCS fit, install it from a copy of the repository:

```
git clone https://github.com/jovian-explorer/corona3d.git
cd corona3d
pip install -e ".[test,fit]"
```

The dependencies are NumPy, SciPy, Matplotlib, PyVista (Sullivan & Kaszynski [2019](https://doi.org/10.21105/joss.01450)), pyhdf, Astropy, SunPy ([2020](https://doi.org/10.3847/1538-4357/ab4f7a)), and PyThea 1.3, which needs Python 3.12 or newer, and pip installs all of them. On macOS computers for which pyhdf provides no wheel, the HDF4 library must be installed first (`brew install hdf4`).

### Command-line options

```
corona3d init                                   # write config.json (the event of 2013-10-25) in the working directory
corona3d fetch                                  # download the MAS/MAST run and the white-light images
corona3d fetch --model                          # only the model run (or --wl: only the images; --cor1a: also COR1-A)
corona3d fetch --ovro-lwa                       # only the OVRO-LWA images of the radio_images entries
corona3d render                                 # the three-panel figure with the default layers
corona3d render --c2 --cor1 --cor2              # with the coronagraph images
corona3d render --no-open --no-closed           # without the open and closed field lines
corona3d render --radio                         # with the radio sources of the NRH images in data/nrh/
corona3d render --view b=300,20 --view c=north  # panel b from Carrington longitude 300, latitude 20; panel c from the north pole
corona3d render --half 3 --wl-time 2013-10-25T15:12 --out figures/wide.png
corona3d --version
```

| Layer | Switch | Default |
|---|---|---|
| electron-density surface coloured by the Alfvén speed | `--surface` / `--no-surface` | on |
| 1.2 R⊙ reference sphere | `--sphere` / `--no-sphere` | on |
| open field lines | `--open` / `--no-open` | on |
| closed field lines | `--closed` / `--no-closed` | on |
| closed field lines through the CME shell (streamer loops) | `--streamer` / `--no-streamer` | on |
| GCS shell of the CME | `--gcs` / `--no-gcs` | on |
| radio sources | `--radio` / `--no-radio` | off |
| LASCO C2 (Earth view) | `--c2` / `--no-c2` | off |
| STEREO-B COR1 (STEREO-B view) | `--cor1` / `--no-cor1` | off |
| STEREO-B COR2 (STEREO-B view) | `--cor2` / `--no-cor2` | off |

Radio sources are best written under `radio_sources` in `config.json` (user guide, Section 5.7), since the shells of Windows treat quotes inside a command-line argument differently from those of macOS and Linux.

From Python, the same figure is drawn with

```python
from corona3d import make_figure
make_figure("figures/corona3d_wl.png", data_root="corona3d_work", c2=True, cor1=True, cor2=True, open=False)
```

which reads `config.json` of the working directory and writes the figure inside it.

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
python -m pytest                          # offline checks (about 40 s)
```

To include the online test, which downloads the public data of 2013 October 25 (about 130 MB) and draws both figures of this README in 5 to 10 minutes, set `CORONA3D_ONLINE=1` before `python -m pytest` (`set CORONA3D_ONLINE=1` in the Windows Command Prompt, `$env:CORONA3D_ONLINE="1"` in PowerShell, `export CORONA3D_ONLINE=1` in macOS and Linux).

The offline checks use synthetic data with known answers. In particular, the field-line tracer must follow the analytic field lines of a dipole with its default step, and the MAS units and the Alfvén speed must reproduce their analytic values. Moreover, a feature placed in synthetic coronagraph images at five roll angles must appear at its true position in three dimensions, inputs that would draw a wrong or blank image must be refused, and halving `--half` must double the size of the Sun in the panel. The full list of checks is given in [docs/MODEL_AUDIT.md](docs/MODEL_AUDIT.md). The workflow `.github/workflows/ci.yml` runs the offline checks on Windows, macOS, and Linux with Python 3.12 for every push to `main` and every pull request into `main`, although, while the repository is private, it runs only when the account has GitHub Actions minutes left.
