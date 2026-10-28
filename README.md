<p align="center"><img src="docs/logo.png" width="180" alt="corona3d logo"></p>

# corona3d

> **corona3d is private** (pre-release 0.1.0). It is ready for use and for testing, although its results can still change between versions, and problems, wrong results, and suggestions are welcome on the [issue tracker](https://github.com/jovian-explorer/corona3d/issues/new/choose). To cite corona3d, please use [CITATION.cff](CITATION.cff), to which a DOI will be added when the repository is made public.

## Purpose

Coronal mass ejections (CMEs) drive shocks through a corona whose density, magnetic field, and Alfvén speed vary strongly with height, latitude, and longitude. The radio bursts that these shocks excite, such as the herringbones of type II bursts, are emitted where the shock meets plasma of a particular density and field geometry. Whether a shock becomes super-Alfvénic, and whether it accelerates electrons along closed streamer loops or along open field, therefore depends on the three-dimensional structure of the corona at the moment the CME passes, which no single image can show. Morosan et al. ([2022](https://doi.org/10.1051/0004-6361/202244432)) combined a magnetohydrodynamic (MHD) model of the corona, a graduated cylindrical shell (GCS) of the CME, and radio-source positions in one three-dimensional reconstruction seen from three viewpoints (their Fig. 7), but such a reconstruction remains laborious to repeat for another event. In this study, we present corona3d, a Python package that follows the method described in Morosan et al. (2022) and builds the reconstruction for an event from its time and flare position, using a public MHD run, public coronagraph images, and radio images of the Nançay Radioheliograph (NRH; Kerdraon & Delouis [1997](https://doi.org/10.1007/BFb0106458)) or of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA; Chhabra et al. [2021](https://doi.org/10.3847/1538-4357/abc94b)). Every layer can be switched on or off, every viewpoint can be chosen, and every layer can be exported as arrays or tables for further analysis. Since the model run of Morosan et al. (2022) is no longer served and their CME shell and radio positions were not published, the figure that corona3d draws follows their approach and the layout of their Fig. 7 without reproducing it exactly.

## Install

corona3d needs Python 3.12 or newer, because PyThea 1.3, which builds the GCS shell, requires it, and it runs on Windows, macOS, and Linux. Once corona3d is released on PyPI, one command installs it with every optional part:

```
pip install "corona3d[all]"
```

Until then, pip installs it from the repository, which is private at present, so the account used by git must have access to it:

```
pip install "corona3d[all] @ git+https://github.com/jovian-explorer/corona3d.git"
```

The plain package is the light core (downloads, event set-up, radio interval, `corona3d doctor`). The extras `render` (PyVista and PyThea, for the figures), `gcs` (PyThea and scikit-image, for the shell fit), and `radio` (h5py, for OVRO-LWA) add the heavy packages, `all` adds the three, and a command that needs a missing one names the extra to install. `corona3d --version` prints the version, and `corona3d doctor` reports the extras, the working directory (by default the platform cache directory, changed with `--home` or `CORONA3D_HOME`), and which archives answer. On macOS with an Intel processor, the HDF4 library must be installed first with `brew install hdf4`, and `conda install -c conda-forge pyhdf` is an alternative on any system. The user guide ([Section 1](docs/USER_GUIDE.md#1-install-and-start)) gives the details.

## Quick start: one command for an event

```
mkdir event_work
cd event_work
corona3d event --time 2024-05-09T17:32 --flare S14W28
```

The command writes the configuration of the event from its time and the flare position (14° south and 28° west of the disc centre), downloads the PSI MAS/MAST run and the white-light images, and, from 2024 April onward, the OVRO-LWA images, chooses the interval of the radio burst from the images, fits the GCS shell to the coronagraph images, draws the figures, and writes `event_report.md`. The report separates what was measured from what was assumed, lists the missing data, and gives the quality of the fit. When the free fit of the shell fails or is poor, the fit with the direction held at the flare site follows, and the report says so. The 2-D images of the NRH cannot be downloaded automatically (see the data sources below), and the report then names the manual step. The download takes a few minutes and the fit several more, whereas a figure takes about 40 s.

The step-by-step commands, for the event of 2013 October 25 studied by Morosan et al. (2022), are

```
mkdir corona3d_work
cd corona3d_work
corona3d init                  # writes config.json of that event
corona3d fetch                 # model run and white-light images (about 130 MB)
corona3d render                # figures/corona3d.png
corona3d render --c2 --cor1 --cor2 --out figures/corona3d_wl.png
```

The figures are written to `figures/` together with the image of each panel (`view_a.png`, `view_b.png`, `view_c.png`). Moreover, the [user guide](docs/USER_GUIDE.md) describes every command, option, and parameter, [docs/METHODS.md](docs/METHODS.md) derives how each element of the reconstruction is computed, [docs/LIMITATIONS.md](docs/LIMITATIONS.md) states where the results cannot be trusted, and [docs/MODEL_AUDIT.md](docs/MODEL_AUDIT.md) lists every method and number with its source and its test.

### Worked example: the CME of 2013 October 25

The default configuration is the event studied by Morosan et al. (2022), in which a CME was launched from the eastern hemisphere at about 15:00 UT and observed from Earth and from STEREO-B, 142° apart, while the NRH recorded herringbones at 173 and 228 MHz. The command `corona3d render` draws this event as seen from Earth, from above the solar north pole, and from STEREO-B. The electron-density surface lies at 9 × 10⁷ cm⁻³, the level that Morosan et al. used for the 173 MHz harmonic emission, and is coloured by the Alfvén speed between 100 and 700 km s⁻¹. The open and closed field lines are traced from 4,500 footpoints, and the closed lines that pass through the CME shell, the loops of the streamer that the CME meets, are drawn in purple.

<p align="center"><img src="docs/corona3d_2013-10-25.png" width="100%" alt="The reconstruction of 2013 October 25 seen from Earth, the solar north pole, and STEREO-B"></p>

With `--c2 --cor1 --cor2`, the LASCO C2 image of 15:24 UT appears in the Earth view and the COR1 and COR2 images of STEREO-B appear in the STEREO-B view, each divided by the first image of the download window (14:35 to 14:39 UT). The shell is then drawn at 15:24 UT, the time of the images.

<p align="center"><img src="docs/corona3d_2013-10-25_wl.png" width="100%" alt="The same reconstruction with LASCO C2 and STEREO-B COR1 and COR2 images"></p>

With `--radio` and the NRH images of the event copied into `data/nrh/`, the sources of the burst are added as coloured blobs (user guide, Section 5.7). The images available for the afternoon of that day are those at 150.9 and 432 MHz, and between 14:59 and 15:03 UT one image per minute gives five blobs at 150.9 MHz and three at 432 MHz, at 1.647 and 1.221 R⊙. These frequencies differ from the 173 and 228 MHz of Morosan et al. (2022), so that the blobs follow their figure in kind but not in frequency.

<p align="center"><img src="docs/corona3d_2013-10-25_radio.png" width="100%" alt="The reconstruction with the NRH radio sources at 150.9 and 432 MHz as coloured blobs"></p>

## What is measured and what is assumed

| Element | Measured from data | Assumed or chosen |
|---|---|---|
| Corona before the CME | none | the steady-state PSI MAS/MAST run of the Carrington rotation; the density multiplied by 2 after Morosan et al. (2022) |
| Density surface and Alfvén speed | none | the density level `ne_level` (9 × 10⁷ cm⁻³); the Alfvén speed from the scaled density |
| Field lines and streamer loops | none | traced in the static model; loops selected as the closed lines through the shell |
| CME shell | six GCS parameters, fitted to the envelopes in LASCO C2 and COR1 or COR2 at the time of the images; apex height and speed extrapolated from the front | the direction held at the flare site when the free fit fails; an apex height given by hand; a constant speed between the frames |
| Coronagraph images | brightness ratio or difference against the first frame of the window | the image is placed in the plane of the sky through Sun centre |
| Radio sources | sky-plane centroids of the NRH or OVRO-LWA images above `min_tb` | the height of the four-fold Newkirk model at the second harmonic, on the side facing the observer; no correction for scattering |
| Observers | Carrington positions of Earth and STEREO from JPL Horizons through SunPy | none |

[docs/LIMITATIONS.md](docs/LIMITATIONS.md) describes the cases in which these inputs fail, namely halo CMEs, data gaps, a small separation between Earth and STEREO, assumed heights, the manual step for the NRH, and the resolution of the model.

## Data sources

The model run and the images are public, and `corona3d fetch` downloads them into the working directory.

| Data | Source and access | Used for |
|---|---|---|
| MAS/MAST thermodynamic MHD run (density and magnetic field; the mesh of the Alfvén-speed file), by default `cr2143-medium/hmi_mast_mas_std_0201` | [PSI model archive](https://www.predsci.com/data/runs/), public (Mikić et al. [1999](https://doi.org/10.1063/1.873474); Lionello et al. [2009](https://doi.org/10.1088/0004-637X/690/1/902); Mikić et al. [2018](https://doi.org/10.1038/s41550-018-0562-5)) | density surface, Alfvén speed, field lines |
| SOHO/LASCO C2 images, level 0.5 | Virtual Solar Observatory through SunPy, public (Brueckner et al. [1995](https://doi.org/10.1007/BF00733434)) | Earth-view overlay, GCS fit |
| STEREO-A or STEREO-B SECCHI COR1 polarisation triplets and COR2 total-brightness images, level 0.5 | Virtual Solar Observatory through SunPy, public (Howard et al. [2008](https://doi.org/10.1007/s11214-008-9341-4)); `--cor1a` adds STEREO-A COR1 | STEREO-view overlay, GCS fit |
| NRH 2-D images (Stokes I, about 10 s) | **On request only.** BASS2000 (Paris Observatory) holds only four snapshots per day at 164 and 327 MHz, movies, and visibility files; the 10 s images (`nrh2_*_i.fts`) are delivered by the Radio Solar Database at Nançay, [https://rsdb.obs-nancay.fr/](https://rsdb.obs-nancay.fr/), after a request there. The user copies the files into `data/nrh/`, and `corona3d fetch --nrh` reports which frequencies they cover | radio layer |
| OVRO-LWA spectral images (Stokes I, 10 s, 12 bands at 34 to 85 MHz), from 2024 April | [OVRO-LWA solar data](https://ovsa.njit.edu/docs/pages/ovro-lwa-solar-data-products/), public bucket on Amazon Web Services (CC BY 4.0); `corona3d fetch --ovro-lwa` | radio layer |
| Other radio-source positions (optional) | supplied by the user | radio layer |

Work that uses these data should acknowledge the PSI model runs, the LASCO and SECCHI teams, the Nançay Radioheliograph, and OVRO-LWA with the references above. All references are collected in [docs/refs.bib](docs/refs.bib).

## Using corona3d from Python

`corona3d.build(home, **layers)` returns the layers of a working directory as NumPy arrays in a `Reconstruction`, and `corona3d export` writes them as `npz`, `csv`, or VTK files together with `metadata.json`, which records the version, the configuration, and the sources of every input (user guide, Section 8).

## Citing

If you use corona3d, please cite it with [CITATION.cff](CITATION.cff), and cite Morosan et al. (2022) for the method that it follows, Thernisien et al. ([2006](https://doi.org/10.1086/508254)) and Thernisien ([2011](https://doi.org/10.1088/0067-0049/194/2/33)) for the GCS model, Kouloumvakos et al. ([2022](https://doi.org/10.3389/fspas.2022.974137)) for PyThea, and the references of the data above.

## Licence

corona3d is free software under the MIT licence (`LICENSE`). Copyright 2026 Keshav Aggarwal.

## Author

Keshav Aggarwal, Student Visitor, Space Physics Laboratory (SPL), Vikram Sarabhai Space Centre (VSSC), Indian Space Research Organisation (ISRO), Thiruvananthapuram.

## Tests

```
pip install -e ".[all,test]"
python -m pytest                          # offline checks (about 1 min)
CORONA3D_ONLINE=1 python -m pytest        # also downloads the data of 2013 October 25
```
