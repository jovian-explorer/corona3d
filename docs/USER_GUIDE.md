# corona3d user guide

This guide describes corona3d as it is in this repository: how to install it and prepare a working directory, how the public data are downloaded, what each layer of the reconstruction shows and how it is switched on or off, how the viewpoint of each panel is chosen, how every parameter of the configuration acts on the figure, how corona3d is used from Python, how far its results can be trusted, and what to do when something does not work. Moreover, the computation behind each layer is described with its sources in [METHODS.md](METHODS.md), every method and number is listed with its source and its test in [MODEL_AUDIT.md](MODEL_AUDIT.md), and the changes between versions are recorded in [CHANGELOG.md](../CHANGELOG.md). Where the results cannot be trusted is stated in [LIMITATIONS.md](LIMITATIONS.md).

## At a glance

| Task | Command |
|---|---|
| Everything for one event, from its time and flare position | `corona3d event --time T --flare S14W28` (with a date only, the bursts of the day are listed and `--burst N` chooses one) |
| The default event of Morosan et al. (2022) | `corona3d init`, `corona3d fetch`, `corona3d render` |
| Download the data only | `corona3d fetch [--model] [--wl] [--cor1a] [--ovro-lwa] [--nrh]` |
| Fit the CME shell | `corona3d gcs --t-ref T --t-out T --flare LON LAT` |
| Choose the radio interval from the images | `corona3d radio-interval` |
| Draw figures | `corona3d render [--<layer>] [--view PANEL=VIEW]` |
| Write the layers as data files | `corona3d export [--format npz, csv, or vtk]` |

The sections below describe each command and its options, Section 7 lists the keys of `config.json`, and Section 8 describes the Python interface and the export formats with their metadata.

Contents

1. [Install and start](#1-install-and-start)
2. [The working directory](#2-the-working-directory)
3. [Downloading the data](#3-downloading-the-data)
4. [Drawing the figure](#4-drawing-the-figure)
5. [The layers](#5-the-layers)
6. [The viewpoints](#6-the-viewpoints)
7. [The configuration file](#7-the-configuration-file)
8. [Using corona3d from Python](#8-using-corona3d-from-python)
9. [Another event or another model run](#9-another-event-or-another-model-run)
10. [How the results were checked](#10-how-the-results-were-checked)
11. [Differences from Morosan et al. (2022) and limits of the results](#11-differences-from-morosan-et-al-2022-and-limits-of-the-results)
12. [Data sources, terms of use, and citation](#12-data-sources-terms-of-use-and-citation)
13. [Troubleshooting](#13-troubleshooting)

---

## 1. Install and start

corona3d runs on Windows, macOS, and Linux and needs Python 3.12 or newer, because PyThea 1.3, which builds the shell of the coronal mass ejection (CME), requires it. Once corona3d is released on PyPI, one command installs it with every optional part:

```
pip install "corona3d[all]"
```

Until then, pip installs it from the GitHub repository, which is private at present, so the account used by git must have access to it:

```
pip install "corona3d[all] @ git+https://github.com/jovian-explorer/corona3d.git"
```

The plain `pip install corona3d` is the light core. It installs NumPy, SciPy, Matplotlib, Astropy, SunPy, pyhdf for the HDF4 files of the magnetohydrodynamic (MHD) model, and platformdirs, and it runs the downloads, `corona3d new`, `corona3d fetch`, `corona3d radio-interval`, and `corona3d doctor`. The heavy packages are optional extras, imported only when a command needs them, and a command that needs a missing one stops with the name of the extra to install:

| Extra | Adds | Needed for |
|---|---|---|
| `render` | PyVista (with the Visualization Toolkit, VTK) and PyThea | `corona3d render`, `corona3d event`, `corona3d export`, `corona3d.build` |
| `gcs` | PyThea and scikit-image | the fit of the graduated cylindrical shell (GCS) of Section 9 |
| `radio` | h5py | the OVRO-LWA images |
| `all` | the three extras above | the complete workflow |
| `test` | pytest and Pillow | the test suite |
| `docs` | MkDocs | building the documentation as a web site |

To work on the code or run the tests, install a copy of the repository in editable mode:

```
git clone https://github.com/jovian-explorer/corona3d.git
cd corona3d
pip install -e ".[all,test]"
```

pyhdf provides wheels for Windows, Linux, and recent macOS on Apple processors. On other Macs, the HDF4 library must be installed first with `brew install hdf4`, while `conda install -c conda-forge pyhdf` is an alternative on any system. Since the renderings are drawn off screen, no window opens, and corona3d therefore also runs on a server, provided that VTK finds an OpenGL implementation (Section 13).

`corona3d doctor` reports the version, the installed core packages and extras, the working directory and where it came from, the package data, and whether the archives that corona3d uses (PSI, the Virtual Solar Observatory, the OVRO-LWA bucket, NOAA NCEI, and JPL Horizons) answer; `corona3d doctor --offline` skips the network check. Nothing in the installed package depends on the layout of the repository.

After the installation, the command `corona3d` is available in the terminal (`python -m corona3d` does the same), and `corona3d --version` prints the version.

## 2. The working directory

Every command of corona3d works inside one working directory, which holds the configuration of the event (`config.json`), the downloaded data (`data/`), and the figures (`figures/`). The directory is created with

```
mkdir corona3d_work
cd corona3d_work
corona3d init
```

`corona3d init` writes `config.json` with the event of 2013 October 25 studied by Morosan et al. (2022). An existing `config.json` is never overwritten, so that an edited configuration is kept when `init` is given again. The option `--event` selects another preset instead: `corona3d init --event 2024-05-09` writes the configuration of the CME and type II burst of 2024 May 9 (Section 9.1), and `corona3d init --help` lists the presets that are installed.

corona3d finds the working directory from the option `--home` given before the command (`corona3d --home D:\corona3d_work render`), then from the environment variable `CORONA3D_HOME`, then the current folder when it holds a `config.json`, and otherwise the platform cache directory (`%LOCALAPPDATA%\corona3d\Cache` on Windows, `~/.cache/corona3d` on Linux, `~/Library/Caches/corona3d` on macOS), so that nothing depends on the folder from which the command is given. `corona3d doctor` prints the directory in use and where it came from. Since one working directory holds one event, an event in the cache directory is replaced with `corona3d new --force`, and another event is better kept in its own directory with `--home`. When the chosen directory has no `config.json`, `fetch`, `gcs`, and `render` stop with the advice to run `corona3d init` first.

After a download and two renderings, the directory contains

```
corona3d_work/
  config.json
  data/
    cr2143-medium_hmi_mast_mas_std_0201/   rho002.hdf, br002.hdf, bt002.hdf, bp002.hdf, va002.hdf (53 MB)
    wl/                                     LASCO C2 and STEREO COR1 and COR2 images (about 75 MB)
    nrh/                                    NRH image files, copied by the user (Section 5.7)
  figures/
    corona3d.png                            the three-panel figure
    view_a.png, view_b.png, view_c.png      the three panels at full resolution
```

The panel images are written next to the figure and are replaced at every rendering into the same folder, so that renderings whose panels should be kept are best given separate folders with `--out`.

## 3. Downloading the data

```
corona3d fetch            # the model run and the white-light images
corona3d fetch --model    # only the model run
corona3d fetch --wl       # only the white-light images
corona3d fetch --cor1a    # as above, and also the STEREO-A COR1 images used by the GCS fit
corona3d fetch --nrh      # no download: reports which NRH files are present (see below)
```

The model run is downloaded file by file from the archive of Predictive Science Inc. (PSI), at `https://www.predsci.com/data/runs/<run>/corona/`. A file that is already present is not downloaded again, a new file is written under its final name only when it is complete, and the partial file of a failed download is removed, so that giving the command again downloads exactly the files that are missing.

The white-light images, namely the C2 images of the Large Angle and Spectrometric Coronagraph (LASCO) and the COR1 and COR2 images of the Solar Terrestrial Relations Observatory (STEREO) spacecraft named by `stereo_sc` in the configuration (B for the event of 2013, A for events after 2014 October), are searched for and downloaded through the Virtual Solar Observatory with SunPy for the time window `wl_window` of the configuration. The search also returns products that corona3d does not draw (polarised C2 images, the COR1 sequences `s5c1` and `n4c1`, and the polarised COR2 sequence `n4c2`), which are kept in `data/wl/` but are not read, except `s5c1` by the GCS fit. After the download, every image is read once, and an image that cannot be read completely, such as one truncated by an interrupted transfer, is deleted. When files fail or are deleted, the command says so and ends with exit code 1, and giving it again downloads the missing files. On an ordinary connection, the whole download (about 130 MB) takes a few minutes.

The 2-D images of the Nançay Radioheliograph (NRH) cannot be downloaded automatically. The BASS2000 archive of the Paris Observatory (`ftp://ftpbass2000.obspm.fr/pub/nancay/`) holds only four snapshots per day at 164 and 327 MHz, together with movies, and it holds the daily visibility files, from which images can be made only with the NRH routines of the IDL SolarSoft package. The 10 s images at the frequencies of the NRH (files `nrh2_FFFF_h60_YYYYMMDD_HHMMSScCC_i.fts`) are supplied by the Radio Solar Database at Nançay (`https://rsdb.obs-nancay.fr/`) after a request there. Therefore `corona3d fetch --nrh` downloads nothing. For every NRH entry of `radio_images`, it lists the frequencies whose files in `data/nrh/` overlap the interval of the entry and those that do not, and, for the missing ones, it repeats the manual step: request the Stokes I images for the frequencies and the interval from the Radio Solar Database, and copy the files into `data/nrh/`. The check is also made by `corona3d fetch` without options.

## 4. Drawing the figure

```
corona3d render
```

draws the three-panel figure into `figures/corona3d.png`. By default, panel a shows the reconstruction from Earth, panel b from above the solar north pole with the direction of the Earth at the bottom, and panel c from the STEREO spacecraft of the event, and Section 6 explains how to change them. The panels are drawn with an orthographic projection, so that every structure appears at its true size in R⊙ whatever its distance along the line of sight. The bottom of each panel states the model run and the time of the CME shell, and, where a coronagraph image is shown, the times of the image and of its base frame. With the default layers, the rendering takes about 40 s, most of which is spent on tracing the field lines. The options of `render` are:

| Option | Effect |
|---|---|
| `--<layer>` / `--no-<layer>` | switch a layer on or off (Section 5) |
| `--view PANEL=VIEW` | the viewpoint of a panel (Section 6); may be given for each panel |
| `--out PATH` | file of the figure, relative to the working directory (default `figures/corona3d.png`) |
| `--half H` | half-height of all three panels in R⊙; by default each panel has the framing of the configuration (about 1.5 R⊙), widened to 6 R⊙ in the panel of C2 or COR2 and to 4 R⊙ in the panel of COR1 alone |
| `--wl-time T` | time of the coronagraph frames, e.g. `2013-10-25T15:12`; the frame nearest to it within `wl_window` is taken for each instrument |
| `--radio-sources JSON` | radio sources for the radio layer (Section 5.7) |

`corona3d render --help` lists every option and, at its end, the names of the layers. Some combinations that are useful in practice are the following:

```
corona3d render --no-open --no-closed --out figures/surface.png          # the density surface, the streamer, and the CME only
corona3d render --no-surface --out figures/field.png                       # the field lines without the surface
corona3d render --c2 --cor1 --cor2 --out figures/corona3d_wl.png           # with the coronagraph images
corona3d render --cor1 --no-open --no-closed --out figures/cor1.png        # COR1 in panel c (widened to 4 R_sun); the other layers stay on
corona3d render --cor2 --half 15 --out figures/cor2_full.png               # the whole COR2 field, out to 15 R_sun, in all panels
corona3d render --c2 --wl-time 2013-10-25T15:36 --out figures/c2_1536.png  # a later C2 frame
corona3d render --view b=300,20 --view c=north --out figures/views.png     # other viewpoints for panels b and c
```

## 5. The layers

Each layer is switched on with `--<name>` and off with `--no-<name>`, and only the layers that are on are computed, so that a figure without field lines is drawn in a few seconds.

### 5.1 `surface`: the density surface coloured by the Alfvén speed (on)

This layer draws the surface on which the electron density of the MHD model equals `ne_level`, by default 9 × 10⁷ cm⁻³, the level that Morosan et al. (2022) used for the 173 MHz harmonic emission (the exact value for 173 MHz is 9.28 × 10⁷ cm⁻³). The colour gives the Alfvén speed on the surface over `va_range` (100 km s⁻¹ in blue to 700 km s⁻¹ in red), on the colour bar below the panels, and this colour bar is drawn only when the layer is on. By default, the Alfvén speed is computed from the doubled density of the surface, whereas `"va_density": "model"` computes it from the density of the model, which raises it by a factor of √2. The blue regions mark where a shock of a few hundred km s⁻¹ becomes strongly super-Alfvénic.

### 5.2 `sphere`: the reference sphere (on)

This layer draws a grey sphere of radius `sphere_r` = 1.2 R⊙, which provides a reference of height, since the density level lies below 1.2 R⊙ wherever the sphere is visible through the surface.

### 5.3 `open` and `closed`: field lines (on)

These layers draw the field lines traced from 4,500 footpoints spread uniformly at 1.03 R⊙, and every line is traced through the whole model, out to 30.5 R⊙. A line that leaves the outer boundary is open and drawn in black, a line whose two ends return to the Sun is closed and drawn in grey, and a line that is neither is counted as undetermined and not drawn, while the numbers of each kind are printed at every rendering. Both kinds are nearly transparent, so that the surface remains visible, and both are drawn only within `line_r_max` = 3 R⊙. The open lines mark where electrons accelerated by the shock can escape into the heliosphere, whereas the closed lines mark where they remain trapped.

### 5.4 `streamer`: the closed field lines through the CME shell (on)

This layer draws in purple the closed field lines that pass through the CME shell. The lines are traced from every vertex of the shell and only the closed ones are kept, so that they are the loops of the streamer that the CME meets, selected without any free parameter. For the default event, 478 of the 943 lines through the shell are closed. When a coronagraph layer is on and the shell is drawn at the time of the images (Section 5.5), these lines still belong to the shell at the time of the reconstruction (15:00 UT).

### 5.5 `gcs`: the CME shell (on)

This layer draws the graduated cylindrical shell of the CME as rings around its axis in magenta. Without coronagraph images, the shell is drawn at the time of the reconstruction (key `gcs`, 15:00 UT). With any coronagraph layer on, it is drawn at the time of the images (key `gcs_wl`, 15:24 UT), so that the shell and the observed CME can be compared directly. The default parameters were fitted to the LASCO C2 and COR1 images with `corona3d gcs` (METHODS.md, Section 7), and the parameters of another fit are entered under `gcs` and `gcs_wl` in `config.json` (Section 7).

### 5.6 `c2`, `cor1`, and `cor2`: coronagraph images (off)

These layers draw the coronagraph images, compared with the first frame of the download window, with LASCO C2 in a view from Earth and COR1 and COR2 of STEREO-B in a view from STEREO-B. A view from elsewhere shows no image, since an image covers the model at its true place only when seen from its own spacecraft. Each image is the frame nearest to `wl_time` divided by the first frame of the window minus one (`wl_mode` = `"ratio"`), or the difference of the two (`"difference"`), after both were normalised by their exposure times, and the ratio removes vignetting and other multiplicative instrumental effects. The result is rotated to solar north up and drawn on the plane through Sun centre perpendicular to the line of sight of the spacecraft. Bright regions mark where the brightness increased since the base frame, in particular the CME and the material that it displaces, and dark regions mark where it decreased. The occulter and the region outside the field of view are transparent: C2 covers 2.2 to 6 R⊙, COR1 1.5 to 4 R⊙, and COR2 4 to 15 R⊙. When one of these layers is on, the panel of its spacecraft is widened (Section 4), and the lines and the surface keep their size in R⊙, so that they appear smaller. Since the default panel reaches 6 R⊙, only the inner part of the COR2 field appears; `--half 15` shows all of it.

The terminal shows which frames were used, for example

```
c2: 2013-10-25T15:24:05 vs base 2013-10-25T14:36:07 (ratio)
cor1: 2013-10-25T15:25:51 vs base 2013-10-25T14:35:51 (ratio)
cor2: 2013-10-25T15:24:51 vs base 2013-10-25T14:39:51 (ratio)
```

A `wl_time` outside `wl_window`, or one whose nearest frame is the base frame itself, stops the rendering with a message, since either would draw a wrong or a blank image. A time zone may be given (`2013-10-25T15:24Z`).

### 5.7 `radio`: radio sources (off)

This layer draws the sources of a radio burst as coloured blobs, in the manner of the blobs in Fig. 7 of Morosan et al. (2022). Its input is radio images of two instruments, which may be used together. The first is the 2-D images of the Nançay Radioheliograph (NRH) at 150 to 450 MHz: the files `nrh2_FFFF_..._i.fts`, one per frequency and observing period, each holding a Stokes I image of 64 × 64 pixels about every 10 s. Since these files are supplied only on request by the Radio Solar Database at Nançay (Section 3), `corona3d fetch` does not download them; they are copied into `data/nrh/` of the working directory, and `corona3d fetch --nrh` checks which frequencies they cover. The second is the spectral images of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA) at 32 to 87 MHz, available from 2024 April onward: one file `ovro-lwa-352.lev1_mfs_10s.<time>Z.image_I.hdf` per 10 s integration, with 12 band images centred at 34.1 to 84.6 MHz. These files are public, in the OVRO-LWA solar data bucket on Amazon Web Services, and `corona3d fetch --ovro-lwa` downloads the band-averaged files of the interval of every OVRO-LWA entry.

The image sets are listed under `radio_images` in `config.json`, one entry per instrument, for example

```json
{"radio_images": [
  {"instrument": "nrh", "folder": "data/nrh", "freqs": [150.9, 228.0], "start": "2024-05-14T14:08:00",
   "end": "2024-05-14T14:14:00", "step_s": 60, "level": 0.5, "min_tb": 1e8, "contours": false},
  {"instrument": "ovro_lwa", "folder": "data/ovro_lwa", "freqs": [52.4, 70.8], "start": "2024-05-14T14:08:00",
   "end": "2024-05-14T14:14:00", "step_s": 60, "level": 0.5, "min_tb": 1e6, "contours": false}]}
```

where the numbers only illustrate the format. For OVRO-LWA, each frequency selects the nearest band within 3 MHz, and each time the nearest file within 2 min.

For each frequency of an entry (`freqs`), the layer takes one image every `step_s` seconds between `start` and `end`, keeps the images whose peak brightness temperature reaches `min_tb`, and measures in each the centroid of the pixels brighter than `level` times the peak. Each centroid becomes a blob at the height that the four-fold Newkirk model gives for that frequency at the second harmonic, on the side of the plane of the sky that faces Earth, and the blobs of one frequency are coloured along one colour map from the first image to the last (`radio_colormaps`: spring, winter, autumn, and summer for the first to the fourth frequency, counted over all entries). With `"contours": true`, the contour of the last image at `level` times its peak is drawn as well, at the same height. The terminal lists every image used, with its time, peak brightness temperature, centroid, and height:

```
radio: NRH 150.9 MHz 14:58:55 UT, peak 1.15e+10 K, centroid (-1.20, -0.01) R_sun, at 1.647 R_sun
```

For the event of 2013 October 25, the default settings take one image per minute between 14:59 and 15:03 UT at 150.9 and 432 MHz, which gives five blobs at 1.647 R⊙ and three at 1.221 R⊙. A frequency or a time that the files do not cover stops the rendering with a message that names the times the files cover. The NRH files of the afternoon of that day end at 15:04:52 UT, and the afternoon files at 228 and 327 MHz were not available, so that the blobs follow the figure of Morosan et al. (2022), drawn at 173 and 228 MHz, in kind but not in frequency. When the layer is on, the purple streamer lines are drawn more faintly, so that the blobs remain visible; a blob that lies inside the density surface, as the 432 MHz blobs do from Earth, remains hidden behind it, since the reconstruction keeps true depths.

Positions measured with other instruments can be added under `radio_sources` in `config.json`, as a list in which every source has its frequency in MHz (`freq`) and its sky-plane position seen from Earth, `x` towards the west limb and `y` towards solar north, in R⊙ or, with `"unit": "arcsec"`, in arcsec. A source may be one position or a series of positions (lists of equal length for `x` and `y`):

```json
{"radio_sources": [{"freq": 173, "x": [-1.19, -1.14, -1.09], "y": [-0.15, -0.12, -0.08]},
                   {"freq": 228, "x": -1100, "y": -190, "unit": "arcsec"}]}
```

The numbers above only illustrate the format and are not measured positions. In macOS and Linux, the same list can also be given on the command line enclosed in single quotes (`--radio-sources '[{"freq": 173, ...}]'`); the shells of Windows treat quotes inside an argument differently, so that on Windows the configuration file is the reliable route. These sources are placed in the same way as the NRH blobs and are coloured by frequency, from light blue at 150 MHz to magenta at 450 MHz and above. Arcseconds are converted with the angular radius of the Sun seen from Earth at `time` (964.8 arcsec on 2013 October 25), and a position farther from Sun centre than its height is kept in the plane of the sky, with a warning (METHODS.md, Section 8).

## 6. The viewpoints

Each panel looks at the reconstruction from a direction given by one entry of `views` in `config.json`, and the default entries are Earth, the solar north pole, and the STEREO spacecraft of the event (STEREO-B for the default event):

```json
{"views": [{"from": "earth",    "frame": [1.513, 0.087, 0.158]},
           {"from": "north",    "frame": [1.554, 0.076, 0.104]},
           {"from": "stereo", "frame": [1.455, -0.181, 0.093]}]}
```

An entry names an observer with `"from"` (`"earth"`, `"stereo"`, which is STEREO-A or STEREO-B according to `stereo_sc`, or `"north"`), or gives the direction from which the Sun is seen as a Carrington longitude and latitude, `{"lon": 300, "lat": 20}`, in degrees. A view from a longitude and latitude has solar north up, and `"roll"` turns the panel counter-clockwise by that angle. The north-pole view looks down the rotation axis with the direction of the Earth at the bottom. Moreover, `"frame"` gives the half-height of the panel and the offset of Sun centre to the right and upwards, in R⊙ (default `[1.6, 0, 0]`), and `"title"` replaces the panel title. To look at the CME from behind, for example, panel b can be set to look from Carrington longitude 124°, opposite to the direction of the CME (Carrington longitude 304°, Stonyhurst longitude −53°).

On the command line, `--view` sets one panel and may be repeated:

```
corona3d render --view b=north          # panel b from above the north pole
corona3d render --view c=300,20         # panel c from Carrington longitude 300, latitude 20
corona3d render --view c=300,20,45      # the same, turned by 45 degrees
corona3d render --view a=stereo --view c=earth
```

The legacy names `stereo_a` and `stereo_b` are accepted by `--view` as `stereo`, provided that they name the spacecraft of `stereo_sc`; a name for the other spacecraft is refused with a message, since it would silently change the observer. The panels are labelled a, b, and c from left to right. A panel set with `--view` to a longitude and latitude uses the framing `[1.6, 0, 0]`, whereas a panel set to the observer it already had keeps its framing.

## 7. The configuration file

`config.json` holds every parameter of the event and of the figure. Only the keys that are present replace the defaults, and nested keys are merged one by one, so that a file containing `{"gcs": {"height": 2.0}}` alone changes only the apex height of the shell. Lists such as `views`, however, are replaced as a whole. The keys are listed below.

| Key | Default | Meaning |
|---|---|---|
| `event` | text | description of the event, for the user only |
| `time` | `2013-10-25T15:00` | time of the CME shell `gcs` and of the angular solar radius used for the radio sources (the observers are set by `earth` and `stereo`) |
| `run` | `cr2143-medium/hmi_mast_mas_std_0201` | PSI MAS/MAST run, as its path in the PSI archive |
| `density_scale` | 2.0 | factor applied to the model density (Morosan et al. 2022) |
| `va_density` | `"scaled"` | density for the Alfvén speed: `"scaled"` (as the surface) or `"model"` |
| `model_dlon` | 0.0 | rotation of the model in Carrington longitude, degrees (METHODS.md, Section 3) |
| `ne_level` | 9 × 10⁷ | electron density of the surface, cm⁻³ |
| `va_range` | [100, 700] | range of the Alfvén-speed colours, km s⁻¹ |
| `sphere_r` | 1.2 | radius of the reference sphere, R⊙ |
| `earth`, `stereo` | [357.372, 4.993], [215.329, −0.984] | Carrington longitude and latitude of Earth and of the STEREO spacecraft, degrees |
| `stereo_sc` | `"B"` | the STEREO spacecraft at `stereo`, `"A"` or `"B"`; it selects the COR1 and COR2 images and the title of the third panel |
| `gcs` | height 1.931, kappa 0.396, alpha_deg 23.46, lon_stonyhurst −53.262, lat_deg 3.146, tilt_deg −146.654 | GCS shell at `time`: apex height (R⊙), aspect ratio, half-angle (deg), direction (Stonyhurst, deg), and tilt (deg; PyThea's tilt minus 90°) |
| `gcs_wl` | time `2013-10-25T15:24`, height 5.538, the same shape | GCS shell at the time of the coronagraph images |
| `gcs_note` | none (set in the preset `2024-05-09`) | free text of the user about the values of `gcs` and `gcs_wl`, for example that a height is assumed and why; it is not used by the figure and is copied into the metadata of `corona3d export` |
| `gcs_fit` | none | provenance of a fit written by `corona3d gcs` (Section 9.3); later fits read it, the figure does not |
| `mesh_rings` | 5, 10, width 3, `#c23ad8` | vertices of the shell along its legs and front, line width, and colour |
| `n_seeds`, `seed_r` | 4500, 1.03 | number and radius of the footpoints of the open and closed lines |
| `line_r_max` | 3.0 | field lines are drawn only within this radius, R⊙ |
| `streamer_max_lines` | 150 | at most this many streamer lines are drawn, an even subset |
| `radio_images` | one NRH entry: folder `data/nrh`, 14:59 to 15:03 UT every 60 s, 150.9 and 432 MHz, level 0.5, min_tb 10⁸ K, no contours | radio image sets for the blobs, NRH or OVRO-LWA (Section 5.7); the optional key `interval_source` is `"default"` when `corona3d new` wrote the interval, `"data"` when `corona3d radio-interval` chose it from the images, and absent when it was given by hand (Section 9.5) |
| `radio_colormaps` | spring, winter, autumn, summer | colour maps of the blob chains, one per frequency |
| `radio_sources` | [] | other radio sources (Section 5.7) |
| `radio_radius` | 0.12 | radius of the source spheres, R⊙ |
| `wl_window` | 14:35 to 15:40 UT | window of the white-light download; its first frame of each instrument is the base |
| `wl_time` | `2013-10-25T15:24` | time of the white-light frames |
| `wl_mode` | `"ratio"` | `"ratio"` (frame / base − 1) or `"difference"` (frame − base) |
| `wl_opacity` | 0.75 | opacity of the coronagraph images |
| `views` | Earth, north pole, STEREO | viewpoint, framing, and title of each panel (Section 6) |
| `projection` | `"orthographic"` | `"orthographic"`, which draws every depth at the scale of the plane of the sky, or `"perspective"`, the camera of the published figure |
| `camera_distance` | 6.65 | distance of the perspective camera at the default framing of each panel, R⊙ |
| `window` | [1200, 918] | size of each panel image, pixels (the aspect of the panels of the published figure) |

Configuration files written by version 0.1.0 and by the commits before this release remain valid: `stereo_b` is read as `stereo` with `stereo_sc` set to `"B"`, a view from `stereo_a` or `stereo_b` is read as a view from `stereo`, and `radio_nrh` is read as an entry of `radio_images` with the instrument `nrh` (its `colormaps` become `radio_colormaps`). The file itself is not rewritten.

A density surface at the plasma level of another frequency is obtained from the harmonic plasma relation, n_e = (f / 2 / 8980 Hz)² cm⁻³, so that, for example, 228 MHz at the harmonic gives 1.61 × 10⁸ cm⁻³:

```json
{"ne_level": 1.61e8}
```

The Python function `corona3d.newkirk_r(f_mhz)` gives the four-fold Newkirk height of a frequency, which indicates where such a surface is expected.

## 8. Using corona3d from Python

```python
from corona3d import LAYERS, make_figure, newkirk_r, reconstruct

# the README figures, written into corona3d_work/figures
make_figure("figures/corona3d.png", data_root="corona3d_work")
make_figure("figures/corona3d_wl.png", data_root="corona3d_work", c2=True, cor1=True, cor2=True)

# parameters changed for one figure only, without editing config.json
make_figure("figures/scale25.png", data_root="corona3d_work", params={"density_scale": 2.5},
            open=False, closed=False)
make_figure("figures/far.png", data_root="corona3d_work",
            params={"views": [{"from": "earth"}, {"lon": 124, "lat": 0}, {"from": "north"}]})

print(LAYERS)            # the names of the layers
print(newkirk_r(173))    # 1.576 R_sun
```

`make_figure` takes the path of the figure (relative to the working directory unless absolute), the working directory as `data_root`, the half-height `half`, a dictionary `params` whose keys replace those of the working directory's `config.json` (nested keys may be partial), and any layer as a keyword set to `True` or `False`, and a misspelt layer raises an error that lists the valid names. The reconstruction itself is built by `reconstruct(params, layers, data_root)`, where `params` and `layers` may hold only the keys that differ from the defaults (`reconstruct({"density_scale": 2.5}, {"open": False}, "corona3d_work")`). To start from the working directory's configuration, pass `corona3d.paths.config("corona3d_work")` as `params`. The function returns a dictionary with the density surface (a PyVista mesh with the Alfvén speed as `va`), the field lines (lists of N × 3 arrays in R⊙, Carrington-Cartesian coordinates, within `line_r_max`), the GCS shell (an array of rings), the coronagraph planes (PyVista grids with the image as `wl`), and the times used, for analyses beyond the figure. `corona3d.geometry.FieldTracer` traces lines through the whole model and returns each with its kind (`open`, `closed`, or `undetermined`).

### 8.1 The layers as arrays

For an analysis without the figure, `corona3d.build` reads the configuration of a working directory, builds the layers with the same code as `corona3d render`, and returns them as NumPy arrays in a `Reconstruction`:

```python
import corona3d

rec = corona3d.build("corona3d_work", radio=True)      # the layer switches are those of render
rec.surface["vertices"]       # (n, 3) iso-surface of ne_level
rec.surface["faces"]          # (m, 3) vertex indices of the triangles
rec.surface["va_kms"]         # (n,) Alfvén speed at each vertex
rec.lines["closed"][0]        # (k, 3) one closed field line; also "open" and "streamer"
rec.gcs["time"]["rings"]      # (n_axis, n_circle, 3) shell at the time of the reconstruction
rec.gcs["image"]["rings"]     # the shell at the time of the coronagraph images (gcs_wl)
rec.radio["lon_stonyhurst_deg"]   # one entry per radio source position
rec.observers["earth"]        # Carrington longitude, latitude, and direction of Earth
```

Every position is in solar radii in Heliographic Carrington Cartesian coordinates, centred on the Sun, with X towards Carrington longitude 0 and Z towards solar north; this frame rotates with the Sun, so that the positions refer to the time of the reconstruction. The Alfvén speed is in km s⁻¹, the electron density `ne_cm3` in cm⁻³, the brightness temperature `tb_k` in K, and angles in degrees. The field lines are cut at `line_r_max`, and the lines of the class `streamer` are the closed lines that pass through the CME shell. The table `rec.radio` has one row for each position of the sources in `radio_sources` (instrument `user`, no time and no brightness temperature) and for each image of the NRH and OVRO-LWA sets, with the columns `instrument`, `freq_mhz`, `time_utc`, the sky-plane position seen from Earth (`x_rsun`, `y_rsun`, west and north, and `x_arcsec`, `y_arcsec`), `tb_k`, the distance from Sun centre at the Newkirk height of the frequency (`r_rsun`) and the height above the photosphere (`height_rsun`, which is `r_rsun` − 1), the Carrington Cartesian position (`x`, `y`, `z`), and the Carrington longitude, latitude, and Stonyhurst longitude (`lon_carrington_deg`, `lat_deg`, `lon_stonyhurst_deg`). The Stonyhurst longitude is the Carrington longitude minus that of Earth at the time of the reconstruction. The layers `c2`, `cor1`, and `cor2` give, for each instrument, the points of the image plane, the image values `wl`, and the shape of the grid in `rec.whitelight`. `build(home, params, **layers)` takes the same `params` and layer switches as `make_figure`, and `corona3d.load_config(home)` returns the configuration alone, with the defaults filled in. The public names are `build`, `Reconstruction`, `load_config`, `export`, `reconstruct`, `make_figure`, and `newkirk_r`.

### 8.2 `corona3d export`: the layers as data files

```
corona3d export                                  # every layer that render draws by default, as one npz file
corona3d export --radio --format csv             # tables, also with the radio sources
corona3d export --format vtk --c2 --out ex_vtk   # files for ParaView, also with the C2 plane
```

The command writes the enabled layers into the folder `--out` (default `export`, relative to the working directory) in one of three formats, and every format is accompanied by `metadata.json`. With `npz`, the file `layers.npz` holds all arrays under the names `surface_vertices`, `surface_faces`, `surface_va_kms`, `surface_ne_cm3`, `lines_xyz`, `lines_id`, `lines_class`, `gcs_time_rings`, `gcs_image_rings`, `radio_<column>`, `<instrument>_points`, `<instrument>_wl`, and `observer_earth` and `observer_stereo` (longitude, latitude, and the three components of the direction). With `csv`, each layer is a table (`surface_vertices.csv`, `surface_faces.csv`, `field_lines.csv` with the columns `line_id`, `class`, `x`, `y`, `z`, `gcs_time.csv`, `gcs_image.csv`, `radio_sources.csv`, `observers.csv`, and `<instrument>_plane.csv`). With `vtk`, the files are written with PyVista: `surface.vtp` with the arrays `va_kms` and `ne_cm3`, `lines_open.vtp`, `lines_closed.vtp`, and `lines_streamer.vtp` with the point array `line_class` (0 open, 1 closed, 2 streamer), `gcs_time.vts` and `gcs_image.vts`, `radio_sources.vtp`, and `<instrument>_plane.vts`. The same function is available as `corona3d.export(home, out, fmt, **layers)`.

The sidecar `metadata.json` records the corona3d version, the time of the export, the SHA-256 of the configuration used together with the configuration itself, the event, its time, and the model run, the coordinate frame and the units, the meaning of the radio columns, the positions of the observers, the times of the shells, and the sources of every input that was used: the PSI MAS/MAST files with their names and sizes and the download address, the density scale and the rotation of the model, the GCS parameters with the note `gcs_note` and the fit record `gcs_fit` if present, the radio image sets with the names of the files read, and the coronagraph frames with their times.

## 9. Another event or another model run

Every element of an event is a key of `config.json` (Section 7), and another event is therefore prepared from a preset (Section 9.1), from its time and flare position (Section 9.2), or by editing the keys directly (Section 9.4). In each case the shell of the CME is then fitted to the coronagraph images (Section 9.3).

### 9.1 Event presets and STEREO-A

`corona3d init --event NAME` writes the configuration of a preset. The preset `2013-10-25` is the default event of Morosan et al. (2022), and `2024-05-09` is the CME and type II burst of 2024 May 9, with the PSI run `cr2284-high`, the positions of Earth and STEREO-A, and the radio images of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA). The GCS shell of the second preset is the fit to the coronagraph images of that event, and a new fit with Section 9.3 replaces it. The preset `2024-05-09` carries the GCS shell fitted to COR1-A and COR2-A at 17:55 UT, and, since the coronagraph images show no rising front at 17:32 UT, an apex height of 1.77 R⊙ at that time that is assumed and not measured; the key `gcs_note` states this (Section 7).

The key `stereo_sc` (`"A"` or `"B"`) states which STEREO spacecraft lies at the position `stereo`. STEREO-B lost contact with Earth in October 2014, and the events after that date therefore use STEREO-A. The third panel, the layers COR1 and COR2, the download of Section 3, and the entry `"from": "stereo"` of `views` all follow `stereo_sc`, so that the same commands work for both spacecraft.

### 9.2 `corona3d new`: a new event from its time and flare position

```
mkdir event_work
cd event_work
corona3d new --time 2024-05-09T17:32 --flare S14W28
```

writes `config.json` for the event. The positions of Earth and of the STEREO spacecraft are computed at `--time` with SunPy and JPL Horizons, which needs a network connection, and the spacecraft is STEREO-B before 2014 October and STEREO-A afterwards. The PSI MAS/MAST run is the one of the Carrington rotation of the event, at medium resolution when it is listed in the PSI archive and otherwise at high resolution. The windows of the coronagraph and radio images are set around the time, the NRH images are set for 150.9 and 228 MHz, and the OVRO-LWA images are added from 2024 April onward. The flare position is given as latitude and then longitude of the heliographic position, so that S14W28 is 14° south and 28° west of the disc centre, and the shell starts above that site with typical parameters (aspect ratio 0.4, half-angle 30°). An existing `config.json` is replaced only when `--force` is given. `corona3d fetch` then downloads the data of the event.

### 9.3 `corona3d gcs`: fitting the CME shell

```
corona3d fetch --wl --cor1a
corona3d gcs --t-ref 2024-05-09T18:52 --t-out 2024-05-09T17:32 --flare 28 -14 --views c2,cor2a --fix-direction
```

fits the GCS shell to the images in `data/wl/` of the working directory in the two steps of METHODS.md, Section 7, which follow the procedure of Morosan et al. (2022). The six GCS parameters are fitted at `--t-ref` to the CME envelopes in the chosen instruments, and the shape found is kept while the apex height at `--t-out`, the time of the reconstruction, is extrapolated from the height-time profile of the CME front in the limb view. The fit needs scikit-image, so that the package must be installed with `pip install ".[fit]"`. The options are the following.

| Option | Effect |
|---|---|
| `--wl DIR` | folder of the images (default `data/wl` of the working directory) |
| `--t-ref T` | time of the six-parameter fit (default: for the event of 2013 October 25, 15:24:30 UT; for any other event the earliest frame after `--t-out` at which the CME is well developed, found from the images as described in METHODS.md, Section 7.2, and recorded in `gcs_fit` as `t_ref_choice`); it also sets `wl_time` of `config.json` to this time and records it in `gcs_fit` (`wl_time`, `wl_time_from` = `"--t-ref"` or `"t_ref chosen automatically"`), so that the coronagraph figure is drawn at the time of the fit, for example at 17:55 UT on 2024 May 9, when LASCO C2 has no frame between 14:48 and 18:52 UT |
| `--t-out T` | time of the reconstruction (default: `time` of the configuration) |
| `--flare LON LAT` | Stonyhurst longitude and latitude of the flare site in degrees, which seeds the envelopes (default: the site of the previous fit, otherwise the direction of the shell in the configuration; S06E69 for the default event) |
| `--views LIST` | instruments, from `c2`, `cor1a`, `cor1b`, `cor2a`, `cor2b` (default `c2,cor1a,cor1b`); the legacy names `stereo_a` and `stereo_b` are accepted and stand for `cor1a,cor2a` and `cor1b,cor2b` |
| `--fix-direction` | keeps the axis of the shell at the flare site and fits only the height, the aspect ratio, the half-angle, and the tilt |
| `--front INSTR` | instrument of the height-time profile (default: the instrument that sees the CME most nearly at the limb) |
| `--running-diff MINUTES` | builds the envelopes from running differences, each frame minus the frame of the same observing programme nearest to MINUTES earlier, instead of the difference from the first frame; useful when earlier ejecta or a slow rise leave no signal against the first frame |
| `--height-out R` and `--height-note TEXT` | give the apex height at `--t-out` in R⊙ by hand, for a CME whose rise at that time cannot be measured; step 2 is skipped, and the reason in `--height-note` (required) is written into `gcs_fit` together with `"height_source": "assumed"` |
| `--drop-poor-views` / `--keep-all-views` | with the default `--drop-poor-views`, a view whose IoU is below 0.5 and below 0.75 times the mean of the others is left out and the fit repeated; the refit is kept unless the mean IoU of the remaining views falls by more than 0.01, and at least one view always remains (rule in METHODS.md, Section 7.3); both fits are written to `gcs_fit.json` (`view_check`) and the views left out to `views_dropped` of `gcs_fit`; `--keep-all-views` fits all views as given |
| `--direction-cone DEG` | with one view, or views whose observers are closer than 30°, the direction is not free but bounded within DEG of the flare site (default 15; 0 fixes it at the flare site). The constraint is recorded in `gcs_fit` (`direction_constraint`) and stated in `event_report.md` as "direction constrained by the source region (single viewpoint)"; METHODS.md, Section 7.4 |
| `--onset T` and `--base-minutes N` | the base frame of the envelopes is the last frame of each observing programme at or before the flare onset `T` minus `N` minutes (defaults: `--t-out` and 10 min), so that the fit does not depend on when the download starts; the choice is recorded in `gcs_fit` (`base`) and explained in METHODS.md, Section 7.2 |
| `--skip-step1` | takes step 1 from the `gcs_fit.json` of the working directory and repeats only step 2 |

The command writes two files into the working directory. `gcs_fit.json` holds the full result, namely the six parameters at both times, the height-time frames of the front, and the speed. In `config.json`, the keys `gcs` (at `--t-out`) and `gcs_wl` (at `--t-ref`) are replaced by the fitted shells, with the tilt of the configuration equal to the PyThea tilt minus 90°, and the new key `gcs_fit` records how they were obtained: the two times, the instruments, the intersection over union (IoU) of the shell with the envelope in each instrument, the instrument of the front, its speed in km s⁻¹, the apex height at `--t-out` with its uncertainty (`"height_source": "measured"`), the number of frames, the flare site, the interval of the running difference if one was used, and the parameters that were held fixed. An assumed height has no speed or uncertainty and is marked as such. The other keys of `config.json` are left unchanged. The fit of the default event reproduces the shell of Section 5.5 (apex height 1.931 ± 0.158 R⊙ at 15:00 UT, speed 1039 km s⁻¹), and its result is kept in `scripts/gcs_fit_2013-10-25.json` under the keys `t_ref` and `t_out`. The fit takes several minutes, since the shell is projected into every instrument for each trial of the optimisation.

### 9.4 The manual route

1. Choose a MAS/MAST run of the Carrington rotation of the event in the [PSI archive](https://www.predsci.com/data/runs/) and enter its path (for example `cr2143-medium/hmi_mast_mas_std_0201`) as `run`. The run must provide `rho002`, `br002`, `bt002`, `bp002`, and `va002` in its `corona/` folder.
2. Enter the time as `time`, keep `model_dlon` at 0, and set the white-light window and time as `wl_window` and `wl_time`.
3. Compute the Carrington coordinates of Earth and of the STEREO spacecraft at that time, enter them as `earth` and `stereo`, and set `stereo_sc`:
   ```python
   from sunpy.coordinates import get_earth, get_horizons_coord, HeliographicCarrington
   t = "2013-10-25T15:00"
   for c in (get_earth(t), get_horizons_coord("STEREO-B", t)):
       h = c.transform_to(HeliographicCarrington(observer="self", obstime=t))
       print(round(h.lon.deg, 3), round(h.lat.deg, 3))
   ```
4. Run `corona3d fetch --cor1a` and fit the shell with `corona3d gcs` (Section 9.3).
5. Copy the NRH image files of the event into `data/nrh/`, set the burst interval and the frequencies of each instrument under `radio_images`, and, for OVRO-LWA, run `corona3d fetch --ovro-lwa`; other positions can be entered under `radio_sources`. Then run `corona3d render --radio`.

### 9.5 `corona3d radio-interval` and `corona3d event`: from a time to the figures

```
corona3d radio-interval [--search MINUTES] [--force]
corona3d event --time 2024-05-09T17:32 --flare S14W28 [--home DIR] [--fit auto|fixed|skip]
               [--height-out R --height-note TEXT] [--iou-min 0.5] [--views LIST] [--t-ref T]
               [--drop-poor-views | --keep-all-views] [--burst N] [--radio-search MINUTES] [--force]
```

`corona3d radio-interval` chooses the interval of the radio blobs from the images themselves. For every frequency of a `radio_images` entry, it takes the brightest pixel of each image within `--search` minutes of the event time (default 10) and keeps the images whose brightest pixel reaches `min_tb`. Images separated by more than three times the cadence, and at least 60 s, belong to different runs, and the run that holds the brightest image is the burst of that frequency. The entry takes the earliest start and the latest end of its frequencies, and its `step_s` is lengthened so that at most about thirteen blobs are drawn per frequency. The command writes `start`, `end`, `step_s`, and `"interval_source": "data"` into `config.json` and prints the interval, the number of images, and the peak brightness temperature of every frequency. It changes only the entries that `corona3d new` wrote (`"interval_source": "default"`); an interval given by hand is kept unless `--force` is given. A frequency without any image at or above `min_tb` is named in the message, and an entry in which no frequency reaches `min_tb` keeps its interval.

When `--time` is a date without a time, for example `corona3d event --time 2024-05-09 --flare S14W28`, the command lists the type II and type III bursts of that day and stops. The list is the daily file of the edited events of the NOAA Space Weather Prediction Center (the observatory entries of type RSP), which the NCEI archive serves without a login at `https://www.ngdc.noaa.gov/stp/space-weather/swpc-products/daily_reports/solar_event_reports/` (the FTP server of SWPC, `ftp://ftp.swpc.noaa.gov/pub/indices/events/`, holds the recent files and is the fallback). Each line gives a number, the burst type, the begin and end time in UT, the frequency range in MHz, the observatory, and, for type II, the shock speed reported. `--burst N` then runs the event with the begin time of burst N as `--time`, for example `--time 2024-05-09 --burst 4` for the type II burst of 17:32 UT. The file is kept in `data/bursts/` of the working directory. The list is made by observers at several stations and is not complete, so that `--time` with a time remains the way to choose any other moment; when the file cannot be fetched or the day has no burst, the command says so and asks for `--time`.

`corona3d event` runs all the steps for one event in the working directory (`--home`, otherwise as in Section 2) and ends with the file `event_report.md`. The steps are the following.

1. `new` writes `config.json` from `--time` and `--flare` (Section 9.2) when it does not exist; a `config.json` of another time stops the command.
2. `fetch` downloads the model run, the white-light images, and, from 2024 April onward, the OVRO-LWA images of the search window of `--radio-search` minutes (default 10), so that the interval can be chosen from them, and it checks the NRH files (Section 3).
3. `radio-interval` chooses the burst interval as described above.
4. `gcs` fits the shell to the coronagraph images of `--views` (default C2, COR1, and COR2 of the spacecraft of `stereo_sc`, each when its images exist) at `--t-ref` (default: chosen from the images, the earliest coronagraph frame after `--time` at which the CME is well developed in most views, that is, an envelope of at least 3 R⊙² whose front lies at least 0.2 R⊙ inside the outer edge of the field, within 150 min of `--time`; when no frame qualifies, `--time` plus 60 min is used and the report says so). The time of the reconstruction, which is also the time of the shell at the radio sources and of the radio interval, is always `--time`, that is, the begin time of the burst when `--burst` is used. With `--fit auto`, the free fit is made on all views, a view whose IoU is far below the others is left out and the free fit repeated (Section 9.3), and the quality of the result is judged on the views that remain. Only when it is poor, or when the free fit failed, for instance because the CME envelope is not found, the fit with the direction held at the flare site (Section 9.3, `--fix-direction`) is made, also with the view rule, and the better of the two fits is kept (a fit that has only one view left, or views whose observers are less than 30° apart, does not free the direction: it stays within `--direction-cone` of the flare site, 15° by default, METHODS.md, Section 7.4): a good fit beats a poor one, otherwise the higher mean IoU of the kept views wins, and a difference below 0.02 goes to the free fit (METHODS.md, Section 7.5). A fit is poor when the mean IoU over the kept views is below `--iou-min` (default 0.5), when no instrument gave an envelope, or when a fitted parameter lies within 0.5 % of its range from a bound of the optimisation. `event_report.md` lists every step, both fits, and the views used in the accepted fit, and `gcs_fit.json` holds them under `decision`; with `--keep-all-views` the view rule is switched off. `--fit fixed` runs only the fixed-direction fit, and `--fit skip` keeps the shell of `config.json`. `--height-out` with `--height-note` gives the apex height by hand as in Section 9.3. `--t-ref` sets the time of the fit, and `wl_time` follows the time of the fit (given or chosen) as in Section 9.3.
5. `render` draws the default figure (`figures/corona3d.png`), the figure with the coronagraph images that exist (`figures/wl/corona3d.png`, as with `--c2`, `--cor1`, and `--cor2`), and the figure with the radio sources (`figures/radio/corona3d.png`, as with `--radio`) when radio files exist. A figure without data is skipped, and the rendering defaults are those of Section 4.
6. `report` writes `event_report.md`.

Every step can be repeated. The file `event_state.json` of the working directory records the downloads and the fit, a fit is repeated only when the flare, the times, or the height given by hand differ from the recorded fit, and a figure is drawn again only when it is older than `config.json`. `--force` repeats the downloads, the fit, and the figures. A partial result ends with exit code 0, and every failed download, fit, or figure is listed under *Warnings* in the report; the exit code is 1 only when no figure could be drawn, or, with `--strict`, when any step failed.

`event_report.md` has six sections. *Measured* lists the fitted shell with its parameters, the apex height when it was extrapolated from the front, and the radio intervals chosen from the images. *Assumed* lists, with the reason of each, what the data do not give: a height given by hand, an axis held at the flare site, the steady-state MAS/MAST corona, the Newkirk placement of the radio sources, and a radio interval that remained the default window. *Missing data* lists the absent model files, coronagraph images, and radio files (with the manual step for the NRH), the frequencies without a burst, and the instruments left out of the coronagraph figure. *Warnings* lists the steps that failed and, when no coronagraph had a frame usable at `wl_time`, the new `wl_time`. *Fit quality* gives the outcome of each fit and the IoU per instrument. *Figures* lists the figures drawn. An instrument with no frame within 30 min of `wl_time`, or whose nearest frame is the base frame itself, is left out of the coronagraph figure and named in the report; for 2024 May 9, for example, LASCO C2 has no frame between 14:48 and 18:52 UT, and its first frame is the base frame of the window. When no instrument has a usable frame, `wl_time` moves to the time of the nearest frame that can be drawn, the figure is drawn with that frame, and the report states the old and the new time, so that the panel is not shown as the time of the fit (`corona3d render` itself still refuses such a time with a message).

## 10. How the results were checked

The offline tests run with `python -m pytest` in a copy of the repository, take about 40 s, and use synthetic data with known answers, from which they verify that

* the field-line tracer reproduces the analytic field lines of a dipole with its default step;
* the MAS units and the Alfvén speed reproduce their analytic values;
* the rotation of the model turns the density surface by the given angle;
* the streamer loops are closed lines through the shell;
* the radio heights are those of Morosan et al. (2022);
* a feature in synthetic coronagraph frames at five roll angles appears at its true position in three dimensions, also when the base frame has another size;
* inputs that would draw a wrong or blank image are refused;
* the frames of each instrument are selected correctly;
* the configuration merges partial keys;
* halving `--half` doubles the size of the Sun in the panel;
* the viewpoints face the expected directions;
* damaged and partial downloads are detected and removed.

The online test (`CORONA3D_ONLINE=1 python -m pytest`) downloads the public data of 2013 October 25 into a temporary working directory and draws both figures of the README. The list of checks, together with the tolerances, is given in METHODS.md, Section 11, and in [MODEL_AUDIT.md](MODEL_AUDIT.md). Moreover, both figures of the README were drawn from a new working directory after `corona3d init` and `corona3d fetch`, with the C2, COR1, and COR2 frames printed in Section 5.6.

## 11. Differences from Morosan et al. (2022) and limits of the results

The MAST run used by Morosan et al. (2022) is no longer served by PSI. We therefore chose the public run of the same Carrington rotation that reproduces the outline of their density surface best (intersection over union 0.87 to 0.93), although its Alfvén speeds and its streamer differ from theirs. In particular, their colours are matched more closely with a density scale of 2.2 to 2.5 and a rotation of 7.5° east, which corona3d does not apply, and the public run does not form their regular fan of streamer loops.

Morosan et al. (2022) do not list their GCS parameters, and the default parameters were therefore fitted again from public white-light images, with the height at 15:00 UT extrapolated from three COR1-B frames at a constant speed of 1039 km s⁻¹ (METHODS.md, Section 7). Moreover, the herringbone positions of Morosan et al. (2022) are not part of corona3d, because they were measured from radio images that are not openly downloadable, and the radio layer instead draws the positions that the user measures, at Newkirk heights and without a correction for scattering. The default figure shows which closed field lines of the public model the shell intersects, but it does not test independently their conclusion that the herringbones originate on such loops, since that test requires the radio positions.

At the same time, the coronagraph images are level-0.5 ratios or differences without photometric calibration, which are suited to the comparison of positions, while brightness or mass estimates would need calibrated images. Each image is placed in the plane of the sky of its own spacecraft, which is exact only from that spacecraft. Finally, the field lines are traced in a static model of the corona before the CME, so that they show the corona that the CME meets, without its deformation by the CME.

## 12. Data sources, terms of use, and citation

The MAS/MAST runs are provided by Predictive Science Inc. at [predsci.com](https://www.predsci.com/data/runs/), and work that uses them should cite Lionello et al. (2009) and Mikić et al. (2018). LASCO is a project of the ESA/NASA Solar and Heliospheric Observatory (SOHO) mission, and its data should be acknowledged as such together with Brueckner et al. (1995). Similarly, the data of the Sun Earth Connection Coronal and Heliospheric Investigation (SECCHI) are produced by the SECCHI consortium of the NASA STEREO mission and should be cited with Howard et al. (2008). The GCS model should be cited with Thernisien (2011), PyThea with Kouloumvakos et al. (2022), and the reconstruction of the event with Morosan et al. (2022), and all references are collected in `docs/refs.bib`.

To cite corona3d itself, use [CITATION.cff](../CITATION.cff), to which a DOI will be added when the repository is made public.

## 13. Troubleshooting

| Problem | Cause and remedy |
|---|---|
| `... is not a corona3d working directory (no config.json)` | The command was given outside a working directory; run `corona3d init` there first, or give `--home` |
| `FileNotFoundError: rho002, ... .hdf missing in ...; run corona3d fetch --model` | The model run has not been downloaded, or not completely; run `corona3d fetch --model` |
| `FileNotFoundError: no c2 frames between ... ; run corona3d fetch --wl` | The white-light images have not been downloaded, or `wl_window` contains no frame of the instrument; run `corona3d fetch --wl` |
| `ValueError: no ... frame within 30 min of wl_time ...` | The images have a gap at `wl_time`; the message names the nearest frame times, and `wl_time` (or `--t-ref`) must be set to one of them |
| `ValueError: wl_time ... lies outside wl_window` or `... is the base frame itself` | Choose a `wl_time` inside the window and after its first frame |
| `corona3d event --time DATE` stops with `the SWPC event list of ... could not be fetched` | Neither the NCEI archive nor the SWPC FTP server answered, or the day is older than 1996 or newer than the last edited file; give `--time` with a time instead |
| `corona3d fetch` reports failed or damaged files and ends with exit code 1 | Give the command again; complete files are skipped and the missing ones downloaded |
| `ImportError` for `pyhdf`, or pip fails to build pyhdf | No wheel for the system; install HDF4 (`brew install hdf4` on macOS) or use `conda install -c conda-forge pyhdf` |
| pip reports that no version of PyThea matches | Python is older than 3.12; install Python 3.12 or newer |
| The rendering fails with an OpenGL or display error on a Linux server | VTK finds no OpenGL implementation; install `libosmesa6` or run under a virtual display (`xvfb-run corona3d render`) |
| The panels are blank or show only the white background | The camera sees no layer, for example after all layers are switched off; check the layers and the half-height |
| `corona3d render: error: unrecognized arguments: --...` | A misspelt option or layer; `corona3d render --help` lists the options and, at its end, the layers |
| The coronagraph image does not appear | C2 appears only in a view from Earth and COR1 and COR2 only in a view from STEREO-B; check `views` |
| `UserWarning: ... radio position(s) lie farther than ...` | A source lies farther from Sun centre in the plane of the sky than its Newkirk height; it is drawn in the plane of the sky; check its position and unit |
| Warnings from Astropy or SunPy about missing metadata of LASCO | The level-0.5 LASCO headers lack the observer position; SunPy then places the observer at Earth, which is correct for SOHO to within the 1 % difference in distance |

Problems, wrong results, and suggestions are welcome on the [issue tracker](https://github.com/jovian-explorer/corona3d/issues/new/choose).
