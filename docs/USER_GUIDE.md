# corona3d user guide

This guide describes corona3d as it is in this repository: how to install it and prepare a working directory, how the public data are downloaded, what each layer of the reconstruction shows and how it is switched on or off, how the viewpoint of each panel is chosen, how every parameter of the configuration acts on the figure, how corona3d is used from Python, how far its results can be trusted, and what to do when something does not work. Moreover, the computation behind each layer is described with its sources in [METHODS.md](../METHODS.md), every method and number is listed with its source and its test in [MODEL_AUDIT.md](MODEL_AUDIT.md), and the changes between versions are recorded in [CHANGELOG.md](../CHANGELOG.md).

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

corona3d runs on Windows, macOS, and Linux and needs Python 3.12 or newer, because PyThea 1.3, which builds the shell of the coronal mass ejection (CME), requires it. Since corona3d is not on PyPI, pip installs it from the GitHub repository, which is private at present, and the account used by git must therefore have access to it:

```
pip install "git+https://github.com/jovian-explorer/corona3d.git"
```

To work on the code, run the tests, or run the fit of the graduated cylindrical shell (GCS) described in Section 9, install it from a copy of the repository in editable mode, together with pytest, Pillow, and scikit-image:

```
git clone https://github.com/jovian-explorer/corona3d.git
cd corona3d
pip install -e ".[test,fit]"
```

pip installs NumPy, SciPy, and Matplotlib for the computation and the figure, PyVista and the Visualization Toolkit (VTK) for the three-dimensional rendering, pyhdf for the HDF4 files of the magnetohydrodynamic (MHD) model, Astropy and SunPy for the coronagraph images and their coordinates, and PyThea for the GCS shell of the CME. pyhdf provides wheels for Windows, Linux, and recent macOS on Apple processors. On other Macs, the HDF4 library must be installed first with `brew install hdf4`, while `conda install -c conda-forge pyhdf` is an alternative on any system. Since the renderings are drawn off screen, no window opens, and corona3d therefore also runs on a server, provided that VTK finds an OpenGL implementation (Section 13).

After the installation, the command `corona3d` is available in the terminal (`python -m corona3d` does the same), and `corona3d --version` prints the version.

## 2. The working directory

Every command of corona3d works inside one working directory, which holds the configuration of the event (`config.json`), the downloaded data (`data/`), and the figures (`figures/`). The directory is created with

```
mkdir corona3d_work
cd corona3d_work
corona3d init
```

`corona3d init` writes `config.json` with the event of 2013 October 25 studied by Morosan et al. (2022). An existing `config.json` is never overwritten, so that an edited configuration is kept when `init` is given again.

corona3d finds the working directory from the option `--home` given before the command (`corona3d --home D:\corona3d_work render`), then from the environment variable `CORONA3D_HOME`, and otherwise takes the current folder. Since a command given in a folder without `config.json` would fill that folder with data, `fetch` and `render` stop there with the advice to run `corona3d init` first.

After a download and two renderings, the directory contains

```
corona3d_work/
  config.json
  data/
    cr2143-medium_hmi_mast_mas_std_0201/   rho002.hdf, br002.hdf, bt002.hdf, bp002.hdf, va002.hdf (53 MB)
    wl/                                     LASCO C2 and STEREO-B COR1 and COR2 images (about 75 MB)
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
```

The model run is downloaded file by file from the archive of Predictive Science Inc. (PSI), at `https://www.predsci.com/data/runs/<run>/corona/`. A file that is already present is not downloaded again, a new file is written under its final name only when it is complete, and the partial file of a failed download is removed, so that giving the command again downloads exactly the files that are missing.

The white-light images, namely the C2 images of the Large Angle and Spectrometric Coronagraph (LASCO) and the COR1 and COR2 images of the Solar Terrestrial Relations Observatory B (STEREO-B), are searched for and downloaded through the Virtual Solar Observatory with SunPy for the time window `wl_window` of the configuration. The search also returns products that corona3d does not draw (polarised C2 images, the COR1 sequences `s5c1B` and `n4c1B`, and the polarised COR2 sequence `n4c2B`), which are kept in `data/wl/` but are not read, except `s5c1B` by the GCS fit. After the download, every image is read once, and an image that cannot be read completely, such as one truncated by an interrupted transfer, is deleted. When files fail or are deleted, the command says so and ends with exit code 1, and giving it again downloads the missing files. On an ordinary connection, the whole download (about 130 MB) takes a few minutes.

## 4. Drawing the figure

```
corona3d render
```

draws the three-panel figure into `figures/corona3d.png`. By default, panel a shows the reconstruction from Earth, panel b from above the solar north pole with the direction of the Earth at the bottom, and panel c from STEREO-B, and Section 6 explains how to change them. The panels are drawn with an orthographic projection, so that every structure appears at its true size in R⊙ whatever its distance along the line of sight. The bottom of each panel states the model run and the time of the CME shell, and, where a coronagraph image is shown, the times of the image and of its base frame. With the default layers, the rendering takes about 40 s, most of which is spent on tracing the field lines. The options of `render` are:

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

This layer draws the graduated cylindrical shell of the CME as rings around its axis in magenta. Without coronagraph images, the shell is drawn at the time of the reconstruction (key `gcs`, 15:00 UT). With any coronagraph layer on, it is drawn at the time of the images (key `gcs_wl`, 15:24 UT), so that the shell and the observed CME can be compared directly. The default parameters were fitted to the LASCO C2 and COR1 images with `scripts/gcs_fit.py` (METHODS.md, Section 7), and the parameters of another fit are entered under `gcs` and `gcs_wl` in `config.json` (Section 7).

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

This layer draws the sources of a radio burst as coloured blobs, in the manner of the blobs in Fig. 7 of Morosan et al. (2022). Its input is radio images of two instruments, which may be used together. The first is the 2-D images of the Nançay Radioheliograph (NRH) at 150 to 450 MHz: the files `nrh2_FFFF_..._i.fts`, one per frequency and observing period, each holding a Stokes I image of 64 × 64 pixels about every 10 s. Since the BASS2000 archive of the Paris Observatory delivers these files only on request, `corona3d fetch` does not download them; they are copied into `data/nrh/` of the working directory. The second is the spectral images of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA) at 32 to 87 MHz, available from 2024 April onward: one file `ovro-lwa-352.lev1_mfs_10s.<time>Z.image_I.hdf` per 10 s integration, with 12 band images centred at 34.1 to 84.6 MHz. These files are public, in the OVRO-LWA solar data bucket on Amazon Web Services, and `corona3d fetch --ovro-lwa` downloads the band-averaged files of the interval of every OVRO-LWA entry.

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

Each panel looks at the reconstruction from a direction given by one entry of `views` in `config.json`, and the default entries are Earth, the solar north pole, and STEREO-B:

```json
{"views": [{"from": "earth",    "frame": [1.513, 0.087, 0.158]},
           {"from": "north",    "frame": [1.554, 0.076, 0.104]},
           {"from": "stereo_b", "frame": [1.455, -0.181, 0.093]}]}
```

An entry names an observer with `"from"` (`"earth"`, `"stereo_b"`, or `"north"`), or gives the direction from which the Sun is seen as a Carrington longitude and latitude, `{"lon": 300, "lat": 20}`, in degrees. A view from a longitude and latitude has solar north up, and `"roll"` turns the panel counter-clockwise by that angle. The north-pole view looks down the rotation axis with the direction of the Earth at the bottom. Moreover, `"frame"` gives the half-height of the panel and the offset of Sun centre to the right and upwards, in R⊙ (default `[1.6, 0, 0]`), and `"title"` replaces the panel title. To look at the CME from behind, for example, panel b can be set to look from Carrington longitude 124°, opposite to the direction of the CME (Carrington longitude 304°, Stonyhurst longitude −53°).

On the command line, `--view` sets one panel and may be repeated:

```
corona3d render --view b=north          # panel b from above the north pole
corona3d render --view c=300,20         # panel c from Carrington longitude 300, latitude 20
corona3d render --view c=300,20,45      # the same, turned by 45 degrees
corona3d render --view a=stereo_b --view c=earth
```

The panels are labelled a, b, and c from left to right. A panel set with `--view` to a longitude and latitude uses the framing `[1.6, 0, 0]`, whereas a panel set to the observer it already had keeps its framing.

## 7. The configuration file

`config.json` holds every parameter of the event and of the figure. Only the keys that are present replace the defaults, and nested keys are merged one by one, so that a file containing `{"gcs": {"height": 2.0}}` alone changes only the apex height of the shell. Lists such as `views`, however, are replaced as a whole. The keys are listed below.

| Key | Default | Meaning |
|---|---|---|
| `event` | text | description of the event, for the user only |
| `time` | `2013-10-25T15:00` | time of the CME shell `gcs` and of the angular solar radius used for the radio sources (the observers are set by `earth` and `stereo_b`) |
| `run` | `cr2143-medium/hmi_mast_mas_std_0201` | PSI MAS/MAST run, as its path in the PSI archive |
| `density_scale` | 2.0 | factor applied to the model density (Morosan et al. 2022) |
| `va_density` | `"scaled"` | density for the Alfvén speed: `"scaled"` (as the surface) or `"model"` |
| `model_dlon` | 0.0 | rotation of the model in Carrington longitude, degrees (METHODS.md, Section 3) |
| `ne_level` | 9 × 10⁷ | electron density of the surface, cm⁻³ |
| `va_range` | [100, 700] | range of the Alfvén-speed colours, km s⁻¹ |
| `sphere_r` | 1.2 | radius of the reference sphere, R⊙ |
| `earth`, `stereo_b` | [357.372, 4.993], [215.329, −0.984] | Carrington longitude and latitude of the observers, degrees |
| `gcs` | height 1.931, kappa 0.396, alpha_deg 23.46, lon_stonyhurst −53.262, lat_deg 3.146, tilt_deg −146.654 | GCS shell at `time`: apex height (R⊙), aspect ratio, half-angle (deg), direction (Stonyhurst, deg), and tilt (deg; PyThea's tilt minus 90°) |
| `gcs_wl` | time `2013-10-25T15:24`, height 5.538, the same shape | GCS shell at the time of the coronagraph images |
| `mesh_rings` | 5, 10, width 3, `#c23ad8` | vertices of the shell along its legs and front, line width, and colour |
| `n_seeds`, `seed_r` | 4500, 1.03 | number and radius of the footpoints of the open and closed lines |
| `line_r_max` | 3.0 | field lines are drawn only within this radius, R⊙ |
| `streamer_max_lines` | 150 | at most this many streamer lines are drawn, an even subset |
| `radio_images` | one NRH entry: folder `data/nrh`, 14:59 to 15:03 UT every 60 s, 150.9 and 432 MHz, level 0.5, min_tb 10⁸ K, no contours | radio image sets for the blobs, NRH or OVRO-LWA (Section 5.7) |
| `radio_colormaps` | spring, winter, autumn, summer | colour maps of the blob chains, one per frequency |
| `radio_sources` | [] | other radio sources (Section 5.7) |
| `radio_radius` | 0.12 | radius of the source spheres, R⊙ |
| `wl_window` | 14:35 to 15:40 UT | window of the white-light download; its first frame of each instrument is the base |
| `wl_time` | `2013-10-25T15:24` | time of the white-light frames |
| `wl_mode` | `"ratio"` | `"ratio"` (frame / base − 1) or `"difference"` (frame − base) |
| `wl_opacity` | 0.75 | opacity of the coronagraph images |
| `views` | Earth, north pole, STEREO-B | viewpoint, framing, and title of each panel (Section 6) |
| `projection` | `"orthographic"` | `"orthographic"`, which draws every depth at the scale of the plane of the sky, or `"perspective"`, the camera of the published figure |
| `camera_distance` | 6.65 | distance of the perspective camera at the default framing of each panel, R⊙ |
| `window` | [1200, 918] | size of each panel image, pixels (the aspect of the panels of the published figure) |

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

## 9. Another event or another model run

corona3d is set up for the event of 2013 October 25, but every element of the event is a key of `config.json`, and another event is therefore prepared in the following steps.

1. Choose a MAS/MAST run of the Carrington rotation of the event in the [PSI archive](https://www.predsci.com/data/runs/) and enter its path (for example `cr2143-medium/hmi_mast_mas_std_0201`) as `run`. The run must provide `rho002`, `br002`, `bt002`, `bp002`, and `va002` in its `corona/` folder.
2. Enter the time as `time`, and keep `model_dlon` at 0.
3. Compute the Carrington coordinates of the observers at that time with SunPy, and enter them as `earth` and `stereo_b`:
   ```python
   from sunpy.coordinates import get_earth, get_horizons_coord, HeliographicCarrington
   t = "2013-10-25T15:00"
   for c in (get_earth(t), get_horizons_coord("STEREO-B", t)):
       h = c.transform_to(HeliographicCarrington(observer="self", obstime=t))
       print(round(h.lon.deg, 3), round(h.lat.deg, 3))
   ```
4. Set the white-light window and time as `wl_window` and `wl_time`, run `corona3d fetch --cor1a`, and fit the shell with
   ```
   python scripts/gcs_fit.py <working directory>/data/wl --out gcs_fit.json
   ```
   from a copy of the repository. The script fits the six GCS parameters to the CME envelopes in C2, COR1-A, and COR1-B at the time `T_REF` set at its top, and it extrapolates the apex height to the time `T_OUT` from the front in the limb view (METHODS.md, Section 7). The flare site `FLARE`, which seeds the envelopes, is also set at its top, and the results `t1500` and `t1524` go under `gcs` and `gcs_wl` (tilt minus 90°).
5. Copy the NRH image files of the event into `data/nrh/`, set the burst interval and the frequencies of each instrument under `radio_images`, and, for OVRO-LWA, run `corona3d fetch --ovro-lwa`; other positions can be entered under `radio_sources`. Then run `corona3d render --radio`.

The third observer and the coronagraph layers COR1 and COR2 are those of STEREO-B, which lost contact with Earth in October 2014. For events after that date, any panel can still be set to look from the position of STEREO-A entered as a longitude and latitude under `views`, but the COR1 and COR2 layers read only STEREO-B images.

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

The online test (`CORONA3D_ONLINE=1 python -m pytest`) downloads the public data of 2013 October 25 into a temporary working directory and draws both figures of the README. The list of checks, together with the tolerances, is given in METHODS.md, Section 10, and in [MODEL_AUDIT.md](MODEL_AUDIT.md). Moreover, both figures of the README were drawn from a new working directory after `corona3d init` and `corona3d fetch`, with the C2, COR1, and COR2 frames printed in Section 5.6.

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
| `ValueError: wl_time ... lies outside wl_window` or `... is the base frame itself` | Choose a `wl_time` inside the window and after its first frame |
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
