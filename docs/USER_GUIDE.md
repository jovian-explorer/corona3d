# corona3d user guide

This guide describes corona3d as it is in this repository: how to install it and prepare a working directory, how the public data are downloaded, what each layer of the scene shows and how it is switched on or off, how every parameter of the configuration acts on the figure, how corona3d is used from Python, how far its results can be trusted, and what to do when something does not work. The computation behind each layer, with its sources, is described in [METHODS.md](../METHODS.md), and the changes between versions are recorded in [CHANGELOG.md](../CHANGELOG.md).

Contents

1. [Install and start](#1-install-and-start)
2. [The working directory](#2-the-working-directory)
3. [Downloading the data](#3-downloading-the-data)
4. [Drawing the scene](#4-drawing-the-scene)
5. [The layers](#5-the-layers)
6. [The configuration file](#6-the-configuration-file)
7. [Using corona3d from Python](#7-using-corona3d-from-python)
8. [Another event or another model run](#8-another-event-or-another-model-run)
9. [How the results were checked](#9-how-the-results-were-checked)
10. [Differences from Morosan et al. (2022) and limits of the results](#10-differences-from-morosan-et-al-2022-and-limits-of-the-results)
11. [Data sources, terms of use, and citation](#11-data-sources-terms-of-use-and-citation)
12. [Troubleshooting](#12-troubleshooting)

---

## 1. Install and start

corona3d needs Python 3.10 or newer and runs on Windows, macOS, and Linux. It is not on PyPI, and pip installs it from the GitHub repository, which is private at present, so that the account used by git must have access to it:

```
pip install "git+https://github.com/jovian-explorer/corona3d.git"
```

To work on the code or to run the tests, install it from a copy of the repository in editable mode, together with pytest:

```
git clone https://github.com/jovian-explorer/corona3d.git
cd corona3d
pip install -e ".[test]"
```

pip installs NumPy, SciPy, and Matplotlib for the computation and the figure, PyVista and the Visualization Toolkit (VTK) for the three-dimensional rendering, pyhdf for the HDF4 files of the MHD model, Astropy and SunPy for the coronagraph images and their coordinates, and PyThea for the GCS shell of the CME. pyhdf provides compiled wheels for Windows, macOS, and Linux; on a system without a wheel, it needs the HDF4 library (`conda install -c conda-forge pyhdf` is the simplest route). The renderings are drawn off screen, so that no window opens and corona3d also runs on a server, provided that VTK finds an OpenGL implementation (Section 12).

After the installation, the command `corona3d` is available in the terminal (`python -m corona3d` does the same), and `corona3d --version` prints the version.

## 2. The working directory

Every command of corona3d works inside one working directory, which holds the configuration of the event (`config.json`), the downloaded data (`data/`), and the figures (`figures/`). The directory is created with

```
mkdir corona3d_work
cd corona3d_work
corona3d init
```

`corona3d init` writes `config.json` with the event of 2013 October 25 studied by Morosan et al. (2022). An existing `config.json` is never overwritten, so that a configuration that has been edited is kept when `init` is given again.

The working directory is found in the following order: the option `--home` given before the command (`corona3d --home D:\corona3d_work render`), the environment variable `CORONA3D_HOME`, and otherwise the current folder. Since a command given in a folder without `config.json` would fill that folder with data, `fetch` and `render` stop there with the advice to run `corona3d init` first.

After a download and two renderings, the directory contains

```
corona3d_work/
  config.json
  data/
    cr2143-medium_hmi_mast_mas_std_0201/   rho002.hdf, br002.hdf, bt002.hdf, bp002.hdf, va002.hdf (53 MB)
    wl/                                     LASCO C2 and STEREO-B COR1 and COR2 images (about 75 MB)
  figures/
    corona3d.png                            the three-panel figure
    view_earth.png, view_north.png, view_stb.png   the three panels at full resolution
```

The panel images are written next to the figure and are replaced at every rendering into the same folder, so that renderings whose panels should be kept are best given separate folders with `--out`.

## 3. Downloading the data

```
corona3d fetch            # the model run and the white-light images
corona3d fetch --model    # only the model run
corona3d fetch --wl       # only the white-light images
```

The model run is downloaded from the archive of Predictive Science Inc. (PSI), file by file, from `https://www.predsci.com/data/runs/<run>/corona/`. A file that is already present is not downloaded again, and a file is written under its final name only when it is complete, so that an interrupted download can be resumed by giving the command again.

The white-light images are searched for and downloaded through the Virtual Solar Observatory with SunPy, for the time window `wl_window` of the configuration: LASCO C2 images, and the COR1 and COR2 images of STEREO-B. The search also returns products that corona3d does not use (polarised C2 images, the COR1 sequence `s5c1B`, and the polarised COR2 sequence `n4c2B`); they are kept in `data/wl/` but are not read. When some files fail, `corona3d fetch --wl` given again downloads the missing ones, and the files already present are skipped by SunPy. On an ordinary connection, the whole download takes a few minutes.

## 4. Drawing the scene

```
corona3d render
```

draws the three-panel figure into `figures/corona3d.png`: panel a shows the scene from Earth, panel b from above the solar north pole with the direction of the Earth at the bottom, and panel c from STEREO-B. With the default layers, the rendering takes 30 to 45 s, most of it the tracing of about 4,500 field lines. The options of `render` are:

| Option | Effect |
|---|---|
| `--<layer>` / `--no-<layer>` | switch a layer on or off (Section 5) |
| `--out PATH` | file of the figure, relative to the working directory (default `figures/corona3d.png`) |
| `--half H` | half-height of all three panels in R⊙; by default the panels have the framing of the published figure (about 1.5 R⊙), widened to 6 R⊙ for C2 or COR2 and to 4 R⊙ for COR1 |
| `--wl-time T` | time of the white-light frames, e.g. `2013-10-25T15:12`; the frame nearest to it is taken for each instrument |
| `--radio-sources JSON` | radio sources for the radio layer (Section 5.7) |

Some combinations that are useful in practice:

```
corona3d render --no-open --no-closed --out figures/surface.png          # the density surface, the streamer, and the CME only
corona3d render --no-surface --out figures/field.png                       # the field lines without the surface
corona3d render --c2 --cor1 --cor2 --out figures/scene_wl.png              # with the coronagraph images
corona3d render --cor1 --no-open --no-closed --out figures/cor1.png        # COR1 alone, 4 R_sun panels
corona3d render --c2 --wl-time 2013-10-25T15:36 --out figures/c2_1536.png  # a later C2 frame
```

## 5. The layers

Each layer is switched on with `--<name>` and off with `--no-<name>`, and only the layers that are on are computed, so that a figure without field lines is drawn in a few seconds.

### 5.1 `surface`: the density surface coloured by the Alfvén speed (on)

The surface on which the electron density of the MHD model equals `ne_level` (9 × 10⁷ cm⁻³ by default). This is the plasma level of 173 MHz emitted at the second harmonic, the height at which the 173 MHz herringbones of the event are emitted. The colour gives the Alfvén speed on the surface, from 100 km s⁻¹ (blue) to 700 km s⁻¹ (red), on the colour bar below the panels. Blue regions are where a shock of a few hundred km s⁻¹ becomes strongly super-Alfvénic; Morosan et al. (2022) found the herringbone sources above such a region, where the CME flank meets the streamer. The colour bar is drawn only when this layer is on.

### 5.2 `sphere`: the reference sphere (on)

A grey sphere of radius `sphere_r` = 1.2 R⊙. Where the sphere is visible through the surface, the density level lies below 1.2 R⊙, so that the sphere is a reference of height.

### 5.3 `open` and `closed`: field lines (on)

Field lines traced from 4,500 footpoints spread uniformly at 1.03 R⊙: lines that reach the outer boundary of the model (3 R⊙) are drawn in black, and lines that return to the Sun in grey, both nearly transparent so that the surface remains visible. Open lines mark where electrons accelerated by the shock can escape into the heliosphere; closed lines mark where they remain trapped.

### 5.4 `streamer`: the streamer loops above the CME flank (on)

The closed loops of the helmet streamer near the CME, in purple. They are traced from a line of 40 seeds at 1.8 R⊙ along the meridian at Stonyhurst longitude −40°, over ±60° in latitude about the point of the CME shell that faces Earth, and only loops whose apex lies above 1.9 R⊙ are kept. These are the loops that the CME flank meets as it expands (METHODS.md, Section 5).

### 5.5 `gcs`: the CME shell (on)

The graduated cylindrical shell of the CME at 15:00 UT, drawn as rings around its axis in magenta. The default parameters are those fitted from public images (METHODS.md, Section 7); the parameters of another fit are entered under `gcs` in `config.json` (Section 6).

### 5.6 `c2`, `cor1`, and `cor2`: coronagraph images (off)

Base-difference images of the coronagraphs: LASCO C2 in the Earth view, and COR1 and COR2 of STEREO-B in the STEREO-B view (the north-pole view shows no image). Each image is the frame nearest to `wl_time` minus the first frame of the downloaded window, normalised by the exposure times, rotated to solar north up, and drawn on the plane through Sun centre perpendicular to the line of sight of its spacecraft. The image therefore covers the model at its true scale and orientation in the view of the spacecraft that recorded it. Bright regions are where the brightness increased since the base frame, in particular the CME and the material that it displaces, and dark regions are where it decreased. The occulter and the region outside the field of view are transparent: C2 covers 2.2 to 6 R⊙, COR1 1.5 to 4 R⊙, and COR2 4 to 15 R⊙. When one of these layers is on, the panel in which the image appears is widened (Section 4), and the lines and the surface keep their size in R⊙, so that they appear smaller. For the default event, the CME front lies on the east limb of the C2 image and on the west limb of the COR2 image, on the same side as the GCS shell in each view.

The terminal shows which frames were used, for example

```
c2: 2013-10-25T15:24:05 - base 2013-10-25T14:36:07
cor1: 2013-10-25T15:25:51 - base 2013-10-25T14:35:51
cor2: 2013-10-25T15:24:51 - base 2013-10-25T14:39:51
```

### 5.7 `radio`: radio sources (off)

Spheres at the positions of radio sources, coloured by frequency (light blue at 150 MHz to magenta at 450 MHz and above). corona3d does not derive the sources from radio images; their positions are given by the user as a list in which every source has its frequency in MHz (`freq`) and its sky-plane position seen from Earth, `x` towards the west limb and `y` towards solar north, in R⊙ or, with `"unit": "arcsec"`, in arcsec. A source may be one position or a series of positions (lists of equal length for `x` and `y`), for example the positions of a centroid at successive times:

```
corona3d render --radio --radio-sources '[{"freq": 173, "x": [-1.19, -1.14, -1.09], "y": [-0.15, -0.12, -0.08]}, {"freq": 228, "x": -1100, "y": -190, "unit": "arcsec"}]'
```

(in macOS and Linux; the Windows shells treat quotes inside an argument differently, so that on Windows the list is best written in `config.json`). Since such lists are long, they are more conveniently written once under `radio_sources` in `config.json` (Section 6), after which `corona3d render --radio` draws them. Each source is placed at the height where the four-fold Newkirk density model gives the plasma frequency of its frequency at the second harmonic (1.576 R⊙ at 173 MHz and 1.449 R⊙ at 228 MHz), on the side of the plane of the sky that faces Earth (METHODS.md, Section 8). The numbers above only illustrate the format; they are not measured positions. For the event of 2013 October 25, the positions are measured from the images of the Nançay Radioheliograph, which the BASS2000 archive of the Paris Observatory delivers on request as visibilities, to be imaged with the NRH software of SolarSoft.

## 6. The configuration file

`config.json` holds every parameter of the event and of the figure. Only the keys that are present replace the defaults, so that a file containing `{"density_scale": 2.5}` alone is a valid configuration. The keys are:

| Key | Default | Meaning |
|---|---|---|
| `event` | text | description of the event, for the user only |
| `time` | `2013-10-25T15:00` | time of the scene (observers and CME shell) |
| `run` | `cr2143-medium/hmi_mast_mas_std_0201` | PSI MAS/MAST run, as its path in the PSI archive |
| `density_scale` | 2.0 | factor applied to the model density (Morosan et al. 2022) |
| `model_dlon` | −8.243 | rotation of the model in Carrington longitude, degrees (METHODS.md, Section 3) |
| `ne_level` | 9 × 10⁷ | electron density of the surface, cm⁻³ |
| `va_range` | [100, 700] | range of the Alfvén-speed colours, km s⁻¹ |
| `sphere_r` | 1.2 | radius of the reference sphere, R⊙ |
| `earth`, `stereo_b` | [357.372, 4.993], [215.329, −0.984] | Carrington longitude and latitude of the observers at `time`, degrees |
| `gcs` | height 1.443, kappa 0.334, alpha_deg 1.071, lon_stonyhurst −53.262, lat_deg 3.146, tilt_deg −146.654 | GCS parameters: apex height (R⊙), aspect ratio, half-angle (deg), direction (Stonyhurst, deg), and tilt (deg; PyThea's tilt minus 90°) |
| `mesh_rings` | 5, 10, width 3, `#c23ad8` | vertices of the shell along its legs and front, line width, and colour |
| `n_seeds`, `seed_r` | 4500, 1.03 | number and radius of the footpoints of the open and closed lines |
| `streamer_flank` | r 1.8, half_span 60, n 40, lon −40 | seed line of the streamer loops: radius, half-span in latitude, number of seeds, and Stonyhurst longitude |
| `streamer_min_apex` | 1.9 | lowest apex of a streamer loop that is drawn, R⊙ |
| `radio_sources` | [] | radio sources (Section 5.7) |
| `radio_radius` | 0.12 | radius of the source spheres, R⊙ |
| `wl_window` | 14:35 to 15:40 UT | window of the white-light download; its first frame of each instrument is the base |
| `wl_time` | `2013-10-25T15:24` | time of the white-light frames |
| `wl_opacity` | 0.75 | opacity of the coronagraph images |
| `camera_distance` | 6.65 | distance of the camera from the focal point at the default framing, R⊙ |
| `frame` | three entries | half-height and offset of Sun centre (right, up) of each panel, R⊙ |
| `window` | [1080, 918] | size of each panel image, pixels |

A density surface at the plasma level of another frequency is obtained from the harmonic plasma relation, n_e = (f / 2 / 8980 Hz)² cm⁻³; for example, 228 MHz at the harmonic gives 1.61 × 10⁸ cm⁻³:

```json
{"ne_level": 1.61e8}
```

The Python function `corona3d.newkirk_r(f_mhz)` gives the four-fold Newkirk height of a frequency, which helps to judge where such a surface is expected.

## 7. Using corona3d from Python

```python
from corona3d import LAYERS, make_figure, newkirk_r

# the README figures
make_figure("figures/scene.png", data_root="corona3d_work")
make_figure("figures/scene_wl.png", data_root="corona3d_work", c2=True, cor1=True, cor2=True)

# parameters changed for one figure only, without editing config.json
make_figure("figures/scale25.png", data_root="corona3d_work", params={"density_scale": 2.5},
            open=False, closed=False)

print(LAYERS)            # the names of the layers
print(newkirk_r(173))    # 1.576 R_sun
```

`make_figure` takes the path of the figure, the working directory as `data_root`, the half-height `half`, a dictionary `params` whose keys replace those of the default configuration, and any layer as a keyword set to `True` or `False`. A misspelt layer raises an error that lists the valid names. The parameters of `config.json` are not read by `make_figure`; to use them from Python, read them with `corona3d.paths.config("corona3d_work")` and pass them as `params`. The scene itself is built by `corona3d.scene.build_scene(params, layers, data_root)`, which returns a dictionary with the density surface (a PyVista mesh with the Alfvén speed as `va`), the field lines (lists of N × 3 arrays in R⊙, Carrington-Cartesian coordinates), the GCS shell (an array of rings), and the coronagraph planes (PyVista grids with the difference image as `wl`), for analyses beyond the figure.

## 8. Another event or another model run

corona3d is set up for the event of 2013 October 25, but every element of the event is a key of `config.json`. For another event:

1. Choose a MAS/MAST run of the Carrington rotation of the event in the [PSI archive](https://www.predsci.com/data/runs/) and enter its path (for example `cr2143-medium/hmi_mast_mas_std_0201`) as `run`. The run must provide `rho002`, `br002`, `bt002`, `bp002`, and `va002` in its `corona/` folder.
2. Enter the time as `time`, and set `model_dlon` to 0 unless there is a reason to rotate the model.
3. Compute the Carrington coordinates of the observers at that time with SunPy, and enter them as `earth` and `stereo_b`:
   ```python
   from sunpy.coordinates import get_earth, get_horizons_coord, HeliographicCarrington
   t = "2013-10-25T15:00"
   for c in (get_earth(t), get_horizons_coord("STEREO-B", t)):
       h = c.transform_to(HeliographicCarrington(observer="self", obstime=t))
       print(round(h.lon.deg, 3), round(h.lat.deg, 3))
   ```
4. Enter the GCS parameters of the CME under `gcs`, the white-light window and time as `wl_window` and `wl_time`, and the radio sources under `radio_sources`.
5. Run `corona3d fetch` and `corona3d render`.

The third view and the coronagraph layers COR1 and COR2 are those of STEREO-B, which lost contact with Earth in October 2014. For events after that date, the third panel can still show the scene from the position given under `stereo_b` (for example the position of STEREO-A entered there), but the COR1 and COR2 layers read only STEREO-B images. The streamer seed line is placed relative to the CME shell and to the meridian under `streamer_flank.lon`, which should be chosen anew for each event from the closed loops of its model.

## 9. How the results were checked

The offline tests run with `python -m pytest` in a copy of the repository and take about half a minute. They check that the field-line tracer reproduces the analytic field lines of a dipole, that the radio heights are those of Morosan et al. (2022), that the axes of each observer are right-handed with x towards the west limb, that a feature in a synthetic coronagraph image at any roll angle is placed at its true position in three dimensions, that COR1 frames are grouped into complete triplets, that the command line behaves as described in Sections 2 and 4, and that the layers that need no data render. The online test (`CORONA3D_ONLINE=1 python -m pytest`) downloads the public data of 2013 October 25 into a temporary working directory and draws both figures of the README. The list of checks, with the tolerances, is in METHODS.md, Section 10.

Both figures of the README were drawn from a new working directory after `corona3d init` and `corona3d fetch`, and the C2, COR1, and COR2 frames used were those printed in Section 5.6.

## 10. Differences from Morosan et al. (2022) and limits of the results

* **Model run.** The MAST run used by Morosan et al. (2022) is no longer served by PSI. The public run of the same Carrington rotation reproduces the outline of their density surface (intersection over union 0.87 to 0.93), but its Alfvén speeds and its streamer differ: their colours are matched more closely with a density scale of 2.2 to 2.5 and a further rotation of 7.5° east, which corona3d does not apply, and the public run does not form their regular fan of streamer loops.
* **CME shell.** Morosan et al. (2022) do not list their GCS parameters. The default parameters were fitted again from public images and agree well with the white-light CME at 15:24 UT but less well with the EUV eruption at 15:00 UT (METHODS.md, Section 7).
* **Radio sources.** The herringbone positions of Morosan et al. (2022) are not part of corona3d, because they were measured from radio images that are not openly downloadable; the radio layer draws the positions that the user measures.
* **Coronagraph images.** The images are level-0.5 base differences without photometric calibration, suited to the comparison of positions, while brightness or mass estimates would need calibrated images. Each image is placed in the plane of the sky of its own spacecraft, which is exact only from that spacecraft.
* **Field lines.** The lines are traced in a static model of the corona before the CME; they show the corona that the CME meets, without its deformation by the CME.

## 11. Data sources, terms of use, and citation

The MAS/MAST runs are provided by Predictive Science Inc. at [predsci.com](https://www.predsci.com/data/runs/); work that uses them should cite Lionello et al. (2009) and Mikić et al. (2018). LASCO is a project of the ESA/NASA SOHO mission, and its data should be acknowledged as such together with Brueckner et al. (1995); SECCHI data are produced by the SECCHI consortium of the NASA STEREO mission and should be cited with Howard et al. (2008). The GCS model should be cited with Thernisien (2011) and PyThea with Kouloumvakos et al. (2022), and the scene of the event with Morosan et al. (2022). The references are collected in `docs/refs.bib`.

To cite corona3d itself, use [CITATION.cff](../CITATION.cff). A DOI will be added when the repository is made public.

## 12. Troubleshooting

| Problem | Cause and remedy |
|---|---|
| `... is not a corona3d working directory (no config.json)` | The command was given outside a working directory; run `corona3d init` there first, or give `--home` |
| `FileNotFoundError: ... rho002.hdf` | The model run has not been downloaded; run `corona3d fetch --model` |
| `FileNotFoundError: no c2 files in ...` | The white-light images have not been downloaded, or the window `wl_window` contains no frame of the instrument; run `corona3d fetch --wl` |
| The download stops or some files fail | Give `corona3d fetch` again; complete files are skipped and the missing ones downloaded |
| `ImportError` for `pyhdf` | No compiled pyhdf for the system; install it with `conda install -c conda-forge pyhdf` |
| The rendering fails with an OpenGL or display error on a Linux server | VTK finds no OpenGL implementation; install `libosmesa6` or run under a virtual display (`xvfb-run corona3d render`) |
| The panels are blank or show only the white background | The camera sees no layer, for example after all layers are switched off, or a very small `--half`; check the layers and the half-height |
| `unknown layers {...}` | A misspelt layer name; the message lists the valid names |
| The coronagraph image appears in the wrong panel | C2 appears only in the Earth view and COR1 and COR2 only in the STEREO-B view, by design |
| Warnings from Astropy or SunPy about missing metadata of LASCO | The level-0.5 LASCO headers lack the observer position; SunPy then places the observer at Earth, which is correct for SOHO to within the 1 % difference in distance |

Problems, wrong results, and suggestions are welcome on the [issue tracker](https://github.com/jovian-explorer/corona3d/issues/new/choose).
