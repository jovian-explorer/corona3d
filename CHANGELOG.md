# Changelog

## 0.1.0 (2026-10-27)

A review of version 0.0.1 found errors in the physics, the code, the documentation, and the packaging; this version corrects all of them (`docs/MODEL_AUDIT.md`, items F1 to F9).

- The model is no longer rotated by −8.243° in Carrington longitude (`model_dlon` = 0); MAS is a steady-state Carrington model, and the observers' longitudes already contain the solar rotation.
- The panels are drawn with an orthographic projection, so that the model, the shell, the radio sources, and the coronagraph images share one scale; the perspective camera had magnified structures in front of the plane of the sky by up to 18 %.
- The CME shell at 15:00 UT keeps the shape of the 15:24 UT white-light fit, and its apex height (1.931 ± 0.158 R⊙) is extrapolated from the COR1-B front (1039 km s⁻¹); the shell at 15:24 UT is drawn when coronagraph images are shown. The fit is reproducible with `scripts/gcs_fit.py`, and its result is kept in `scripts/gcs_fit_2013-10-25.json`.
- The streamer loops are the closed field lines through the CME shell, selected without free parameters, in place of a tuned seed line.
- Field lines are traced through the whole model (to 30.5 R⊙) with a step of 0.01 r and classed as open, closed, or undetermined; they are drawn within 3 R⊙ (`line_r_max`).
- The viewpoint of every panel can be chosen (`views` in `config.json`, `--view` on the command line): Earth, STEREO-B, the north pole, or any Carrington longitude and latitude with a roll angle; the panel images are now `view_a.png`, `view_b.png`, and `view_c.png`.
- `va_density` chooses the density from which the Alfvén speed is computed; radio positions in arcsec use the angular solar radius at the time of the event, and a position that cannot lie at its height produces a warning.
- Coronagraph images are compared with the first frame of `wl_window` as a ratio by default (`wl_mode`), which removes vignetting; inputs that would draw a wrong or blank image are refused, COR1 triplets are grouped by time, frames of different size are resampled, and times with a time zone are accepted. Each panel states the times of the model, the shell, and the images.
- Nested configuration keys are merged one by one; `make_figure` reads the working directory's `config.json` and writes inside it; `reconstruct` accepts partial parameters.
- Downloads: damaged and partial files are removed, and `corona3d fetch` ends with exit code 1 when files fail; `--cor1a` adds the STEREO-A COR1 images for the GCS fit.
- Python 3.12 or newer is required (PyThea 1.3); dependencies have lower bounds and PyThea is pinned to 1.3; the version has one source; a release workflow; CI installs HDF4 on macOS and runs the authorship tests.
- The command line no longer loads VTK for `--version` and `init`, and figures are drawn without a GUI backend.
- Offline tests of the MAS units and the Alfvén speed, the model rotation, the streamer selection, the undetermined lines, the radio warning, the coronagraph frames (five roll angles, ratio and difference, mixed sizes, refused inputs, C2 headers, COR1 triplets across a minute), the configuration merge, the panel scale, the viewpoints, and the downloads; 40 offline tests in all.
- README, METHODS.md, and the user guide corrected against the code, and `docs/MODEL_AUDIT.md` added.

## 0.0.1 (2026-10-27)

First pre-release, in a private repository.

- Three-dimensional reconstruction of a coronal mass ejection and the surrounding corona in the layout of Fig. 7 of Morosan et al. (2022), seen from Earth, from above the solar north pole, and from STEREO-B, built from public data only, with every layer switchable from the command line (`--<layer>` / `--no-<layer>`) and from Python (`make_figure(..., layer=True)`).
- Layers: the electron-density surface of a PSI MAS/MAST run at a chosen plasma level, coloured by the Alfvén speed computed from the scaled density; open and closed field lines; the streamer loops above the CME flank, from a seed line; the 1.2 R⊙ reference sphere; the GCS shell of the CME (PyThea); radio sources placed at their four-fold Newkirk harmonic heights from sky-plane positions given by the user; and base-difference images of LASCO C2 and of STEREO-B COR1 (total brightness of the polarisation triplets) and COR2, drawn in the plane of the sky of their spacecraft.
- A working directory with its own `config.json` (`corona3d init`), downloads of the model run from the PSI archive and of the coronagraph images through the Virtual Solar Observatory (`corona3d fetch`), and figures with their full-resolution panels (`corona3d render`).
- Field-line tracer written for corona3d (fourth-order Runge-Kutta along the unit field vector), since the tracer of psipy did not follow these files correctly.
- GCS parameters of the event fitted again from public images, as the paper does not list them (METHODS.md, Section 7).
- Offline tests: the tracer against the analytic field lines of a dipole, the radio heights of the paper, the observer axes and the deprojection, the placement of a feature from synthetic coronagraph images at three roll angles, the grouping of COR1 triplets, the command line, and a rendering without data; an online test that downloads the data of 2013 October 25 and draws both figures of the README (`CORONA3D_ONLINE=1`).
- README, user guide (`docs/USER_GUIDE.md`), and METHODS.md, with every reference checked against its registry record (`docs/refs.bib`).
