# Changelog

## 0.0.1 (2026-10-27)

First pre-release, in a private repository.

- Three-dimensional reconstruction of a coronal mass ejection and the surrounding corona in the layout of Fig. 7 of Morosan et al. (2022), seen from Earth, from above the solar north pole, and from STEREO-B, built from public data only, with every layer switchable from the command line (`--<layer>` / `--no-<layer>`) and from Python (`make_figure(..., layer=True)`).
- Layers: the electron-density surface of a PSI MAS/MAST run at a chosen plasma level, coloured by the Alfvén speed computed from the scaled density; open and closed field lines; the streamer loops above the CME flank, from a seed line; the 1.2 R⊙ reference sphere; the GCS shell of the CME (PyThea); radio sources placed at their four-fold Newkirk harmonic heights from sky-plane positions given by the user; and base-difference images of LASCO C2 and of STEREO-B COR1 (total brightness of the polarisation triplets) and COR2, drawn in the plane of the sky of their spacecraft.
- A working directory with its own `config.json` (`corona3d init`), downloads of the model run from the PSI archive and of the coronagraph images through the Virtual Solar Observatory (`corona3d fetch`), and figures with their full-resolution panels (`corona3d render`).
- Field-line tracer written for corona3d (fourth-order Runge-Kutta along the unit field vector), since the tracer of psipy did not follow these files correctly.
- GCS parameters of the event fitted again from public images, as the paper does not list them (METHODS.md, Section 7).
- Offline tests: the tracer against the analytic field lines of a dipole, the radio heights of the paper, the observer axes and the deprojection, the placement of a feature from synthetic coronagraph images at three roll angles, the grouping of COR1 triplets, the command line, and a rendering without data; an online test that downloads the data of 2013 October 25 and draws both figures of the README (`CORONA3D_ONLINE=1`).
- README, user guide (`docs/USER_GUIDE.md`), and METHODS.md, with every reference checked against its registry record (`docs/refs.bib`).
