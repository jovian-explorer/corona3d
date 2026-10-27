# Methods of corona3d: a three-dimensional reconstruction of a coronal mass ejection and the corona from public data

This document describes how corona3d computes every element of its reconstruction, which values it uses, and where each value comes from. Every number below is a key of `config.json` that the user can change. [docs/MODEL_AUDIT.md](docs/MODEL_AUDIT.md) lists each method and value with its source and its test. The references are collected in `docs/refs.bib`, and every reference with a DOI was checked against its registry record with Argus.

## 1. Overview

The interaction of a coronal mass ejection (CME) with the surrounding corona depends on quantities that vary in all three dimensions: the electron density sets the height at which a given radio frequency is emitted, the magnetic field and the density together set the Alfvén speed that the shock must exceed, and the topology of the field decides whether accelerated electrons escape along open field or remain on closed loops. Morosan et al. (2022) brought these quantities together in one three-dimensional reconstruction of the CME of 2013 October 25, with a magnetohydrodynamic (MHD) model of the corona, a graduated cylindrical shell (GCS) of the CME, and the positions of the herringbone radio sources, seen from Earth, from above the solar north pole, and from STEREO-B. corona3d rebuilds this reconstruction from a public MHD run and public coronagraph images, and adds the white-light images of LASCO C2 and STEREO-B COR1 and COR2. Every element is a layer that can be switched on or off (Sections 3 to 9).

## 2. Data

* **Coronal model.** The thermodynamic MHD model MAS/MAST of Predictive Science Inc. (PSI; Lionello et al. 2009; Mikić et al. 2018), run `cr2143-medium/hmi_mast_mas_std_0201` (Carrington rotation 2143, boundary condition from an SDO/HMI magnetogram). It is read from the PSI archive as HDF4 files of the mass density (`rho002`) and of the three components of the magnetic field (`br002`, `bt002`, `bp002`); of the Alfvén-speed file (`va002`), only the axes of its mesh are used. The run spans 1.0 to 30.5 R⊙. The run used by Morosan et al. (2022) is no longer served by PSI, and among the twelve public runs of Carrington rotations 2141 to 2144 this one reproduces the outline of their density surface best (intersection over union 0.87 to 0.93 in the three views). This choice was therefore made by comparison with the published figure.
* **White-light images.** SOHO/LASCO C2 (Brueckner et al. 1995) images of 1024 × 1024 pixels taken with the clear polariser and the orange filter, and STEREO-B/SECCHI (Howard et al. 2008) COR1 polarisation triplets (sequence `s4c1B`, at 0°, 120°, and 240°) and COR2 total-brightness images (sequence `d4c2B`), all at level 0.5. They are downloaded from the Virtual Solar Observatory with SunPy (The SunPy Community et al. 2020) for the window `wl_window` (default 14:35 to 15:40 UT); `--cor1a` adds the COR1 images of STEREO-A, which only the GCS fit uses.
* **CME geometry.** The six GCS parameters of the CME at 15:24 and 15:00 UT, fitted to the white-light images with `scripts/gcs_fit.py` (Section 7).
* **Radio sources.** Sky-plane positions supplied by the user (Section 8).

## 3. The MHD model and the Alfvén speed

The four fields are interpolated linearly onto the cell-centre mesh of the Alfvén-speed file over the whole run. The MAS code units are converted with ρ₀ = 10⁸ cm⁻³ (number density, one proton mass per particle) and B₀ = 2.2068914 G. Following Morosan et al. (2022), who multiplied the model density by 2 after Wang et al. (2017) to obtain more accurate emission heights, the electron density is n_e = 10⁸ ρ × `density_scale` cm⁻³, with `density_scale` = 2. The Alfvén speed,

v_A = |B| / (4π n m_p)^{1/2},

is computed by default from the scaled density (`va_density` = `"scaled"`), so that it is consistent with the density surface on which it is drawn, which lowers it by a factor of √2 relative to the model density. Morosan et al. (2022) do not state which density they used for the Alfvén speed, and `va_density` = `"model"` computes it from the unscaled density instead.

The model is placed in the Carrington frame of the observers without rotation (`model_dlon` = 0). MAS is a steady-state model in Carrington coordinates, and the Carrington longitudes of the observers at 15:00 UT already contain the solar rotation since the start of the day, so that no further rotation is justified physically. `model_dlon` turns the model about the rotation axis by the given angle (positive towards increasing Carrington longitude) for comparisons with other time conventions; −8.243°, the change of the Earth's Carrington longitude between 00:00 and 15:00 UT, places the model as for observers at the start of the day.

## 4. The density surface (layer `surface`)

The scaled density is placed on a structured grid in Carrington-Cartesian coordinates (x towards Carrington longitude 0°, z along the rotation axis), and the iso-surface at `ne_level` = 9 × 10⁷ cm⁻³ is extracted with the contour filter of the Visualization Toolkit (VTK), through PyVista (Sullivan & Kaszynski 2019). This is the density that Morosan et al. (2022) used for the plasma level of 173 MHz emitted at the second harmonic; the exact value for 173 MHz is 9.28 × 10⁷ cm⁻³, and 9 × 10⁷ cm⁻³ corresponds to 170.4 MHz. The Alfvén speed is interpolated onto the surface and drawn with the `Spectral_r` colour map over `va_range` (100 to 700 km s⁻¹). The surface is not smoothed; only its normals are computed, so that the shading is smooth while the shape remains as contoured.

## 5. Field lines (layers `open`, `closed`, and `streamer`)

The field lines are integrated with a fourth-order Runge-Kutta scheme along the unit vector of B, interpolated trilinearly in (φ, θ, r), with a step of 0.01 r (0.01 R⊙ at the solar surface and 0.3 R⊙ at 30 R⊙) and at most 1500 steps in each direction from the seed. A line ends when it leaves the model below its inner boundary (1.0 R⊙) or beyond its outer boundary (30.5 R⊙). It is classed as open when one of its ends leaves the outer boundary, as closed when both ends return to the inner boundary, and as undetermined otherwise (a line that runs out of steps or stalls where B vanishes). Undetermined lines are counted and not drawn; for the default event there are none. Lines are traced through the whole model for their classification, so that a loop that closes above 3 R⊙ near the cusp of a streamer is classed as closed, but they are drawn only within `line_r_max` = 3 R⊙. The tracer was written for corona3d because the tracer of psipy (streamtracer) did not follow these files correctly in our tests: its lines were not tangent to B, with a median |cos(line, B)| of 0.83, against 1.000 for the present tracer. The offline tests check the tracer against the analytic field lines of a dipole (Section 10).

* **Open and closed lines.** The lines are seeded at `n_seeds` = 4500 points distributed uniformly over a sphere of radius `seed_r` = 1.03 R⊙ (a Fibonacci lattice). Open lines are drawn in black and closed lines in grey, both at low opacity, so that the density surface remains visible through them.
* **Streamer loops.** Morosan et al. (2022) concluded that the herringbones originate where the CME flank meets the closed loops of a streamer. corona3d selects these loops without any free parameter: lines are traced through every vertex of the GCS shell above 1.02 R⊙, and the closed ones are drawn in purple. They are therefore the closed field lines that the CME shell intersects at the time of the reconstruction. For the default event, 478 of the 943 lines through the shell are closed and 465 are open, and the program prints these numbers at every rendering.

## 6. Observers, viewpoints, and projection

The observers are given by their Carrington longitude and latitude at 15:00 UT on 2013 October 25: Earth at (357.372°, 4.993°) and STEREO-B at (215.329°, −0.984°), 142° apart. These values were computed with SunPy from the JPL Horizons ephemerides (user guide, Section 8). For each observer, corona3d uses the heliocentric-Cartesian axes of the observer's sky: x towards the west limb, y towards solar north (projected onto the plane of the sky), and z towards the observer.

Each of the three panels is given by an entry of `views`: an observer (`"from": "earth"`, `"stereo_b"`, or `"north"`, the last looking down the rotation axis with the direction of the Earth at the bottom of the panel), or any viewing direction given by a Carrington longitude and latitude, with solar north up and an optional roll angle. The panels are drawn with an orthographic (parallel) projection. Every depth along the line of sight is therefore drawn at the scale of the plane of the sky, which keeps the model, the shell, the radio sources, and the coronagraph images at one common scale in R⊙. A perspective camera would magnify structures in front of the plane of the sky by up to 18 % for a camera at 6.65 R⊙. The half-height and the offset of Sun centre of each panel (`frame` in each view) follow the layout of the published figure. When coronagraph images are shown, the panel of their spacecraft is widened to 6 R⊙ (C2 or COR2) or 4 R⊙ (COR1 alone); `--half` sets one half-height for all panels.

## 7. The CME shell (layer `gcs`)

The shell is the GCS model of Thernisien (2011), built with PyThea (Kouloumvakos et al. 2022), from the apex height, the aspect ratio κ, the half-angle α between the legs, the Stonyhurst longitude and latitude of the axis, and the tilt. corona3d takes the tilt in its own convention, and PyThea receives it plus 90°. With half the default number of vertices along the axis and around it (`mesh_rings`: 5 straight and 10 circular), the shell is drawn as rings around its axis in magenta.

Morosan et al. (2022) fitted the shell to white-light images at 15:24 UT and scaled it back in time to the radio observations at 15:00 UT, but they do not list its parameters. The default values were fitted again with `scripts/gcs_fit.py` in two steps that follow this procedure.

* **Step 1 (15:24 UT).** All six parameters are fitted by differential evolution to the CME envelopes in the LASCO C2, COR1-A, and COR1-B frames nearest to 15:24:30 UT. Each envelope is the region of the base difference brighter than three robust standard deviations above its median, connected to the radial line above the flare site (S06E69) and filled to its convex hull. The fit maximises the mean intersection over union (IoU) of the projected shell with the three envelopes. It gives an apex height of 5.538 R⊙, κ = 0.396, α = 23.46°, a direction of Stonyhurst longitude −53.262° and latitude 3.146°, and a tilt of −56.654° in the PyThea convention, with IoU values of 0.76 (C2), 0.85 (COR1-A), and 0.74 (COR1-B), each frame being compared with the first frame of its own observing programme.
* **Step 2 (15:00 UT).** The shape (κ, α, direction, and tilt) is kept, and the apex height at 15:00 UT is extrapolated from the height-time profile of the CME front.
  * The front is measured in COR1-B, the instrument whose line of sight is most nearly perpendicular to the fitted direction (projection factor 1.00), as the largest distance of the envelope from Sun centre. It is measured in every frame of the two COR1-B programmes (`s4c1B` and `s5c1B`, each differenced with its own pre-event frame) in which the envelope covers at least 1 R⊙², up to the first frame in which the front reaches 0.2 R⊙ inside the 4 R⊙ edge of the field of view.
  * The front lies at 2.646, 3.192, and 3.542 R⊙ at 15:08:21, 15:13:21, and 15:18:21 UT. A straight line through these three points gives a speed of 1039 km s⁻¹ and an apex height at 15:00 UT of 1.931 ± 0.158 R⊙.

The shell at 15:00 UT (layer `gcs`, key `gcs`) is drawn when no coronagraph image is shown; when one is, the shell at 15:24 UT (key `gcs_wl`) is drawn instead, so that the shell and the images belong to the same time, and each panel states the time of its shell. The apex at 15:00 UT lies above both radio heights of Section 8. The extrapolation rests on three frames and assumes a constant speed between 15:00 and 15:18 UT. Since an accelerating CME moves more slowly early in its rise, the true height at 15:00 UT may lie below the extrapolated value. The fitted values and the measured heights are kept in `scripts/gcs_fit_2013-10-25.json`.

## 8. Radio sources (layer `radio`)

The sources are given as sky-plane positions seen from Earth (x west, y north; in R⊙, or in arcsec with `"unit": "arcsec"`, converted with the angular radius of the Sun seen from Earth at `time`, 964.8 arcsec at 15:00 UT on 2013 October 25), each with its frequency. Each source is placed at the radial distance where the four-fold Newkirk (1961) density model,

n_e(r) = 4 × 4.2 × 10⁴ × 10^{4.32/r} cm⁻³,

has the plasma frequency f_p = 8980 n_e^{1/2} Hz equal to half the observed frequency (emission at the second harmonic, as assumed by Morosan et al. 2022). This gives 1.576 R⊙ at 173 MHz and 1.449 R⊙ at 228 MHz. The point is taken on the line of sight of the observed position, on the side of the plane of the sky that faces the observer, so that x and y are kept exactly. A position that lies farther from Sun centre than its height cannot be placed at that height; it is kept in the plane of the sky, and a warning is printed. The heights follow the Newkirk model and not the density of the MHD model, so that a source and the density surface locate the same plasma frequency at different heights wherever the two density models differ. Scattering of radio waves in the corona, which displaces observed sources outward, is not corrected. The spheres are coloured by frequency.

## 9. Coronagraph images (layers `c2`, `cor1`, and `cor2`)

For each instrument, the frame nearest to `wl_time` (default 15:24 UT) within `wl_window` and the first frame of that window are read and normalised by their exposure times. For COR1, the total brightness of each triplet is formed as tB = (2/3)(I₀ + I₁₂₀ + I₂₄₀), where the three exposures of a triplet are those taken within 60 s. Each frame is binned by block averages to about 256 × 256 pixels, and a frame of another size is resampled onto the grid of the base frame. The comparison is a ratio by default (`wl_mode` = `"ratio"`), I / I_base − 1, which removes multiplicative instrumental effects such as vignetting and the stray light that scales with the coronal brightness; `"difference"` gives I − I_base instead. A `wl_time` outside the window, or a frame that would be the base frame itself, stops the rendering with a message, since either would draw a wrong or blank image. The result is rotated to solar north up with bilinear interpolation. The position of each pixel in the plane of the sky follows from the world coordinate system of the image and the angular radius of the Sun seen by the spacecraft,

(X, Y) = (pixel − Sun centre) × (plate scale / angular solar radius) R⊙.

The image is drawn on the plane through Sun centre perpendicular to the line of sight of its spacecraft, at the point X e_x + Y e_y of the observer's axes (Section 6). Seen from that spacecraft, in the orthographic projection, it therefore covers the model at its true scale and orientation. Pixels inside the occulter and outside the field of view are transparent: LASCO C2 from 2.2 to 6 R⊙, COR1 from 1.5 to 4 R⊙, and COR2 from 4 to 15 R⊙, so that COR2 begins where COR1 ends and the two coplanar images do not overlap. The grey scale spans ± the 98th percentile of the absolute values, and the opacity is `wl_opacity` = 0.75. LASCO C2 is drawn only in a view from Earth, and COR1 and COR2 only in a view from STEREO-B.

## 10. Validation

The offline tests (`python -m pytest`) use synthetic data with known answers:

* **Field-line tracer.** The tracer is checked against a dipole field, for which a line obeys r = L sin²θ.
  * With the default step and step count, the line from 1.05 R⊙ at colatitude 45° is closed. Its apex lies at 2.1 R⊙ within 0.02 R⊙, it returns to the inner boundary at latitude −46.4°, and L is constant along it to 1.5 %.
  * The line from colatitude 8° is open.
  * A line given only 20 steps is undetermined.
* **MAS units and Alfvén speed.** On a synthetic run with uniform density and a dipole field, the density, |B|, and v_A equal their analytic values in code units times 10⁸ cm⁻³ and 2.2068914 G. v_A from the model density is √2 times v_A from the doubled density.
* **Model rotation.** A density enhancement at Carrington longitude 0° appears at longitude `model_dlon` (0° and 40°) on the density surface.
* **Streamer loops.** The loops selected through the shell are closed, with both feet on the inner boundary, and pass within 0.05 R⊙ of a vertex of the shell.
* **Radio sources.** The radio heights are 1.576 R⊙ at 173 MHz and 1.449 R⊙ at 228 MHz. A position farther from Sun centre than its height produces a warning.
* **Observer axes.** The axes are orthonormal and right-handed, with x towards the west limb, and the deprojected sources lie at their radius with their sky-plane coordinates unchanged.
* **Coronagraph planes.**
  * A feature placed at (5, 2) R⊙ west and north in synthetic COR2 frames, at roll angles of 0°, 6.6°, 90°, 175.3°, and 271°, is found in three dimensions at its true position within 0.15 R⊙ for both the ratio and the difference. It is also found there when the base frame has half the size of the frame shown, and the occulter is transparent.
  * Frames smaller than 256 pixels are not binned.
  * A `wl_time` outside the window, the base frame itself, and an empty window are refused.
  * A time with a time zone is accepted.
* **Frame selection.** COR1 frames are grouped into triplets by time, including a triplet across a minute boundary, and incomplete triplets are dropped. C2 frames are selected from their DETECTOR, POLAR, and NAXIS1 keywords, and files without them are skipped.
* **Command line and configuration.**
  * The command line refuses to work outside a working directory, and `corona3d init` leaves an edited configuration unchanged.
  * Partial nested keys of `config.json` are merged with the defaults, and an unknown layer is refused from Python.
  * Halving `--half` doubles the radius of the Sun in the panel to within 3 %, and a zero half-height is refused.
  * The viewpoints (observers, Carrington directions, roll, north pole) give the expected directions and screen orientation.
  * The layers that need no data render.
* **Downloads.** The model download skips files that are present, writes a file under its final name only when it is complete, and removes the partial file of a failed download. A truncated image is detected, deleted, and counted as failed, so that the next `corona3d fetch --wl` downloads it again.

With `CORONA3D_ONLINE=1`, the online test also downloads the public data of 2013 October 25 and draws both figures of the README.

## 11. Limitations

* The public MAS/MAST run is not the run of Morosan et al. (2022), so that the density surface agrees with theirs in outline but not in every detail. Their Alfvén-speed colours are matched more closely with a density scale of 2.2 to 2.5 and the model rotated by 7.5° east, which indicates an older model version or a different magnetogram date. corona3d keeps the scale of 2 stated in their paper and applies no rotation.
* The streamer loops of the public run do not form the regular fan of loops shown by Morosan et al. (2022). corona3d draws the closed lines that the shell intersects, a selection that does not depend on their figure, so that the default figure does not test their conclusion independently of the model.
* The height of the shell at 15:00 UT rests on three COR1-B frames and on a constant speed (Section 7).
* The radio sources are not derived from the Nançay Radioheliograph images by corona3d, since these are delivered by the BASS2000 archive on request as visibilities that must be imaged with SolarSoft.
* The coronagraph images are level 0.5 and are not calibrated photometrically. The ratio removes multiplicative instrumental effects, but it is a relative change and is used here for the comparison of positions only.
* The images are drawn in the plane of the sky of their spacecraft, which is exact only when seen from that spacecraft. From any other viewpoint, the plane shows only where the image lies, without the three-dimensional structure of the CME.
* The field lines are traced in a static model of the corona before the CME, so that they show the corona that the CME meets, without its deformation by the CME.

## References

Brueckner, G. E., Howard, R. A., Koomen, M. J., et al. 1995, Solar Physics, 162, 357. https://doi.org/10.1007/BF00733434

Howard, R. A., Moses, J. D., Vourlidas, A., et al. 2008, Space Science Reviews, 136, 67. https://doi.org/10.1007/s11214-008-9341-4

Kouloumvakos, A., Rodriguez-Garcia, L., Gieseler, J., et al. 2022, Frontiers in Astronomy and Space Sciences, 9, 974137. https://doi.org/10.3389/fspas.2022.974137

Lionello, R., Linker, J. A., & Mikić, Z. 2009, The Astrophysical Journal, 690, 902. https://doi.org/10.1088/0004-637X/690/1/902

Mikić, Z., Downs, C., Linker, J. A., et al. 2018, Nature Astronomy, 2, 913. https://doi.org/10.1038/s41550-018-0562-5

Morosan, D. E., Pomoell, J., Kumari, A., Vainio, R., & Kilpua, E. K. J. 2022, Astronomy & Astrophysics, 668, A15. https://doi.org/10.1051/0004-6361/202244432

Newkirk, G., Jr. 1961, The Astrophysical Journal, 133, 983. https://doi.org/10.1086/147104

Sullivan, C. B., & Kaszynski, A. A. 2019, Journal of Open Source Software, 4, 1450. https://doi.org/10.21105/joss.01450

The SunPy Community, Barnes, W. T., Bobra, M. G., et al. 2020, The Astrophysical Journal, 890, 68. https://doi.org/10.3847/1538-4357/ab4f7a

Thernisien, A. 2011, The Astrophysical Journal Supplement Series, 194, 33. https://doi.org/10.1088/0067-0049/194/2/33

Wang, T., Reginald, N. L., Davila, J. M., St. Cyr, O. C., & Thompson, W. T. 2017, Solar Physics, 292, 97. https://doi.org/10.1007/s11207-017-1130-3
