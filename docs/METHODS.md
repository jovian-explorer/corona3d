# Methods of corona3d

This document describes how corona3d computes every element of its reconstruction, which values it uses, and where each value comes from. Every number below that is written in `code font` is a key of `config.json`, which the user can change. Each method and value is also listed with its source and its test in [MODEL_AUDIT.md](MODEL_AUDIT.md), the assumptions that limit the results are collected in [LIMITATIONS.md](LIMITATIONS.md), and the references are in [refs.bib](refs.bib); every reference with a digital object identifier (DOI) was checked against its registry record with Argus.

## 1. Overview

The interaction of a coronal mass ejection (CME) with the surrounding corona depends on quantities that vary in all three dimensions. The electron density sets the height at which a given radio frequency is emitted, the magnetic field and the density together set the Alfvén speed that a shock must exceed, and the topology of the field decides whether accelerated electrons escape along open field or remain on closed loops. Morosan et al. (2022) brought these quantities together in one three-dimensional reconstruction of the CME of 2013 October 25, which combines a magnetohydrodynamic (MHD) model of the corona, a graduated cylindrical shell (GCS) of the CME, and the positions of herringbone radio sources, seen from Earth, from above the solar north pole, and from the Solar Terrestrial Relations Observatory B (STEREO-B). corona3d follows the method described in Morosan et al. (2022), builds the reconstruction from a public MHD run and public coronagraph images, and adds the white-light images of the Large Angle and Spectrometric Coronagraph (LASCO) C2 and of the STEREO coronagraphs COR1 and COR2. Every element is a layer that can be switched on or off (Sections 3 to 9), and Sections 10 and 11 describe the fit of the shell and the checks of the code.

All positions are expressed in solar radii (R⊙) in one of two frames (Section 6). The model, the field lines, the shell, and the radio sources are placed in Heliographic Carrington Cartesian coordinates, whereas an image of an observer is described by the sky-plane axes of that observer.

## 2. Data

The coronal model is the thermodynamic MHD model MAS/MAST of Predictive Science Inc. (PSI; Mikić et al. 1999; Lionello et al. 2009; Mikić et al. 2018). The default run, `cr2143-medium/hmi_mast_mas_std_0201`, covers Carrington rotation 2143 with a boundary condition from a magnetogram of the Helioseismic and Magnetic Imager on the Solar Dynamics Observatory (SDO/HMI). It is read from the PSI archive as HDF4 files of the mass density (`rho002`) and of the three components of the magnetic field (`br002`, `bt002`, `bp002`), whereas of the Alfvén-speed file (`va002`) only the axes of its mesh are used. The run spans 1.0 to 30.5 R⊙. Since the run used by Morosan et al. (2022) is no longer served by PSI, we chose, among the twelve public runs of Carrington rotations 2141 to 2144, the one that reproduces the outline of their density surface best (intersection over union 0.87 to 0.93 in the three views). This choice was therefore made by comparison with the published figure. For a new event, `corona3d new` selects the run of the Carrington rotation of the event, at medium resolution when it is listed in the PSI archive and otherwise at high resolution.

The white-light images are those of the C2 coronagraph of LASCO on the Solar and Heliospheric Observatory (SOHO; Brueckner et al. 1995), with 1024 × 1024 pixels, together with those of the COR1 and COR2 coronagraphs of the Sun Earth Connection Coronal and Heliospheric Investigation (SECCHI; Howard et al. 2008) on STEREO-A or STEREO-B. COR1 provides polarisation triplets (at 0°, 120°, and 240°) and COR2 total-brightness images, all at level 0.5. These images are downloaded from the Virtual Solar Observatory with SunPy (The SunPy Community et al. 2020) for the window `wl_window`, and `--cor1a` adds the COR1 images of STEREO-A for events observed by STEREO-B, which only the GCS fit uses.

The CME geometry is given by the six GCS parameters, fitted to the white-light images with `corona3d gcs` (Section 10). The radio sources are measured in the 2-D images of the Nançay Radioheliograph (NRH; Kerdraon & Delouis 1997) or of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA; Chhabra et al. 2021), or they are sky-plane positions supplied by the user (Section 8).

## 3. The MHD model and the Alfvén speed

The four fields are interpolated linearly onto the cell-centre mesh of the Alfvén-speed file over the whole run. The MAS code units are converted with a number density $\rho_0 = 10^{8}\ \mathrm{cm^{-3}}$ (one proton mass per particle) and a magnetic field $B_0 = 2.2068914\ \mathrm{G}$. Following Morosan et al. (2022), who multiplied the model density by 2 after Wang et al. (2017) to obtain more accurate emission heights, the electron density is

$$ n_e = 10^{8}\,\rho\,s \quad [\mathrm{cm^{-3}}], \qquad s = \texttt{density\_scale} = 2, $$

where $\rho$ is the density in code units. The Alfvén speed,

$$ v_A = \frac{|\mathbf{B}|}{\sqrt{4\pi\, n\, m_p}}, $$

with $m_p$ the proton mass, is computed by default from the scaled density, $n = n_e$ (`va_density` = `"scaled"`), so that it is consistent with the density surface on which it is drawn; this lowers it by a factor of $\sqrt{2}$ relative to the model density. Morosan et al. (2022) do not state which density they used for the Alfvén speed, and `va_density` = `"model"` computes it from the unscaled density instead.

The model is placed in the Carrington frame of the observers without rotation (`model_dlon` = 0). MAS is a steady-state model in Carrington coordinates, and the Carrington longitudes of the observers at the time of the event already contain the solar rotation, so that no further rotation is justified physically. `model_dlon` turns the model about the rotation axis by the given angle (positive towards increasing Carrington longitude) for comparisons with other time conventions.

## 4. The density surface (layer `surface`)

The scaled density is placed on a structured grid in Carrington Cartesian coordinates, and the iso-surface at `ne_level` = $9\times10^{7}\ \mathrm{cm^{-3}}$ is extracted with the contour filter of the Visualization Toolkit (VTK), through PyVista (Sullivan & Kaszynski 2019). This is the density that Morosan et al. (2022) used for the plasma level of 173 MHz emitted at the second harmonic. The plasma frequency and the density are related by

$$ f_p = 8980\sqrt{n_e}\ \mathrm{Hz}, \qquad n_e = \left(\frac{f}{2\times 8980\ \mathrm{Hz}}\right)^{2}, $$

where the factor 2 in the second relation refers to emission at the second harmonic, $f = 2 f_p$. The exact value for 173 MHz is $9.28\times10^{7}\ \mathrm{cm^{-3}}$, and $9\times10^{7}\ \mathrm{cm^{-3}}$ corresponds to 170.4 MHz. The Alfvén speed is interpolated onto the surface and drawn with the `Spectral_r` colour map over `va_range` (100 to 700 km s⁻¹). The surface is not smoothed, and only its normals are computed, so that the shading is smooth while the shape remains as contoured.

## 5. Field lines (layers `open`, `closed`, and `streamer`)

### 5.1 Runge-Kutta tracing

A field line is the solution of

$$ \frac{d\mathbf{x}}{ds} = \frac{\mathbf{B}(\mathbf{x})}{|\mathbf{B}(\mathbf{x})|}, $$

which corona3d integrates with the classical fourth-order Runge-Kutta scheme,

$$ \mathbf{x}_{n+1} = \mathbf{x}_n + \tfrac{1}{6}\,(\mathbf{k}_1 + 2\mathbf{k}_2 + 2\mathbf{k}_3 + \mathbf{k}_4), $$

where $\mathbf{k}_i$ are the usual stage increments of the step $\Delta s = \delta\, r$ with $\delta = 0.01$ and $r$ the distance from Sun centre. The step is thus 0.01 R⊙ at the solar surface and 0.3 R⊙ at 30 R⊙. The field is interpolated trilinearly in the spherical coordinates $(\phi, \theta, r)$ of the model mesh. At most 1500 steps are taken in each direction from the seed, and a line ends when it leaves the model below its inner boundary (1.0 R⊙) or beyond its outer boundary (30.5 R⊙).

A line is classed as open when one of its ends leaves the outer boundary, as closed when both ends return to the inner boundary, and as undetermined otherwise, for example when it runs out of steps or stalls where $\mathbf{B}$ vanishes. Undetermined lines are counted and not drawn. Lines are traced through the whole model for their classification, so that a loop that closes above 3 R⊙ near the cusp of a streamer is classed as closed, but they are drawn only within `line_r_max` = 3 R⊙. We wrote the tracer for corona3d because the tracer of psipy (streamtracer) did not follow these files correctly in our tests, since its lines were not tangent to $\mathbf{B}$, with a median $|\cos(\text{line},\mathbf{B})|$ of 0.83 against 1.000 for the present tracer. The offline tests check the tracer against the analytic field lines of a dipole (Section 11).

### 5.2 Seeds and streamer loops

The open and closed lines are seeded at `n_seeds` = 4500 points distributed uniformly over a sphere of radius `seed_r` = 1.03 R⊙ (a Fibonacci lattice). Open lines are drawn in black and closed lines in grey, both at low opacity, so that the density surface remains visible through them.

Morosan et al. (2022) concluded that the herringbones originate where the CME flank meets the closed loops of a streamer. corona3d selects these streamer loops without any free parameter, by tracing lines through every vertex of the GCS shell above 1.02 R⊙ and drawing the closed ones in purple. These loops are the closed field lines that the CME shell intersects at the time of the reconstruction. For the default event, 478 of the 943 lines through the shell are closed and 465 are open, and the program prints these numbers at every rendering. At most `streamer_max_lines` = 150 of them are drawn, as an even subset.

## 6. Coordinate frames, viewpoints, and projection

### 6.1 Frames

Model, field lines, shell, and radio sources are expressed in Heliographic Carrington Cartesian coordinates in R⊙, centred on the Sun, with the $x$ axis towards Carrington longitude 0° and the $z$ axis along the rotation axis towards solar north. A point at Carrington longitude $\ell$, latitude $b$, and distance $r$ is

$$ \mathbf{x} = r\,(\cos b\cos\ell,\ \cos b\sin\ell,\ \sin b). $$

This frame rotates with the Sun, so that the positions refer to the time of the reconstruction. The Stonyhurst longitude of a point is its Carrington longitude minus that of Earth at that time.

The observers are given by their Carrington longitude and latitude (`earth` and `stereo`). For the event of 2013 October 25 at 15:00 UT, Earth is at (357.372°, 4.993°) and STEREO-B at (215.329°, −0.984°), 142° apart. These values were computed with SunPy from the Horizons ephemerides of the Jet Propulsion Laboratory (JPL), and `corona3d new` computes them in the same way for a new event. For each observer with the unit vector $\hat{\mathbf{z}}$ towards the observer, corona3d uses the sky-plane axes

$$ \hat{\mathbf{x}} \ \text{towards the west limb},\qquad \hat{\mathbf{y}} \ \text{towards solar north (projected on the sky plane)},\qquad \hat{\mathbf{z}} = \hat{\mathbf{x}}\times\hat{\mathbf{y}}, $$

an orthonormal, right-handed set.

### 6.2 Viewpoints and projection

Each of the three panels is given by an entry of `views`: an observer (`"from": "earth"`, `"stereo"`, or `"north"`, the last looking down the rotation axis with the direction of the Earth at the bottom of the panel), or any viewing direction given by a Carrington longitude and latitude, with solar north up and an optional roll angle. The panels are drawn with an orthographic (parallel) projection, so that every depth along the line of sight is drawn at the scale of the plane of the sky, which keeps the model, the shell, the radio sources, and the coronagraph images at one common scale in R⊙. In contrast, a perspective camera at 6.65 R⊙ would magnify structures in front of the plane of the sky by up to 18 %. Such a camera is nevertheless available with `"projection": "perspective"` at `camera_distance` = 6.65 R⊙, so that a figure can be drawn with the perspective of the published Fig. 7 for comparison. The half-height and the offset of Sun centre of each panel (`frame` in each view) follow the layout of the published figure. When coronagraph images are shown, the panel of their spacecraft is widened to 6 R⊙ (C2 or COR2) or 4 R⊙ (COR1 alone), while `--half` sets one half-height for all panels.

### 6.3 Plane-of-sky projection

A point $\mathbf{x}$ of the reconstruction appears in the panel of an observer at the sky-plane position

$$ X = \mathbf{x}\cdot\hat{\mathbf{x}}, \qquad Y = \mathbf{x}\cdot\hat{\mathbf{y}}, $$

and its depth is $Z = \mathbf{x}\cdot\hat{\mathbf{z}}$, which an orthographic projection discards. Conversely, a source measured only in the plane of the sky has no depth, and Section 8 describes how its depth is assigned.

## 7. The CME shell (layer `gcs`)

### 7.1 Geometry

The shell is the GCS model of Thernisien et al. (2006) and Thernisien (2011), built with PyThea (Kouloumvakos et al. 2022). Its parameters are the apex height $h$ (the distance of the apex from Sun centre), the aspect ratio $\kappa$, the half-angle $\alpha$ between the legs, the Stonyhurst longitude and latitude of the axis, and the tilt. In the notation of PyThea, the height of the cone formed by the legs is

$$ h_c = h\,\frac{(1-\kappa)\cos\alpha}{1+\sin\alpha}, $$

the radius of the circular cross-section at the apex is

$$ a_{\rm apex} = \kappa\,\frac{h_c/\cos\alpha + h_c\tan\alpha}{1-\kappa^{2}}, $$

and the centre of the apex circle lies at the distance $r_c = h - a_{\rm apex}$ from Sun centre. The cross-section radius grows with distance along the axis at the rate $\kappa$ (the aspect ratio of the model). corona3d takes the tilt in its own convention, $\texttt{tilt\_deg} = \texttt{tilt}_{\rm PyThea} - 90^\circ$, and passes it to PyThea plus 90°. With half the default number of vertices along the axis and around it (`mesh_rings`: 5 straight and 10 circular), the shell is drawn as rings around its axis in magenta.

### 7.2 Fit to the coronagraph images

Morosan et al. (2022) fitted the shell to white-light images at 15:24 UT and scaled it back in time to the radio observations at 15:00 UT, but they do not list its parameters. The command `corona3d gcs` therefore fits the six parameters in two steps that follow this procedure.

In the first step (at the time `--t-ref`), all six parameters are fitted by differential evolution to the CME envelopes in the chosen instruments, in the frame nearest to `--t-ref`. Each envelope is the region of the base difference (or, with `--running-diff`, of the running difference) brighter than three robust standard deviations above its median, connected to the radial line above the flare site and filled to its convex hull. The shell is projected into each instrument with PyThea, and the objective is the mean intersection over union (IoU),

$$ \mathrm{IoU} = \frac{|S\cap E|}{|S\cup E|}, $$

of the projected silhouette $S$ with the envelope $E$ in each instrument, evaluated on the usable pixels (outside the occulter and inside the field of view); each frame is compared with the first frame of its own observing programme. With `--fix-direction`, the longitude and latitude of the axis are held at the flare site and only the other four parameters are fitted.

In the second step (at the time `--t-out`, the time of the reconstruction), the shape ($\kappa$, $\alpha$, direction, and tilt) is kept, and the apex height is extrapolated from the height-time profile of the CME front. The front is measured in the instrument whose line of sight is most nearly perpendicular to the fitted direction, as the largest distance of the envelope from Sun centre, divided by the projection factor of that direction, in the uninterrupted run of frames of the observing programmes in which the envelope covers at least 2 R⊙² (a frame below that area ends the run before it, so that small envelopes before the eruption are not used). A straight line through the measured heights gives the speed and the height at `--t-out`, with the standard error of the fitted intercept as its uncertainty. Alternatively, `--height-out` gives the height by hand, and it is then recorded as assumed.

For the event of 2013 October 25 this gives, at 15:24 UT, an apex height of 5.538 R⊙, $\kappa$ = 0.396, $\alpha$ = 23.46°, a direction of Stonyhurst longitude −53.262° and latitude 3.146°, and a tilt of −56.654° in the PyThea convention, with IoU values of 0.76 (C2), 0.85 (COR1-A), and 0.74 (COR1-B). The front in COR1-B (projection factor 1.00) lies at 2.646, 3.192, and 3.542 R⊙ at 15:08:21, 15:13:21, and 15:18:21 UT, and a straight line through these points gives a speed of 1039 km s⁻¹ and an apex height at 15:00 UT of 1.931 ± 0.158 R⊙. The fitted values are kept in `scripts/gcs_fit_2013-10-25.json`.

The base frame is the earliest frame of each observing programme that the folder holds, so that the envelope areas, and with them the height-time points, depend on how early the download starts. For the COR1-B frames of 2013 October 25, the folder used for `scripts/gcs_fit_2013-10-25.json` starts at 14:43 UT and gives an envelope of 2.66 R⊙² at 15:08:21, whereas a folder that starts at 14:35 UT (the window of `corona3d fetch`) gives 2.16 R⊙²; the same code, applied to both folders, reproduces both values. The change of the base frame by eight minutes alters the front height at that frame by 0.03 R⊙; with the different run of frames that results (four frames from 15:08 instead of three), the extrapolated apex height at 15:00 UT is 1.80 R⊙ at 1126 km s⁻¹ instead of 1.93 R⊙ at 1039 km s⁻¹, which agrees within the uncertainty of 0.16 R⊙. Where the result must not depend on the download window, `--running-diff` removes the dependence on the first frame.

### 7.3 Views that disagree with the others

A single instrument with a poor envelope can pull the fit away from the others, as for 2024 May 9, where the fit to COR1-A and COR2-A gave IoU values of 0.548 and 0.295. After the first fit, `corona3d gcs` and `corona3d event` therefore apply the following rule (`gcs.check_views`), which is active by default and is switched off with `--keep-all-views`.

1. A view is poor when its IoU is below 0.5 (`--iou-min` in `corona3d event`) and below 0.75 times the mean IoU of the other views. Of several poor views the one with the lowest IoU is taken, and a single view is never poor, since there is nothing to compare it with.
2. The six (or four) parameters are fitted again without that view.
3. The refit is kept only when the mean IoU over the remaining views rises by at least 0.05 compared with the fit to all views evaluated on the same views. Otherwise the poor view did not bias the shell, the fit to all views is kept, and the refit is only reported.
4. When the refit is kept, the rule is applied again to the remaining views, so that views are removed one at a time, the worst first. The number of views never falls below one.

Both fits are written to the key `view_check` of `gcs_fit.json` (the instruments, the IoU of each, the mean, and whether the refit was kept) and to the section *Fit quality* of `event_report.md`, and the instruments left out are listed in `views_dropped` of the record `gcs_fit` of `config.json`. The IoU and the quality criteria of Section 7.4 then refer to the views that were kept, and the front of the height-time step is chosen among them. The limits of the rule are that the thresholds are conventions and not statistics, that a fit to one view leaves the direction and the extent along the line of sight poorly constrained (a fixed-direction fit is more reliable then), and that the criterion tests only the agreement of the shell with the envelopes, not whether the dropped envelope was the correct one.

### 7.4 Quality criteria

A fit is judged poor, by `gcs.quality`, when the mean IoU over the instruments is below `--iou-min` (default 0.5), when no instrument gave an envelope, or when a fitted parameter lies within 0.5 % of its range from a bound of the optimisation. `corona3d event` uses these criteria to decide whether the fixed-direction fit replaces the free one, and the report states both outcomes.

### 7.5 Shells that are drawn

The shell at the time of the reconstruction (layer `gcs`, key `gcs`) is drawn when no coronagraph image is shown. When an image is shown, the shell at the time of the images (key `gcs_wl`) is drawn instead, so that the shell and the images belong to the same time, and each panel states the time of its shell. The extrapolation of the height is based on a few frames and assumes a constant speed. Since an accelerating CME moves more slowly early in its rise, the true height at the earlier time may lie below the extrapolated value.

## 8. Radio sources (layer `radio`)

### 8.1 Centroids

The main input of the layer is the 2-D images of the NRH (Kerdraon & Delouis 1997), the Stokes I brightness temperature of 64 × 64 pixels at one frequency about every 10 s, rotated to solar north up, with Sun centre at the reference pixel and one solar radius spanning a number of pixels given in the header. For each frequency, one image is taken every `step_s` seconds between `start` and `end`, and an image is kept when its peak brightness temperature reaches `min_tb` ($10^{8}$ K). In each kept image, the centroid is the brightness-weighted mean position of the pixels brighter than `level` (0.5) times the peak. The uncertainty of a centroid is of the order of one pixel, 1/16 R⊙, in the plane of the sky. The spectral images of OVRO-LWA (Chhabra et al. 2021), available from 2024 April onward at 32 to 87 MHz, are treated in the same way. Each file holds one 10 s integration with 12 band images, stored at reduced size and restored here to the full grid of 256 × 256 pixels by cubic interpolation. The helioprojective header (reference pixel, plate scale, reference coordinates, and angular solar radius) gives the position of every pixel in R⊙, and for each requested frequency the nearest band within 3 MHz is used. The centroids are drawn as a chain of blobs coloured from the first to the last image, as Morosan et al. (2022) drew their sources at successive times.

The interval of the burst can be chosen from the images themselves with `corona3d radio-interval`: for each frequency, the images within the search window whose brightest pixel reaches `min_tb` are grouped into runs separated by gaps of more than $\max(60\ \mathrm{s},\,3\times\text{median cadence})$, and the run that holds the brightest image is taken as the burst.

### 8.2 Newkirk heights

Other radio sources, and the centroids, are given as sky-plane positions seen from Earth ($x$ west, $y$ north; in R⊙, or in arcsec with `"unit": "arcsec"`, converted with the angular radius of the Sun seen from Earth at `time`, 964.8 arcsec at 15:00 UT on 2013 October 25), each with its frequency. Each source is placed at the radial distance where the four-fold Newkirk (1961) density model,

$$ n_e(r) = 4 \times 4.2\times10^{4}\times 10^{4.32/r}\ \mathrm{cm^{-3}}, $$

with $r$ in R⊙, has a plasma frequency $f_p = 8980\sqrt{n_e}$ Hz equal to half the observed frequency (emission at the second harmonic, as assumed by Morosan et al. 2022). Solving for $r$ gives

$$ r(f) = \frac{4.32}{\log_{10}\!\big[n_e(f)\,/\,(4\times 4.2\times 10^{4})\big]}, \qquad n_e(f) = \left(\frac{f}{2\times 8980\ \mathrm{Hz}}\right)^{2}, $$

which is 1.576 R⊙ at 173 MHz and 1.449 R⊙ at 228 MHz (the Python function `corona3d.newkirk_r` evaluates it). At OVRO-LWA frequencies, the heights lie between 2.04 R⊙ (84.6 MHz) and 3.25 R⊙ (34 MHz).

### 8.3 Deprojection

The point is taken on the line of sight of the observed position, on the side of the plane of the sky that faces the observer, that is,

$$ Z = +\sqrt{r^{2} - X^{2} - Y^{2}}, $$

so that $X$ and $Y$ are kept exactly. A position that lies farther from Sun centre than its height, $X^{2}+Y^{2}>r^{2}$, cannot be placed at that height, and it is therefore kept in the plane of the sky with a printed warning. The heights follow the Newkirk model and not the density of the MHD model, so that a source and the density surface locate the same plasma frequency at different heights wherever the two density models differ. The scattering of radio waves in the corona, which displaces observed sources outward, is not corrected. The sources are drawn as spheres coloured by frequency.

## 9. Coronagraph images (layers `c2`, `cor1`, and `cor2`)

For each instrument, the frame nearest to `wl_time` within `wl_window` and the first frame of that window are read and normalised by their exposure times. For COR1, the total brightness of each triplet is formed as $tB = \tfrac{2}{3}(I_0 + I_{120} + I_{240})$, where the three exposures of a triplet are those taken within 60 s. Each frame is binned by block averages to about 256 × 256 pixels, and a frame of another size is resampled onto the grid of the base frame. The comparison is a ratio by default (`wl_mode` = `"ratio"`),

$$ \frac{I}{I_{\rm base}} - 1, $$

which removes multiplicative instrumental effects such as vignetting and the stray light that scales with the coronal brightness, whereas `"difference"` gives $I - I_{\rm base}$. A `wl_time` outside the window, or a frame that would be the base frame itself, stops the rendering with a message, since either would draw a wrong or blank image. When the nearest frame to `wl_time` lies more than 30 min away, which happens at a data gap, the rendering stops as well, and the message names the nearest frame times before and after `wl_time`; `corona3d gcs --t-ref` and `corona3d event --t-ref` then set `wl_time` to the chosen time. The result is rotated to solar north up with bilinear interpolation. The position of each pixel in the plane of the sky follows from the world coordinate system of the image and the angular radius of the Sun seen by the spacecraft,

$$ (X, Y) = (\text{pixel} - \text{Sun centre}) \times \frac{\text{plate scale}}{\text{angular solar radius}}\ \mathrm{R_\odot}. $$

The image is drawn on the plane through Sun centre perpendicular to the line of sight of its spacecraft, at the point $X\hat{\mathbf{x}} + Y\hat{\mathbf{y}}$ of the observer's axes (Section 6), and seen from that spacecraft in the orthographic projection it therefore covers the model at its true scale and orientation. Pixels inside the occulter and outside the field of view are transparent: LASCO C2 from 2.2 to 6 R⊙, COR1 from 1.5 to 4 R⊙, and COR2 from 4 to 15 R⊙, so that COR2 begins where COR1 ends and the two coplanar images do not overlap. The grey scale spans plus and minus the 98th percentile of the absolute values, and the opacity is `wl_opacity` = 0.75. LASCO C2 is drawn only in a view from Earth, and COR1 and COR2 only in a view from the STEREO spacecraft of the event.

## 10. From a time to a reconstruction

For an event given only by its time and flare position, `corona3d event` runs the steps in the order of the earlier sections: it writes the configuration (observers, model run, and image windows, Section 2 and 6), downloads the data, chooses the radio interval (Section 8.1), fits the shell (Section 7.2), and draws the figures. The free fit comes first, and the fit with the direction held at the flare site follows when the free fit fails or is judged poor (Section 7.4). A halo CME, a data gap, or a small angular separation between Earth and the STEREO spacecraft limits the free fit, and [LIMITATIONS.md](LIMITATIONS.md) describes these cases. The file `event_report.md` separates what was measured from what was assumed, lists the missing data, and states the fit quality.

## 11. Validation

The offline tests (`python -m pytest`) use synthetic data with known answers:

* **Field-line tracer.** The tracer is checked against a dipole field, for which a line obeys $r = L\sin^{2}\theta$.
  * With the default step and step count, the line from 1.05 R⊙ at colatitude 45° is closed. Its apex lies at 2.1 R⊙ within 0.02 R⊙, it returns to the inner boundary at latitude −46.4°, and $L$ is constant along it to 1.5 %.
  * The line from colatitude 8° is open.
  * A line given only 20 steps is undetermined.
* **MAS units and Alfvén speed.** On a synthetic run with uniform density and a dipole field, the density, $|\mathbf{B}|$, and $v_A$ equal their analytic values in code units times $10^{8}\ \mathrm{cm^{-3}}$ and 2.2068914 G. $v_A$ from the model density is $\sqrt{2}$ times $v_A$ from the doubled density.
* **Model rotation.** A density enhancement at Carrington longitude 0° appears at longitude `model_dlon` (0° and 40°) on the density surface.
* **Streamer loops.** The loops selected through the shell are closed, with both feet on the inner boundary, and pass within 0.05 R⊙ of a vertex of the shell.
* **Radio sources.** The radio heights are 1.576 R⊙ at 173 MHz and 1.449 R⊙ at 228 MHz. A position farther from Sun centre than its height produces a warning.
* **Observer axes.** The axes are orthonormal and right-handed, with $\hat{\mathbf{x}}$ towards the west limb, and the deprojected sources lie at their radius with their sky-plane coordinates unchanged.
* **Coronagraph planes.**
  * A feature placed at (5, 2) R⊙ west and north in synthetic COR2 frames, at roll angles of 0°, 6.6°, 90°, 175.3°, and 271°, is found in three dimensions at its true position within 0.15 R⊙ for both the ratio and the difference. It is also found there when the base frame has half the size of the frame shown, and the occulter is transparent.
  * Frames smaller than 256 pixels are not binned.
  * A `wl_time` outside the window, the base frame itself, an empty window, and a gap of more than 30 min are refused, and the message names the nearest frame times.
  * A time with a time zone is accepted.
* **Frame selection.** COR1 frames are grouped into triplets by time, including a triplet across a minute boundary, and incomplete triplets are dropped. C2 frames are selected from their DETECTOR, POLAR, and NAXIS1 keywords, and files without them are skipped.
* **Event pipeline and fit.** The order of the steps, the skipping of finished work, the fallback to the fixed-direction fit, the report, the interval chosen from synthetic NRH images, the height-time functions on a synthetic rise with a known speed, and the effect of `--t-ref` on `wl_time` are tested.
* **Command line and configuration.**
  * The command line refuses to work outside a working directory, and `corona3d init` leaves an edited configuration unchanged.
  * Partial nested keys of `config.json` are merged with the defaults, and an unknown layer is refused from Python.
  * Halving `--half` doubles the radius of the Sun in the panel to within 3 %, and a zero half-height is refused.
  * The viewpoints (observers, Carrington directions, roll, north pole, and the legacy names `stereo_a` and `stereo_b`) give the expected directions and screen orientation.
  * The layers that need no data render.
* **Downloads.** The model download skips files that are present, writes a file under its final name only when it is complete, and removes the partial file of a failed download. A truncated image is detected, deleted, and counted as failed, so that the next `corona3d fetch --wl` downloads it again.

With `CORONA3D_ONLINE=1`, the online test also downloads the public data of 2013 October 25 and draws both figures of the README.

## References

Brueckner, G. E., Howard, R. A., Koomen, M. J., et al. 1995, Solar Physics, 162, 357. https://doi.org/10.1007/BF00733434

Chhabra, S., Gary, D. E., Hallinan, G., et al. 2021, The Astrophysical Journal, 906, 132. https://doi.org/10.3847/1538-4357/abc94b

Howard, R. A., Moses, J. D., Vourlidas, A., et al. 2008, Space Science Reviews, 136, 67. https://doi.org/10.1007/s11214-008-9341-4

Kerdraon, A., & Delouis, J.-M. 1997, in Coronal Physics from Radio and Space Observations, Lecture Notes in Physics, 483, 192. https://doi.org/10.1007/BFb0106458

Kouloumvakos, A., Rodriguez-Garcia, L., Gieseler, J., et al. 2022, Frontiers in Astronomy and Space Sciences, 9, 974137. https://doi.org/10.3389/fspas.2022.974137

Lionello, R., Linker, J. A., & Mikić, Z. 2009, The Astrophysical Journal, 690, 902. https://doi.org/10.1088/0004-637X/690/1/902

Mikić, Z., Linker, J. A., Schnack, D. D., Lionello, R., & Tarditi, A. 1999, Physics of Plasmas, 6, 2217. https://doi.org/10.1063/1.873474

Mikić, Z., Downs, C., Linker, J. A., et al. 2018, Nature Astronomy, 2, 913. https://doi.org/10.1038/s41550-018-0562-5

Morosan, D. E., Pomoell, J., Kumari, A., Vainio, R., & Kilpua, E. K. J. 2022, Astronomy & Astrophysics, 668, A15. https://doi.org/10.1051/0004-6361/202244432

Newkirk, G., Jr. 1961, The Astrophysical Journal, 133, 983. https://doi.org/10.1086/147104

Sullivan, C. B., & Kaszynski, A. A. 2019, Journal of Open Source Software, 4, 1450. https://doi.org/10.21105/joss.01450

The SunPy Community, Barnes, W. T., Bobra, M. G., et al. 2020, The Astrophysical Journal, 890, 68. https://doi.org/10.3847/1538-4357/ab4f7a

Thernisien, A. F. R., Howard, R. A., & Vourlidas, A. 2006, The Astrophysical Journal, 652, 763. https://doi.org/10.1086/508254

Thernisien, A. 2011, The Astrophysical Journal Supplement Series, 194, 33. https://doi.org/10.1088/0067-0049/194/2/33

Wang, T., Reginald, N. L., Davila, J. M., St. Cyr, O. C., & Thompson, W. T. 2017, Solar Physics, 292, 97. https://doi.org/10.1007/s11207-017-1130-3
