# Methods of corona3d

This document describes how corona3d computes every element of its reconstruction, which values it uses, and where each value comes from. Every number below that is written in `code font` is a key of `config.json`, which the user can change. Each method and value is also listed with its source and its test in [MODEL_AUDIT.md](MODEL_AUDIT.md), the assumptions that limit the results are collected in [LIMITATIONS.md](LIMITATIONS.md), and the references are in [refs.bib](refs.bib); every reference with a digital object identifier (DOI) was checked against its registry record with Argus.

## 1. Overview

The interaction of a coronal mass ejection (CME) with the surrounding corona depends on quantities that vary in all three dimensions. The electron density sets the height at which a given radio frequency is emitted, the magnetic field and the density together set the Alfvén speed that a shock must exceed, and the topology of the field decides whether accelerated electrons escape along open field or remain on closed loops. Morosan et al. (2022) brought these quantities together in one three-dimensional reconstruction of the CME of 2013 October 25, which combines a magnetohydrodynamic (MHD) model of the corona, a graduated cylindrical shell (GCS) of the CME, and the positions of herringbone radio sources, seen from Earth, from above the solar north pole, and from the Solar Terrestrial Relations Observatory B (STEREO-B). Such a reconstruction, however, has so far been made by hand for single events, and its steps are rarely available in a form that can be repeated for another event. corona3d follows the method described in Morosan et al. (2022) and makes each of its steps a reproducible computation from public data for any event, namely the MHD model of Predictive Science Inc. (PSI) rendered in three dimensions, a GCS shell fitted to coronagraph images, the white-light images of the Large Angle and Spectrometric Coronagraph (LASCO) C2 and of the STEREO coronagraphs COR1 and COR2, and radio sources measured in radio images and placed at the heights of a coronal density model on the side of the plane of the sky that faces the observer. Every element is a layer that can be switched on or off (Sections 3 to 9), and the shell is fitted as described in Section 7. Moreover, Section 10 describes how these steps are run for an event given only by its time and flare position, and Section 11 describes the checks of the code and the validation against the event of 2013 October 25, which is the benchmark of corona3d.

All positions are expressed in solar radii (R⊙) in one of two frames (Section 6). The model, the field lines, the shell, and the radio sources are placed in Heliographic Carrington Cartesian coordinates, whereas an image of an observer is described by the sky-plane axes of that observer.

## 2. Data

The coronal model is the thermodynamic MHD model MAS/MAST of PSI (Mikić et al. 1999; Lionello et al. 2009; Mikić et al. 2018). The default run, `cr2143-medium/hmi_mast_mas_std_0201`, is that of the benchmark event and covers Carrington rotation 2143 with a boundary condition from a magnetogram of the Helioseismic and Magnetic Imager on the Solar Dynamics Observatory (SDO/HMI). It is read from the PSI archive as HDF4 files of the mass density (`rho002`) and of the three components of the magnetic field (`br002`, `bt002`, `bp002`), whereas of the Alfvén-speed file (`va002`) only the axes of its mesh are used. The run spans 1.0 to 30.5 R⊙. Since the run used by Morosan et al. (2022) is no longer served by PSI, we chose, among the twelve public runs of Carrington rotations 2141 to 2144, the one that reproduces the outline of their density surface best (intersection over union 0.87 to 0.93 in the three views). This choice was therefore made by comparison with the published figure, and it applies to the benchmark alone. For a new event, `corona3d new` selects the run of the Carrington rotation of the event, at medium resolution when it is listed in the PSI archive and otherwise at high resolution.

The white-light images are those of the C2 coronagraph of LASCO on the Solar and Heliospheric Observatory (SOHO; Brueckner et al. 1995), with 1024 × 1024 pixels, together with those of the COR1 and COR2 coronagraphs of the Sun Earth Connection Coronal and Heliospheric Investigation (SECCHI; Howard et al. 2008) on STEREO-A or STEREO-B. COR1 provides polarisation triplets (at 0°, 120°, and 240°) and COR2 total-brightness images, all at level 0.5. These images are downloaded from the Virtual Solar Observatory with SunPy (The SunPy Community et al. 2020) for the window `wl_window`, and `--cor1a` adds the COR1 images of STEREO-A for events observed by STEREO-B, which only the GCS fit uses.

The CME geometry is given by the six GCS parameters, fitted to the white-light images with `corona3d gcs` (Section 7). The radio sources are measured in the 2-D images of the Nançay Radioheliograph (NRH; Kerdraon & Delouis 1997) or of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA; Chhabra et al. 2021), or they are sky-plane positions supplied by the user (Section 8).

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

Each of the three panels is given by an entry of `views`: an observer (`"from": "earth"`, `"stereo"`, or `"north"`, the last looking down the rotation axis with the direction of the Earth at the bottom of the panel), or any viewing direction given by a Carrington longitude and latitude, with solar north up and an optional roll angle. The panels are drawn with an orthographic (parallel) projection, so that every depth along the line of sight is drawn at the scale of the plane of the sky, which keeps the model, the shell, the radio sources, and the coronagraph images at one common scale in R⊙. In contrast, a perspective camera at 6.65 R⊙ would magnify structures in front of the plane of the sky by up to 18 %. Such a camera is nevertheless available with `"projection": "perspective"` at `camera_distance` = 6.65 R⊙, so that a figure can be drawn with the perspective of the published Fig. 7 for comparison. In the configuration of the benchmark event, the half-height and the offset of Sun centre of each panel (`frame` in each view) follow the layout of the published figure, and for any other event they are set freely. When coronagraph images are shown, the panel of their spacecraft is widened to 6 R⊙ (C2 or COR2) or 4 R⊙ (COR1 alone), while `--half` sets one half-height for all panels.

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

In the first step (at the time `--t-ref`), all six parameters are fitted by differential evolution to the CME envelopes in the chosen instruments, in the frame nearest to `--t-ref`. Each envelope is the region of the base difference (or, with `--running-diff`, of the running difference) brighter than three robust standard deviations above its median, connected to the radial line above the flare site and filled to its convex hull. Before the threshold, the pixels of the base frame that are more than three times brighter than the median of their annulus of 0.1 R⊙ are masked (`gcs.SAT_FACTOR`), because saturated or stray-light pixels near the occulter, such as the band south of the COR1-A occulter on 2024 May 9, leave a flat region in every difference that interrupts the envelope. The difference is then smoothed over the usable pixels with a Gaussian of 0.06 R⊙ (`gcs.ENV_SMOOTH`), so that the smoothing has the same physical scale in every instrument, and the thresholded region is opened with a radius of 0.05 R⊙ (`gcs.ENV_OPEN`), and bright regions smaller than 0.1 R⊙² (`gcs.ENV_MIN_AREA`) are left out, so that noise specks connected to the CME do not widen its convex hull. For the event of 2013 October 25, a fit repeated with these envelopes gives a mean IoU of 0.825 instead of 0.804 and nearly the same shell (apex height 5.37 R⊙ at 15:24 UT instead of 5.25 R⊙, and 1.83 R⊙ at 15:00 UT instead of 1.94 R⊙, within the uncertainty of 0.6 R⊙), and the default parameters of the 2013 figure are not changed. For 2024 May 9, the COR1-A envelope becomes developed at 17:51 UT instead of 18:06 UT, and the fit to COR1-A reaches an IoU of 0.89 instead of 0.59, with the axis at latitude −23° in the south, as the flare site S14W28, instead of +1°. The apex of 5.1 R⊙ still lies outside the field of COR1-A, so that this fit remains poor. The shell is projected into each instrument with PyThea, and the objective is the mean intersection over union (IoU),

$$ \mathrm{IoU} = \frac{|S\cap E|}{|S\cup E|}, $$

of the projected silhouette $S$ with the envelope $E$ in each instrument, evaluated on the usable pixels (outside the occulter and inside the field of view); each frame is compared with the base frame of its own observing programme (below). With `--fix-direction`, the longitude and latitude of the axis are held at the flare site and only the other four parameters are fitted.

The bounds of the optimisation are 0.1 to 0.8 for $\kappa$, 0 to 60° for $\alpha$, ±90° for the tilt, and the direction within ±40° in longitude and ±35° in latitude of the flare site (or within the cone of Section 7.4). The apex height lies between 1.5 R⊙ and an upper bound that depends on the views being fitted: it is 1.5 times the largest outer field-of-view radius among them (`gcs.HEIGHT_FOV_FACTOR`, `gcs.height_limit`), but never more than 15 R⊙. The bound is 6 R⊙ for COR1 alone, 9 R⊙ when C2 is among the views, and 15 R⊙ for COR2. The reason is that an apex far outside the field is not seen in any image, so that the fit can match the envelope only with the legs of a thin shell and the height is not measured. On 2024 May 9 the fit to COR1-A alone, whose field ends at 4 R⊙, had reached an apex height of 14.6 R⊙ with a half-angle of 1.7° before this limit was introduced. The bound is evaluated for the views of each fit, so that it applies to the fit to all views, to every refit of the view rule (Section 7.3), to a fit to one view, and to the fixed-direction fit, and the value that was used is written to `height_bound` of `gcs_fit.json` and of the record `gcs_fit` in `config.json` and named in `event_report.md`.

In the second step (at the time `--t-out`, the time of the reconstruction), the shape ($\kappa$, $\alpha$, direction, and tilt) is kept, and the apex height is extrapolated from the height-time profile of the CME front. The front is measured in the instrument whose line of sight is most nearly perpendicular to the fitted direction, as the largest distance of the envelope from Sun centre, divided by the projection factor of that direction, in the uninterrupted run of frames of the observing programmes in which the envelope covers at least 2 R⊙² (a frame below that area ends the run before it, so that small envelopes before the eruption are not used). A straight line through the measured heights gives the speed and the height at `--t-out`, with the standard error of the fitted intercept as its uncertainty. Alternatively, `--height-out` gives the height by hand, and it is then recorded as assumed.

For the event of 2013 October 25 this gives, at 15:24 UT, an apex height of 5.538 R⊙, $\kappa$ = 0.396, $\alpha$ = 23.46°, a direction of Stonyhurst longitude −53.262° and latitude 3.146°, and a tilt of −56.654° in the PyThea convention, with IoU values of 0.76 (C2), 0.85 (COR1-A), and 0.74 (COR1-B). The front in COR1-B (projection factor 1.00) lies at 2.646, 3.192, and 3.542 R⊙ at 15:08:21, 15:13:21, and 15:18:21 UT, and a straight line through these points gives a speed of 1039 km s⁻¹ and an apex height at 15:00 UT of 1.931 ± 0.158 R⊙. The fitted values are kept in `scripts/gcs_fit_2013-10-25.json`.

When `--t-ref` is not given, the time of the fit is chosen from the images. For every frame of each instrument within 150 min after the time of the reconstruction, the envelope is formed as above, and its area and the largest distance of its pixels from Sun centre are measured. At a candidate time, an instrument counts as developed when its frame within 6 min of that time has an envelope of at least 3 $R_\odot^2$ whose front lies at most 0.9 of the outer radius of its field of view from Sun centre (and at least 0.2 $R_\odot$ inside the edge). The front is the largest distance of the envelope from Sun centre, the same measurement that `front_profile` uses in step 2. The limit of 0.9 keeps the apex of the shell inside the field of the view from which t_ref is chosen, because a front near the edge of the field leaves the apex height unconstrained by the image. The earliest candidate time with the largest number of instruments that satisfy both conditions is taken (`gcs.pick_t_ref`, `gcs.FRONT_FRAC`). When no time satisfies the limit of 0.9, the earlier rule (developed envelope, front 0.2 $R_\odot$ inside the edge) is kept, `front_inside_fov` is false in `t_ref_choice`, and the report says that the apex height is then poorly constrained. When no frame qualifies at all, the time of the reconstruction plus 60 min is used, and the report says that it was not chosen from the images. A time given with `--t-ref` is never changed. The chosen time is never earlier than the time of the reconstruction and always follows the base frame of the view, since frames at or before the base frame are not measured. The height-time front of step 2 does not depend on t_ref, and `wl_time` follows t_ref. The time of the reconstruction itself, to which the shell is scaled and at which the radio sources are placed, is always the time given with `--time`, which is the begin time of the burst when `--burst` is used. The choice is recorded in the key `t_ref_choice` of `gcs_fit`.

The base frame is fixed by a pre-event window and not by the start of the download. For each instrument and observing programme, the base is the last frame at or before the cutoff time

$$ t_{\rm base} = t_{\rm onset} - \Delta t, $$

where $t_{\rm onset}$ is the flare onset given with `--onset` (default: the time of the reconstruction, `--t-out`) and $\Delta t$ is given with `--base-minutes` (default 10 min). Since the base depends only on $t_{\rm onset}$ and $\Delta t$, a download that starts earlier or later gives the same envelopes, as long as it contains the frame before the cutoff. When a programme has no frame that early, its earliest frame is used instead, the command prints a warning, and the programme is listed in `fallback_earliest`. The record `base` of `gcs_fit` in `config.json` (and of `gcs_fit.json`) holds the onset, $\Delta t$, the cutoff, the time of the base frame of every instrument and programme, and the fallbacks, so that the fit can be repeated. Frames at or before the base frame are not used in the height-time profile. The rule replaces the earlier choice of the earliest frame in the folder, with which the envelope areas, and with them the height-time points, depended on how early the download started: for the COR1-B frames of 2013 October 25, a folder that starts at 14:43 UT gave an envelope of 2.66 R⊙² at 15:08:21, whereas the folder of `corona3d fetch` (from 14:35 UT) gave 2.16 R⊙², and the extrapolated apex height at 15:00 UT changed from 1.93 R⊙ at 1039 km s⁻¹ to 1.80 R⊙ at 1126 km s⁻¹, within the uncertainty of 0.16 R⊙. The base frame should lie before the flare onset and before any earlier ejecta that are visible in the field of view, and `--running-diff` removes the dependence on a base frame altogether. The rule is tested offline on synthetic frame lists (`tests/test_gcs_fit.py`).

The coronagraph figure can also be drawn at the time of the burst (`wl_at` burst, `reconstruction.burst_shell`). The shell is then the one at the reconstruction time, with the apex height of that time, which is called assumed whenever the fit or the user (`wl_height`) did not measure it, and each coronagraph contributes its own frame nearest to the burst that follows its base frame (`whitelight.nearest_frame`). A frame farther than 30 min from the burst is used when it is the nearest one, and the panel and the event report state its distance in minutes, because a shell at the burst time over an image taken much later is then not a simultaneous comparison.

### 7.3 Views that disagree with the others

A single instrument with a poor envelope can pull the fit away from the others, as for 2024 May 9, where the fit to COR1-A and COR2-A gave IoU values of 0.548 and 0.295. After the first fit, `corona3d gcs` and `corona3d event` therefore apply the following rule (`gcs.check_views`), which is active by default and is switched off with `--keep-all-views`.

1. A view is poor when its IoU is below 0.5 (`--iou-min` in `corona3d event`) and below 0.75 times the mean IoU of the other views. Of several poor views the one with the lowest IoU is taken, and a single view is never poor, since there is nothing to compare it with.
2. The six (or four) parameters are fitted again without that view.
3. The refit is kept unless the mean IoU over the remaining views falls by more than 0.01 compared with the fit to all views evaluated on the same views. A view that the shell cannot match, for instance COR2-A with an IoU of 0 on 2024 May 9, costs the other views nothing (COR1-A gave 0.9009 with it and 0.9017 without it), but it halves the mean IoU of the fit and would make an acceptable fit appear poor, and for this reason a gain is not required. When the mean IoU of the remaining views falls, the poor view did constrain the shell, the fit to all views is kept, and the refit is only reported.
4. When the refit is kept, the rule is applied again to the remaining views, so that views are removed one at a time, the worst first. The number of views never falls below one.

When the observers of the views lie within 30° of each other, as for COR1-A and COR2-A of the same spacecraft on 2024 May 9, the fit to all views is not allowed to decide which view survives. Such a pair constrains the direction only weakly, so the view with the lower IoU in the joint fit is not necessarily the poorer one: on 2024 May 9 the joint fit gave 0.28 for COR1-A and 0.62 for COR2-A, whereas COR1-A fitted alone reached 0.90 and COR2-A alone reached 0.80 with the aspect ratio at its lower bound. In this case every view, or every subset obtained by leaving out one view when more than two are present, is fitted separately, and the refits are compared in three steps. A fit with no parameter within 0.5 % of its range from a bound (a direction at the edge of the cone does not count) is preferred to one with such a parameter. Between fits of the same kind the higher mean IoU decides. Refits whose mean IoU differs by less than 0.02 are tied, and the tie goes to the refit that keeps more views and then to the view nearer the Sun, that is, with the smaller field of view. The chosen refit is then subjected to the same test against the joint fit as in point 3, and every candidate, with its IoU and the parameters at a bound, is written to `refits` of `view_check`.

Views from widely separated observers, such as C2 and the two COR1 instruments on 2013 October 25, follow the rule above unchanged. If the fit that is finally kept still has the aspect ratio $\kappa$ or the half-angle $\alpha$ at a bound, `event_report.md` names the parameter and the fit stays poor, since the bounds are not widened. Both fits are written to the key `view_check` of `gcs_fit.json` (the instruments, the IoU of each, the mean, and whether the refit was kept) and to the section *Fit quality* of `event_report.md`. The instruments left out are listed in `views_dropped` of the record `gcs_fit` of `config.json`. The IoU and the quality criteria of Section 7.5 then refer to the views that were kept, and the front of the height-time step is chosen among them. The thresholds of the rule are conventions and not statistics. A fit to one view leaves the direction and the extent along the line of sight poorly constrained, and the direction is then tied to the source region (Section 7.4). Moreover, the criterion tests only the agreement of the shell with the envelopes and not whether the dropped envelope was the correct one.

### 7.4 Direction of the shell with one viewpoint

A single coronagraph, or two coronagraphs whose observers are close to each other, constrains the apex height, the aspect ratio, and the half-angle through the outline of the envelope, but not the direction of the axis along the line of sight, because a shell tilted towards or away from the observer projects to nearly the same outline. In the practice of the GCS model, the direction is therefore tied to the source region when only one viewpoint is available (Thernisien 2011). Following this practice, corona3d does not free the Stonyhurst longitude and latitude of the axis when the fit has only one view, or when the largest angle between the observers of the views that are fitted is below 30°. The angle is computed from the positions of the spacecraft at the base frame (`gcs.observer_separation`, `gcs.direction_limit`). The direction is then bounded within a cone around the flare site, $|\Delta\mathrm{lon}| \le c$ and $|\Delta\mathrm{lat}| \le c$, where $c$ is set with `--direction-cone DEG` (default 15°; with 0 the direction is fixed at the flare site, as with `--fix-direction`, and a cone wider than the free bounds of ±40° in longitude and ±35° in latitude changes nothing). The rule is applied to the views that remain after the view rule of Section 7.3, so that a fit that starts with two views and ends with one is also constrained. A fit with three views at wide separation, such as C2 and COR1 of both STEREO spacecraft on 2013 October 25, remains free.

A direction at the edge of the cone is a prior from the source region and not a limit of the data, so that the check of the bounds in Section 7.5 skips the longitude and latitude in this case. The constraint is written to `direction` of `gcs_fit.json` and to `direction_constraint` of the record `gcs_fit` of `config.json` (the views, their separation, the cone, and the reason), and `event_report.md` states "direction constrained by the source region (single viewpoint)" under *Assumed* and in *Fit quality*. The direction of such a shell is assumed and not measured, and the freedom that remains is that of the cone.

### 7.5 Quality criteria

A fit is judged poor, by `gcs.quality`, when any of the following holds.

- The mean IoU over the views that were kept is below `--iou-min` (default 0.5), or no instrument gave an envelope.
- A fitted parameter lies within 0.5 % of its range from a bound of the optimisation. The upper bound of the apex height is the one that was used in that fit (Section 7.2), which is recorded in `height_bound` of `gcs_fit.json`, so that a height at the bound of the fit is reported as a bound hit even though the bound is below 15 R⊙.
- The half-angle is below 5° (`gcs.ALPHA_MIN`). A shell this narrow is a degenerate needle, whose outline in one view can be matched by its legs while its apex remains unconstrained, and a high IoU does not make such a fit physical.
- The apex height at the time of the fit lies outside the outer field of view of every view that was fitted, that is, beyond 4 R⊙ for COR1, 6 R⊙ for C2, and 15 R⊙ for COR2, so that no image shows the front of the shell.

With `--fit auto`, `corona3d event` decides between the free and the fixed-direction fit in the following order (`pipeline.fit_step`, `gcs.choose_fit`).

1. The free fit is made on all views.
2. The view rule of Section 7.3 is applied to it, and the free fit is repeated without the poor views.
3. The quality of the result is judged with the criteria above, on the views that remain.
4. Only when this result is poor, or when the free fit failed, the fit with the direction held at the flare site is made, again with the view rule.
5. The better of the two fits is kept: a good fit beats a poor one; between two fits of the same quality the one with the higher mean IoU over its kept views wins, and a difference below 0.02 goes to the free fit, since the fixed-direction fit is a special case of it. A free fit that is good is therefore never replaced, and the fixed-direction fit is then not run.

The order matters because the quality of a fit that still contains a view the shell cannot match is judged too low. On 2024 May 9 the free fit gave IoU values of 0.90 (COR1-A) and 0 (COR2-A), a mean of 0.45, and the fixed-direction fit with both views a mean of 0.38 only; after the view rule, the free fit on COR1-A alone is good, as judged by its mean IoU. The steps, both fits, the views used, and the criterion are written to the key `decision` of `gcs_fit.json` and to the section *Fit quality* of `event_report.md`, and the list of fitted views in the report is that of the views used in the accepted fit (`views_used`), so that an instrument without an envelope is never listed.

### 7.6 Shells that are drawn

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

The image is drawn on the plane through Sun centre perpendicular to the line of sight of its spacecraft, at the point $X\hat{\mathbf{x}} + Y\hat{\mathbf{y}}$ of the observer's axes (Section 6), and seen from that spacecraft in the orthographic projection it therefore covers the model at its true scale and orientation. Pixels inside the occulter and outside the field of view are transparent: LASCO C2 from 2.2 to 6 R⊙, COR1 from 1.5 to 4 R⊙, and COR2 from 4 to 15 R⊙, so that COR2 begins where COR1 ends and the two coplanar images do not overlap. The grey scale spans plus and minus the 98th percentile of the absolute values, and the opacity is `opacity.wl` = 0.75, which can also be made to vary with the pixel value (USER_GUIDE.md, Section 4). LASCO C2 is drawn only in a view from Earth, and COR1 and COR2 only in a view from the STEREO spacecraft of the event.

## 10. From a time to a reconstruction

For an event given only by its time and flare position, `corona3d event` runs the steps in the order of the earlier sections: it writes the configuration (observers, model run, and image windows, Section 2 and 6), downloads the data, chooses the radio interval (Section 8.1), fits the shell (Section 7.2), and draws the figures. The free fit comes first, and the fit with the direction held at the flare site follows when the free fit fails or is judged poor (Section 7.5). A halo CME, a data gap, or a small angular separation between Earth and the STEREO spacecraft limits the free fit, and [LIMITATIONS.md](LIMITATIONS.md) describes these cases. The file `event_report.md` separates what was measured from what was assumed, lists the missing data, and states the fit quality.

## 11. Validation

corona3d is validated in two ways. The offline tests check each step of the code against synthetic data with known answers, and the event of 2013 October 25 checks the whole reconstruction against the published reconstruction of Morosan et al. (2022), from which the run, the shell, and the radio heights can be compared. The offline tests (`python -m pytest`) use the following synthetic data:

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

With `CORONA3D_ONLINE=1`, the online test also downloads the public data of 2013 October 25 and draws both figures of the README. The comparison of this benchmark with the published figure, together with the differences that remain because some of its inputs are not public, is described in Section 2 above and in Section 11 of [USER_GUIDE.md](USER_GUIDE.md).

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
