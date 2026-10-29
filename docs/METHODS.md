# Methods of corona3d

This document describes how corona3d computes every element of its reconstruction, which values it uses, and where each value comes from. Every number below that is written in `code font` is a key of `config.json`, which the user can change. Each method and value is also listed with its source and its test in [MODEL_AUDIT.md](MODEL_AUDIT.md), the assumptions that limit the results are collected in [LIMITATIONS.md](LIMITATIONS.md), and the references are in [refs.bib](refs.bib); every reference with a digital object identifier (DOI) was checked against its registry record with Argus.

## 1. Overview

The interaction of a coronal mass ejection (CME) with the surrounding corona depends on quantities that vary in all three dimensions. The electron density sets the height at which a given radio frequency is emitted, the magnetic field and the density together set the Alfvén speed that a shock must exceed, and the topology of the field decides whether accelerated electrons escape along open field or remain on closed loops. Morosan et al. (2022) brought these quantities together in a three-dimensional reconstruction of the CME of 2013 October 25, which combines a magnetohydrodynamic (MHD) model of the corona, a graduated cylindrical shell (GCS) of the CME, and the positions of herringbone radio sources, seen from Earth, from above the solar north pole, and from the Solar Terrestrial Relations Observatory B (STEREO-B). Such reconstructions, however, have so far been made by hand for single events, with steps and inputs that are rarely available in a form that can be repeated for another event, so that the relation between the shock, the coronal field, and the radio sources remains poorly quantified for most events. corona3d is a package that we developed to make every step of such a reconstruction an objective and reproducible computation from public data for any event, namely the MHD model of Predictive Science Inc. (PSI) rendered in three dimensions, a GCS shell fitted to coronagraph images, the white-light images of the Large Angle and Spectrometric Coronagraph (LASCO) C2 and of the STEREO coronagraphs COR1 and COR2, and radio sources measured in radio images and placed at the heights of a coronal density model on the side of the plane of the sky that faces the observer. Every element is a layer that can be switched on or off (Sections 3 to 9), and the shell is fitted as described in Section 7. Moreover, Section 10 describes how these steps are run for an event given only by its time and flare position, and Section 11 describes the checks of the code and the validation on the event of 2013 October 25, for which the result of corona3d is compared with the reconstruction of Morosan et al. (2022). Beyond the figure, Sections 12 to 15 describe the physical analysis of the reconstruction: the Mach numbers and the shock geometry of the CME shell in the coronal model, the speeds of the shell and of the radio sources, the range of the radio heights over several density models, the circular polarisation of the radio sources, and the apex height measured in EUV images at the time of the burst.

All positions are expressed in solar radii (R⊙) in one of two frames (Section 6). The model, the field lines, the shell, and the radio sources are placed in Heliographic Carrington Cartesian coordinates, whereas an image of an observer is described by the sky-plane axes of that observer.

## 2. Data

The coronal model is the thermodynamic MHD model MAS/MAST of PSI (Mikić et al. 1999; Lionello et al. 2009; Mikić et al. 2018). The default run, `cr2143-medium/hmi_mast_mas_std_0201`, is that of the validation event and covers Carrington rotation 2143 with a boundary condition from a magnetogram of the Helioseismic and Magnetic Imager on the Solar Dynamics Observatory (SDO/HMI). It is read from the PSI archive as HDF4 files of the mass density (`rho002`) and of the three components of the magnetic field (`br002`, `bt002`, `bp002`), whereas of the Alfvén-speed file (`va002`) only the axes of its mesh are used. The run spans 1.0 to 30.5 R⊙. Since the run used by Morosan et al. (2022) is no longer served by PSI, we chose, among the twelve public runs of Carrington rotations 2141 to 2144, the one whose density surface agrees best with the outline of theirs (intersection over union 0.87 to 0.93 in the three views). This choice was therefore made by comparison with the published figure, and it applies to the validation event alone. For a new event, `corona3d new` selects the run of the Carrington rotation of the event, at medium resolution when it is listed in the PSI archive and otherwise at high resolution.

The white-light images are those of the C2 coronagraph of LASCO on the Solar and Heliospheric Observatory (SOHO; Brueckner et al. 1995), with 1024 × 1024 pixels, together with those of the COR1 and COR2 coronagraphs of the Sun Earth Connection Coronal and Heliospheric Investigation (SECCHI; Howard et al. 2008) on STEREO-A or STEREO-B. COR1 provides polarisation triplets (at 0°, 120°, and 240°) and COR2 total-brightness images, all at level 0.5. These images are downloaded from the Virtual Solar Observatory with SunPy (The SunPy Community et al. 2020) for the window `wl_window`, and `--cor1a` adds the COR1 images of STEREO-A for events observed by STEREO-B, which only the GCS fit uses.

The CME geometry is given by the six GCS parameters, fitted to the white-light images with `corona3d gcs` (Section 7). The radio sources are measured in the 2-D images of the Nançay Radioheliograph (NRH; Kerdraon & Delouis 1997) or of the Owens Valley Radio Observatory Long Wavelength Array (OVRO-LWA; Chhabra et al. 2021), or they are sky-plane positions supplied by the user (Section 8).

## 3. The MHD model and the Alfvén speed

The four fields are interpolated linearly onto the cell-centre mesh of the Alfvén-speed file over the whole run. The MAS code units are converted with a number density $\rho_0 = 10^{8}\ \mathrm{cm^{-3}}$ (one proton mass per particle) and a magnetic field $B_0 = 2.2068914\ \mathrm{G}$. The model density is multiplied by 2, the factor that Morosan et al. (2022) applied after Wang et al. (2017) to obtain more accurate emission heights, so that the electron density is

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

A line is classed as open when one of its ends leaves the outer boundary, as closed when both ends return to the inner boundary, and as undetermined otherwise, for example when it runs out of steps or stalls where $\mathbf{B}$ vanishes. Undetermined lines are counted and not drawn. Lines are traced through the whole model for their classification, so that a loop that closes above 3 R⊙ near the cusp of a streamer is classed as closed, but they are drawn only within `line_r_max` = 3 R⊙. corona3d uses its own tracer, since in our tests on these files the lines of the psipy tracer (streamtracer) were not tangent to $\mathbf{B}$, with a median $|\cos(\text{line},\mathbf{B})|$ of 0.83 against 1.000 for the tracer of corona3d. The offline tests check the tracer against the analytic field lines of a dipole (Section 11).

### 5.2 Seeds and streamer loops

The open and closed lines are seeded at `n_seeds` = 4500 points distributed uniformly over a sphere of radius `seed_r` = 1.03 R⊙ (a Fibonacci lattice). Open lines are drawn in black and closed lines in grey, both at low opacity, so that the density surface remains visible through them.

Morosan et al. (2022) concluded that the herringbones originate where the CME flank meets the closed loops of a streamer. The earlier reconstruction does not describe how its streamer loops were selected, whereas corona3d selects them without any free parameter, by tracing lines through every vertex of the GCS shell above 1.02 R⊙ and drawing the closed ones in purple. These loops are the closed field lines that the CME shell intersects at the time of the reconstruction. For the default event, 478 of the 943 lines through the shell are closed and 465 are open, and the program prints these numbers at every rendering. At most `streamer_max_lines` = 150 of them are drawn, as an even subset.

## 6. Coordinate frames, viewpoints, and projection

### 6.1 Frames

Model, field lines, shell, and radio sources are expressed in Heliographic Carrington Cartesian coordinates in R⊙, centred on the Sun, with the $x$ axis towards Carrington longitude 0° and the $z$ axis along the rotation axis towards solar north. A point at Carrington longitude $\ell$, latitude $b$, and distance $r$ is

$$ \mathbf{x} = r\,(\cos b\cos\ell,\ \cos b\sin\ell,\ \sin b). $$

This frame rotates with the Sun, so that the positions refer to the time of the reconstruction. The Stonyhurst longitude of a point is its Carrington longitude minus that of Earth at that time.

The observers are given by their Carrington longitude and latitude (`earth` and `stereo`). For the event of 2013 October 25 at 15:00 UT, Earth is at (357.372°, 4.993°) and STEREO-B at (215.329°, −0.984°), 142° apart. These values were computed with SunPy from the Horizons ephemerides of the Jet Propulsion Laboratory (JPL), and `corona3d new` computes them in the same way for a new event. For each observer with the unit vector $\hat{\mathbf{z}}$ towards the observer, corona3d uses the sky-plane axes

$$ \hat{\mathbf{x}} \ \text{towards the west limb},\qquad \hat{\mathbf{y}} \ \text{towards solar north (projected on the sky plane)},\qquad \hat{\mathbf{z}} = \hat{\mathbf{x}}\times\hat{\mathbf{y}}, $$

an orthonormal, right-handed set.

### 6.2 Viewpoints and projection

Each of the three panels is given by an entry of `views`: an observer (`"from": "earth"`, `"stereo"`, or `"north"`, the last looking down the rotation axis with the direction of the Earth at the bottom of the panel), or any viewing direction given by a Carrington longitude and latitude, with solar north up and an optional roll angle. The panels are drawn with an orthographic (parallel) projection, so that every depth along the line of sight is drawn at the scale of the plane of the sky, which keeps the model, the shell, the radio sources, and the coronagraph images at one common scale in R⊙. In contrast, a perspective camera at 6.65 R⊙ would magnify structures in front of the plane of the sky by up to 18 %. Such a camera is nevertheless available with `"projection": "perspective"` at `camera_distance` = 6.65 R⊙, so that a figure can be drawn with a perspective comparable to that of the published Fig. 7. In the configuration of the validation event, the half-height and the offset of Sun centre of each panel (`frame` in each view) match the panels of the published figure, so that the two reconstructions can be compared directly, and for any other event they are set freely. When coronagraph images are shown, the panel of their spacecraft is widened to 6 R⊙ (C2 or COR2) or 4 R⊙ (COR1 alone), while `--half` sets one half-height for all panels.

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

A radio burst usually precedes the time at which the CME is well developed in the coronagraph images, so that the shell cannot be fitted at the time of the burst itself. For the event of 2013 October 25, Morosan et al. (2022) fitted their shell at 15:24 UT and reconstructed it at earlier times from extreme-ultraviolet images, without listing its parameters or a criterion of the fit. The command `corona3d gcs` instead fits the six parameters automatically in two steps, with an objective measure of the agreement with the images and an uncertainty of the height at the time of the burst.

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

The main input of the layer is the 2-D images of the NRH (Kerdraon & Delouis 1997), the Stokes I brightness temperature of 64 × 64 pixels at one frequency about every 10 s, rotated to solar north up, with Sun centre at the reference pixel and one solar radius spanning a number of pixels given in the header. For each frequency, one image is taken every `step_s` seconds between `start` and `end`, and an image is kept when its peak brightness temperature reaches `min_tb` ($10^{8}$ K). In each kept image, the centroid is the brightness-weighted mean position of the pixels brighter than `level` (0.5) times the peak. The uncertainty of a centroid is of the order of one pixel, 1/16 R⊙, in the plane of the sky. The spectral images of OVRO-LWA (Chhabra et al. 2021), available from 2024 April onward at 32 to 87 MHz, are treated in the same way. Each file holds one 10 s integration with 12 band images, stored at reduced size and restored here to the full grid of 256 × 256 pixels by cubic interpolation. The helioprojective header (reference pixel, plate scale, reference coordinates, and angular solar radius) gives the position of every pixel in R⊙, and for each requested frequency the nearest band within 3 MHz is used. The centroids are drawn as a chain of blobs coloured from the first to the last image, so that the motion of a source appears in the reconstruction, in the same way as the herringbone sources of Morosan et al. (2022) were coloured by time.

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

so that $X$ and $Y$ are kept exactly. A position that lies farther from Sun centre than its height, $X^{2}+Y^{2}>r^{2}$, cannot be placed at that height, and it is therefore kept in the plane of the sky with a printed warning. The heights follow the Newkirk model and not the density of the MHD model, so that a source and the density surface locate the same plasma frequency at different heights wherever the two density models differ. The scattering of radio waves in the corona, which displaces observed sources outward, is not corrected. The sources are drawn as spheres coloured by frequency. The heights of the same positions in other density models, including the MHD density itself, and their range are given by `corona3d shock` (Section 13).

## 9. Coronagraph images (layers `c2`, `cor1`, and `cor2`)

For each instrument, the frame nearest to `wl_time` within `wl_window` and the first frame of that window are read and normalised by their exposure times. For COR1, the total brightness of each triplet is formed as $tB = \tfrac{2}{3}(I_0 + I_{120} + I_{240})$, where the three exposures of a triplet are those taken within 60 s. Each frame is binned by block averages to about 256 × 256 pixels, and a frame of another size is resampled onto the grid of the base frame. The comparison is a ratio by default (`wl_mode` = `"ratio"`),

$$ \frac{I}{I_{\rm base}} - 1, $$

which removes multiplicative instrumental effects such as vignetting and the stray light that scales with the coronal brightness, whereas `"difference"` gives $I - I_{\rm base}$. A `wl_time` outside the window, or a frame that would be the base frame itself, stops the rendering with a message, since either would draw a wrong or blank image. When the nearest frame to `wl_time` lies more than 30 min away, which happens at a data gap, the rendering stops as well, and the message names the nearest frame times before and after `wl_time`; `corona3d gcs --t-ref` and `corona3d event --t-ref` then set `wl_time` to the chosen time. The result is rotated to solar north up with bilinear interpolation. The position of each pixel in the plane of the sky follows from the world coordinate system of the image and the angular radius of the Sun seen by the spacecraft,

$$ (X, Y) = (\text{pixel} - \text{Sun centre}) \times \frac{\text{plate scale}}{\text{angular solar radius}}\ \mathrm{R_\odot}. $$

The image is drawn on the plane through Sun centre perpendicular to the line of sight of its spacecraft, at the point $X\hat{\mathbf{x}} + Y\hat{\mathbf{y}}$ of the observer's axes (Section 6), and seen from that spacecraft in the orthographic projection it therefore covers the model at its true scale and orientation. Pixels inside the occulter and outside the field of view are transparent: LASCO C2 from 2.2 to 6 R⊙, COR1 from 1.5 to 4 R⊙, and COR2 from 4 to 15 R⊙, so that COR2 begins where COR1 ends and the two coplanar images do not overlap. The grey scale spans plus and minus the 98th percentile of the absolute values, and the opacity is `opacity.wl` = 0.75, which can also be made to vary with the pixel value (USER_GUIDE.md, Section 4). LASCO C2 is drawn only in a view from Earth, and COR1 and COR2 only in a view from the STEREO spacecraft of the event.

## 10. From a time to a reconstruction

For an event given only by its time and flare position, `corona3d event` runs the steps in the order of the earlier sections: it writes the configuration (observers, model run, and image windows, Section 2 and 6), downloads the data, chooses the radio interval (Section 8.1), fits the shell (Section 7.2), and draws the figures. The free fit comes first, and the fit with the direction held at the flare site follows when the free fit fails or is judged poor (Section 7.5). A halo CME, a data gap, or a small angular separation between Earth and the STEREO spacecraft limits the free fit, and [LIMITATIONS.md](LIMITATIONS.md) describes these cases. The file `event_report.md` separates what was measured from what was assumed, lists the missing data, and states the fit quality.

## 11. Validation

corona3d is validated in two ways. The offline tests check each step of the code against synthetic data with known answers, and the event of 2013 October 25 compares the whole reconstruction with the independent reconstruction of Morosan et al. (2022), with which the density surface, the shell, and the radio heights can be compared. The offline tests (`python -m pytest`) use the following synthetic data:

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

With `CORONA3D_ONLINE=1`, the online test also downloads the public data of 2013 October 25 and draws both figures of the README. The comparison of this validation case with the published reconstruction, together with the differences that remain because some of its inputs are not public, is described in Section 2 above and in Section 11 of [USER_GUIDE.md](USER_GUIDE.md).

## 12. Shock diagnostics of the CME shell (`corona3d shock`)

Whether a CME drives a shock, and how strong that shock becomes, depends on the speed of the CME front relative to the characteristic speeds of the corona ahead of it, and these speeds vary strongly over the surface of the CME, since the front meets open and closed field, dense streamers and tenuous coronal holes at the same time. The command `corona3d shock` therefore evaluates the Mach numbers and the shock geometry at every point of the GCS shell and at the radio sources. The steady MAS/MAST run contains no CME, so that the plasma it gives at a point of the shell is the undisturbed (upstream) corona into which the shell advances. The shell is used as the surface of the shock, which assumes that the shock lies close to the CME front at these low heights; a stand-off distance is not added.

### 12.1 Upstream plasma

Besides the density and the field (Section 3), the command reads the temperature (`t002`) and the three components of the flow (`vr002`, `vt002`, `vp002`) of the same run, which it downloads when they are missing (`corona3d fetch --model-extra`). The MAS code units of these quantities follow from those of the field and the density: the speed unit is $v_0 = B_0/\sqrt{4\pi n_0 m_p} = 481.37\ \mathrm{km\,s^{-1}}$, and the temperature unit is $T_0 = m_p v_0^2/k_B = B_0^2/(4\pi n_0 k_B) = 2.807\times10^{7}$ K, with $n_0 = 10^{8}\ \mathrm{cm^{-3}}$ and $B_0 = 2.2068914$ G. On the run of the validation event the pressure, the density, and the temperature of the files satisfy $p = 2\rho T$ in code units to within $10^{-7}$, which shows that the run holds electrons and protons at one temperature, and the unit gives about 1.2 MK at 1.5 R⊙. The Alfvén speed is computed as in Section 3, with the density chosen by `va_density`, and the sound speed is

$$ c_s = \sqrt{\frac{\gamma\, p}{\rho}} = \sqrt{\frac{2\gamma\, k_B T}{m_p}}, \qquad \gamma = 5/3, $$

since $p = 2 n k_B T$ and $\rho = n m_p$. The sound speed does not depend on `density_scale`, whereas the Alfvén speed, and with it the Alfvén Mach number, changes by $\sqrt{2}$ between the two choices of `va_density`. The magnetic field is converted from its spherical components to Carrington Cartesian components at each point, and the model rotation `model_dlon` is applied to positions and vectors alike.

### 12.2 Normal, shock angle, and fast-magnetosonic speed

The outward normal $\hat{\mathbf{n}}$ of the shell is computed at every vertex from the tangents along its axis and around its cross-section, and it points away from the centre of the cross-section. The angle between the upstream field and the normal,

$$ \theta_{Bn} = \arccos\frac{|\mathbf{B}\cdot\hat{\mathbf{n}}|}{|\mathbf{B}|}, \qquad 0^\circ \le \theta_{Bn} \le 90^\circ, $$

separates quasi-perpendicular ($\theta_{Bn} > 45^\circ$) from quasi-parallel regions of the shock (Treumann 2009), and the threshold is the key `shock.quasi_perp_deg`. The speed of a fast-magnetosonic wave that propagates along the normal is (Priest 2014)

$$ v_f^2 = \tfrac12\left[v_A^2 + c_s^2 + \sqrt{(v_A^2 + c_s^2)^2 - 4 v_A^2 c_s^2 \cos^2\theta_{Bn}}\right], $$

which equals $\sqrt{v_A^2 + c_s^2}$ at a perpendicular shock and the larger of $v_A$ and $c_s$ at a parallel one.

### 12.3 Speed of the shell

With a fixed aspect ratio and half-angle, every length of the GCS shell is proportional to its apex height $h$, and its legs start at Sun centre, so that the shell expands self-similarly about Sun centre (Thernisien 2011). A point $\mathbf{x}$ of the shell therefore moves with the velocity $(\dot h/h)\,\mathbf{x}$, and its speed along the normal, relative to the upstream flow $\mathbf{u}$ of the run, is

$$ V_n = \frac{\dot h}{h}\,(\mathbf{x}\cdot\hat{\mathbf{n}}) - \mathbf{u}\cdot\hat{\mathbf{n}}. $$

The apex speed $\dot h$ is that of the height-time fit of Section 7.2 (`speed_km_s` of `gcs_fit`), or a speed given with `--speed` and `--speed-note`, which is then recorded as assumed. The normal speed equals $\dot h$ at the apex, and it decreases towards the flanks in proportion to $\mathbf{x}\cdot\hat{\mathbf{n}}$, so that the flank speeds are set by the geometry of the fitted shell and by the assumption of self-similar expansion; a CME that expands laterally faster than this, as is often observed early in an eruption, has flank speeds above those of the model. The half-width of the shell, its largest distance from the axis, grows at $\dot h\,w/h$, and the summary gives this lateral speed. The side of the shell that faces the Sun has $V_n \le 0$, does not advance into the corona, and is left without Mach numbers.

### 12.4 Mach numbers and summaries

The Alfvén, fast-magnetosonic, and sonic Mach numbers are

$$ M_A = \frac{V_n}{v_A}, \qquad M_f = \frac{V_n}{v_f}, \qquad M_s = \frac{V_n}{c_s}, $$

at every advancing point. The points are grouped by the angle between their normal and the axis of the CME, the nose within `shock.nose_deg` = 30° and the flanks beyond `shock.flank_deg` = 60°, and by the field line through each point, traced with the tracer of Section 5.1 and classed as open, closed, or undetermined. The summary gives the apex values and, for the nose, the flanks, and the points on closed and on open field, the medians and extremes of $V_n$, $v_A$, $v_f$, $\theta_{Bn}$, $M_A$, and $M_f$, together with the fraction of quasi-perpendicular points. The analysis uses a finer shell than the figure (`shock.mesh`), and points below `shock.r_min` = 1.02 R⊙ are left out.

At each radio source (Section 8) the upstream plasma is sampled at the deprojected position of the source, and the shell, scaled self-similarly to the time of the source with the apex speed, gives the point nearest to the source, its normal, and its normal speed, which does not change with the scaling. The table records the distance of the source from the shell and whether the source lies ahead of or behind the front ($\hat{\mathbf{n}}\cdot(\mathbf{q}-\mathbf{x}) < 0$ for a source $\mathbf{q}$ behind it), so that the Mach number of a source inside the shell is not mistaken for that of the shock ahead of it.

The command writes `shock_shell.csv` (one row per point of the shell), `shock_radio.csv` (one row per radio position), `radio_tracks.csv`, `shock.npz`, `shock_summary.json`, and `shock_summary.md` to the folder `analysis/`, and the section *Shock diagnostics* of `event_report.md` when the report exists. `corona3d render --shell-color mach_f` (or `mach_a`, `theta_bn`, `vn`, `va`) colours the drawn shell by the same quantities, with `--shell-clim` for the colour limits, and the points that do not advance are drawn in grey.

### 12.5 Speeds of the radio sources

For every instrument and frequency, a straight line is fitted to the three-dimensional positions of the centroids against time, and the length of its velocity vector, with the standard error propagated from the fit, is the speed of the source in three dimensions; the same fit to the sky-plane positions gives the projected speed. Since every centroid of one frequency is placed at the same height of the density model, the radial component of this motion is zero by construction, and the speed describes the motion of the emission site over the plasma level, not the radial speed of an electron beam or of the shock.

### 12.6 Validation on the event of 2013 October 25

The offline tests check the speeds against their analytic limits, the normals against those of a tube, the apex speed of a PyThea shell, the self-similar scaling of its half-width, the sampling of a synthetic run with a dipole field, a uniform temperature, and a radial flow, and the Mach numbers against $V_n/v_A$ point by point (`tests/test_shock.py`). For the validation event, the shell of the stored fit (apex height 1.931 R⊙ at 15:00 UT, 1039 km s⁻¹) gives an Alfvén speed of 212 km s⁻¹ and $M_A$ = 4.8 at the apex, where the field is open and $\theta_{Bn}$ = 4°, whereas the flanks advance at a median of 228 km s⁻¹ and reach a median $M_A$ of only 1.3 (quasi-perpendicular fraction 0.81). With the apex height measured in the EUVI images (Section 15), 1.50 ± 0.04 R⊙ at 15:00 UT, the same shape gives $M_A$ = 6.6 at the apex and a median $M_A$ of 3.2 and $M_f$ of 2.2 on the closed field that the shell meets. The values at the radio sources are listed in Section 11 of USER_GUIDE.md, together with those of the earlier reconstruction.

## 13. Height range of the radio sources

The height of a radio source follows from a density model, and the corona can be denser or more tenuous than any one model, so that the height carries an uncertainty that the radio image alone cannot remove. `corona3d shock` therefore computes the height of every radio position in each model of `shock.height_models` and gives their range. The models are the $K$-fold Newkirk (1961) model (`newkirk1`, `newkirk2`, `newkirk4`, and `newkirk6`, Section 8.2), the equatorial background corona of Saito et al. (1977),

$$ n_e(r) = 1.36\times10^{6}\, r^{-2.14} + 1.68\times10^{8}\, r^{-6.13}\ \mathrm{cm^{-3}}, $$

solved numerically for the plasma level, and the scaled density of the MAS/MAST run itself (`mas`). For the last, the density is sampled along the line of sight of the observed position in steps of 0.005 R⊙, from 10 R⊙ towards the plane of the sky (or towards the solar disc for a position on the disc), and the height is the first point at which the density reaches the plasma level, interpolated in $\log n_e$; a line of sight that does not reach it on the side facing the observer gives no MAS height. All models use emission at the second harmonic, as the radio layer does. The columns `r_<model>`, `r_min`, and `r_max` of `shock_radio.csv` hold the heights and their range, and `x_min` to `z_max` the end points of the range on the line of sight, which `corona3d render --radio --radio-height-range` draws as dark segments through the sources. For the 150.9 MHz sources of 2013 October 25 the heights are 1.15 R⊙ (Saito), 1.34 to 1.77 R⊙ (Newkirk ×1 to ×6, 1.647 R⊙ for ×4), and 1.13 to 1.30 R⊙ in the MAS density, so that the run places these sources lower than the four-fold Newkirk model by 0.35 to 0.52 R⊙.

## 14. Circular polarisation of the radio sources

The degree of circular polarisation of a burst indicates whether it is emitted at the fundamental or at the harmonic of the plasma frequency, since fundamental emission is usually strongly polarised and harmonic emission weakly. The NRH files delivered by the Radio Solar Database hold the Stokes V images of each integration in the column `STOKESV` of the same file as Stokes I, and a separate file `nrh2_..._v.fts` beside the `_i` file is also read. For every NRH position of the radio table, the image of the same frame gives

$$ \frac{V}{I}\bigg|_{\rm source} = \frac{\sum V}{\sum I}, \qquad \frac{V}{I}\bigg|_{\rm peak}, $$

where the sums run over the pixels whose Stokes I reaches `level` times the peak, the same pixels that define the centroid (Section 8.1), and the second value is taken at the brightest pixel. The sign is that of the NRH files. A degree of polarisation cannot exceed one in magnitude, and a value $|V/I| > 1$ shows that the Stokes V image holds artefacts larger than the source; such a row is kept in `shock_radio.csv` with `dcp_valid` false and left out of the summary. The files of 2013 October 25 at hand hold Stokes V, and at 150.9 MHz the degree is between −0.02 and +0.04 at the five positions, which indicates harmonic emission, whereas at 432 MHz the Stokes V image exceeds the Stokes I peak by a factor of up to 4.6 and two of the three values are unphysical.

## 15. EUV images and the apex height at early times

A radio burst often begins before the CME is visible in the coronagraphs, so that the apex height at the time of the burst is extrapolated from later frames (Section 7.2), whereas the extreme-ultraviolet (EUV) images of the low corona show the CME and the EUV front at that time. `corona3d fetch --euv` downloads the images of the Atmospheric Imaging Assembly (AIA; Lemen et al. 2012) of the Solar Dynamics Observatory at 193 Å and 211 Å and of the Extreme Ultraviolet Imager (EUVI; Wuelser et al. 2004) of the STEREO spacecraft of the event at 195 Å from the Virtual Solar Observatory, within `euv_window` (by default from 20 min before to 15 min after the event time) and one image every `euv_sample_min` = 5 min.

With `corona3d gcs --euv` (or `corona3d event --euv`), the apex height at t_out is measured in these images instead of being extrapolated, while the shape of the shell (aspect ratio, half-angle, direction, and tilt) is kept from the white-light fit. For every frame of the first EUV view within 10 min of t_out (AIA 193 Å when present, with the frames of the other views within 90 s), the envelope is formed as for the coronagraphs, in the annulus above the limb from 1.02 R⊙ to 1.3 R⊙ (AIA) or 1.7 R⊙ (EUVI), so that the on-disc part of the EUV wave is left out, and the apex height is chosen on a grid of 0.01 R⊙ from 1.05 to 3 R⊙ so that the mean IoU of the projected shell with the envelopes is largest. A frame whose best height lies beyond the outer edge of the EUV field, or at an end of the grid, is not used, because its envelope then constrains only the flanks of the shell. With three or more usable frames, a straight line through their heights gives the height at t_out with the standard error of its intercept and the early speed of the apex. With one or two, the frame nearest to t_out is used and moved to t_out with the speed of the white-light front, and the half-width of the heights whose mean IoU lies within 0.02 of the best is its uncertainty. The EUV height replaces the extrapolated one in `gcs` of `config.json`, the extrapolated height is kept as `height_t_out_wl`, and the record `euv` of `gcs_fit.json` lists every frame with its height, IoU, and the reason why it was not used. Without `--euv`, no EUV image is read and the fit is unchanged. For 2013 October 25, the EUVI-B frame at 15:01:21 UT gives 1.63 R⊙ (IoU 0.61, 1.60 to 1.68 R⊙), and moved to 15:00 UT at 1108 km s⁻¹ it gives 1.50 ± 0.04 R⊙, below the 1.83 ± 0.61 R⊙ extrapolated from the COR1-B front of the same download, as expected for a CME that accelerates early in its rise; the frame at 15:06:22 UT, whose best height of 2.98 R⊙ lies beyond the EUVI field, was not used. The AIA images could not be obtained from the Virtual Solar Observatory when this was written, so that the AIA path is tested only offline.

## References

Brueckner, G. E., Howard, R. A., Koomen, M. J., et al. 1995, Solar Physics, 162, 357. https://doi.org/10.1007/BF00733434

Chhabra, S., Gary, D. E., Hallinan, G., et al. 2021, The Astrophysical Journal, 906, 132. https://doi.org/10.3847/1538-4357/abc94b

Howard, R. A., Moses, J. D., Vourlidas, A., et al. 2008, Space Science Reviews, 136, 67. https://doi.org/10.1007/s11214-008-9341-4

Kerdraon, A., & Delouis, J.-M. 1997, in Coronal Physics from Radio and Space Observations, Lecture Notes in Physics, 483, 192. https://doi.org/10.1007/BFb0106458

Kouloumvakos, A., Rodriguez-Garcia, L., Gieseler, J., et al. 2022, Frontiers in Astronomy and Space Sciences, 9, 974137. https://doi.org/10.3389/fspas.2022.974137

Lemen, J. R., Title, A. M., Akin, D. J., et al. 2012, Solar Physics, 275, 17. https://doi.org/10.1007/s11207-011-9776-8

Lionello, R., Linker, J. A., & Mikić, Z. 2009, The Astrophysical Journal, 690, 902. https://doi.org/10.1088/0004-637X/690/1/902

Mikić, Z., Linker, J. A., Schnack, D. D., Lionello, R., & Tarditi, A. 1999, Physics of Plasmas, 6, 2217. https://doi.org/10.1063/1.873474

Mikić, Z., Downs, C., Linker, J. A., et al. 2018, Nature Astronomy, 2, 913. https://doi.org/10.1038/s41550-018-0562-5

Morosan, D. E., Pomoell, J., Kumari, A., Vainio, R., & Kilpua, E. K. J. 2022, Astronomy & Astrophysics, 668, A15. https://doi.org/10.1051/0004-6361/202244432

Newkirk, G., Jr. 1961, The Astrophysical Journal, 133, 983. https://doi.org/10.1086/147104

Priest, E. 2014, Magnetohydrodynamics of the Sun (Cambridge: Cambridge University Press). https://doi.org/10.1017/CBO9781139020732

Saito, K., Poland, A. I., & Munro, R. H. 1977, Solar Physics, 55, 121. https://doi.org/10.1007/BF00150879

Sullivan, C. B., & Kaszynski, A. A. 2019, Journal of Open Source Software, 4, 1450. https://doi.org/10.21105/joss.01450

The SunPy Community, Barnes, W. T., Bobra, M. G., et al. 2020, The Astrophysical Journal, 890, 68. https://doi.org/10.3847/1538-4357/ab4f7a

Thernisien, A. F. R., Howard, R. A., & Vourlidas, A. 2006, The Astrophysical Journal, 652, 763. https://doi.org/10.1086/508254

Thernisien, A. 2011, The Astrophysical Journal Supplement Series, 194, 33. https://doi.org/10.1088/0067-0049/194/2/33

Treumann, R. A. 2009, The Astronomy and Astrophysics Review, 17, 409. https://doi.org/10.1007/s00159-009-0024-2

Wang, T., Reginald, N. L., Davila, J. M., St. Cyr, O. C., & Thompson, W. T. 2017, Solar Physics, 292, 97. https://doi.org/10.1007/s11207-017-1130-3

Wuelser, J.-P., Lemen, J. R., Tarbell, T. D., et al. 2004, in Telescopes and Instrumentation for Solar Astrophysics, Proc. SPIE, 5171, 111. https://doi.org/10.1117/12.506877
