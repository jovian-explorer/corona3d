# Methods of corona3d: a three-dimensional coronal scene of a coronal mass ejection from public data

This document describes how corona3d computes every element of its scene, which values it uses, and where each value comes from. The references are collected in `docs/refs.bib`, and every reference with a DOI was checked against its registry record with Argus.

## 1. Overview

The interaction of a coronal mass ejection (CME) with the surrounding corona depends on quantities that vary in all three dimensions: the electron density sets the height at which a given radio frequency is emitted, the magnetic field and the density together set the Alfvén speed that the shock must exceed, and the topology of the field decides whether accelerated electrons escape along open field or remain on closed loops. Morosan et al. (2022) brought these quantities together in one scene for the CME of 2013 October 25, with a magnetohydrodynamic (MHD) model of the corona, a graduated cylindrical shell (GCS) of the CME, and the positions of the herringbone radio sources, seen from Earth, from above the solar north pole, and from STEREO-B. corona3d rebuilds this scene from public data and adds the white-light coronagraph images of LASCO C2 and STEREO-B COR1 and COR2. Every element is a layer that can be switched on or off (Sections 3 to 9), and every number below is a key of `config.json` that the user can change.

## 2. Data

* **Coronal model.** The thermodynamic MHD model MAS/MAST of Predictive Science Inc. (PSI; Lionello et al. 2009; Mikić et al. 2018), run `cr2143-medium/hmi_mast_mas_std_0201` (Carrington rotation 2143, boundary condition from an SDO/HMI magnetogram), read from the PSI archive as HDF4 files of the mass density (`rho002`), the three components of the magnetic field (`br002`, `bt002`, `bp002`), and the Alfvén speed (`va002`). The run used by Morosan et al. (2022) is no longer served by PSI. Among the twelve public runs of Carrington rotations 2141 to 2144, this one reproduces the outline of their density surface best, with an intersection over union of 0.87 to 0.93 in the three views.
* **White-light images.** SOHO/LASCO C2 (Brueckner et al. 1995) total-brightness images with the clear filter, and STEREO-B/SECCHI (Howard et al. 2008) COR1 polarisation triplets (sequence `s4c1B`, at 0°, 120°, and 240°) and COR2 total-brightness images (sequence `d4c2B`), all at level 0.5, downloaded from the Virtual Solar Observatory with SunPy (The SunPy Community et al. 2020) for the window `wl_window` (default 14:35 to 15:40 UT).
* **CME geometry.** The six GCS parameters of the CME. They are not part of the data downloads; the values of the default configuration were fitted from public images as described in Section 7.
* **Radio sources.** Sky-plane positions supplied by the user (Section 8).

## 3. The MHD model and the Alfvén speed

The four fields are interpolated linearly onto the cell-centre mesh of the Alfvén speed file and cut at 3 R⊙. The MAS code units are converted with ρ₀ = 10⁸ cm⁻³ (number density, one proton mass per particle) and B₀ = 2.2068914 G. Following Morosan et al. (2022), who multiplied the model density by 2 after Wang et al. (2017) to obtain more accurate emission heights, the electron density is n_e = 10⁸ ρ × `density_scale` cm⁻³, with `density_scale` = 2. The Alfvén speed is then computed from the scaled density,

v_A = |B| / (4π n_e m_p)^{1/2},

so that it is consistent with the density surface on which it is drawn.

The model is rotated about the solar rotation axis by `model_dlon` = −8.243°, which is the change of the Earth's Carrington longitude between 00:00 and 15:00 UT on 2013 October 25 (5.615° and 357.372°). With this rotation, the model is placed as it would be for observers at the start of the day, while the observers themselves are kept at their positions of 15:00 UT. Setting `model_dlon` to 0 places the model and the observers at the same time.

## 4. The density surface (layer `surface`)

The scaled density is placed on a structured grid in Carrington-Cartesian coordinates (x towards Carrington longitude 0°, z along the rotation axis), and the iso-surface at `ne_level` = 9 × 10⁷ cm⁻³ is extracted with the contour filter of the Visualization Toolkit (VTK), through PyVista (Sullivan & Kaszynski 2019). This density is the plasma level of 173 MHz emitted at the second harmonic, the lower of the two herringbone frequencies of Morosan et al. (2022). The Alfvén speed is interpolated onto the surface and drawn with the `Spectral_r` colour map between 100 and 700 km s⁻¹ (`va_range`). The surface is not smoothed; only its normals are computed, so that the shading is smooth while the shape remains as contoured.

## 5. Field lines (layers `open`, `closed`, and `streamer`)

The field lines are integrated with a fourth-order Runge-Kutta scheme along the unit vector of B, interpolated trilinearly in (φ, θ, r), with a step of 0.01 R⊙ and at most 800 steps in each direction from the seed. A line ends when it leaves the model between the inner boundary and 3 R⊙. It is classed as open when one of its ends lies within 0.05 R⊙ of the outer boundary, and as closed when both ends return to the inner boundary. The tracer was written for corona3d because the tracer of psipy (streamtracer) did not follow these files correctly in our tests: its lines were not tangent to B, with a median |cos(line, B)| of 0.83, against 1.000 for the present tracer. The offline tests check the tracer against the analytic field lines of a dipole (Section 10).

* **Open and closed lines.** The lines are seeded at `n_seeds` = 4500 points distributed uniformly over a sphere of radius `seed_r` = 1.03 R⊙ (a Fibonacci lattice). Open lines are drawn in black and closed lines in grey, both at low opacity, so that the density surface remains visible through them.
* **Streamer loops.** Morosan et al. (2022) found that the herringbones originate where the CME flank meets the closed loops of a streamer. corona3d seeds these loops on a line of `streamer_flank.n` = 40 points at `streamer_flank.r` = 1.8 R⊙, along the meridian at Stonyhurst longitude `streamer_flank.lon` = −40° (the meridian of the highest closed loops of the model near the CME). The line spans ±60° (`streamer_flank.half_span`) in latitude about the latitude of the point of the GCS shell that faces Earth most directly. Of the lines traced from these seeds, the closed lines whose apex exceeds `streamer_min_apex` = 1.9 R⊙ are drawn in purple. In this model, they are the arcades of the helmet streamer that overlie the CME, above the low loops at its foot.

## 6. Observers and viewing geometry

The observers are given by their Carrington longitude and latitude at 15:00 UT on 2013 October 25: Earth at (357.372°, 4.993°) and STEREO-B at (215.329°, −0.984°), 142° apart, computed with SunPy from the JPL Horizons ephemerides. For each observer, corona3d uses the heliocentric-Cartesian axes of the observer's sky: x towards the west limb, y towards solar north (projected onto the plane of the sky), and z towards the observer. The north-pole view looks down the rotation axis with the direction of the Earth at the bottom of the panel.

The panels are perspective renderings, with one camera at `camera_distance` = 6.65 R⊙ from the focal point for all three views and a view angle of 2 arctan(h/D), where h is the half-height of the panel. The half-height and the offset of Sun centre in each panel (`frame`) follow the layout of the published figure. When coronagraph images are shown, or when `--half` is given, the half-height grows, and the camera moves back in proportion so that the perspective does not change.

## 7. The CME shell (layer `gcs`)

The shell is the GCS model of Thernisien (2011), built with PyThea (Kouloumvakos et al. 2022), from the apex height, the aspect ratio κ, the half-angle α between the legs, the Stonyhurst longitude and latitude of the axis, and the tilt. corona3d takes the tilt in its own convention, and PyThea receives it plus 90°. With half the default number of vertices along the axis and around it (`mesh_rings`: 5 straight and 10 circular), the shell is drawn as rings around its axis in magenta.

Morosan et al. (2022) fitted the shell to white-light images at 15:24 UT and scaled it back to 15:00 UT, but they do not list its parameters. The default values were therefore fitted again from public images in the same two steps. At 15:24 UT, all six parameters were fitted by differential evolution to the CME envelopes in LASCO C2 and COR1-A and -B (base differences, with the significant regions connected to the radial line above the flare site at S06E69 filled to their convex hull). The fit gives a direction of Stonyhurst longitude −53.3° and latitude 3.1°, a tilt of −56.7° in the PyThea convention, and an intersection over union of 0.76, 0.90, and 0.79 in the three images. At 15:00 UT, the direction and the tilt were kept, and the height, κ, and α were fitted to the eruption envelopes in SDO/AIA 193 Å and STEREO/EUVI-A and -B 195 Å, giving a height of 1.443 R⊙, κ = 0.334, and α = 1.07°. The agreement at 15:00 UT is weaker (intersection over union 0.63 for AIA, 0.34 for EUVI-B, and 0 for EUVI-A), because the low EUV front of the eruption does not have the shape of the white-light CME 24 minutes later. The fitting scripts are not part of corona3d; users who have their own GCS fit enter its parameters under `gcs` in `config.json`.

## 8. Radio sources (layer `radio`)

The sources are given as sky-plane positions seen from Earth (x west, y north; in R⊙, or in arcsec with `"unit": "arcsec"`, converted with 959.6 arcsec per R⊙), each with its frequency. Each source is placed at the radial distance where the four-fold Newkirk (1961) density model,

n_e(r) = 4 × 4.2 × 10⁴ × 10^{4.32/r} cm⁻³,

has the plasma frequency f_p = 8980 n_e^{1/2} Hz equal to half the observed frequency (emission at the second harmonic). This gives 1.576 R⊙ at 173 MHz and 1.449 R⊙ at 228 MHz, the heights of Morosan et al. (2022). The point is the one at this distance on the line of sight of the observed position, on the side of the plane of the sky that faces the observer, so that x and y are kept exactly. The spheres are coloured by frequency.

## 9. Coronagraph images (layers `c2`, `cor1`, and `cor2`)

For each instrument, the frame nearest to `wl_time` (default 15:24 UT) and the first frame of the downloaded window are read, normalised by their exposure times, and binned by block averages to about 256 × 256 pixels. For COR1, the total brightness of each triplet is formed as tB = (2/3)(I₀ + I₁₂₀ + I₂₄₀) before the binning. The first frame is subtracted from the chosen one (a base difference, which removes the static stray light and the F corona of the level-0.5 images), and the difference is rotated to solar north up with bilinear interpolation. The position of each pixel in the plane of the sky follows from the world coordinate system of the image and the angular radius of the Sun seen by the spacecraft,

(X, Y) = (pixel − Sun centre) × (plate scale / angular solar radius) R⊙.

The image is drawn on the plane through Sun centre perpendicular to the line of sight of its spacecraft, at the point X e_x + Y e_y of the observer's axes (Section 6), so that from that spacecraft it covers the model at its true scale and orientation. Pixels inside the occulter and outside the field of view are transparent: LASCO C2 from 2.2 to 6 R⊙, COR1 from 1.5 to 4 R⊙, and COR2 from 4 to 15 R⊙, so that COR2 begins where COR1 ends and the two coplanar images do not overlap. The grey scale spans ± the 98th percentile of the absolute difference, and the opacity is `wl_opacity` = 0.75. LASCO C2 is drawn in the Earth view and COR1 and COR2 in the STEREO-B view; the north-pole view shows no image.

## 10. Validation

The offline tests (`python -m pytest`) check the following:

* The tracer follows a dipole field, for which a line obeys r = L sin²θ. From 1.05 R⊙ at colatitude 45°, the traced line is closed, its apex lies at 2.1 R⊙ within 0.02 R⊙, it returns to the inner boundary at latitude −46.4°, and L is constant along it to 1.5 %. From colatitude 8°, the line is open.
* The radio heights are 1.576 R⊙ at 173 MHz and 1.449 R⊙ at 228 MHz.
* The observer axes are orthonormal and right-handed, with x towards the west limb, and the deprojected sources lie at their radius with their sky-plane coordinates unchanged.
* A feature placed at (5, 2) R⊙ west and north in a synthetic COR2 image, at roll angles of 0°, 6.6°, and 175.3°, is found in three dimensions at its true position within 0.15 R⊙, and the occulter is transparent.
* COR1 frames are grouped into complete triplets only.
* The command line refuses to work outside a working directory, `corona3d init` leaves an edited configuration unchanged, an unknown layer is refused, and the layers that need no data render.

With `CORONA3D_ONLINE=1`, the online test also downloads the public data of 2013 October 25 and draws both figures of the README.

## 11. Limitations

* The public MAS/MAST run is not the run of Morosan et al. (2022), so that the density surface agrees with theirs in outline but not in every detail. Their Alfvén-speed colours are matched more closely with a density scale of 2.2 to 2.5 and the model rotated a further 7.5° east, which indicates an older model version or a different magnetogram date. corona3d keeps the scale of 2 stated in their paper.
* The streamer loops of the public run do not form the regular fan of loops shown by Morosan et al. (2022); corona3d draws the loops of the public model by the same method, not a copy of their shape.
* The GCS fit at 15:00 UT is weak in the EUV images (Section 7).
* The radio sources are not derived from the Nançay Radioheliograph images by corona3d, since these are delivered by the BASS2000 archive on request as visibilities that must be imaged with SolarSoft.
* The coronagraph images are level 0.5 and are not calibrated photometrically; base differences show where the brightness changed and are used for comparison of positions only.
* The images are drawn in the plane of the sky of their spacecraft, which is exact only when seen from that spacecraft; from any other viewpoint, the plane shows only where the image lies, without the three-dimensional structure of the CME.

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
