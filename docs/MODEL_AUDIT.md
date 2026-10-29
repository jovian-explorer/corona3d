# Model audit: every method, model, and number of corona3d

corona3d reconstructs the corona, a coronal mass ejection (CME), and its radio sources in three dimensions for any event, and its results can be trusted only as far as each method and number behind them can be traced. This file therefore lists every place where a method or a number of corona3d determines what the reconstruction shows. For each item, it gives the location in the code, the value used, the literature source when there is one, and the offline test that checks it. Every reference with a digital object identifier (DOI) was checked against its registry record with Argus and is listed in `docs/refs.bib`. Where a value depends on the event, the value given is that of the benchmark event of 2013 October 25 (Morosan et al. 2022), which is the default configuration, whereas another event takes such values from its own `config.json`.

The status of each item is given by one of the following words.
* **sourced**: the method or value follows the cited work.
* **choice**: a setting of corona3d, documented in METHODS.md but taken from no publication.
* **fitted**: derived from the data by `corona3d gcs`.
* **fixed**: wrong in version 0.0.1 and corrected in 0.1.0, with a test.

## Fixed in version 0.1.0

### F1. The model was rotated by −8.243° without a physical reason
* **Where**: `default_config.json` (`model_dlon`), `reconstruction.py` (`reconstruct`).
* **Problem**: MAS is a steady-state model in Carrington coordinates, and the observers' Carrington longitudes at 15:00 UT already contain the solar rotation, so that rotating the model as well offset it from the GCS shell by 8.2° in longitude.
* **Fix**: `model_dlon` = 0 by default, while the rotation remains available as an option.
* **Test**: `test_model.py::test_model_rotation_turns_the_surface_eastward_or_westward`.

### F2. The perspective camera magnified structures in front of the coronagraph planes
* **Where**: `reconstruction.py` (`render`).
* **Problem**: with the camera 6.65 R⊙ from the Sun, a structure 1 R⊙ in front of the plane of the sky was drawn 18 % larger than the coronagraph image, so that the model and the image did not share one scale.
* **Fix**: the panels are drawn with an orthographic projection, and `camera_distance` was removed.
* **Test**: `test_cli.py::test_half_sets_the_scale_of_the_panels`.

### F3. The CME shell at 15:00 UT was collapsed and lay below the radio heights
* **Where**: `default_config.json` (`gcs`), `corona3d gcs`.
* **Problem**: the shell had been refitted at 15:00 UT to EUV fronts, giving a half-angle of 1.07° (the two legs collapsed into one), an intersection over union of 0 in EUVI-A, and an apex at 1.443 R⊙, below both radio heights.
* **Fix**: the shape of the 15:24 UT white-light fit is kept, and the apex height at 15:00 UT, 1.931 ± 0.158 R⊙, is extrapolated from the COR1-B front, as Morosan et al. (2022) scaled their shell back in time (METHODS.md, Section 7).
* **Test**: there is no offline test, but the fit is reproducible with `corona3d gcs` and its result is kept in `scripts/gcs_fit_2013-10-25.json`.

### F4. The streamer loops were selected by tuned parameters
* **Where**: `reconstruction.py` (`reconstruct`).
* **Problem**: the seed line (meridian −40°, 1.8 R⊙, ±60°) and the apex limit (1.9 R⊙) had been chosen so that the loops overlie the CME, which built the result of the figure into its selection.
* **Fix**: the streamer loops are the closed field lines through the vertices of the shell, with no free parameter.
* **Test**: `test_model.py::test_streamer_lines_are_the_closed_lines_threading_the_shell`.

### F5. Field lines were classed as open at 3 R⊙, and unfinished lines as closed
* **Where**: `geometry.py` (`FieldTracer`), `mast_model.py` (`load`).
* **Problem**: the model was cut at 3 R⊙, so that loops closing above it near the streamer cusp counted as open, and any line that was not open, including one that ran out of steps, counted as closed.
* **Fix**: lines are traced through the whole run (to 30.5 R⊙) with a step of 0.01 r and classed as open, closed, or undetermined. Undetermined lines are counted and not drawn, and the other lines are drawn within `line_r_max` = 3 R⊙.
* **Tests**: `test_geometry.py::test_tracer_follows_dipole_field_lines`, `test_a_line_that_runs_out_of_steps_is_undetermined`.

### F6. Radio positions in arcsec used the solar radius at 1 AU, and clipped positions passed silently
* **Where**: `reconstruction.py` (`radio_points`), `geometry.py` (`deproject`).
* **Fix**: arcsec are converted with the angular radius of the Sun seen from Earth at `time` (964.8 arcsec on 2013 October 25, instead of 959.6), and a position that cannot lie at its height produces a warning.
* **Test**: `test_geometry.py::test_deproject_warns_when_a_source_cannot_lie_at_its_height`.

### F7. Coronagraph frames
* **Where**: `whitelight.py`.
* **Problems**:
  * The base was the earliest file in the folder, whatever `wl_window` said.
  * A `wl_time` before every frame drew a blank image.
  * COR1 triplets were grouped by the minute in the file name.
  * Images smaller than 256 pixels failed to bin.
  * A time with a time zone raised a TypeError.
  * Level-0.5 differences kept the vignetting.
* **Fixes**:
  * The base is the first frame within `wl_window`, and a `wl_time` outside it, or one whose nearest frame is the base itself, is refused.
  * Triplets are grouped by time (60 s), and the binning factor is at least 1.
  * Frames of another size are resampled to the base frame.
  * Times are parsed with Astropy.
  * The default comparison is a ratio, which removes multiplicative instrumental effects.
* **Tests**: all of `test_whitelight.py`.

### F8. The configuration, the downloads, and the command line
* **Where**: `paths.py`, `fetch.py`, `cli.py`, `reconstruction.py`.
* **Fixes**:
  * Nested keys are merged one by one.
  * Interrupted or damaged downloads are removed and counted, and `fetch` ends with exit code 1 when files fail.
  * `make_figure` reads the working directory's `config.json` and writes inside it.
  * `reconstruct` accepts partial parameters and layers.
  * A zero `--half` is refused, and the colour-bar ticks follow `va_range`.
  * The command line no longer loads VTK for `--version` and `init`.
  * Figures are drawn without a GUI backend.
* **Tests**: `test_cli.py`.

### F9. The documentation disagreed with the code
In version 0.0.1, twenty-one statements of the README, METHODS.md, and the user guide disagreed with the code. These statements concerned the error messages, the meaning of `time`, the base frame, the output paths of `make_figure`, the C2 filter, the use of `va002`, the download sizes, the COR2 framing, the 173 MHz density level, and the descriptions of the tests, and all of them were corrected against the code of version 0.1.0.

## Methods and values

| Item | Where | Value | Status | Source | Test |
|---|---|---|---|---|---|
| MAS density unit | `mast_model.py` | 10⁸ cm⁻³ | sourced | PSI MAS documentation; Lionello et al. (2009) | `test_units_and_alfven_speed` |
| MAS field unit | `mast_model.py` | 2.2068914 G | sourced | PSI MAS documentation | `test_units_and_alfven_speed` |
| Density scaling | `density_scale` | 2 | sourced | Morosan et al. (2022) after Wang et al. (2017) | `test_units_and_alfven_speed` |
| Alfvén speed | `mast_model.py` | B / (4π n m_p)^{1/2}, scaled density | sourced (formula); choice (density) | standard MHD | `test_units_and_alfven_speed` |
| Density level | `ne_level` | 9 × 10⁷ cm⁻³ | sourced | Morosan et al. (2022) | none |
| Model run | `run` | `cr2143-medium/hmi_mast_mas_std_0201` | choice | best outline match of 12 public runs | none |
| Model rotation | `model_dlon` | 0 | fixed (F1) | METHODS.md, Section 3 | `test_model_rotation_...` |
| Tracer | `geometry.py` | RK4, step 0.01 r, 1500 steps each way | choice | | `test_tracer_follows_dipole_field_lines` |
| Open / closed / undetermined | `geometry.py` | outer boundary / both ends at inner boundary / otherwise | fixed (F5) | | `test_a_line_that_runs_out_...` |
| Footpoints | `n_seeds`, `seed_r` | 4500 at 1.03 R⊙ | choice | | `test_fibonacci_sphere` |
| Drawing radius of lines | `line_r_max` | 3 R⊙ | choice | | none |
| Streamer loops | `reconstruction.py` | closed lines through the shell vertices above 1.02 R⊙ | fixed (F4) | | `test_streamer_lines_...` |
| GCS model | `geometry.py` (`pythea_shell`) | PyThea, tilt + 90° | sourced | Thernisien (2011); Kouloumvakos et al. (2022) | `test_render_of_the_layers_that_need_no_data` |
| GCS at 15:24 UT | `gcs_wl` | 5.538 R⊙, κ 0.396, α 23.46°, −53.262°, 3.146°, tilt −56.654° (PyThea) | fitted | `corona3d gcs`, step 1 | none |
| GCS at 15:00 UT | `gcs` | apex 1.931 ± 0.158 R⊙, shape as at 15:24 UT | fitted (F3) | `corona3d gcs`, step 2 | none |
| Radio heights | `newkirk_r` | 4 × Newkirk, harmonic; 1.576 and 1.449 R⊙ | sourced | Newkirk (1961); Morosan et al. (2022) | `test_newkirk_heights_of_the_paper` |
| Plasma frequency | `newkirk_r` | 8980 n_e^{1/2} Hz | sourced | standard | `test_newkirk_heights_of_the_paper` |
| Arcsec to R⊙ | `radio_points` | angular solar radius at `time` | fixed (F6) | SunPy | none |
| Observer axes | `observer_basis` | x west, y north, z to observer | sourced | heliocentric-Cartesian convention of SunPy | `test_observer_basis_is_right_handed_with_x_west` |
| Projection | `render` | orthographic | fixed (F2) | | `test_half_sets_the_scale_of_the_panels` |
| Viewpoints | `views`, `view_geometry` | observers or Carrington directions, roll | choice | | `test_views_can_be_set_for_every_panel` |
| COR1 total brightness | `whitelight.py` | (2/3)(I₀ + I₁₂₀ + I₂₄₀) | sourced | polariser triplet relation | `test_cor1_triplets_...` |
| Coronagraph comparison | `wl_mode` | ratio (default) or difference | fixed (F7) | | `test_plane_places_a_feature_in_3d` |
| Fields of view | `INSTR` | C2 2.2-6, COR1 1.5-4, COR2 4-15 R⊙ | choice | instrument fields, trimmed at the occulter edges | `test_plane_places_a_feature_in_3d` |
| Binning | `NBIN` | about 256 pixels across | choice | | `test_small_images_are_not_binned_to_nothing` |
| Panel framing | `views[].frame` | from the published figure | choice | layout of Morosan et al. (2022), Fig. 7 | none |
| NRH image geometry | `nrh.py` | Sun centre at CRPIX (from 1), SOLAR_R pixels per R⊙, north up | sourced | NRH2 file headers (ROT_P = 1); Kerdraon & Delouis (1997) | `test_reader_places_a_source_at_its_sky_position` |
| NRH source centroid | `nrh.source` | brightness-weighted, pixels above 0.5 of the peak | choice | | `test_reader_places_a_source_at_its_sky_position` |
| Radio blobs | `radio_images` | one image per 60 s, peak ≥ 10⁸ K, coloured along a colour map per frequency | choice | the manner of Morosan et al. (2022), Fig. 7 | `test_radio_layer_draws_a_chain_of_blobs_per_frequency` |
| Streamer lines drawn | `streamer_max_lines` | at most 150, an even subset; fainter with the radio layer | choice | | none |
| Perspective option | `projection`, `camera_distance` | 6.65 R⊙ at the default framing | choice | the camera fitted to the published Fig. 7 | `test_perspective_projection_magnifies_the_near_side` |
| OVRO-LWA image geometry | `ovro_lwa.py` | helioprojective header (CRPIX, CDELT, CRVAL in arcsec, RSUN_OBS); band images restored to the full grid by cubic interpolation | sourced | file headers; lwa-solar-util `recover_fits_from_h5` | `test_reader_restores_the_grid_and_places_the_source` |
| OVRO-LWA selection | `ovro_lwa.image_at` | nearest band within 3 MHz, nearest file within 2 min | choice | | `test_reader_restores_the_grid_and_places_the_source` |
