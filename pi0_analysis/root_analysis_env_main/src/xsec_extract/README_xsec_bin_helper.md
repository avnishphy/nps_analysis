# LT bin selection notebook

`xsec_bin_helper.ipynb` compares data coverage and proposes bins. It does not
extract cross sections. Run it in the Hall C Python environment:

```bash
cd /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main
PATH=/group/nps/singhav/software/python/bin:$PATH \
  /group/nps/singhav/software/python/bin/python -m jupyterlab \
  src/xsec_extract/xsec_bin_helper.ipynb
```

Choose Python 3 and Run All Cells. The first code cell verifies the interpreter.
If using an existing Jupyter server, register and select an explicit kernel:

```bash
/group/nps/singhav/software/python/bin/python -m ipykernel install --user \
  --name nps-lt --display-name 'NPS LT (Hall C Python)'
```

Default inputs are the two ROOT files specified in the configuration cell.
For one setting, change `KINEMATICS` and `REFERENCE_KINEMATIC`, or launch with
`NPS_BIN_SETTINGS=KinC_x36_4` preceding the launch command. A comma-separated
list selects several settings. Each setting needs a `DATA_FILES` entry and
exactly one LH2 row in `config/nps_simulation_kinematics*.csv`. Its beam energy
determines the diagnostic epsilon; the table also displays the spectrometer
configuration. No configuration is inferred from the filename alone.

## Interaction

- Click anywhere inside either upper axis to place a corner. Coordinates come
  from the SVG screen transform, not a data point or heatmap cell. All four
  marker IDs persist on both axes. Rapid clicks queue until acknowledged.
- Corners may be entered in any order, but must form four distinct vertices of
  a convex polygon. Their labels retain insertion order. After corner 4 the
  target stays at 4. To revise another corner, choose its number first.
- Select a number and click a new position, or enter xB/Q2 and press **Set
  corner**. **Undo** restores a previous edit, including a move or clear.
- Outer edge fields commit together with **Apply outer limits**. The button
  reads the visible numbers directly, including a field still focused when you
  click; its status confirms the applied selection or explains invalid limits.
  Counts/modes apply immediately. Manual edges are JSON arrays; xB has one array per Q2
  bin and may have different bin counts in different Q2 bins. Array lengths
  govern manual bin counts; the xB count box governs automatic modes.
- **Map min** is the per-setting threshold for coloring a coverage-map cell.
  **Sparse <** flags analysis bins using a separate count threshold. The
  t-prime/phi panel is a marginal sum unless a Q2/xB slice is selected. The
  expandable table lists every full Q2/xB/t-prime/phi bin and every setting,
  including zeros, with sparse bins first. Counts are unweighted entries,
  not statistical precision or effective weighted event counts.
- The copy-configuration panel preserves settings across kernel restarts.
  After editing, rerun the downstream edge review and SIMC coordinate cells.
  Incomplete/invalid selections block downstream use; changing the dashboard
  also invalidates earlier representative-coordinate/export results.

The six-panel figure gives W/Q2 and xB/Q2 the full upper row, with the
other views below. It distinguishes exclusive setting coverage with
blue/orange (additional colors for additional settings), shared subsets with
gray, and all-setting overlap with purple. Pale coverage is outside the
proposed selection. Solid black is the diamond, dashed gray is the outer
window, and dotted lines are bins. The accepted region is their intersection.
The overlap is a finite-resolution occupancy diagnostic, not an acceptance
correction or proof of shared support throughout a cell. The W/Q2 coverage map
uses a 120 x 100 grid, while xB/Q2 remains 60 x 50. **Map min** starts at one
event per setting so the finer W/Q2 cells remain visible; adjust it in the
dashboard to suppress isolated cells. These display settings do not change the
event cut or analysis bin edges.

## Physics and exports

The only polygon mask is in `(xB, Q2)` and is applied identically to every
setting. In the W view, each polygon edge is sampled before transforming
`W = sqrt(Mp**2 + Q2*(1/xB - 1))`; the edges are generally curved. The module
checks this relation against the ROOT W branch before drawing (maximum
allowed difference 0.0005 GeV). Retaining the finite physical missing-mass
sample allows the outer t-prime cut to expand without rereading ROOT.

Uniform, count-quantile, positive-weight yield-quantile, and manual edges use
only selected data. Quantiles pool settings using the existing per-setting
display weights, without introducing a new LT luminosity normalization.
Yield quantiles discard nonpositive/nonfinite weights, as in the original
helper. Phi is uniform on `[0, 2*pi]` after wrapping. All bin interiors are
left-closed/right-open; only the final bin includes the upper outer edge.

Static export is available in the dashboard, or in a later cell:

```python
dashboard.export_figure(REPO_ROOT / 'output/lt_bin_selection/selected')
```

This writes PDF and SVG from the same Matplotlib figure, without controls or
hover elements. PDF embeds fonts; SVG uses vector glyph outlines so physics
symbols do not depend on fonts installed in the browser. Axes and boundaries
are vector graphics; occupancy maps are rasterized at 300 dpi in the export.
The user-specified filenames are replaced if they already exist.

After the interactive plots, rerun **Review edges and diamond for the xsec
configs**. Copy its JSON fields, including `diamond_xb_q2_vertices`, into each
selected `xsec_config/xsec_config_*.json` preset. Four `[xB, Q2]` pairs apply
one common diamond cut to data and SIMC in both extractors; `null` leaves it
off. The W/Q2 view is a projection of these same corners. Validate each edited
preset with `python3 src/xsec_extract/generate_xsec_config.py <preset.json>
/tmp/xsec_config_check.h` before extraction. The generator requires xB edge
rows with matching counts and outer limits.

Scientific JSON export retains the original conditions: one setting, its
configured SIMC file, and explicitly manual Q2/xB edges. SIMC only determines
representative generated Q2i/xBi values weighted by `full_weight` after
reconstructed-bin/exclusive/diamond selection. Geometric t-prime/phi centers
and the existing JSON schema are unchanged. Revision checks reject stale
coordinates. Static figures need neither SIMC nor manual edges.

## Validation commands

The following were run against the full supplied ROOT files, without running
production analysis. The scripts execute copies, preserving the source
notebook. They write executed notebooks, figures, and JSON results under
`output/lt_bin_validation`.

```bash
/group/nps/singhav/software/python/bin/python \
  src/xsec_extract/validate_xsec_bin_helper.py --case both
/group/nps/singhav/software/python/bin/python \
  src/xsec_extract/validate_xsec_bin_helper.py --case single
```

The tests independently compare the polygon with convex half-plane cuts,
audit every full analysis bin and overlap map, check downstream arrays and
representative selection, test all modes, expanded limits, manual conditional
bins, three-setting coverage, and invalid/incomplete/stale-state guards.

For real frontend tests, start a separate authenticated localhost server:

```bash
PATH=/group/nps/singhav/software/python/bin:$PATH \
  /group/nps/singhav/software/python/bin/python -m jupyterlab \
  --no-browser --ip=127.0.0.1 --port=8893 --ServerApp.port_retries=0 \
  --IdentityProvider.token=lt-bin-local-validation \
  --ServerApp.root_dir="$PWD"
```

In another shell, the isolated browser setup used in this session was:

```bash
/group/nps/singhav/software/python/bin/python -m pip install \
  --target /tmp/lt-bin-browser playwright
npm install --prefix /tmp/lt-browser-npm @sparticuz/chromium@131.0.1 \
  --no-audit --no-fund
node -e "require('/tmp/lt-browser-npm/node_modules/@sparticuz/chromium').executablePath().then(p=>console.log(p))"
```

On this host the last command returned `/scratch/singhav/chromium`. Use the
printed path if yours differs. The ordinary Playwright CDN download failed
with `ECONNRESET`; the npm-packaged Chromium ran successfully. No packages
were installed into the shared Python environment.

```bash
PYTHONPATH=/tmp/lt-bin-browser /group/nps/singhav/software/python/bin/python \
  src/xsec_extract/validate_xsec_bin_frontend.py --case both \
  --token lt-bin-local-validation --chromium /scratch/singhav/chromium
PYTHONPATH=/tmp/lt-bin-browser /group/nps/singhav/software/python/bin/python \
  src/xsec_extract/validate_xsec_bin_frontend.py --case single \
  --token lt-bin-local-validation --chromium /scratch/singhav/chromium
```

These open copies of the actual notebook in JupyterLab and run all cells.
All selection edits use browser pointer/form events. A separate kernel
channel reads state and independently audits the resulting data. Tests
cover both 2D entry views, all eight visible corner labels, corner-4 moves,
undo, clear, rapid clicks, numeric coordinates, modes, bins, and outer limits.
Screenshots and browser exports are in `frontend_both/` and `frontend_single/`.
Browser coverage is JupyterLab with headless Chromium; VS Code, classic
Notebook, Firefox, Safari, and touch devices were not verified.

## Preservation and diagnosis

The original untracked notebook and module, including outputs, are preserved
byte-for-byte in `backups/xsec_bin_20260927_145011/`. The edited notebook keeps
its cell IDs and unrelated analysis code; stale outputs/widget state were
cleared. Existing edits elsewhere were left in place.

The previous source coupled editing to Plotly heatmap/scatter callbacks,
which report trace positions rather than arbitrary axis coordinates. Undo
removed a slot instead of restoring a move. On an empty/invalid proposal,
refresh returned before replacing the previous successful arrays and plots.
The W outline connected transformed endpoints with straight lines, and
downstream xB manual overrides could differ from the dashboard. These
mechanisms were replaced rather than patched with another marker-size change.
The earlier user's exact fourth-corner failure was not reproduced in their
original frontend; the replacement's behavior was tested in JupyterLab.
