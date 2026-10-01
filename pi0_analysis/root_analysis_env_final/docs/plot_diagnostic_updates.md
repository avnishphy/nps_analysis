# Plot diagnostics without changing production selection

## Scope

The producer creates all new histograms from actual events, clusters or pairs.
The collector only lays out existing images, retaining both PNG and PDF files.
No masses, weights, acceptance cuts, pair ranking, efficiency calculations,
cross sections or existing selector branches are changed. The missing-mass
correction investigation is deferred.

Existing plot filenames are retained. Four-panel producer canvases receive a
2-by-2 block in the collector; single panels receive one slot. Only exterior
pure-white image margins are trimmed, with padding. No subplot is inferred,
redrawn, stretched disproportionately or reconstructed from an image.

## Before and after populations

- HMS: events after the existing good-event filter; the after histogram passes the
  complete existing HMS selection. Unavailable branches stay unavailable.
- NPS: packed clusters after HMS and the existing multiplicity prerequisite;
  after histogram passes cluster energy/position/time and dead-block cuts, before
  the existing top-four truncation. These are cluster counts, not event yields.
- The `nclust` panel inside `cut_debug_run<RUN>.pdf` draws the existing
  `h_nclusters_run<RUN>` histogram:
  clusters per event after HMS, before the two-cluster prerequisite, packing,
  and NPS cluster cuts. It includes zero/one-cluster events and retains the
  unit-width integer bins from 0 through 1080 (the detector block count), with
  flow bins retained. Multiplicity is taken before the 20-cluster processing
  cap; the displayed range follows occupied bins.
  Its red overlay (`h_nclusters_run<RUN>_after`) counts input clusters passing
  energy, position, timing and dead-block cuts in each of the same HMS events, including
  clusters beyond the processing cap. Zero surviving clusters fill zero;
  no two-cluster gate or top-four truncation is applied to this diagnostic.
  Legend totals sum the actual clusters before/after cuts across those events,
  independently of display limits or flow bins. Histogram entries still count
  the same HMS events, including zero-cluster events. The view includes zero
  with independent axes: red after-cuts multiplicity/events on the bottom/left,
  blue before-cuts multiplicity/events on the top/right. Each x range follows
  its own population's occupied bins with padding and includes zero. Each y
  axis uses its own population's visible maximum.
  Only a disposable display copy is mapped between axis coordinates; stored
  histogram counts/errors and cluster totals are unchanged. There is no separate `nclust` PNG.
- Pairs: all unordered pairs from the existing eligible/top-four cluster list
  versus the pair selected by the unchanged producer. The two-cluster path
  does not impose the same timing check as the higher-multiplicity path.
  Timing lines are therefore labeled as a reference, not a universal gate.
- Existing corrected-mass cut diagnostic: same candidate population, before
  and after its existing gate. Its mass definition and unweighted counts are
  unchanged; it does not represent the separate positive-weight requirement.
- Physics/mass comparisons: the same positive pi0-weight sample and actual
  stored selectors. Ellipse/MCD curves are overlaid on existing canvases.
  Combined diagnostic histograms use event flags, not bin-center membership.

1D before/after histograms remain overlaid, with legends at the top right.
2D cut diagnostics show two variables per page: each row has before on the
left and after on the right, with identical spatial axes. NPS cluster-position
maps use each population's own color maximum so the after-cut structure is
visible; other paired maps keep shared color limits. Each panel
contains only one 2D histogram. Mass canvases likewise keep separate before/after
density maps; analytical cut/reference lines remain visible, and the selectors
are compared in the 1D projections.

New per-run ROOT histograms append `_after`, `_before_selection`,
`_ellipse_compare`, or `_mcd_compare` to the corresponding existing histogram
name. Cluster x/y diagnostics, including before/after histograms and the 2D
position map, use 2 cm bins: 34 along x and 40 along y. They are filled directly
from clusters during analysis. Other existing binning is preserved, except that the older
`h_clustXY` diagnostic map expands from x=[-30,30], y=[-36,36] to
x=[-34,34], y=[-40,40] cm at its original 2 cm bin width. Its extra bins are
filled directly by the producer; coordinates are not stretched or reconstructed.
New `h_mass_compare_*` histograms and `mass_cut_comparison` canvas live in the
existing `mass_cut_run<RUN>.root` file.

## Ellipse qualification and fallback

Numerical solver validity alone is not adequate: run 2463 previously reported
a valid ellipse built from four bins, accepting 0.0303725 percent of the
weighted population inside the fit window, while MCD accepted 46.3516 percent.

The plotting qualification requires the existing solver validity, at least
the configured `auto_min_core_bins` (8) occupied fit bins, at least the configured
`auto_min_core_total_fraction` (0.005) fit population, and finite positive
covariance determinant. These are the fitter's existing qualification limits,
not a new threshold chosen to obtain a desired accepted yield.

Failure is printed in the log and on the existing comparison canvas. Its
after-selection panel uses the existing `is_exclusive` de-correlation selector.
The original ellipse and MCD selectors remain visible and unchanged, including
an explicit failed-ellipse label. A failure does not silently relabel a stored
ellipse flag as de-correlation. This qualification catches documented collapsed
fits; it is not proof of signal purity or full cross-kinematics validation.

Per-run status is saved as `diagnostic_mass_cut_status` in both diagnostics
and mass-cut ROOT outputs. Parameter CSV/debug text fields
`diagnostic_ellipse_valid` and `diagnostic_fallback_decorrelation` distinguish
plot qualification from the unchanged solver's `ellipse_valid` field.

The dashed reference is

```
Mmiss = 0.938 - 31.95 * (Minv - 0.135)  # both masses in GeV
```

It is an overlay, not a new mass correction or selection. No new band width
or fitted slope is silently chosen. Production adoption of a new band requires
the separate signal/background and cross-run validation discussed with the user.

## Plotting workflow

The only user-facing plotting command remains:

```bash
python3 scripts/collect_kinematic_plots.py --kin KinC_x60_4a
```

Its existing input/output defaults and options are unchanged. No range
preparation command, CSV setup, or extra plotting script is required.

During normal analysis, the producer automatically selects the view from each
pre-cut histogram's occupied bins, adds a 5 percent margin (at least two bins),
and includes configured cut boundaries. All overlaid selections use that view.
NPS cluster x and y always use [-34,34] and [-40,40] cm, respectively, in both
1D histograms and 2D position maps, regardless of occupancy or cut boundaries.
Paired 2D panels use the same pre-cut spatial view, even if no entries pass.
An empty NPS after-cut map uses a 0-to-1 color scale. Cluster-energy and corrected
missing-mass cut-debug panels start at exactly 0.4 GeV; the energy frame may
clip a displayed bin without changing its contents. Legend counts and retained
fractions still use the full populations, including entries below 0.4 GeV.
The cluster overview energy panel and combined missing-mass overlay also start
at 0.4 GeV. Mass-fit windows and physics selections are unchanged.
"Outside view" and related flow-count annotations are removed from plots;
the actual flow bins and numerical mass-window metadata are retained.
Ranges may differ between runs; combined plots use their combined pre-cut
population. Full histogram contents and flow bins remain in ROOT, with
`diagnostic_range_source` recording how the view was chosen.

The normal analysis workflow supplies the histograms and source plots. The
collector assembles whatever outputs already exist, including older outputs;
it never invents missing before/after histograms or launches production jobs.
New before/after content appears after the normal analysis has produced it.

The `nclust` histogram is cleaned up only through `cut_debug_plots`. Its earlier
addition to that list left it in the cluster cleanup array as well, causing a
second deletion after outputs were written. Removing that duplicate ownership
fixes the exit-time segmentation violation without changing any event processing.

## Validation commands

Python regressions:

```bash
PYTHONPYCACHEPREFIX=/tmp/nps-pycache MPLCONFIGDIR=/tmp/nps-mpl-cache \
  OPENBLAS_NUM_THREADS=1 python3 -m unittest discover -s tests -p 'test_*py'
```

ROOT macro syntax and diagnostic smoke tests (Hall C environment required):

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ -std=c++17 -fsyntax-only -include tests/root_compile_prelude.h `root-config --cflags` src/analysis/nps_analysis_main.C'
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; g++ -std=c++17 `root-config --cflags` tests/test_nps_plot_diagnostics.cpp -o /tmp/nps_plot_diagnostics_test `root-config --libs` && /tmp/nps_plot_diagnostics_test /tmp'
```

The test prelude supplies the standard-library names already exposed in ROOT's
interactive environment; legacy physics helper headers are not modified.
The smoke test verifies unchanged coordinates/weights, flow-bin preservation,
actual selector histogram integrals, and visible failure/fallback outputs.
It also checks fixed detector bounds, one 2D histogram per panel, shared
before/after spatial axes, independent NPS color limits, and the empty-after case.

Checks completed during implementation: eight Python regressions, ROOT macro
syntax check, compiled ROOT diagnostic smoke test, and visual inspection of a
collector preview using existing run-4224 images. Full production was not rerun.

Cleanup regression: a complete run 5237 with `--gevnum-cut yes` passed with zero
failed jobs; its summary CSV matches the output written before the cleanup fix.
The original cleanup blocks reproduce a segmentation violation; the fixed blocks
exit successfully. To check the full producer in an isolated output directory:

```bash
csh -c 'source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash src/analysis/run_parallel_nps_analysis_main.sh --kin KinC_x50_0a --target LH2 --run 5237 --jobs 1 --no-combine --gevnum-cut yes --output-base /tmp/nps_cleanup_validation'
```
