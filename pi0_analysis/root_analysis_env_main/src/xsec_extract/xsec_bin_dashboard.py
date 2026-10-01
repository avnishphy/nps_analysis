"""Data-only LT bin selection with an exact-coordinate SVG editor.

Diagnostic tool: no SIMC input or cross-section calculation. Matplotlib supplies
both the Jupyter display and publication exports. Counts are never sampled.
"""

from html import escape
from io import StringIO
from pathlib import Path
import json

import anywidget
import ipywidgets as widgets
import matplotlib as mpl
from matplotlib.figure import Figure
from matplotlib.path import Path as PolygonPath
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np
import pandas as pd
import traitlets
from IPython.display import display

COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
MODE_OPTIONS = ("uniform", "count_quantile", "yield_quantile", "manual")
LABELS = {
    "tprime": r"$t'\;[\mathrm{GeV}^2]$",
    "Q2": r"$Q^2\;[\mathrm{GeV}^2]$",
    "xB": r"$x_B$",
    "W": r"$W\;[\mathrm{GeV}]$",
}


def _quantile_edges(values, weights, limits, count, mode, manual=None):
    lo, hi = limits
    if not np.isfinite([lo, hi]).all() or lo >= hi or count < 1:
        raise ValueError("Outer limits must increase; bin count must be positive")
    if mode == "manual":
        if manual is None:
            raise ValueError("Manual mode requires an edge array")
        edges = np.asarray(manual, dtype=float)
    elif mode == "uniform":
        edges = np.linspace(lo, hi, count + 1)
    elif mode in ("count_quantile", "yield_quantile"):
        values, weights = np.asarray(values, float), np.asarray(weights, float)
        keep = np.isfinite(values) & (values >= lo) & (values <= hi)
        if mode == "yield_quantile":
            keep &= np.isfinite(weights) & (weights > 0)
        values, weights = values[keep], weights[keep]
        if not len(values):
            raise ValueError(
                "No selected finite data (positive weights for yield quantiles)"
            )
        probabilities = np.arange(1, count) / count
        if mode == "count_quantile":
            interior = np.quantile(values, probabilities)
        else:
            order = np.argsort(values)
            values, weights = values[order], weights[order]
            cumulative = (np.cumsum(weights) - 0.5 * weights) / weights.sum()
            interior = np.interp(probabilities, cumulative, values)
        edges = np.r_[lo, interior, hi]
    else:
        raise ValueError(f"Unknown edge mode: {mode}")
    if (
        edges.ndim != 1
        or len(edges) < 2
        or not np.isfinite(edges).all()
        or np.any(np.diff(edges) <= 0)
    ):
        raise ValueError(f"Invalid or collapsed {mode} edges: {edges}")
    if not np.allclose([edges[0], edges[-1]], limits, rtol=0, atol=1e-10):
        raise ValueError(f"Manual edges must span the current limits {limits}")
    return edges


def ordered_polygon(vertices):
    """Convex boundary in CCW order; marker IDs retain insertion order."""
    if vertices is None or any(v is None for v in vertices):
        return None
    v = np.asarray(vertices, float)
    if (
        v.shape != (4, 2)
        or not np.isfinite(v).all()
        or (v <= 0).any()
        or (v[:, 0] >= 1).any()
    ):
        raise ValueError("Use four physical corners: 0 < xB < 1 and Q2 > 0")
    centered = v - v.mean(axis=0)
    v = v[np.argsort(np.arctan2(centered[:, 1], centered[:, 0]))]
    d = np.roll(v, -1, axis=0) - v
    cross = d[:, 0] * np.roll(d[:, 1], -1) - d[:, 1] * np.roll(d[:, 0], -1)
    if np.any(cross <= 1e-12):
        raise ValueError(
            "Four distinct convex corners required; move the interior/collinear corner"
        )
    return v


def polygon_mask(frame, polygon):
    if polygon is None:
        return np.ones(len(frame), dtype=bool)
    return PolygonPath(np.vstack([polygon, polygon[0]]), closed=True).contains_points(
        frame[["xB", "Q2"]].to_numpy(), radius=1e-9
    )


def in_bin(values, edges, index):
    return (values >= edges[index]) & (
        values <= edges[index + 1]
        if index == len(edges) - 2
        else values < edges[index + 1]
    )


class SelectionCanvas(anywidget.AnyWidget):
    _esm = Path(__file__).with_name("xsec_bin_canvas.js")
    svg = traitlets.Unicode("").tag(sync=True)
    geometry = traitlets.List([]).tag(sync=True)


class OuterLimitApply(anywidget.AnyWidget):
    _esm = Path(__file__).with_name("xsec_bin_limits.js")


class BinDashboard:
    def __init__(
        self,
        settings,
        coverage_frames,
        proton_mass,
        limits,
        modes,
        counts,
        vertices=None,
        manual_edges=None,
        on_change=None,
        overlap_bins=(60, 50),
        w_q2_bins=None,
        min_events_per_cell=3,
        max_scatter=4000,
    ):
        if (
            not settings
            or len(settings) != len(coverage_frames)
            or len(settings) > len(COLORS)
        ):
            raise ValueError(
                "Supply one data frame per setting (1 to 6 distinct settings)"
            )
        if len(set(settings)) != len(settings) or any(f.empty for f in coverage_frames):
            raise ValueError(
                "Setting names must be distinct; coverage samples must be nonempty"
            )
        self.settings, self.coverage_frames = list(settings), coverage_frames
        self.proton_mass, self.on_change = proton_mass, on_change
        self.overlap_bins = overlap_bins
        self.w_q2_bins = overlap_bins if w_q2_bins is None else w_q2_bins
        self.manual_edges = manual_edges or {}
        self.vertices = (
            [None] * 4 if vertices is None else [tuple(map(float, v)) for v in vertices]
        )
        if len(self.vertices) != 4:
            raise ValueError("Exactly four corner slots required")
        self._history, self._busy = [], False
        self.revision, self.valid, self.error = 0, False, ""
        self.selected_frames, self.edges = [], {}
        self.last_click = "Choose four corners on either upper plot."
        self.projection_checks = []
        for name, frame in zip(settings, coverage_frames):
            delta = np.abs(self._w(frame.xB, frame.Q2) - frame.W)
            if not np.isfinite(delta).all() or delta.max() > 5e-4:
                raise ValueError(
                    f"{name}: W branch disagrees with xB/Q2 projection (max {delta.max():.6g} GeV)"
                )
            self.projection_checks.append(
                {"setting": name, "max_abs_W_delta_GeV": float(delta.max())}
            )

        self.canvas = SelectionCanvas()
        self.canvas.on_msg(self._canvas_message)
        self.status, self.feedback = widgets.HTML(), widgets.HTML()
        self.corner = widgets.Dropdown(
            description="Corner",
            options=[(str(i), i - 1) for i in range(1, 5)],
            layout=widgets.Layout(width="135px"),
        )
        self.corner_x = widgets.FloatText(
            description="xB", value=0.36, layout=widgets.Layout(width="165px")
        )
        self.corner_q = widgets.FloatText(
            description="Q² [GeV²]", value=4.0, layout=widgets.Layout(width="190px")
        )
        self.apply_corner = widgets.Button(
            description="Set corner", layout=widgets.Layout(width="100px")
        )
        self.undo = widgets.Button(
            description="Undo", layout=widgets.Layout(width="80px")
        )
        self.reset = widgets.Button(
            description="Clear", layout=widgets.Layout(width="80px")
        )
        self.corner.observe(self._corner_changed, names="value")
        self.apply_corner.on_click(
            lambda _: self.place_corner(self.corner_x.value, self.corner_q.value)
        )
        self.undo.on_click(self._undo)
        self.reset.on_click(self._reset)
        rows, self.limit_fields, self.manual_fields = [], {}, {}
        for key, short, prefix in (
            ("tprime", "t'", "t"),
            ("Q2", "Q²", "q"),
            ("xB", "xB", "x"),
        ):
            # Numeric outer edges commit together and can expand beyond the initial cut.
            control = widgets.FloatRangeSlider(value=limits[key], min=-1e6, max=1e6)
            setattr(self, prefix + "_limits", control)
            mode = widgets.Dropdown(
                options=MODE_OPTIONS,
                value=(
                    "manual"
                    if key == "xB" and self.manual_edges.get(key) is not None
                    else modes[key]
                ),
                description=f"{short} mode",
                layout=widgets.Layout(width="245px"),
            )
            count = widgets.BoundedIntText(
                value=counts[key],
                min=1,
                max=max(20, counts[key]),
                description="Bins",
                layout=widgets.Layout(width="125px"),
            )
            setattr(self, prefix + "_mode", mode)
            setattr(self, "n_" + prefix, count)
            fields = [
                widgets.FloatText(
                    value=v,
                    description=f"{short} {bound}",
                    layout=widgets.Layout(width="180px"),
                )
                for v, bound in zip(limits[key], ("min", "max"))
            ]
            self.limit_fields[key] = fields
            for field, bound in zip(fields, ("min", "max")):
                field.add_class(f"outer-limit-{prefix}-{bound}")
            rows.append(widgets.HBox([*fields, mode, count]))
            self.manual_fields[key] = widgets.Text(
                value=json.dumps(self.manual_edges.get(key)),
                description=short,
                layout=widgets.Layout(width="95%"),
            )
            for c in (control, mode, count):
                c.observe(self._control_changed, names="value")
        self.apply_limits = OuterLimitApply()
        self.apply_limits.on_msg(self._apply_limits_message)
        self.apply_manual = widgets.Button(
            description="Apply manual arrays", layout=widgets.Layout(width="180px")
        )
        self.apply_manual.on_click(self._apply_manual)
        self.n_phi = widgets.BoundedIntText(
            value=counts["phi"],
            min=1,
            max=72,
            description="φ bins",
            layout=widgets.Layout(width="140px"),
        )
        self.cell_min = widgets.BoundedIntText(
            value=min_events_per_cell,
            min=1,
            max=10000,
            description="Map min",
            layout=widgets.Layout(width="150px"),
        )
        self.sparse_min = widgets.BoundedIntText(
            value=10,
            min=1,
            max=100000,
            description="Sparse <",
            layout=widgets.Layout(width="160px"),
        )
        self.slice = widgets.Dropdown(
            description="t'–φ slice",
            options=[("All Q²/xB", -1)],
            layout=widgets.Layout(width="280px"),
        )
        for c in (self.n_phi, self.cell_min, self.sparse_min, self.slice):
            c.observe(self._control_changed, names="value")
        self.copy_code = widgets.Textarea(
            layout=widgets.Layout(width="100%", height="175px")
        )
        self.occupancy_table = widgets.HTML(
            layout=widgets.Layout(max_height="260px", overflow="auto")
        )
        self.export_path = widgets.Text(
            value="lt_bin_selection",
            description="Figure path",
            layout=widgets.Layout(width="440px"),
        )
        self.export_button = widgets.Button(
            description="Export PDF + SVG", layout=widgets.Layout(width="175px")
        )
        self.export_message = widgets.HTML()
        self.export_button.on_click(self._export_clicked)
        details = widgets.Accordion(
            children=[
                widgets.VBox(
                    [
                        widgets.HTML(
                            "JSON arrays; xB: one array per Q² bin. Array lengths determine manual bin counts."
                        ),
                        *self.manual_fields.values(),
                        self.apply_manual,
                    ]
                ),
                self.occupancy_table,
                self.copy_code,
                widgets.VBox(
                    [
                        widgets.HBox([self.export_path, self.export_button]),
                        self.export_message,
                    ]
                ),
            ],
            selected_index=None,
        )
        for i, title in enumerate(
            (
                "Manual edges",
                "Every proposed bin: per-setting counts (sparse first)",
                "Copy configuration for a future run",
                "Static figure export",
            )
        ):
            details.set_title(i, title)
        self.widget = widgets.VBox(
            [
                widgets.HTML(
                    "<b>LT bin selection · data only</b> — click to place; select a corner to move it. Numbers persist after closure. Numeric fields + Set corner perform the same operation."
                ),
                widgets.HBox(
                    [
                        self.corner,
                        self.corner_x,
                        self.corner_q,
                        self.apply_corner,
                        self.undo,
                        self.reset,
                    ]
                ),
                self.feedback,
                self.status,
                self.canvas,
                *rows,
                widgets.HBox(
                    [
                        self.apply_limits,
                        self.n_phi,
                        self.cell_min,
                        self.sparse_min,
                        self.slice,
                    ]
                ),
                details,
            ],
            layout=widgets.Layout(width="100%", max_width="1120px"),
        )
        self.widget.add_class("lt-dashboard")
        self.refresh()

    def _w(self, xb, q2):
        return np.sqrt(self.proton_mass**2 + np.asarray(q2) * (1 / np.asarray(xb) - 1))

    def _polygon(self):
        return ordered_polygon(self.vertices)

    _polygon_mask = staticmethod(polygon_mask)

    def _corner_changed(self, _):
        vertex = self.vertices[self.corner.value]
        if vertex is not None:
            self.corner_x.value, self.corner_q.value = vertex

    def place_corner(self, xb, q2):
        if not np.isfinite([xb, q2]).all() or not 0 < xb < 1 or q2 <= 0:
            self.feedback.value = (
                "<b>Corner rejected: require 0 &lt; xB &lt; 1 and Q² &gt; 0.</b>"
            )
            return
        self._history.append((self.vertices.copy(), self.corner.value))
        index = self.corner.value
        self.vertices[index] = (float(xb), float(q2))
        self.last_click = f"Corner {index+1} accepted: xB={xb:.6f}, Q²={q2:.6f} GeV²."
        # After closure stay on corner 4 instead of silently wrapping to corner 1.
        if any(v is None for v in self.vertices):
            self.corner.value = self.vertices.index(None)
        self._corner_changed(None)
        self.refresh()

    def _canvas_message(self, _, content, buffers):
        if content.get("kind") != "place":
            return
        try:
            x, y = float(content["x"]), float(content["y"])
            xb, q2 = (
                (x, y)
                if content["view"] == "xB"
                else (x / (y * y - self.proton_mass**2 + x), x)
            )
            self.place_corner(xb, q2)
        except (ValueError, KeyError, ZeroDivisionError) as exc:
            self.feedback.value = f"<b>Corner rejected: {escape(str(exc))}</b>"
        finally:
            self.canvas.send(
                {
                    "kind": "ack",
                    "message": self.feedback.value.replace("<b>", "").replace(
                        "</b>", ""
                    ),
                }
            )

    def _undo(self, _=None):
        if self._history:
            self.vertices, self.corner.value = self._history.pop()
            self.last_click = "Undid last corner edit (including moves or clear)."
            self._corner_changed(None)
            self.refresh()

    def _reset(self, _=None):
        self._history.append((self.vertices.copy(), self.corner.value))
        self.vertices, self.corner.value = [None] * 4, 0
        self.last_click = "Diamond cleared; outer limits remain active."
        self.refresh()

    def _control_changed(self, _):
        if not self._busy:
            self.refresh()

    def _apply_limits_message(self, _, content, buffers):
        if content.get("kind") != "apply_limits":
            return
        try:
            applied = self._apply_limits(limits=content.get("limits"))
        except Exception:
            self.apply_limits.send({
                "kind": "limits_ack", "message": "Apply failed; see the kernel error."
            })
            raise
        if applied:
            message = (
                self.last_click if self.valid else
                f"Outer limits applied; selection invalid: {self.error}"
            )
        else:
            message = "Limits rejected; check the values and the dashboard feedback."
        self.apply_limits.send({"kind": "limits_ack", "message": message})

    def _apply_limits(self, _=None, limits=None):
        if limits is None:
            limits = {
                k: tuple(f.value for f in fields)
                for k, fields in self.limit_fields.items()
            }
        try:
            if not isinstance(limits, dict) or set(limits) != set(self.limit_fields):
                raise ValueError("Supply tprime, Q2 and xB limits")
            limits = {
                key: tuple(float(value) for value in limits[key])
                for key in ("tprime", "Q2", "xB")
            }
            if any(len(pair) != 2 for pair in limits.values()):
                raise ValueError("Each limit needs a minimum and maximum")
        except (TypeError, ValueError, KeyError):
            self.feedback.value = "<b>Limits rejected: enter two numeric edges for t', Q² and xB.</b>"
            return False
        if (
            any(not np.isfinite(v).all() or v[0] >= v[1] for v in limits.values())
            or not 0 < limits["xB"][0] < limits["xB"][1] < 1
            or limits["Q2"][0] <= 0
            or any(abs(value) > 1e6 for pair in limits.values() for value in pair)
        ):
            self.feedback.value = "<b>Limits rejected: require increasing finite edges within ±1e6, 0 &lt; xB &lt; 1, Q² &gt; 0.</b>"
            return False
        self._busy = True
        try:
            for key, c in (
                ("tprime", self.t_limits),
                ("Q2", self.q_limits),
                ("xB", self.x_limits),
            ):
                c.value = limits[key]
                for field, value in zip(self.limit_fields[key], limits[key]):
                    field.value = value
        finally:
            self._busy = False
        self.last_click = "Outer limits applied."
        self.refresh()
        if not self.valid:
            self.feedback.value = (
                f"<b>Outer limits applied; selection invalid: {escape(self.error)}.</b>"
            )
        return True

    def _apply_manual(self, _=None):
        try:
            parsed = {k: json.loads(f.value) for k, f in self.manual_fields.items()}
            for key in ("tprime", "Q2"):
                if parsed[key] is not None and np.asarray(parsed[key]).ndim != 1:
                    raise ValueError(f"{key} requires a flat array")
            if parsed["xB"] is not None and (
                not isinstance(parsed["xB"], list)
                or any(np.asarray(v).ndim != 1 for v in parsed["xB"])
            ):
                raise ValueError("xB requires a list of edge arrays")
        except (ValueError, TypeError) as exc:
            self.feedback.value = f"<b>Manual input rejected: {escape(str(exc))}</b>"
            return
        self.manual_edges = parsed
        self._busy = True
        try:
            for key, c in (
                ("tprime", self.t_mode),
                ("Q2", self.q_mode),
                ("xB", self.x_mode),
            ):
                if parsed[key] is not None:
                    c.value = "manual"
        finally:
            self._busy = False
        self.last_click = "Manual arrays applied."
        self.refresh()

    def require_ready(self, revision=None):
        if not self.valid:
            raise RuntimeError(f"Current dashboard selection is invalid: {self.error}")
        if 0 < sum(v is not None for v in self.vertices) < 4:
            raise RuntimeError(
                "Finish all four corners or Clear before using downstream cells"
            )
        if revision is not None and revision != self.revision:
            raise RuntimeError(
                "Dashboard changed: rerun edge review and SIMC coordinate cells before export"
            )

    def refresh(self):
        if self._busy:
            return
        self._busy = True
        self.revision += 1
        self.valid, self.error, self.edges = False, "", {}
        try:
            self.limits = {
                "tprime": tuple(self.t_limits.value),
                "Q2": tuple(self.q_limits.value),
                "xB": tuple(self.x_limits.value),
            }
            # Keep coverage outside Q2/xB cuts visible. No downsampling in counts.
            self.current_coverage = [
                f.loc[f.tprime.between(*self.limits["tprime"])]
                for f in self.coverage_frames
            ]
            polygon = self._polygon()
            self.selected_frames = [
                f.loc[
                    f.Q2.between(*self.limits["Q2"])
                    & f.xB.between(*self.limits["xB"])
                    & polygon_mask(f, polygon)
                ].copy()
                for f in self.current_coverage
            ]
            pooled = pd.concat(self.selected_frames, ignore_index=True)
            if pooled.empty:
                raise ValueError(
                    "No selected data; move corners or expand outer limits"
                )
            e = {}
            for key, count, mode in (
                ("tprime", self.n_t, self.t_mode),
                ("Q2", self.n_q, self.q_mode),
            ):
                e[key] = _quantile_edges(
                    pooled[key],
                    pooled.data_weight,
                    self.limits[key],
                    count.value,
                    mode.value,
                    self.manual_edges.get(key),
                )
                if mode.value == "manual":
                    count.max = max(count.max, len(e[key]) - 1)
                    count.value = len(e[key]) - 1
            manual_x = self.manual_edges.get("xB")
            if self.x_mode.value == "manual" and (
                manual_x is None or len(manual_x) != len(e["Q2"]) - 1
            ):
                raise ValueError("Manual xB: provide one edge array per current Q² bin")
            e["xB_by_Q2"] = []
            for iq in range(len(e["Q2"]) - 1):
                subset = pooled.loc[in_bin(pooled.Q2, e["Q2"], iq)]
                manual = manual_x[iq] if self.x_mode.value == "manual" else None
                e["xB_by_Q2"].append(
                    _quantile_edges(
                        subset.xB,
                        subset.data_weight,
                        self.limits["xB"],
                        self.n_x.value,
                        self.x_mode.value,
                        manual,
                    )
                )
            e["phi"] = np.linspace(0, 2 * np.pi, self.n_phi.value + 1)
            self.edges, self.valid = e, True
        except (ValueError, TypeError, IndexError) as exc:
            self.error = str(exc)
            self.selected_frames = [f.iloc[:0].copy() for f in self.coverage_frames]
        try:
            for count, mode in (
                (self.n_t, self.t_mode),
                (self.n_q, self.q_mode),
                (self.n_x, self.x_mode),
            ):
                count.disabled = mode.value == "manual"
            self.feedback.value = f"<b>{escape(self.last_click)}</b>"
            self._compute_occupancy()
            self._update_status()
            self._render()
            self._save_configuration()
            if self.on_change is not None:
                self.on_change(self)
        finally:
            self._busy = False

    def _compute_occupancy(self):
        self.occupancy_rows, self.occupancy_matrices = [], []
        if not self.valid:
            self.occupancy_table.value = "No valid bins."
            return
        e = self.edges
        choices, slices = [("All Q²/xB", -1)], []
        for iq, xb in enumerate(e["xB_by_Q2"]):
            for ix in range(len(xb) - 1):
                choices.append((f"Q² bin {iq+1}, xB bin {ix+1}", len(slices)))
                slices.append((iq, ix))
        old_slice = self.slice.value
        self.slice.options = choices
        self.slice.value = old_slice if old_slice in [v for _, v in choices] else -1
        for setting, frame in zip(self.settings, self.selected_frames):
            displayed = frame
            for iq, ix in slices:
                subset = frame.loc[
                    in_bin(frame.Q2, e["Q2"], iq)
                    & in_bin(frame.xB, e["xB_by_Q2"][iq], ix)
                ]
                h = np.histogram2d(
                    subset.tprime, subset.phi_wrapped, bins=(e["tprime"], e["phi"])
                )[0].astype(int)
                for it, ip in np.ndindex(h.shape):
                    self.occupancy_rows.append(
                        {
                            "Q2 bin": iq + 1,
                            "xB bin": ix + 1,
                            "t' bin": it + 1,
                            "phi bin": ip + 1,
                            "setting": setting,
                            "events": int(h[it, ip]),
                            "sparse": bool(h[it, ip] < self.sparse_min.value),
                        }
                    )
                if self.slice.value != -1 and (iq, ix) == slices[self.slice.value]:
                    displayed = subset
            self.occupancy_matrices.append(
                np.histogram2d(
                    displayed.tprime,
                    displayed.phi_wrapped,
                    bins=(e["tprime"], e["phi"]),
                )[0].astype(int)
            )
        self.occupancy_frame = pd.DataFrame(self.occupancy_rows)
        self.occupancy_table.value = self.occupancy_frame.sort_values(
            ["sparse", "events"], ascending=[False, True]
        ).to_html(index=False)

    def _update_status(self):
        if not self.valid:
            self.status.value = f"<b style='color:#a00'>INVALID — {escape(self.error)}. Downstream use/export blocked.</b>"
            return
        counts = " · ".join(
            f"<span style='color:{COLORS[i]}'>{escape(s)}: <b>{len(f):,}</b> selected / {len(c):,} coverage</span>"
            for i, (s, f, c) in enumerate(
                zip(self.settings, self.selected_frames, self.current_coverage)
            )
        )
        n = sum(v is not None for v in self.vertices)
        state = (
            "active"
            if n == 4
            else (
                "none (outer limits only)"
                if n == 0
                else f"DRAFT {n}/4; outer limits only, downstream blocked"
            )
        )
        sparse = self.occupancy_frame.groupby("setting", sort=False).sparse.sum()
        self.status.value = (
            f"{counts}<br>Diamond: <b>{state}</b>. Sparse full bins (&lt;{self.sparse_min.value} events): "
            + "; ".join(f"{escape(s)}: {int(sparse[s])}" for s in self.settings)
            + ". Per-setting counts below; t′–φ sums over Q²/xB unless a slice is selected."
        )

    def _boundary(self, polygon):
        return np.vstack(
            [
                np.linspace(a, b, 65, endpoint=False)
                for a, b in zip(polygon, np.roll(polygon, -1, axis=0))
            ]
            + [polygon[:1]]
        )

    def _render(self):
        with mpl.rc_context(
            {
                "font.family": "DejaVu Sans",
                "font.size": 9,
                "axes.labelsize": 10,
                "axes.titlesize": 10,
                "axes.titleweight": "normal",
                "axes.grid": False,
                "svg.fonttype": "path",
                "pdf.fonttype": 42,
            }
        ):
            # Give the two interactive coverage maps a full row and taller axes.
            fig = Figure(figsize=(12.5, 10.2), dpi=130)
            axes = fig.subplots(
                3, 2, gridspec_kw={"height_ratios": (1.55, 1.0, 1.0)}
            ).ravel()
            fig.subplots_adjust(
                left=0.065, right=0.95, bottom=0.065, top=0.88,
                wspace=0.30, hspace=0.42,
            )
            self.map_stats = []
            for imap, (ax, xcol, ycol) in enumerate(
                ((axes[0], "Q2", "W"), (axes[1], "xB", "Q2"))
            ):
                all_x = np.concatenate([f[xcol] for f in self.coverage_frames])
                all_y = np.concatenate([f[ycol] for f in self.coverage_frames])
                markers = [(i, v) for i, v in enumerate(self.vertices) if v is not None]
                vx = [v[1] if xcol == "Q2" else v[0] for _, v in markers]
                vy = [self._w(*v) if xcol == "Q2" else v[1] for _, v in markers]
                qlo, qhi = self.limits["Q2"]
                xlo, xhi = self.limits["xB"]
                outer = self._boundary(
                    np.array([[xlo, qlo], [xhi, qlo], [xhi, qhi], [xlo, qhi]])
                )
                ox = outer[:, 1] if xcol == "Q2" else outer[:, 0]
                oy = self._w(outer[:, 0], outer[:, 1]) if xcol == "Q2" else outer[:, 1]
                xlim = [
                    min(np.quantile(all_x, 0.001), *ox, *vx),
                    max(np.quantile(all_x, 0.999), *ox, *vx),
                ]
                ylim = [
                    min(np.quantile(all_y, 0.001), *oy, *vy),
                    max(np.quantile(all_y, 0.999), *oy, *vy),
                ]
                for lim in (xlim, ylim):
                    pad = max(lim[1] - lim[0], 0.01) * 0.055
                    lim[0] -= pad
                    lim[1] += pad
                map_bins = self.w_q2_bins if imap == 0 else self.overlap_bins
                xe, ye = np.linspace(*xlim, map_bins[0] + 1), np.linspace(
                    *ylim, map_bins[1] + 1
                )
                hs = np.stack(
                    [
                        np.histogram2d(f[xcol], f[ycol], bins=(xe, ye))[0]
                        for f in self.current_coverage
                    ]
                )
                occupied = hs >= self.cell_min.value
                category = np.zeros(hs.shape[1:], dtype=int)
                for i in range(len(self.settings)):
                    category[occupied[i]] = i + 1
                degree = occupied.sum(axis=0)
                category[degree > 1] = len(self.settings) + 1
                if len(self.settings) > 1:
                    category[degree == len(self.settings)] = len(self.settings) + 2
                cmap = ListedColormap(
                    ["#fafafa", *COLORS[: len(self.settings)], "#999999", "#8064a2"]
                )
                ax.pcolormesh(
                    xe,
                    ye,
                    category.T,
                    cmap=cmap,
                    norm=BoundaryNorm(np.arange(cmap.N + 1) - 0.5, cmap.N),
                    rasterized=True,
                )
                self.map_stats.append(
                    {"counts": hs, "category": category, "x_edges": xe, "y_edges": ye}
                )
                ax.plot(ox, oy, color="#555555", ls="--", lw=1.2)
                try:
                    polygon = self._polygon()
                except ValueError:
                    polygon = None
                cx, cy = np.meshgrid(
                    (xe[:-1] + xe[1:]) / 2, (ye[:-1] + ye[1:]) / 2, indexing="ij"
                )
                cq = cx if xcol == "Q2" else cy
                cxb = cq / (cy**2 - self.proton_mass**2 + cq) if xcol == "Q2" else cx
                inside = (cq >= qlo) & (cq <= qhi) & (cxb >= xlo) & (cxb <= xhi)
                if polygon is not None:
                    inside &= polygon_mask(
                        pd.DataFrame({"xB": cxb.ravel(), "Q2": cq.ravel()}), polygon
                    ).reshape(cx.shape)
                fade = np.ones((*inside.T.shape, 4))
                fade[:, :, 3] = np.where(inside.T, 0, 0.65)
                ax.imshow(
                    fade,
                    extent=[*xlim, *ylim],
                    origin="lower",
                    aspect="auto",
                    interpolation="nearest",
                    zorder=2,
                )
                if polygon is not None:
                    line = self._boundary(polygon)
                    px = line[:, 1] if xcol == "Q2" else line[:, 0]
                    py = self._w(line[:, 0], line[:, 1]) if xcol == "Q2" else line[:, 1]
                    ax.plot(px, py, color="white", lw=3.8, zorder=4)
                    ax.plot(px, py, color="#151515", lw=1.7, zorder=5)
                if self.valid:
                    for edge in self.edges["Q2"][1:-1]:
                        (ax.axvline if xcol == "Q2" else ax.axhline)(
                            edge, color="#333333", lw=0.7, ls=":"
                        )
                    for iq, xe_bin in enumerate(self.edges["xB_by_Q2"]):
                        qq = np.linspace(*self.edges["Q2"][iq : iq + 2], 80)
                        for x in xe_bin[1:-1]:
                            ax.plot(
                                qq if xcol == "Q2" else np.full_like(qq, x),
                                self._w(x, qq) if xcol == "Q2" else qq,
                                color="#333333",
                                lw=0.7,
                                ls=":",
                            )
                for (i, _), x, y in zip(markers, vx, vy):
                    text = ax.text(
                        x,
                        y,
                        str(i + 1),
                        ha="center",
                        va="center",
                        fontsize=9,
                        zorder=10,
                        bbox=dict(
                            boxstyle="circle,pad=.28", fc="white", ec="#111111", lw=1.3
                        ),
                        clip_on=False,
                    )
                    text.set_gid(f"corner-{imap}-{i+1}")
                ax.set(
                    xlim=xlim,
                    ylim=ylim,
                    xlabel=LABELS[xcol],
                    ylabel=LABELS[ycol],
                    title=(
                        r"(a) $W$–$Q^2$ coverage"
                        if imap == 0
                        else r"(b) $x_B$–$Q^2$ coverage"
                    ),
                )
            for ax, key, tag in zip(
                axes[3:], ("tprime", "Q2", "xB"), ("(d)", "(e)", "(f)")
            ):
                if self.valid:
                    bins = np.linspace(*self.limits[key], 45)
                    for i, f in enumerate(self.selected_frames):
                        ax.stairs(
                            np.histogram(f[key], bins)[0],
                            bins,
                            color=COLORS[i],
                            lw=1.4,
                            linestyle="-" if i % 2 == 0 else "--",
                        )
                    edges = (
                        self.edges[key]
                        if key != "xB"
                        else np.unique(np.concatenate(self.edges["xB_by_Q2"]))
                    )
                    for edge in edges:
                        ax.axvline(edge, color="#444444", ls=":", lw=0.8)
                ax.set(
                    xlabel=LABELS[key],
                    ylabel="Events / display interval",
                    title=f"{tag} Selected data",
                )
                ax.grid(axis="y", color="#eeeeee", lw=0.5)
            ax = axes[2]
            if self.valid:
                matrix = np.minimum.reduce(self.occupancy_matrices)
                mesh = ax.pcolormesh(
                    np.degrees(self.edges["phi"]),
                    self.edges["tprime"],
                    matrix,
                    cmap="cividis",
                    shading="flat",
                    rasterized=True,
                )
                fig.colorbar(mesh, ax=ax, pad=0.02, fraction=0.045).set_label(
                    "Min. events / setting", fontsize=8
                )
                for it, ip in zip(*np.where(matrix < self.sparse_min.value)):
                    ax.plot(
                        np.degrees(self.edges["phi"][ip : ip + 2]).mean(),
                        self.edges["tprime"][it : it + 2].mean(),
                        marker="x",
                        color="#c14b25",
                        ms=3,
                        mew=0.7,
                    )
            ax.set(
                xlabel=r"$\phi\;[\mathrm{deg}]$",
                ylabel=LABELS["tprime"],
                title="(c) t′–φ: " + self.slice.label,
            )
            for ax in axes:
                ax.tick_params(labelsize=8, direction="out")
                ax.spines[["top", "right"]].set_visible(False)
            handles = []
            for i, (s, f) in enumerate(zip(self.settings, self.selected_frames)):
                eps = f.epsilon.median() if "epsilon" in f and len(f) else np.nan
                handles.append(
                    Patch(
                        facecolor=COLORS[i], label=f"{s}  (ε≈{eps:.3f}; N={len(f):,})"
                    )
                )
            if len(self.settings) > 2:
                handles.append(
                    Patch(facecolor="#999999", label="Some settings overlap")
                )
            if len(self.settings) > 1:
                handles.append(Patch(facecolor="#8064a2", label="All settings overlap"))
            handles.extend(
                [
                    Line2D([], [], color="#111111", label="Diamond boundary"),
                    Line2D([], [], color="#555555", ls="--", label="Outer limits"),
                    Line2D([], [], color="#333333", ls=":", label="Bin edges"),
                    Line2D(
                        [],
                        [],
                        color="#c14b25",
                        marker="x",
                        ls="",
                        label="Sparse in ≥1 setting",
                    ),
                ]
            )
            fig.legend(
                handles=handles,
                loc="upper center",
                bbox_to_anchor=(0.51, 1.0),
                ncol=3,
                frameon=False,
                fontsize=8,
            )
            fig.text(
                0.065,
                0.015,
                f"Map: ≥{self.cell_min.value} events / setting. Pale = outside selection (diamond ∩ limits). Counts unweighted; data-only edges.",
                fontsize=8,
            )
            if not self.valid:
                fig.text(
                    0.5,
                    0.49,
                    "INVALID SELECTION — adjust controls",
                    color="#a00",
                    ha="center",
                    weight="bold",
                )
            fig.canvas.draw()
            width, height = fig.get_size_inches() * 72
            geometry = []
            for name, ax in zip(("W", "xB"), axes[:2]):
                p = ax.get_position()
                geometry.append(
                    dict(
                        name=name,
                        left=p.x0 * width,
                        top=(1 - p.y1) * height,
                        width=p.width * width,
                        height=p.height * height,
                        xlim=list(ax.get_xlim()),
                        ylim=list(ax.get_ylim()),
                    )
                )
            stream = StringIO()
            fig.savefig(stream, format="svg")
            self.figure, self.axes = fig, axes
            self.canvas.geometry = geometry
            self.canvas.svg = stream.getvalue()

    def _save_configuration(self):
        vertices = self.vertices if all(v is not None for v in self.vertices) else None
        lines = [
            f"DIAMOND_XB_Q2_VERTICES = {vertices!r}",
            f"TPRIME_LIMITS = {self.t_limits.value!r}; Q2_LIMITS = {self.q_limits.value!r}; XB_LIMITS = {self.x_limits.value!r}",
            f"N_TPRIME = {self.n_t.value}; N_Q2 = {self.n_q.value}; N_XB_PER_Q2 = {self.n_x.value}; N_PHI = {self.n_phi.value}",
            f"TPRIME_EDGE_MODE = {self.t_mode.value!r}; Q2_EDGE_MODE = {self.q_mode.value!r}; XB_EDGE_MODE = {self.x_mode.value!r}",
            f"MANUAL_TPRIME_EDGES = {self.manual_edges.get('tprime')!r}",
            f"MANUAL_Q2_EDGES = {self.manual_edges.get('Q2')!r}",
            f"XB_EDGES_BY_Q2 = {self.manual_edges.get('xB') if self.x_mode.value=='manual' else None!r}",
        ]
        self.copy_code.value = "\n".join(lines)

    def export_figure(self, stem):
        self.require_ready()
        stem = Path(stem)
        stem.parent.mkdir(parents=True, exist_ok=True)
        paths = [Path(str(stem) + suffix) for suffix in (".pdf", ".svg")]
        with mpl.rc_context({"svg.fonttype": "path", "pdf.fonttype": 42}):
            for p in paths:
                self.figure.savefig(p, dpi=300)
        return paths

    def _export_clicked(self, _):
        try:
            self.export_message.value = "Saved: " + ", ".join(
                escape(str(p.resolve()))
                for p in self.export_figure(self.export_path.value)
            )
        except (RuntimeError, OSError) as exc:
            self.export_message.value = escape(str(exc))

    def show(self):
        display(self.widget)
