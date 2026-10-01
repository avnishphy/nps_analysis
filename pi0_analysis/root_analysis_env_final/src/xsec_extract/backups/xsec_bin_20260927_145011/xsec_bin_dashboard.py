"""Interactive, data-only bin selection dashboard for xsec_bin_helper.ipynb.

This is a diagnostic UI. The notebook keeps the final bin contract and SIMC
representative-coordinate calculation.
"""

from html import escape

import ipywidgets as widgets
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from IPython.display import display
from matplotlib.path import Path as PolygonPath


COLORS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9")
MODE_OPTIONS = ("uniform", "count_quantile", "yield_quantile", "manual")


def _quantile_edges(values, weights, limits, count, mode, manual=None):
    lo, hi = limits
    if mode == "manual":
        if manual is None:
            raise ValueError("Manual mode needs edges in the configuration cell")
        edges = np.asarray(manual, dtype=float)
    elif mode == "uniform":
        edges = np.linspace(lo, hi, count + 1)
    else:
        values = np.asarray(values, dtype=float)
        keep = np.isfinite(values) & (values >= lo) & (values <= hi)
        if mode == "yield_quantile":
            weights = np.asarray(weights, dtype=float)
            keep &= np.isfinite(weights) & (weights > 0)
        values = values[keep]
        if len(values) == 0:
            raise ValueError("No selected events for quantile edges")
        probabilities = np.arange(1, count) / count
        if mode == "count_quantile":
            interior = np.quantile(values, probabilities)
        elif mode == "yield_quantile":
            weights = weights[keep]
            order = np.argsort(values)
            values, weights = values[order], weights[order]
            cumulative = (np.cumsum(weights) - 0.5 * weights) / weights.sum()
            interior = np.interp(probabilities, cumulative, values)
        else:
            raise ValueError(f"Unknown edge mode {mode}")
        edges = np.r_[lo, interior, hi]
    if len(edges) < 2 or not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise ValueError(f"Invalid or collapsed {mode} edges: {edges}")
    if not np.allclose([edges[0], edges[-1]], limits):
        raise ValueError(f"Manual edges must span the current limits {limits}")
    return edges


def _shape_lines(edges, direction):
    shapes = []
    for edge in edges:
        if direction == "vertical":
            shapes.append(dict(type="line", x0=float(edge), x1=float(edge), y0=0, y1=1,
                               yref="paper", line=dict(color="#202020", width=1, dash="dash")))
        else:
            shapes.append(dict(type="line", y0=float(edge), y1=float(edge), x0=0, x1=1,
                               xref="paper", line=dict(color="#202020", width=1, dash="dash")))
    return shapes


class BinDashboard:
    """Linked FigureWidgets. Click either 2D view to place four diamond corners."""

    def __init__(self, settings, coverage_frames, proton_mass, limits, modes, counts,
                 vertices=None, manual_edges=None, on_change=None,
                 overlap_bins=(100, 90), min_events_per_cell=3, max_scatter=4000):
        if len(settings) > len(COLORS):
            raise ValueError("Add more setting colors to the dashboard")
        self.settings = list(settings)
        self.coverage_frames = coverage_frames
        self.proton_mass = proton_mass
        self.on_change = on_change
        self.manual_edges = manual_edges or {}
        self.overlap_bins = overlap_bins
        self.min_events_per_cell = min_events_per_cell
        self.max_scatter = max_scatter
        self.vertices = [None] * 4 if vertices is None else [tuple(map(float, v)) for v in vertices]
        self.last_click = None
        if len(self.vertices) != 4:
            raise ValueError("The dashboard diamond needs exactly four (xB, Q2) corners")
        self.selected_frames = []
        self.edges = {}
        self._busy = False

        q_values = np.concatenate([frame["Q2"].to_numpy() for frame in coverage_frames])
        x_values = np.concatenate([frame["xB"].to_numpy() for frame in coverage_frames])
        q_bounds = np.quantile(q_values[np.isfinite(q_values)], [0.001, 0.999])
        x_bounds = np.quantile(x_values[np.isfinite(x_values)], [0.001, 0.999])
        self.t_limits = widgets.FloatRangeSlider(description="t'", value=limits["tprime"],
            min=limits["tprime"][0], max=limits["tprime"][1], step=0.005,
            readout_format=".3f", continuous_update=False, layout=widgets.Layout(width="330px"))
        self.q_limits = widgets.FloatRangeSlider(description="Q²", value=limits["Q2"],
            min=min(float(q_bounds[0]), limits["Q2"][0]), max=max(float(q_bounds[1]), limits["Q2"][1]),
            step=0.01, readout_format=".2f", continuous_update=False,
            layout=widgets.Layout(width="330px"))
        self.x_limits = widgets.FloatRangeSlider(description="xB", value=limits["xB"],
            min=min(float(x_bounds[0]), limits["xB"][0]), max=max(float(x_bounds[1]), limits["xB"][1]),
            step=0.002, readout_format=".3f", continuous_update=False,
            layout=widgets.Layout(width="330px"))
        self.t_mode = widgets.Dropdown(description="t' mode", options=MODE_OPTIONS, value=modes["tprime"],
            layout=widgets.Layout(width="210px"))
        self.q_mode = widgets.Dropdown(description="Q² mode", options=MODE_OPTIONS, value=modes["Q2"],
            layout=widgets.Layout(width="210px"))
        self.x_mode = widgets.Dropdown(description="xB mode", options=MODE_OPTIONS, value=modes["xB"],
            layout=widgets.Layout(width="210px"))
        self.n_t = widgets.IntSlider(description="t' bins", value=counts["tprime"], min=1, max=8,
            layout=widgets.Layout(width="220px"), continuous_update=False)
        self.n_q = widgets.IntSlider(description="Q² bins", value=counts["Q2"], min=1, max=8,
            layout=widgets.Layout(width="220px"), continuous_update=False)
        self.n_x = widgets.IntSlider(description="xB bins", value=counts["xB"], min=1, max=8,
            layout=widgets.Layout(width="220px"), continuous_update=False)
        self.n_phi = widgets.IntSlider(description="φ bins", value=counts["phi"], min=4, max=36,
            layout=widgets.Layout(width="190px"), continuous_update=False)
        self.cell_min = widgets.IntSlider(description="Cell min", value=min_events_per_cell,
            min=1, max=20, layout=widgets.Layout(width="190px"), continuous_update=False)
        self.show_dots = widgets.Checkbox(value=False, description="Show event dots",
            layout=widgets.Layout(width="155px"))
        self.show_cell_hover = widgets.Checkbox(value=False, description="Cell counts on hover",
            layout=widgets.Layout(width="185px"))
        self.corner = widgets.Dropdown(description="Corner", options=[(str(i), i - 1) for i in range(1, 5)],
            value=0, layout=widgets.Layout(width="150px"))
        self.undo = widgets.Button(description="Undo corner", icon="undo", layout=widgets.Layout(width="125px"))
        self.reset = widgets.Button(description="Clear diamond", icon="eraser", layout=widgets.Layout(width="145px"))
        self.status = widgets.HTML()
        self.copy_code = widgets.Textarea(description="Save:", layout=widgets.Layout(width="100%", height="110px"))

        self.maps = [self._map_figure("Q2", "W", "W vs Q²"), self._map_figure("xB", "Q2", "Q² vs xB")]
        self.hists = [self._hist_figure("tprime", "t' [GeV²]"),
                      self._hist_figure("Q2", "Q² [GeV²]"),
                      self._hist_figure("xB", "xB")]
        self.occupancy = self._occupancy_figure()
        figures = [*self.maps, *self.hists, self.occupancy]
        grid = widgets.GridBox(figures, layout=widgets.Layout(
            grid_template_columns="repeat(3, 345px)", grid_gap="2px", width="1040px"))
        controls = widgets.VBox([
            widgets.GridBox([
                widgets.VBox([self.t_limits, self.t_mode, self.n_t]),
                widgets.VBox([self.q_limits, self.q_mode, self.n_q]),
                widgets.VBox([self.x_limits, self.x_mode, self.n_x]),
            ], layout=widgets.Layout(grid_template_columns="repeat(3, 345px)", width="1040px")),
            widgets.HBox([self.corner, self.undo, self.reset, self.n_phi, self.cell_min]),
            widgets.HBox([self.show_dots, self.show_cell_hover]),
        ])
        self.widget = widgets.VBox([controls, self.status, grid, self.copy_code])

        for control in (self.t_limits, self.q_limits, self.x_limits,
                        self.t_mode, self.q_mode, self.x_mode,
                        self.n_t, self.n_q, self.n_x, self.n_phi,
                        self.cell_min, self.show_dots, self.show_cell_hover):
            control.observe(self._control_changed, names="value")
        self.undo.on_click(self._undo)
        self.reset.on_click(self._reset)
        self.refresh()

    def _w(self, xb, q2):
        return np.sqrt(self.proton_mass**2 + np.asarray(q2) * (1.0 / np.asarray(xb) - 1.0))

    def _map_figure(self, x_col, y_col, title):
        x_all = np.concatenate([frame[x_col].to_numpy() for frame in self.coverage_frames])
        y_all = np.concatenate([frame[y_col].to_numpy() for frame in self.coverage_frames])
        x_low, x_high = np.nanquantile(x_all, [0.001, 0.999])
        y_low, y_high = np.nanquantile(y_all, [0.001, 0.999])
        x_pad, y_pad = 0.04 * (x_high - x_low), 0.04 * (y_high - y_low)
        x_edges = np.linspace(x_low - x_pad, x_high + x_pad, self.overlap_bins[0] + 1)
        y_edges = np.linspace(y_low - y_pad, y_high + y_pad, self.overlap_bins[1] + 1)
        colorscale = [[0, "#ffffff"], [0.2499, "#ffffff"],
            [0.25, "#1f77b4"], [0.4999, "#1f77b4"],
            [0.5, "#e66100"], [0.7499, "#e66100"],
            [0.75, "#7b1fa2"], [1, "#7b1fa2"]]
        xc, yc = 0.5 * (x_edges[:-1] + x_edges[1:]), 0.5 * (y_edges[:-1] + y_edges[1:])
        cell_hover = "<br>".join(
            f"{name}: %{{customdata[{i}]}}" for i, name in enumerate(self.settings))
        fig = go.FigureWidget()
        fig.add_heatmap(x=xc, y=yc, z=np.zeros((len(yc), len(xc))), zmin=0, zmax=3,
            colorscale=colorscale, showscale=False, opacity=1.0, hoverinfo="none")
        for setting, color in zip(self.settings, COLORS):
            fig.add_scattergl(x=[], y=[], mode="markers", name=setting,
                marker=dict(color=color, opacity=0.13, size=2), visible=False,
                hoverinfo="none")
        polygon_index = len(fig.data)
        fig.add_scatter(x=[], y=[], mode="lines", name="diamond boundary",
            line=dict(color="#111111", width=2), fill="toself",
            fillcolor="rgba(50,50,50,0.0)", hoverinfo="none")
        vertex_index = len(fig.data)
        fig.add_scatter(x=[], y=[], mode="markers+text", name="diamond corners",
            marker=dict(color="#fde725", size=18, line=dict(color="#111111", width=2)),
            textposition="middle center", textfont=dict(color="#111111", size=10),
            hoverinfo="none")
        # Heatmap cells receive clicks throughout the plotted range.
        for trace in fig.data:
            trace.on_click(lambda trace, points, selector, kind=x_col: self._plot_click(kind, points))
        fig.update_layout(title=dict(text=title, font=dict(size=13)), width=340, height=245,
            margin=dict(l=45, r=5, t=32, b=35), template="plotly_white", showlegend=False,
            dragmode="pan", uirevision="diamond-dashboard")
        fig.update_xaxes(title_text="Q² [GeV²]" if x_col == "Q2" else "xB", title_font_size=11)
        fig.update_yaxes(title_text="W [GeV]" if y_col == "W" else "Q² [GeV²]", title_font_size=11)
        fig._x_edges, fig._y_edges = x_edges, y_edges
        fig._cell_hover_template = "%{x:.3f}, %{y:.3f}<br>" + cell_hover + "<extra>selected events</extra>"
        fig._selection_trace_indices = list(range(1, 1 + len(self.settings)))
        fig._polygon_trace_index = polygon_index
        fig._vertex_trace_index = vertex_index
        return fig

    def _hist_figure(self, column, title):
        fig = go.FigureWidget()
        for setting, color in zip(self.settings, COLORS):
            fig.add_histogram(x=[], name=setting, marker_color=color, opacity=0.55,
                              nbinsx=55, histnorm="probability density")
        fig.update_layout(title=dict(text=title, font=dict(size=13)), width=340, height=245,
            margin=dict(l=45, r=5, t=32, b=35), barmode="overlay", template="plotly_white",
            showlegend=False, xaxis_title=title, yaxis_title="density")
        return fig

    def _occupancy_figure(self):
        fig = go.FigureWidget(go.Heatmap(z=[[0]], colorscale="Viridis", showscale=False))
        fig.update_layout(title=dict(text="Least events per setting in t'–φ", font=dict(size=13)),
            width=340, height=245, margin=dict(l=45, r=5, t=32, b=35),
            template="plotly_white", xaxis_title="φ [deg]", yaxis_title="t' [GeV²]")
        return fig

    def _plot_click(self, x_col, points):
        if len(points.xs) == 0 or len(points.ys) == 0:
            return
        x, y = float(points.xs[0]), float(points.ys[0])
        if x_col == "xB":
            xb, q2 = x, y
        else:
            q2 = x
            xb = q2 / (y*y - self.proton_mass**2 + q2)
        if not (np.isfinite(xb) and np.isfinite(q2) and xb > 0 and q2 > 0):
            self.status.value = "<b>Click inside the physical xB/Q² region.</b>"
            return
        index = self.corner.value
        self.vertices[index] = (xb, q2)
        self.last_click = f"Corner {index + 1} placed at xB={xb:.3f}, Q²={q2:.3f}"
        self.corner.value = (index + 1) % 4
        self.refresh()

    def _undo(self, _):
        index = (self.corner.value - 1) % 4
        self.vertices[index] = None
        self.last_click = f"Corner {index + 1} removed"
        self.corner.value = index
        self.refresh()

    def _reset(self, _):
        self.vertices = [None] * 4
        self.last_click = "Diamond cleared"
        self.corner.value = 0
        self.refresh()

    def _control_changed(self, _):
        self.refresh()

    def _polygon(self):
        if any(vertex is None for vertex in self.vertices):
            return None
        vertices = np.asarray(self.vertices, dtype=float)
        center = vertices.mean(axis=0)
        angle = np.arctan2(vertices[:, 1] - center[1], vertices[:, 0] - center[0])
        return vertices[np.argsort(angle)]

    def _polygon_mask(self, frame, polygon):
        if polygon is None:
            return np.ones(len(frame), dtype=bool)
        path = PolygonPath(np.vstack([polygon, polygon[0]]), closed=True)
        return path.contains_points(frame[["xB", "Q2"]].to_numpy(), radius=1e-9)

    def refresh(self):
        if self._busy:
            return
        self._busy = True
        try:
            t_limits = tuple(map(float, self.t_limits.value))
            q_limits = tuple(map(float, self.q_limits.value))
            x_limits = tuple(map(float, self.x_limits.value))
            polygon = self._polygon()
            selected = []
            for frame in self.coverage_frames:
                mask = (frame["tprime"].between(*t_limits, inclusive="both") &
                        frame["Q2"].between(*q_limits, inclusive="both") &
                        frame["xB"].between(*x_limits, inclusive="both") &
                        self._polygon_mask(frame, polygon))
                selected.append(frame.loc[mask].copy())
            pooled = pd.concat(selected, ignore_index=True)
            if pooled.empty:
                raise ValueError("The current cuts select no data")
            t_edges = _quantile_edges(pooled["tprime"], pooled["data_weight"],
                t_limits, self.n_t.value, self.t_mode.value, self.manual_edges.get("tprime"))
            q_edges = _quantile_edges(pooled["Q2"], pooled["data_weight"],
                q_limits, self.n_q.value, self.q_mode.value, self.manual_edges.get("Q2"))
            x_edges_by_q = []
            for iq in range(len(q_edges) - 1):
                upper = pooled["Q2"] <= q_edges[iq + 1] if iq == len(q_edges) - 2 else pooled["Q2"] < q_edges[iq + 1]
                subset = pooled.loc[(pooled["Q2"] >= q_edges[iq]) & upper]
                manual = self.manual_edges.get("xB")
                manual = None if manual is None else manual[iq]
                x_edges_by_q.append(_quantile_edges(subset["xB"], subset["data_weight"],
                    x_limits, self.n_x.value, self.x_mode.value, manual))
            phi_edges = np.linspace(0, 2 * np.pi, self.n_phi.value + 1)
        except ValueError as exc:
            self.status.value = f"<b>Selection needs adjustment:</b> {escape(str(exc))}"
            self._busy = False
            return

        self.selected_frames = selected
        self.edges = {"tprime": t_edges, "Q2": q_edges, "xB_by_Q2": x_edges_by_q, "phi": phi_edges}
        shared_cells = []
        for fig, x_col, y_col in zip(self.maps, ("Q2", "xB"), ("W", "Q2")):
            cell_counts = [np.histogram2d(frame[x_col], frame[y_col],
                bins=(fig._x_edges, fig._y_edges))[0].astype(int) for frame in selected]
            occupied = np.stack(cell_counts) >= self.cell_min.value
            degree = occupied.sum(axis=0)
            if len(self.settings) == 1:
                category = (degree > 0).astype(int)
            elif len(self.settings) == 2:
                category = occupied[0].astype(int) + 2 * occupied[1].astype(int)
            else:
                category = np.where(degree == 0, 0,
                    np.where(degree == len(self.settings), 3,
                             np.where(degree == 1, 1, 2)))
            shared_cells.append(int(np.count_nonzero(degree == len(self.settings))))
            with fig.batch_update():
                fig.data[0].z = category.T
                fig.data[0].customdata = np.stack(cell_counts, axis=-1).transpose(1, 0, 2)
                fig.data[0].hovertemplate = fig._cell_hover_template if self.show_cell_hover.value else None
                fig.data[0].hoverinfo = "all" if self.show_cell_hover.value else "none"
                for index, frame in enumerate(selected):
                    if len(frame) > self.max_scatter:
                        frame = frame.sample(self.max_scatter, random_state=314159)
                    trace = fig.data[fig._selection_trace_indices[index]]
                    trace.x = frame[x_col].to_numpy()
                    trace.y = frame[y_col].to_numpy()
                    trace.visible = self.show_dots.value
                # Preserve click numbers on markers; sorting is only for the polygon path.
                markers = [(i, v) for i, v in enumerate(self.vertices) if v is not None]
                vertex_trace = fig.data[fig._vertex_trace_index]
                vertex_xb = np.asarray([v[0] for _, v in markers])
                vertex_q2 = np.asarray([v[1] for _, v in markers])
                vertex_trace.x = vertex_q2 if x_col == "Q2" else vertex_xb
                vertex_trace.y = self._w(vertex_xb, vertex_q2) if x_col == "Q2" else vertex_q2
                vertex_trace.text = [str(i + 1) for i, _ in markers]
                boundary = polygon if polygon is not None else np.asarray([v for _, v in markers])
                if polygon is not None:
                    boundary = np.vstack([boundary, boundary[0]])
                polygon_trace = fig.data[fig._polygon_trace_index]
                polygon_trace.x = boundary[:, 1] if len(boundary) else []
                polygon_trace.y = (self._w(boundary[:, 0], boundary[:, 1]) if x_col == "Q2" else boundary[:, 1]) if len(boundary) else []
                if x_col != "Q2" and len(boundary):
                    polygon_trace.x = boundary[:, 0]
                polygon_trace.fill = "toself" if polygon is not None else None
                if x_col == "Q2":
                    fig.layout.shapes = tuple(_shape_lines(q_edges[1:-1], "vertical"))
                else:
                    shapes = _shape_lines(q_edges[1:-1], "horizontal")
                    for iq, edges in enumerate(x_edges_by_q):
                        for edge in edges[1:-1]:
                            shapes.append(dict(type="line", x0=float(edge), x1=float(edge),
                                y0=float(q_edges[iq]), y1=float(q_edges[iq + 1]),
                                line=dict(color="#202020", width=1, dash="dash")))
                    fig.layout.shapes = tuple(shapes)
        for fig, col, edges in zip(self.hists, ("tprime", "Q2", "xB"),
                                   (t_edges, q_edges, sorted(set(np.concatenate(x_edges_by_q))))):
            with fig.batch_update():
                for trace, frame in zip(fig.data, selected):
                    trace.x = frame[col].to_numpy()
                fig.layout.shapes = tuple(_shape_lines(edges, "vertical"))
        matrices = [np.histogram2d(frame["tprime"], frame["phi_wrapped"],
                    bins=(t_edges, phi_edges))[0] for frame in selected]
        matrix = np.minimum.reduce(matrices)
        with self.occupancy.batch_update():
            self.occupancy.data[0].z = matrix
            self.occupancy.data[0].x = np.degrees((phi_edges[:-1] + phi_edges[1:]) / 2)
            self.occupancy.data[0].y = (t_edges[:-1] + t_edges[1:]) / 2

        counts = ", ".join(
            f"<span style='color:{COLORS[i]}'><b>{escape(name)}: {len(frame):,}</b></span>"
            for i, (name, frame) in enumerate(zip(self.settings, selected)))
        cut_status = "active" if polygon is not None else f"{sum(v is not None for v in self.vertices)}/4 corners placed"
        if len(self.settings) == 1:
            coverage_key = "<b style='color:#1f77b4'>Blue = selected coverage.</b>"
        elif len(self.settings) == 2:
            coverage_key = ("<b style='color:#1f77b4'>Blue = first only</b>; "
                "<b style='color:#e66100'>orange = second only</b>; "
                "<b style='color:#7b1fa2'>purple = both</b>.")
        else:
            coverage_key = "Blue = one setting; orange = some; purple = all."
        overlap_text = (f" &nbsp; | &nbsp; <b>Shared cells:</b> W–Q² {shared_cells[0]}, "
                        f"xB–Q² {shared_cells[1]}" if len(self.settings) > 1 else "")
        self.status.value = (f"<b>Selected events:</b> {counts} &nbsp; | &nbsp; "
                             f"<b>Diamond:</b> {cut_status} &nbsp; | &nbsp; "
                             f"<b>Smallest common t'–φ cell:</b> {int(matrix.min())} events"
                             + overlap_text + ". " + coverage_key
                             + (f" &nbsp; | &nbsp; <b>{escape(self.last_click)}</b>" if self.last_click else ""))
        vertices_text = "None" if polygon is None else repr([tuple(map(float, v)) for v in polygon])
        self.copy_code.value = (f"DIAMOND_XB_Q2_VERTICES = {vertices_text}\n"
            f"TPRIME_LIMITS = {t_limits}; Q2_LIMITS = {q_limits}; XB_LIMITS = {x_limits}\n"
            f"N_TPRIME = {self.n_t.value}; N_Q2 = {self.n_q.value}; N_XB_PER_Q2 = {self.n_x.value}; N_PHI = {self.n_phi.value}\n"
            f"TPRIME_EDGE_MODE = {self.t_mode.value!r}; Q2_EDGE_MODE = {self.q_mode.value!r}; XB_EDGE_MODE = {self.x_mode.value!r}")
        if self.on_change is not None:
            self.on_change(self)
        self._busy = False

    def show(self):
        display(self.widget)
