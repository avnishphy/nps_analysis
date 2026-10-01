#!/usr/bin/env python3
"""Independent NumPy check of Eq. 5.30-5.31 and data-only point covariance."""
import argparse
import csv
from pathlib import Path

import numpy as np


def read(path):
    with path.open() as source:
        return list(csv.DictReader(source))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("--target-contam", type=float, default=.584)
    parser.add_argument("--target-contam-err", type=float, default=.014)
    args = parser.parse_args()
    directory = args.directory
    points = read(directory / "experimental_points.csv")
    rows = read(directory / "migration_reco_rows.csv")
    parameters = read(directory / "migration_parameters.csv")
    nr, np_ = len(rows), len(parameters)
    assert nr == len(points)
    design = np.zeros((nr, np_))
    for entry in read(directory / "migration_design.csv"):
        design[int(entry["reco_row"]), int(entry["parameter_index"])] = float(entry["response"])
    x = np.array([float(entry["value"]) for entry in parameters])
    data = np.array([float(entry["data"]) for entry in rows])
    data_variance = np.array([float(entry["data_variance"]) for entry in rows])
    fit_index = np.array([int(entry["fit_index"]) for entry in rows])
    selected = np.flatnonzero(fit_index >= 0)
    selected = selected[np.argsort(fit_index[selected])]
    variance = np.array([float(rows[r]["variance_used"]) for r in selected])
    finite_mc = not np.allclose(variance, data_variance[selected], rtol=1e-12, atol=0)
    # Independent normal-equation calculation, used only for validation.
    information = (design[selected].T / variance) @ design[selected]
    covariance = np.linalg.inv(information)
    gain = np.zeros((np_, nr))
    gain[:, selected] = covariance @ (design[selected].T / variance)
    active = {int(entry["truth_block"]): int(entry["active_block_index"])
              for entry in parameters}
    phi_cells, phi_covariance, phi_epsilon_weight = {}, {}, {}
    reco_moments = {name: np.zeros((nr, np_)) for name in ("q2", "xb", "tprime")}
    for entry in read(directory / "truth_phi_response_cells.csv"):
        key = (int(entry["reco_row"]), int(entry["truth_cell"]))
        phi_cells[key] = np.array([float(entry["basis_" + n])
                                   for n in ("U", "LT", "TT")])
        names = ("U", "LT", "TT")
        outer = np.zeros((4, 4))
        outer[:3, :3] = [[float(entry[f"cov_{i}_{j}"])
                         for j in names] for i in names]
        outer[:3, 3] = [float(entry[f"cov_{i}_epsilon_weight"]) for i in names]
        outer[3, :3] = outer[:3, 3]
        outer[3, 3] = float(entry["cov_epsilon_weight_epsilon_weight"])
        phi_covariance[key] = outer
        phi_epsilon_weight[key] = float(entry["epsilon_weight"])
        block = int(entry["truth_block"])
        if block in active:
            a = active[block]
            for name, matrix in reco_moments.items():
                matrix[key[0], 3*a:3*a+3] += [float(entry[f"reco_{name}_{n}"])
                                                for n in names]
    nphi = 1 + max(int(point["reco_phi_bin"]) for point in points)
    reconstructed = np.zeros_like(design)
    for (r, v), basis in phi_cells.items():
        block = v // nphi
        if block in active:
            reconstructed[r, 3 * active[block]:3 * active[block] + 3] += basis
    assert np.allclose(reconstructed, design, rtol=1e-10, atol=1e-20)
    prediction = design @ x
    gradient = np.zeros((nr, nr))
    parameter_gradient = np.zeros((nr, np_))
    point_f, point_m, point_n = np.zeros(nr), np.zeros(nr), np.zeros(nr)
    valid = []
    for entry in points:
        r = int(entry["reco_row"])
        if entry["status"] != "ok":
            continue
        v = int(entry["truth_cell"])
        block = v // nphi
        assert (r, v) in phi_cells
        d = np.zeros(np_)
        d[3 * active[block]:3 * active[block] + 3] = phi_cells[r, v]
        m = d @ x
        f = float(entry["sigma_fit_reference"])
        n = data[r] - (prediction[r] - m)
        assert np.isclose(m, float(entry["diagonal_contribution"]), rtol=1e-10)
        assert np.isclose(n, float(entry["subtracted_yield"]), rtol=1e-10)
        assert np.isclose(n / m, float(entry["correction"]), rtol=1e-10)
        assert np.isclose(f * n / m, float(entry["sigma_exp"]), rtol=1e-10)
        for name, matrix in reco_moments.items():
            expected_mean = (matrix[r] @ x) / prediction[r]
            assert np.isclose(expected_mean, float(entry[f"predicted_reco_{name}"]),
                              rtol=1e-10)
        epsilon = float(entry["epsilon_reference"])
        phi = float(entry["phi_reference"])
        total_weight = sum(2*np.pi*basis[0] for (t, cell), basis in phi_cells.items()
                           if cell == v)
        total_epsilon_weight = sum(value for (t, cell), value in phi_epsilon_weight.items()
                                   if cell == v)
        assert np.isclose(total_epsilon_weight/total_weight, epsilon, rtol=1e-10)
        g = np.zeros(np_)
        g[3 * active[block]:3 * active[block] + 3] = (
            np.array([1, np.sqrt(2 * epsilon * (1 + epsilon)) * np.cos(phi),
                      epsilon * np.cos(2 * phi)]) / (2 * np.pi))
        assert np.isclose(g @ x, f, rtol=1e-10)
        point_f[r], point_m[r], point_n[r] = f, m, n
        parameter_gradient[r] = (n / m) * g - (f / m) * (design[r] - d) - (
            f * n / m**2) * d
        e = np.zeros(nr)
        e[r] = 1
        # Different quotient-rule form from the C++ implementation:
        # E=f*N/m, N=y_r-(D_r-d)Xhat, Xhat=K y.
        derivative_n = e - (design[r] - d) @ gain
        derivative_m = d @ gain
        derivative_f = g @ gain
        gradient[r] = (f / m) * derivative_n + (n / m) * derivative_f - (
            f * n / m**2) * derivative_m
        valid.append(r)
    assert valid
    expected = (gradient * data_variance) @ gradient.T
    exported_cov = np.full((nr, nr), np.nan)
    exported_corr = np.full((nr, nr), np.nan)
    exported_target = np.full((nr, nr), np.nan)
    for entry in read(directory / "experimental_point_covariance.csv"):
        exported_cov[int(entry["reco_row_i"]), int(entry["reco_row_j"])] = float(
            entry["statistical_covariance"])
        exported_target[int(entry["reco_row_i"]), int(entry["reco_row_j"])] = float(
            entry["target_correlated_covariance"])
    for entry in read(directory / "experimental_point_correlation.csv"):
        exported_corr[int(entry["reco_row_i"]), int(entry["reco_row_j"])] = float(
            entry["statistical_correlation"])
    checks = valid if not finite_mc else sorted(set(valid[:3] + valid[len(valid)//2:len(valid)//2+3] + valid[-3:]))
    if finite_mc:
        h_c = parameter_gradient[checks] @ covariance
        mc_cov = np.zeros((len(checks), len(checks)))
        truth_weight = {}
        for (t, v), basis in phi_cells.items():
            truth_weight[v] = truth_weight.get(v, 0.) + 2 * np.pi * basis[0]
        for (t, v), outer in phi_covariance.items():
            block = v // nphi
            if block not in active:
                continue
            a = active[block]
            term = np.zeros((len(checks), 4))
            fi = fit_index[t]
            if fi >= 0:
                residual = data[t] - prediction[t]
                term[:, :3] += h_c[:, 3*a:3*a+3] * (residual / variance[fi])
                term[:, :3] -= np.outer((h_c @ design[t]) / variance[fi], x[3*a:3*a+3])
            for i, r in enumerate(checks):
                if t == r:
                    direct = (-point_f[r] * point_n[r] / point_m[r]**2
                              if v == r else -point_f[r] / point_m[r])
                    term[i, :3] += direct * x[3*a:3*a+3]
                if v == r:
                    eps = float(points[r]["epsilon_reference"])
                    phi = float(points[r]["phi_reference"])
                    k = np.sqrt(2 * eps * (1 + eps))
                    df_deps = (((1 + 2*eps) / k) * np.cos(phi) * x[3*a+1]
                               + np.cos(2*phi) * x[3*a+2]) / (2*np.pi)
                    ge = (point_n[r] / point_m[r]) * df_deps / truth_weight[v]
                    term[i, 0] -= ge * 2*np.pi*eps
                    term[i, 3] += ge
            mc_cov += term @ outer @ term.T
        expected[np.ix_(checks, checks)] += mc_cov
    sub = np.ix_(checks, checks)
    error = np.linalg.norm(expected[sub] - exported_cov[sub]) / np.linalg.norm(expected[sub])
    correlation = expected[sub] / np.sqrt(np.outer(np.diag(expected)[checks],
                                                    np.diag(expected)[checks]))
    corr_error = np.linalg.norm(correlation - exported_corr[sub]) / np.linalg.norm(correlation)
    assert error < 1e-7, error
    assert corr_error < 1e-7, corr_error
    point_values = np.array([float(points[r]["sigma_exp"]) for r in checks])
    target = np.outer(point_values, point_values) * (
        args.target_contam_err / args.target_contam) ** 2
    assert np.allclose(target, exported_target[sub], rtol=1e-10, atol=1e-30)
    print(f"PASS mode={'finite-mc' if finite_mc else 'data'} points={len(valid)}/{nr} "
          f"covariance_checks={len(checks)} covariance_relative_error={error:.3g} "
          f"correlation_relative_error={corr_error:.3g}")


if __name__ == "__main__":
    main()
