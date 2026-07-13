#!/usr/bin/env python3
"""Run SIMC/Geant4 wrapper commands from canonical simulation kinematics CSV.

This wrapper is analysis-side orchestration only. It does not modify SIMC or Geant4 internals.
It reads per-kin values and offsets from config/nps_simulation_kinematics.csv and executes
explicit user-provided command templates.
"""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List


REQUIRED_COLUMNS = [
    "new_kin",
    "kin_old",
    "ebeam_gev",
    "hms_p_gev",
    "hms_theta_deg",
    "nps_theta_deg",
    "nps_target_distance_cm",
    "simc_spec_e_p_offset",
    "simc_spec_e_theta_offset",
    "simc_hms_pointing_offset",
    "geant4_nps_x_offset_cm",
    "geant4_nps_y_offset_cm",
    "geant4_nps_z_offset_cm",
]


@dataclass(frozen=True)
class KinRow:
    values: Dict[str, str]

    @property
    def kin_old(self) -> str:
        return self.values["kin_old"].strip()

    @property
    def kin_tag(self) -> str:
        tag = self.kin_old
        if tag.lower().startswith("kinc_"):
            tag = tag[5:]
        return tag.lower()


def read_rows(sim_config: Path) -> List[KinRow]:
    if not sim_config.exists():
        raise FileNotFoundError(f"Simulation config not found: {sim_config}")

    with sim_config.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise RuntimeError("Simulation config CSV has no header")

        headers = {h.strip().lower(): h for h in reader.fieldnames}
        missing = [c for c in REQUIRED_COLUMNS if c not in headers]
        if missing:
            raise RuntimeError(
                "Simulation config CSV missing required columns: " + ", ".join(missing)
            )

        rows: List[KinRow] = []
        for raw in reader:
            normalized = {
                key: (raw.get(headers[key], "") or "").strip()
                for key in headers.keys()
            }
            kin_old = normalized.get("kin_old", "")
            if not kin_old:
                continue
            rows.append(KinRow(values=normalized))

    if not rows:
        raise RuntimeError("Simulation config CSV has no usable rows")
    return rows


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def run_command(command: str, dry_run: bool) -> None:
    print(f"[sim-wrapper] command: {command}")
    if dry_run:
        return
    subprocess.run(command, shell=True, check=True)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run SIMC/Geant4 wrappers from canonical simulation CSV"
    )
    parser.add_argument(
        "--sim-config",
        default="config/nps_simulation_kinematics.csv",
        help="Path to canonical simulation CSV",
    )
    parser.add_argument(
        "--simc-infile-dir",
        default="config/simc_infiles",
        help="Directory containing generated SIMC infiles",
    )
    parser.add_argument(
        "--simc-file-prefix",
        default="nps_excl_pi0_",
        help="SIMC infile prefix",
    )
    parser.add_argument(
        "--output-base",
        default="output",
        help="Base output directory for simulation artifacts",
    )
    parser.add_argument(
        "--kin",
        action="append",
        default=[],
        help="Filter Kin_old (repeatable)",
    )

    parser.add_argument(
        "--run-simc",
        action="store_true",
        help="Execute SIMC command template",
    )
    parser.add_argument(
        "--run-geant4",
        action="store_true",
        help="Execute Geant4 command template",
    )
    parser.add_argument(
        "--simc-cmd-template",
        default="",
        help="SIMC command template",
    )
    parser.add_argument(
        "--geant4-cmd-template",
        default="",
        help="Geant4 command template",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands and manifests only",
    )
    parser.add_argument(
        "--manifest",
        default="",
        help="Optional manifest CSV path (default: <output-base>/simulation_chain_manifest.csv)",
    )
    return parser.parse_args()


def required_placeholders(template: str, required: Iterable[str], label: str) -> None:
    missing = [name for name in required if "{" + name + "}" not in template]
    if missing:
        raise RuntimeError(
            f"{label} template missing required placeholders: {', '.join(missing)}"
        )


def main() -> int:
    args = parse_args()

    if not args.run_simc and not args.run_geant4:
        raise RuntimeError("Select at least one stage: --run-simc and/or --run-geant4")

    if args.run_simc and not args.simc_cmd_template.strip():
        raise RuntimeError("--run-simc requires --simc-cmd-template")
    if args.run_geant4 and not args.geant4_cmd_template.strip():
        raise RuntimeError("--run-geant4 requires --geant4-cmd-template")

    if args.run_simc:
        required_placeholders(
            args.simc_cmd_template,
            ["simc_infile", "simc_output_file", "kin_old"],
            "SIMC",
        )
    if args.run_geant4:
        required_placeholders(
            args.geant4_cmd_template,
            ["kin_old", "geant4_output_dir"],
            "Geant4",
        )

    repo_root = Path(__file__).resolve().parents[1]
    sim_config = (repo_root / args.sim_config).resolve()
    simc_infile_dir = (repo_root / args.simc_infile_dir).resolve()
    output_base = (repo_root / args.output_base).resolve()
    manifest_path = (
        (repo_root / args.manifest).resolve()
        if args.manifest
        else output_base / "simulation_chain_manifest.csv"
    )

    rows = read_rows(sim_config)
    kin_filter = {k.strip() for k in args.kin if k.strip()}

    selected = [r for r in rows if not kin_filter or r.kin_old in kin_filter]
    if not selected:
        raise RuntimeError("No matching kinematics after --kin filtering")

    manifest_rows: List[Dict[str, str]] = []

    for row in selected:
        kin_old = row.kin_old
        kin_tag = row.kin_tag
        simc_infile = simc_infile_dir / f"{args.simc_file_prefix}{kin_tag}.inp"
        if not simc_infile.exists():
            raise FileNotFoundError(
                f"Required SIMC infile missing for {kin_old}: {simc_infile}"
            )

        output_dir = output_base / kin_tag / "simulation"
        simc_output_dir = output_dir / "simc"
        geant4_output_dir = output_dir / "geant4"
        simc_output_file = simc_output_dir / f"simc_{kin_tag}.root"

        simc_output_dir.mkdir(parents=True, exist_ok=True)
        geant4_output_dir.mkdir(parents=True, exist_ok=True)

        fmt = dict(row.values)
        fmt.update(
            {
                "kin_old": kin_old,
                "kin_tag": kin_tag,
                "simc_infile": shlex.quote(str(simc_infile)),
                "output_dir": shlex.quote(str(output_dir)),
                "simc_output_dir": shlex.quote(str(simc_output_dir)),
                "simc_output_file": shlex.quote(str(simc_output_file)),
                "geant4_output_dir": shlex.quote(str(geant4_output_dir)),
            }
        )

        simc_cmd = ""
        geant4_cmd = ""
        if args.run_simc:
            simc_cmd = args.simc_cmd_template.format(**fmt)
            run_command(simc_cmd, args.dry_run)
        if args.run_geant4:
            geant4_cmd = args.geant4_cmd_template.format(**fmt)
            run_command(geant4_cmd, args.dry_run)

        manifest_rows.append(
            {
                "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                "kin_old": kin_old,
                "new_kin": row.values.get("new_kin", ""),
                "sim_config": str(sim_config),
                "simc_infile": str(simc_infile),
                "simc_command": simc_cmd,
                "geant4_command": geant4_cmd,
                "simc_output_file": str(simc_output_file),
                "geant4_output_dir": str(geant4_output_dir),
                "simc_spec_e_p_offset": row.values.get("simc_spec_e_p_offset", ""),
                "simc_spec_e_theta_offset": row.values.get("simc_spec_e_theta_offset", ""),
                "simc_hms_pointing_offset": row.values.get("simc_hms_pointing_offset", ""),
                "geant4_nps_x_offset_cm": row.values.get("geant4_nps_x_offset_cm", ""),
                "geant4_nps_y_offset_cm": row.values.get("geant4_nps_y_offset_cm", ""),
                "geant4_nps_z_offset_cm": row.values.get("geant4_nps_z_offset_cm", ""),
            }
        )

    ensure_parent(manifest_path)
    with manifest_path.open("w", newline="") as f:
        fieldnames = list(manifest_rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"[sim-wrapper] wrote manifest: {manifest_path}")
    print(f"[sim-wrapper] processed kinematics: {len(manifest_rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
