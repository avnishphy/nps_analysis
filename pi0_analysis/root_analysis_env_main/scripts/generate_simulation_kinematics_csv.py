#!/usr/bin/env python3
"""Generate canonical simulation kinematics table from master run config.

This tool creates config/nps_simulation_kinematics.csv with one row per Kin_old.
Offsets are initialized to zero by design for this phase.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional


@dataclass
class KinRow:
    new_kin: str
    kin_old: str
    source_run: int
    ebeam_gev: float
    hms_p_gev: float
    hms_theta_deg: float
    nps_theta_deg: float
    nps_target_distance_cm: float


def normalize_header(name: str) -> str:
    return name.strip().lower()


def parse_float(value: str) -> Optional[float]:
    try:
        return float(value.strip())
    except Exception:
        return None


def parse_int(value: str) -> Optional[int]:
    try:
        return int(value.strip())
    except Exception:
        return None


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate simulation kinematics CSV")
    parser.add_argument("--config", required=True, help="Path to nps_dvcs_all_kins_main.csv")
    parser.add_argument("--output", required=True, help="Path to output simulation kinematics CSV")
    parser.add_argument(
        "--types",
        default="production,Production",
        help="Comma-separated Type values to include (default: production,Production)",
    )
    parser.add_argument(
        "--selection-rule",
        default="per Kin_old choose minimum run_number among allowed Type rows",
        help="Human-readable provenance selection rule",
    )
    args = parser.parse_args()

    cfg_path = Path(args.config)
    out_path = Path(args.output)

    if not cfg_path.exists():
        raise FileNotFoundError(f"Config CSV not found: {cfg_path}")

    allowed_types = {t.strip().lower() for t in args.types.split(",") if t.strip()}
    if not allowed_types:
        raise ValueError("Allowed type set is empty; provide at least one type")

    with cfg_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise ValueError("Config CSV has no header")

        # Build normalized-name map for robust access against whitespace-padded headers.
        original_by_norm: Dict[str, str] = {
            normalize_header(name): name for name in reader.fieldnames
        }

        required_norm = {
            "new_kin",
            "kin_old",
            "run_number",
            "type",
            "beamenergy",
            "hms_p",
            "hms_thet",
            "nps_thet",
            "nps_target_distance",
        }
        missing = [key for key in sorted(required_norm) if key not in original_by_norm]
        if missing:
            raise ValueError(f"Missing required config columns: {', '.join(missing)}")

        by_kin: Dict[str, KinRow] = {}

        for row in reader:
            def get(norm_key: str) -> str:
                return (row.get(original_by_norm[norm_key], "") or "").strip()

            kin_old = get("kin_old")
            if not kin_old:
                continue

            type_value = get("type").lower()
            if type_value not in allowed_types:
                continue

            run = parse_int(get("run_number"))
            ebeam = parse_float(get("beamenergy"))
            hms_p = parse_float(get("hms_p"))
            hms_thet = parse_float(get("hms_thet"))
            nps_thet = parse_float(get("nps_thet"))
            nps_dist = parse_float(get("nps_target_distance"))
            new_kin = get("new_kin")

            if (
                run is None
                or ebeam is None
                or hms_p is None
                or hms_thet is None
                or nps_thet is None
                or nps_dist is None
            ):
                continue

            candidate = KinRow(
                new_kin=new_kin,
                kin_old=kin_old,
                source_run=run,
                ebeam_gev=ebeam,
                hms_p_gev=hms_p,
                hms_theta_deg=hms_thet,
                nps_theta_deg=nps_thet,
                nps_target_distance_cm=nps_dist,
            )

            current = by_kin.get(kin_old)
            if current is None or candidate.source_run < current.source_run:
                by_kin[kin_old] = candidate

    if not by_kin:
        raise RuntimeError("No rows selected for simulation kinematics table")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    output_columns: List[str] = [
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
        "source_run",
        "source_simc_infile",
        "selection_rule",
        "last_generated_utc",
    ]

    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=output_columns)
        writer.writeheader()

        for kin_old in sorted(by_kin.keys()):
            row = by_kin[kin_old]
            writer.writerow(
                {
                    "new_kin": row.new_kin,
                    "kin_old": row.kin_old,
                    "ebeam_gev": f"{row.ebeam_gev:.6f}",
                    "hms_p_gev": f"{row.hms_p_gev:.6f}",
                    "hms_theta_deg": f"{row.hms_theta_deg:.6f}",
                    "nps_theta_deg": f"{row.nps_theta_deg:.6f}",
                    "nps_target_distance_cm": f"{row.nps_target_distance_cm:.6f}",
                    "simc_spec_e_p_offset": "0.0",
                    "simc_spec_e_theta_offset": "0.0",
                    "simc_hms_pointing_offset": "0.0",
                    "geant4_nps_x_offset_cm": "0.0",
                    "geant4_nps_y_offset_cm": "0.0",
                    "geant4_nps_z_offset_cm": "0.0",
                    "source_run": str(row.source_run),
                    "source_simc_infile": "",
                    "selection_rule": args.selection_rule,
                    "last_generated_utc": now_utc,
                }
            )

    print(f"Wrote {len(by_kin)} kinematic rows to {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
