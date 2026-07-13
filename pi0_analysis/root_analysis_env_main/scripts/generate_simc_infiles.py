#!/usr/bin/env python3
"""Generate per-kin SIMC input files from simulation kinematics CSV.

This script updates key parameters in a template infile for each Kin_old row:
- Ebeam (MeV)
- spec%e%P (MeV/c)
- spec%e%theta (deg)
- selected spect_offset fields from offset columns

Current phase policy: offset columns default to zero in CSV and can be edited later.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List


def kin_old_to_tag(kin_old: str) -> str:
    tag = kin_old.strip()
    if tag.lower().startswith("kinc_"):
        tag = tag[5:]
    return tag.lower()


def replace_assignment(lines: List[str], key: str, value: str) -> bool:
    # Preserve existing inline comment and spacing as much as possible.
    pattern = re.compile(rf'^(\s*{re.escape(key)}\s*=\s*)([^;\n]*)(.*)$')
    for i, line in enumerate(lines):
        m = pattern.match(line)
        if not m:
            continue
        prefix, _old, suffix = m.groups()
        replacement = f"{prefix}{value}{suffix}"
        if line.endswith("\n") and not replacement.endswith("\n"):
            replacement += "\n"
        lines[i] = replacement
        return True
    return False


def ensure_present(lines: List[str], key: str, value: str) -> None:
    if not replace_assignment(lines, key, value):
        raise RuntimeError(f"Template missing required key: {key}")


def main() -> int:
    parser = argparse.ArgumentParser(description="Generate SIMC infiles from simulation kinematics CSV")
    parser.add_argument("--sim-config", required=True, help="Path to config/nps_simulation_kinematics.csv")
    parser.add_argument("--template", required=True, help="Path to SIMC template infile")
    parser.add_argument("--out-dir", required=True, help="Output directory for generated infiles")
    parser.add_argument("--file-prefix", default="nps_excl_pi0_", help="Generated filename prefix")
    parser.add_argument("--kin", action="append", default=[], help="Optional Kin_old filter (repeatable)")
    args = parser.parse_args()

    sim_cfg = Path(args.sim_config)
    template = Path(args.template)
    out_dir = Path(args.out_dir)

    if not sim_cfg.exists():
        raise FileNotFoundError(f"Simulation config not found: {sim_cfg}")
    if not template.exists():
        raise FileNotFoundError(f"Template infile not found: {template}")

    requested_kins = {k.strip() for k in args.kin if k.strip()}

    with sim_cfg.open("r", newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None:
            raise RuntimeError("Simulation config CSV has no header")

        norm_to_header: Dict[str, str] = {
            h.strip().lower(): h for h in reader.fieldnames
        }

        required = [
            "kin_old",
            "ebeam_gev",
            "hms_p_gev",
            "hms_theta_deg",
            "simc_spec_e_p_offset",
            "simc_spec_e_theta_offset",
            "simc_hms_pointing_offset",
            "geant4_nps_x_offset_cm",
            "geant4_nps_y_offset_cm",
            "geant4_nps_z_offset_cm",
        ]
        missing = [k for k in required if k not in norm_to_header]
        if missing:
            raise RuntimeError(
                "Simulation config CSV missing required columns: " + ", ".join(missing)
            )

        rows = list(reader)

    if not rows:
        raise RuntimeError("Simulation config CSV has no rows")

    out_dir.mkdir(parents=True, exist_ok=True)

    template_lines = template.read_text().splitlines(keepends=True)

    generated = 0
    for row in rows:
        def get(norm_key: str) -> str:
            return (row.get(norm_to_header[norm_key], "") or "").strip()

        kin_old = get("kin_old")
        if not kin_old:
            continue
        if requested_kins and kin_old not in requested_kins:
            continue

        ebeam_gev = float(get("ebeam_gev") or "0")
        hms_p_gev = float(get("hms_p_gev") or "0")
        hms_theta_deg = float(get("hms_theta_deg") or "0")

        # Offsets remain explicit and user-editable in simulation CSV.
        p_offset_gev = float(get("simc_spec_e_p_offset") or "0")
        theta_offset_deg = float(get("simc_spec_e_theta_offset") or "0")
        hms_pointing_offset = float(get("simc_hms_pointing_offset") or "0")
        geant4_x = float(get("geant4_nps_x_offset_cm") or "0")
        geant4_y = float(get("geant4_nps_y_offset_cm") or "0")
        geant4_z = float(get("geant4_nps_z_offset_cm") or "0")

        ebeam_mev = 1000.0 * ebeam_gev
        spec_p_mev = 1000.0 * (hms_p_gev + p_offset_gev)
        spec_theta = hms_theta_deg + theta_offset_deg

        lines = list(template_lines)
        ensure_present(lines, "Ebeam", f"{ebeam_mev:.1f}")
        ensure_present(lines, "spec%e%P", f"{spec_p_mev:.2f}")
        ensure_present(lines, "spec%e%theta", f"{spec_theta:.2f}")

        # Keep these offsets explicit in generated SIMC infiles.
        ensure_present(lines, "spec%e%offset%xptar", f"{hms_pointing_offset:.6f}")
        ensure_present(lines, "spec%p%offset%x", f"{geant4_x:.6f}")
        ensure_present(lines, "spec%p%offset%y", f"{geant4_y:.6f}")
        ensure_present(lines, "spec%p%offset%z", f"{geant4_z:.6f}")

        tag = kin_old_to_tag(kin_old)
        out_file = out_dir / f"{args.file_prefix}{tag}.inp"
        out_file.write_text("".join(lines))
        generated += 1

    if generated == 0:
        raise RuntimeError("No SIMC infiles generated; check --kin filters and input CSV rows")

    print(f"Generated {generated} SIMC infile(s) in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
