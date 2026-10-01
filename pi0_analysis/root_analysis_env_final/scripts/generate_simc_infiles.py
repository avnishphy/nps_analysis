#!/usr/bin/env python3
"""Generate canonical simulation kinematics and per-kinematic SIMC infiles."""
# python3 scripts/generate_simc_infiles.py --config config/nps_dvcs_all_kins_main.csv --sim-config config/nps_simulation_kinematics.csv --template /u/group/nps/singhav/simc_gfortran_updated/infiles/nps_excl_pi0_x60_4b.inp
from __future__ import annotations

import argparse
import configparser
import csv
import math
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence


REQUIRED_COLUMNS = [
    "kin_old", "ebeam_gev", "hms_p_gev", "hms_theta_deg", "nps_theta_deg",
    "nps_target_distance_cm", "simc_spec_p_p_gev",
    "simc_ngen", "simc_target_x_offset_cm", "simc_target_y_offset_cm",
    "simc_target_z_offset_cm",
    "simc_spec_e_p_offset_gev", "simc_spec_e_theta_offset_deg",
    "simc_hms_xptar_offset_mrad", "simc_spec_p_x_offset_cm",
    "simc_spec_p_y_offset_cm", "simc_spec_p_z_offset_cm",
]


@dataclass(frozen=True)
class KinRow:
    new_kin: str
    kin_old: str
    target: str
    source_run: int
    ebeam_gev: float
    hms_p_gev: float
    hms_theta_deg: float
    nps_theta_deg: float
    nps_target_distance_cm: float


OUTPUT_COLUMNS: List[str] = [
    "new_kin", "kin_old", "target", "ebeam_gev", "hms_p_gev",
    "hms_theta_deg", "nps_theta_deg", "nps_target_distance_cm",
    "simc_spec_p_p_gev", "simc_ngen", "simc_target_x_offset_cm",
    "simc_target_y_offset_cm", "simc_target_z_offset_cm",
    "simc_spec_e_p_offset_gev", "simc_spec_e_theta_offset_deg",
    "simc_hms_xptar_offset_mrad", "simc_spec_p_x_offset_cm",
    "simc_spec_p_y_offset_cm", "simc_spec_p_z_offset_cm",
    "geant4_nps_x_offset_cm", "geant4_nps_y_offset_cm",
    "geant4_nps_z_offset_cm", "source_run", "source_simc_infile",
    "selection_rule", "last_generated_utc",
]

PRESERVED_COLUMNS: Mapping[str, Sequence[str]] = {
    "simc_spec_p_p_gev": ("simc_spec_p_p_gev",),
    "simc_ngen": ("simc_ngen",),
    "simc_target_x_offset_cm": ("simc_target_x_offset_cm",),
    "simc_target_y_offset_cm": ("simc_target_y_offset_cm",),
    "simc_target_z_offset_cm": ("simc_target_z_offset_cm",),
    "simc_spec_e_p_offset_gev": ("simc_spec_e_p_offset_gev", "simc_spec_e_p_offset"),
    "simc_spec_e_theta_offset_deg": ("simc_spec_e_theta_offset_deg", "simc_spec_e_theta_offset"),
    "simc_hms_xptar_offset_mrad": ("simc_hms_xptar_offset_mrad", "simc_hms_pointing_offset"),
    "simc_spec_p_x_offset_cm": ("simc_spec_p_x_offset_cm",),
    "simc_spec_p_y_offset_cm": ("simc_spec_p_y_offset_cm",),
    "simc_spec_p_z_offset_cm": ("simc_spec_p_z_offset_cm",),
    "geant4_nps_x_offset_cm": ("geant4_nps_x_offset_cm",),
    "geant4_nps_y_offset_cm": ("geant4_nps_y_offset_cm",),
    "geant4_nps_z_offset_cm": ("geant4_nps_z_offset_cm",),
    "source_simc_infile": ("source_simc_infile",),
}

ZERO_COLUMNS = {
    "simc_spec_e_p_offset_gev", "simc_spec_e_theta_offset_deg",
    "simc_hms_xptar_offset_mrad", "simc_spec_p_x_offset_cm",
    "simc_spec_p_y_offset_cm", "simc_spec_p_z_offset_cm",
    "geant4_nps_x_offset_cm", "geant4_nps_y_offset_cm",
    "geant4_nps_z_offset_cm",
}

CHANNELS = {
    "exclusive": {"prefix": "nps_excl_pi0_", "doing_semi": "0", "which_pion": "0"},
    "sidis": {"prefix": "nps_sidis_pi0_", "doing_semi": "1", "which_pion": "0"},
    "delta": {"prefix": "nps_delta_pi0_", "doing_semi": "0", "which_pion": "2"},
}

# Physics-method references copied exactly below:
#   config/acceptance_cuts.conf [global]
#     SIMC generation acceptance; angular slopes are converted rad -> mrad.
#   /u/group/nps/bosted/simc_new/makesimcinfiles.f
#     Egamma_gen_max/Eloss/Coulomb settings.
#   /w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/calc_mispointing.cxx
#     HMS x/y mispointing parameterization (angle in degrees, offsets in cm).
BOSTED_GENERATOR = "/u/group/nps/bosted/simc_new/makesimcinfiles.f"
DEFAULT_ACCEPTANCE_CONFIG = Path(__file__).resolve().parents[1] / "config/acceptance_cuts.conf"
MISPOINTING_REFERENCE = (
    "/w/hallc-scshelf2102/nps/singhav/nps_analysis/elastic_analysis/"
    "calc_mispointing.cxx"
)

GLOBAL_ACCEPTANCE_KEYS = (
    "hms_sim.delta_min", "hms_sim.delta_max",
    "hms_sim.gtr_th_min", "hms_sim.gtr_th_max",
    "hms_sim.gtr_ph_min", "hms_sim.gtr_ph_max",
    "nps_cluster.energy_min", "nps_cluster.x_min", "nps_cluster.x_max",
    "nps_cluster.y_min", "nps_cluster.y_max",
)


def normalize_header(name: str) -> str:
    return name.strip().lower()


def parse_float(value: str) -> Optional[float]:
    try:
        parsed = float(value.strip())
    except (AttributeError, TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def parse_int(value: str) -> Optional[int]:
    try:
        return int(value.strip())
    except (AttributeError, TypeError, ValueError):
        return None


def parse_set(csv_text: str) -> set[str]:
    return {value.strip().lower() for value in csv_text.split(",") if value.strip()}


def is_physical(row: KinRow) -> bool:
    return (
        4.0 <= row.ebeam_gev <= 12.0
        and 0.2 <= row.hms_p_gev < row.ebeam_gev
        and 5.0 <= row.hms_theta_deg <= 50.0
        and 5.0 <= row.nps_theta_deg <= 40.0
        and row.nps_target_distance_cm >= 250.0
    )


def read_existing_values(path: Path) -> Dict[str, Dict[str, str]]:
    """Read user-maintained values, including names from the legacy schema."""
    if not path.exists():
        return {}
    with path.open("r", newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        if reader.fieldnames is None:
            return {}
        headers = {normalize_header(name): name for name in reader.fieldnames}
        result: Dict[str, Dict[str, str]] = {}
        for raw in reader:
            kin_header = headers.get("kin_old")
            kin_old = (raw.get(kin_header, "") or "").strip() if kin_header else ""
            if not kin_old:
                continue
            values: Dict[str, str] = {}
            for canonical, aliases in PRESERVED_COLUMNS.items():
                for alias in aliases:
                    header = headers.get(alias)
                    value = (raw.get(header, "") or "").strip() if header else ""
                    if value:
                        values[canonical] = value
                        break
            result[kin_old] = values
        return result


def repair_surplus_comment_fields(
    raw: Dict[Optional[str], object], fieldnames: Sequence[str]
) -> Dict[Optional[str], object]:
    """Rejoin unquoted commas known to belong to the master CSV comments field."""
    extras = raw.get(None)
    if not isinstance(extras, list) or not extras:
        return raw
    comment_index = next(
        (i for i, name in enumerate(fieldnames) if normalize_header(name) == "comments"),
        -1,
    )
    if comment_index < 0:
        raise ValueError("Cannot repair ragged CSV row without comments column")
    values = [str(raw.get(name, "") or "") for name in fieldnames] + [str(v) for v in extras]
    surplus = len(extras)
    repaired_values = (
        values[:comment_index]
        + [",".join(values[comment_index : comment_index + surplus + 1])]
        + values[comment_index + surplus + 1 :]
    )
    if len(repaired_values) != len(fieldnames):
        raise ValueError("Could not repair ragged master CSV row")
    return dict(zip(fieldnames, repaired_values))


def generate_simulation_csv(
    args: argparse.Namespace, cfg_path: Path, out_path: Path
) -> tuple[int, int, int]:
    """Build the canonical simulation CSV from the master experiment table."""
    if not cfg_path.exists():
        raise FileNotFoundError(f"Config CSV not found: {cfg_path}")
    allowed_types = parse_set(args.types)
    allowed_targets = parse_set(args.targets)
    allowed_goodness = parse_set(args.goodness_codes)
    if not allowed_types or not allowed_targets or not allowed_goodness:
        raise ValueError("--types, --targets, and --goodness-codes cannot be empty")
    if len(allowed_targets) != 1:
        raise ValueError("Generate one target at a time; provide exactly one --targets value")
    if args.simc_ngen <= 0:
        raise ValueError("--simc-ngen must be positive")
    target_offsets = {
        "simc_target_x_offset_cm": args.simc_target_x_offset_cm,
        "simc_target_y_offset_cm": args.simc_target_y_offset_cm,
        "simc_target_z_offset_cm": args.simc_target_z_offset_cm,
    }
    if any(parse_float(value) is None for value in target_offsets.values()):
        raise ValueError("SIMC target offsets must be finite numbers")

    existing = {} if args.reset_offsets else read_existing_values(out_path)
    by_kin: Dict[str, KinRow] = {}
    eligible_kins: set[str] = set()
    repaired_rows = 0
    invalid_rows = 0
    with cfg_path.open("r", newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        if reader.fieldnames is None:
            raise ValueError("Config CSV has no header")
        original_by_norm = {normalize_header(name): name for name in reader.fieldnames}
        required = {
            "new_kin", "kin_old", "run_number", "type", "target", "goodness_code",
            "beamenergy", "hms_p", "hms_thet", "nps_thet", "nps_target_distance",
        }
        missing = sorted(required - original_by_norm.keys())
        if missing:
            raise ValueError(f"Missing required config columns: {', '.join(missing)}")
        for raw in reader:
            if None in raw:
                raw = repair_surplus_comment_fields(raw, reader.fieldnames)
                repaired_rows += 1

            def get(key: str) -> str:
                return (raw.get(original_by_norm[key], "") or "").strip()

            kin_old = get("kin_old")
            if not kin_old or get("type").lower() not in allowed_types:
                continue
            eligible_kins.add(kin_old)
            if get("target").lower() not in allowed_targets:
                continue
            if get("goodness_code").lower() not in allowed_goodness:
                continue
            run = parse_int(get("run_number"))
            values = [
                parse_float(get("beamenergy")), parse_float(get("hms_p")),
                parse_float(get("hms_thet")), parse_float(get("nps_thet")),
                parse_float(get("nps_target_distance")),
            ]
            if run is None or any(value is None for value in values):
                invalid_rows += 1
                continue
            candidate = KinRow(
                new_kin=get("new_kin"), kin_old=kin_old, target=get("target"),
                source_run=run, ebeam_gev=float(values[0]), hms_p_gev=float(values[1]),
                hms_theta_deg=float(values[2]), nps_theta_deg=float(values[3]),
                nps_target_distance_cm=float(values[4]),
            )
            if not is_physical(candidate):
                invalid_rows += 1
                continue
            current = by_kin.get(kin_old)
            if current is None or candidate.source_run < current.source_run:
                by_kin[kin_old] = candidate

    missing_kins = sorted(eligible_kins - by_kin.keys())
    if missing_kins:
        raise RuntimeError("No valid selected row for: " + ", ".join(missing_kins))
    if not by_kin:
        raise RuntimeError("No rows selected for simulation kinematics table")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    selection_rule = args.selection_rule or (
        "per Kin_old choose minimum valid run_number; "
        f"Type={args.types}; target={args.targets}; goodness_code={args.goodness_codes}"
    )
    with out_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
        writer.writeheader()
        for kin_old in sorted(by_kin):
            row = by_kin[kin_old]
            saved = existing.get(kin_old, {})
            simc_p = saved.get("simc_spec_p_p_gev", "")
            if parse_float(simc_p) is None:
                simc_p = f"{row.ebeam_gev - row.hms_p_gev:.6f}"
            output: Dict[str, str] = {
                "new_kin": row.new_kin, "kin_old": row.kin_old, "target": row.target,
                "ebeam_gev": f"{row.ebeam_gev:.6f}",
                "hms_p_gev": f"{row.hms_p_gev:.6f}",
                "hms_theta_deg": f"{row.hms_theta_deg:.6f}",
                "nps_theta_deg": f"{row.nps_theta_deg:.6f}",
                "nps_target_distance_cm": f"{row.nps_target_distance_cm:.6f}",
                "simc_spec_p_p_gev": simc_p, "source_run": str(row.source_run),
                "simc_ngen": saved.get("simc_ngen", str(args.simc_ngen)),
                "source_simc_infile": (
                    args.source_simc_infile
                    or args.template
                    or saved.get("source_simc_infile", "")
                ),
                "selection_rule": selection_rule, "last_generated_utc": now_utc,
            }
            for column in ZERO_COLUMNS:
                output[column] = saved.get(column, "0.0")
            for column, default in target_offsets.items():
                output[column] = saved.get(column, default)
            writer.writerow(output)
    return len(by_kin), repaired_rows, invalid_rows


def kin_old_to_tag(kin_old: str) -> str:
    tag = kin_old.strip()
    if tag.lower().startswith("kinc_"):
        tag = tag[5:]
    tag = tag.lower()
    if not re.fullmatch(r"[a-z0-9][a-z0-9_-]*", tag):
        raise ValueError(f"Unsafe Kin_old output tag: {kin_old!r}")
    return tag


def parse_number(row: Dict[str, str], key: str, kin_old: str) -> float:
    value = (row.get(key, "") or "").strip()
    if not value:
        raise ValueError(f"{kin_old}: empty required value {key}")
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"{kin_old}: invalid {key}={value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{kin_old}: non-finite {key}={value!r}")
    return parsed


def hms_mispointing_cm(theta_deg: float) -> tuple[float, float]:
    """Return (x, y) from calc_mispointing.cxx; input degrees, output cm."""
    theta = abs(theta_deg)
    theta_x = min(theta, 50.0)
    theta_y = min(theta, 40.0)
    offset_x = 0.1 * (2.37 - 0.086 * theta_x + 0.0012 * theta_x**2)
    offset_y = 0.1 * (0.52 - 0.012 * theta_y + 0.002 * theta_y**2)
    return offset_x, offset_y


def load_global_acceptance(path: Path) -> Dict[str, float]:
    """Read the authoritative [global] analysis/simulation acceptance values."""
    parser = configparser.ConfigParser(
        interpolation=None, inline_comment_prefixes=("#", ";")
    )
    with path.open() as stream:
        parser.read_file(stream)
    if "global" not in parser:
        raise RuntimeError(f"Acceptance config has no [global] section: {path}")
    values: Dict[str, float] = {}
    for key in GLOBAL_ACCEPTANCE_KEYS:
        if key not in parser["global"]:
            raise RuntimeError(f"Acceptance config [global] missing {key}: {path}")
        try:
            value = float(parser["global"][key])
        except ValueError as exc:
            raise ValueError(f"Acceptance config [global] has invalid {key}: {path}") from exc
        if not math.isfinite(value):
            raise ValueError(f"Acceptance config [global] has non-finite {key}: {path}")
        values[key] = value
    for stem in ("hms_sim.delta", "hms_sim.gtr_th", "hms_sim.gtr_ph", "nps_cluster.x", "nps_cluster.y"):
        if values[f"{stem}_min"] >= values[f"{stem}_max"]:
            raise ValueError(f"Acceptance config [global] has invalid {stem} range: {path}")
    if values["nps_cluster.energy_min"] <= 0:
        raise ValueError(f"Acceptance config [global] requires positive nps_cluster.energy_min: {path}")
    return values


def replace_assignment(
    lines: List[str], key: str, value: str, comment: Optional[str] = None
) -> None:
    pattern = re.compile(rf"^(\s*{re.escape(key)}\s*=\s*)([^;\n]*)(.*)$")
    matches = [(index, pattern.match(line)) for index, line in enumerate(lines)]
    matches = [(index, match) for index, match in matches if match]
    if len(matches) != 1:
        raise RuntimeError(f"Template must contain exactly one assignment for {key}; found {len(matches)}")
    index, match = matches[0]
    assert match is not None
    prefix, _old, suffix = match.groups()
    if comment is not None:
        suffix = f";  {comment}"
    replacement = f"{prefix}{value}{suffix}"
    if lines[index].endswith("\n") and not replacement.endswith("\n"):
        replacement += "\n"
    lines[index] = replacement


def resolve_template(csv_path: Path, global_template: Optional[Path], row: Dict[str, str]) -> Path:
    configured = (row.get("source_simc_infile", "") or "").strip()
    candidate = Path(configured).expanduser() if configured else global_template
    if candidate is None:
        raise RuntimeError(
            f"{row['kin_old']}: no source_simc_infile and no --template supplied"
        )
    if not candidate.is_absolute():
        candidate = csv_path.parent / candidate
    candidate = candidate.resolve()
    if not candidate.exists():
        raise FileNotFoundError(f"{row['kin_old']}: SIMC template not found: {candidate}")
    return candidate


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate canonical simulation CSV or reviewed-CSV SIMC infiles"
    )
    parser.add_argument(
        "--gen_infile", "--gen-infile", action="store_true",
        help="Generate infiles from an already reviewed simulation CSV",
    )
    parser.add_argument("--config", help="Master experiment CSV; enables CSV generation")
    parser.add_argument("--sim-config", required=True)
    parser.add_argument("--types", default="production,Production")
    parser.add_argument("--targets", default="LH2")
    parser.add_argument("--goodness-codes", default="0")
    parser.add_argument("--selection-rule", default="")
    parser.add_argument("--source-simc-infile", default="")
    parser.add_argument("--simc-ngen", type=int, default=1000000)
    parser.add_argument("--simc-target-x-offset-cm", default="0.131")
    parser.add_argument("--simc-target-y-offset-cm", default="0.005")
    parser.add_argument("--simc-target-z-offset-cm", default="0.0")
    parser.add_argument("--reset-offsets", action="store_true")
    parser.add_argument("--acceptance-config", default=str(DEFAULT_ACCEPTANCE_CONFIG))
    parser.add_argument("--template", help="Fallback template when row provenance is empty")
    parser.add_argument("--out-dir", help="Generate SIMC infiles into this directory")
    parser.add_argument(
        "--channel",
        choices=["exclusive", "sidis", "delta", "all"],
        default="exclusive",
    )
    parser.add_argument("--file-prefix", default="", help="Override prefix for one channel")
    parser.add_argument("--kin", action="append", default=[])
    args = parser.parse_args()

    sim_cfg = Path(args.sim_config).resolve()
    global_template = Path(args.template).expanduser().resolve() if args.template else None
    if not args.gen_infile:
        if not args.config:
            parser.error("CSV generation requires --config")
        if args.out_dir:
            parser.error("--out-dir requires --gen_infile")
        row_count, repaired_rows, invalid_rows = generate_simulation_csv(
            args, Path(args.config).expanduser().resolve(), sim_cfg
        )
        print(
            f"Simulation CSV ready for review: {sim_cfg}\n"
            f"Rows: {row_count}; repaired comments: {repaired_rows}; "
            f"skipped invalid: {invalid_rows}.\n"
            "No SIMC infiles generated. Verify/edit CSV, then run:\n"
            f"  python3 {Path(__file__).resolve()} --gen_infile "
            f"--sim-config {sim_cfg} --out-dir config/simc_infiles --channel all\n"
        )
        return 0
    if args.config:
        parser.error("--config cannot be used with --gen_infile; review the existing CSV first")
    if not args.out_dir:
        parser.error("--gen_infile requires --out-dir")

    acceptance_cfg = Path(args.acceptance_config).expanduser().resolve()
    out_dir = Path(args.out_dir)
    if not sim_cfg.exists():
        raise FileNotFoundError(f"Simulation config not found: {sim_cfg}")
    if not acceptance_cfg.exists():
        raise FileNotFoundError(f"Acceptance config not found: {acceptance_cfg}")
    if global_template is not None and not global_template.exists():
        raise FileNotFoundError(f"Template infile not found: {global_template}")
    acceptance = load_global_acceptance(acceptance_cfg)

    requested_kins = {kin.strip() for kin in args.kin if kin.strip()}
    selected_channels = list(CHANNELS) if args.channel == "all" else [args.channel]
    with sim_cfg.open("r", newline="") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        if reader.fieldnames is None:
            raise RuntimeError("Simulation config CSV has no header")
        headers = {name.strip().lower(): name for name in reader.fieldnames}
        missing = [key for key in REQUIRED_COLUMNS if key not in headers]
        if missing:
            raise RuntimeError("Simulation config missing columns: " + ", ".join(missing))
        rows: List[Dict[str, str]] = []
        for raw in reader:
            if None in raw:
                raise RuntimeError("Simulation config contains ragged CSV row")
            normalized = {
                key: (raw.get(original, "") or "").strip()
                for key, original in headers.items()
            }
            if normalized.get("kin_old"):
                rows.append(normalized)

    if not rows:
        raise RuntimeError("Simulation config CSV has no usable rows")
    template_cache: Dict[Path, List[str]] = {}
    output_paths: set[Path] = set()
    pending_outputs: Dict[Path, str] = {}

    for row in rows:
        kin_old = row["kin_old"]
        if requested_kins and kin_old not in requested_kins:
            continue
        tag = kin_old_to_tag(kin_old)
        ebeam = parse_number(row, "ebeam_gev", kin_old)
        hms_p = parse_number(row, "hms_p_gev", kin_old)
        hms_theta = parse_number(row, "hms_theta_deg", kin_old)
        nps_theta = parse_number(row, "nps_theta_deg", kin_old)
        nps_distance = parse_number(row, "nps_target_distance_cm", kin_old)
        pion_p = parse_number(row, "simc_spec_p_p_gev", kin_old)
        ngen_value = parse_number(row, "simc_ngen", kin_old)
        ngen = int(ngen_value)
        if ngen <= 0 or ngen != ngen_value:
            raise ValueError(f"{kin_old}: simc_ngen must be a positive integer")
        target_x = parse_number(row, "simc_target_x_offset_cm", kin_old)
        target_y = parse_number(row, "simc_target_y_offset_cm", kin_old)
        target_z = parse_number(row, "simc_target_z_offset_cm", kin_old)
        hms_p += parse_number(row, "simc_spec_e_p_offset_gev", kin_old)
        hms_theta += parse_number(row, "simc_spec_e_theta_offset_deg", kin_old)
        hms_xptar = parse_number(row, "simc_hms_xptar_offset_mrad", kin_old)
        simc_p_x = parse_number(row, "simc_spec_p_x_offset_cm", kin_old)
        simc_p_y = parse_number(row, "simc_spec_p_y_offset_cm", kin_old)
        simc_p_z = parse_number(row, "simc_spec_p_z_offset_cm", kin_old)
        if not (
            4.0 <= ebeam <= 12.0 and 0.2 <= hms_p < ebeam
            and 5.0 <= hms_theta <= 50.0 and 5.0 <= nps_theta <= 40.0
            and nps_distance >= 250.0 and 0.1 <= pion_p <= 12.0
        ):
            raise ValueError(f"{kin_old}: effective SIMC kinematics outside physical bounds")

        nps_energy_min = acceptance["nps_cluster.energy_min"]
        if pion_p <= nps_energy_min:
            raise ValueError(f"{kin_old}: NPS central momentum must exceed global energy minimum")
        nps_delta_min = 100.0 * (nps_energy_min / pion_p - 1.0)
        # SIMC xptar is vertical and yptar horizontal; the NPS analysis names
        # those detector coordinates y and x, respectively. Slopes are x/z.
        nps_xptar_min = 1000.0 * acceptance["nps_cluster.y_min"] / nps_distance
        nps_xptar_max = 1000.0 * acceptance["nps_cluster.y_max"] / nps_distance
        nps_yptar_min = 1000.0 * acceptance["nps_cluster.x_min"] / nps_distance
        nps_yptar_max = 1000.0 * acceptance["nps_cluster.x_max"] / nps_distance
        # Bosted uses the 15%-below-central HMS momentum edge, converts GeV
        # to MeV, and writes f7.0: (Ebeam - 0.85*|Phms|)*1000.
        egamma_gen_max = 1000.0 * (ebeam - 0.85 * abs(hms_p))
        hms_offset_x, hms_offset_y = hms_mispointing_cm(hms_theta)

        template = resolve_template(sim_cfg, global_template, row)
        if template not in template_cache:
            template_cache[template] = template.read_text().splitlines(keepends=True)
        common_replacements = {
            "ngen": str(ngen),
            "Ebeam": f"{1000.0 * ebeam:.1f}",
            "spec%e%P": f"{1000.0 * hms_p:.2f}",
            "spec%e%theta": f"{hms_theta:.2f}",
            "spec%p%P": f"{1000.0 * pion_p:.2f}",
            "spec%p%theta": f"{nps_theta:.2f}",
            "SPedge%e%delta%min": f"{acceptance['hms_sim.delta_min']:.6g}",
            "SPedge%e%delta%max": f"{acceptance['hms_sim.delta_max']:.6g}",
            "SPedge%e%yptar%min": f"{1000.0 * acceptance['hms_sim.gtr_ph_min']:.6g}",
            "SPedge%e%yptar%max": f"{1000.0 * acceptance['hms_sim.gtr_ph_max']:.6g}",
            "SPedge%e%xptar%min": f"{1000.0 * acceptance['hms_sim.gtr_th_min']:.6g}",
            "SPedge%e%xptar%max": f"{1000.0 * acceptance['hms_sim.gtr_th_max']:.6g}",
            "SPedge%p%delta%min": f"{nps_delta_min:.2f}",
            "SPedge%p%delta%max": "5.0",
            "SPedge%p%yptar%min": f"{nps_yptar_min:.2f}",
            "SPedge%p%yptar%max": f"{nps_yptar_max:.2f}",
            "SPedge%p%xptar%min": f"{nps_xptar_min:.2f}",
            "SPedge%p%xptar%max": f"{nps_xptar_max:.2f}",
            "drift_to_cal": f"{nps_distance:.2f}",
            "targ%xoffset": f"{target_x:.6f}",
            "targ%yoffset": f"{target_y:.6f}",
            "targ%zoffset": f"{target_z:.6f}",
            "spec%e%offset%x": f"{hms_offset_x:.6f}",
            "spec%e%offset%y": f"{hms_offset_y:.6f}",
            "spec%e%offset%xptar": f"{hms_xptar:.6f}",
            "spec%p%offset%x": f"{simc_p_x:.6f}",
            "spec%p%offset%y": f"{simc_p_y:.6f}",
            "spec%p%offset%z": f"{simc_p_z:.6f}",
            "correct_Eloss": "1",
            "using_Coulomb": "0",
            "Egamma_gen_max": f"{egamma_gen_max:.0f}",
        }
        for channel in selected_channels:
            channel_cfg = CHANNELS[channel]
            prefix = (
                args.file_prefix
                if args.file_prefix and (len(selected_channels) == 1 or channel == "exclusive")
                else channel_cfg["prefix"]
            )
            out_file = out_dir / f"{prefix}{tag}.inp"
            if out_file in output_paths:
                raise RuntimeError(f"Duplicate output path generated: {out_file}")
            output_paths.add(out_file)
            lines = list(template_cache[template])
            for key, value in common_replacements.items():
                replace_assignment(lines, key, value)
            replace_assignment(lines, "doing_semi", channel_cfg["doing_semi"])
            replace_assignment(
                lines,
                "which_pion",
                channel_cfg["which_pion"],
                "0=p->pi+, 1=n->pi-, 2=p->pi0 Delta+, 3=n->pi0 Delta0",
            )
            provenance = (
                f"; generated for {kin_old}; channel: {channel}; "
                f"source template: {template}\n"
                f"; acceptance reference: {acceptance_cfg} [global]\n"
                f"; Eloss/Coulomb/Egamma reference: {BOSTED_GENERATOR}\n"
                f"; HMS x/y offset reference: {MISPOINTING_REFERENCE}\n"
            )
            pending_outputs[out_file] = provenance + "".join(lines)

    if not pending_outputs:
        raise RuntimeError("No SIMC infiles generated; check --kin filters")
    out_dir.mkdir(parents=True, exist_ok=True)
    for out_file, content in pending_outputs.items():
        out_file.write_text(content)
    print(f"Generated {len(pending_outputs)} SIMC infile(s) in {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
