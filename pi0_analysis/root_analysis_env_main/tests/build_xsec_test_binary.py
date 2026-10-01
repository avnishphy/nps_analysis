#!/usr/bin/env python3
"""Build a ratio extractor with explicit test edges in a temporary config copy."""
import argparse
import math
from pathlib import Path
import re
import shutil
import subprocess

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out-dir", required=True, type=Path)
parser.add_argument("--phi-bins", required=True, type=int)
parser.add_argument("--t-edges", required=True)
parser.add_argument("--q2-edges", required=True)
parser.add_argument("--xb-edges", required=True)
args = parser.parse_args()

if args.phi_bins < 3:
    parser.error("--phi-bins must be at least 3")
source_dir = Path(__file__).resolve().parents[1] / "src/xsec_extract"
args.out_dir.mkdir(parents=True, exist_ok=True)
source = args.out_dir / "excl_xsec_pi0_analysis_simc_model.C"
config = args.out_dir / "xsec_config.h"
shutil.copyfile(source_dir / source.name, source)
text = (source_dir / config.name).read_text()


def set_vector(name, values):
    global text
    pattern = rf"std::vector<double> {name} = \{{.*?\}};"
    replacement = f"std::vector<double> {name} = {{{values}}};"
    text, count = re.subn(pattern, replacement, text, count=1, flags=re.S)
    if count != 1:
        raise RuntimeError(f"could not replace {name} in test config")


def parse_edges(raw):
    values = [float(x) for x in raw.split(",")]
    if len(values) < 2 or any(not math.isfinite(x) for x in values):
        parser.error("edge lists need at least two finite values")
    return ", ".join(map(repr, values))


set_vector("phi_bin_edges", ", ".join(
    f"2*TMath::Pi()*{i}/{args.phi_bins}" for i in range(args.phi_bins + 1)
))
set_vector("t_bin_edges", parse_edges(args.t_edges))
set_vector("q2_bin_edges", parse_edges(args.q2_edges))
xb = parse_edges(args.xb_edges)
pattern = r"std::vector<std::vector<double>> xb_bin_edges_by_q2 = \{.*?\};"
rows = ", ".join("{" + xb + "}" for _ in range(len(args.q2_edges.split(",")) - 1))
text, count = re.subn(pattern, f"std::vector<std::vector<double>> xb_bin_edges_by_q2 = {{{rows}}};",
                      text, count=1, flags=re.S)
if count != 1:
    raise RuntimeError("could not replace xB rows in test config")
config.write_text(text)

flags = subprocess.check_output(["root-config", "--cflags", "--libs"], text=True).split()
binary = args.out_dir / f"extractor_{args.phi_bins}"
subprocess.run(["g++", "-O2", "-std=c++17", str(source), f"-I{source_dir}",
                *flags, "-lMinuit2", "-o", str(binary)], check=True)
print(binary)
