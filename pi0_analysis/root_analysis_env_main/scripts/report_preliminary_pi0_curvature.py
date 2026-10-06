#!/usr/bin/env python3
"""Compatibility entry point for the dynamic seven-parameter M0 report.

The former KinC_x36_4-only implementation hard-coded four exterior nuisance
triplets. Keeping it callable would reintroduce fitted Q2/xB coordinates, so
all report generation now uses the campaign metadata-driven implementation.
"""
from report_preliminary_pi0_curvature_dynamic import main


if __name__ == "__main__":
    main()
