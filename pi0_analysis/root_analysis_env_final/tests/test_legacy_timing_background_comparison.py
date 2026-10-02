#!/usr/bin/env python3
"""Unit checks for the independent legacy timing-box comparator."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src" / "background_fit"))

from compare_legacy_timing_background import (  # noqa: E402
    REGION_BITS,
    legacy_estimate_from_region_masks,
)


def main() -> None:
    masks = np.concatenate((
        np.repeat(REGION_BITS["prompt"], 10),
        np.repeat(REGION_BITS["diagonal"], 6),
        np.repeat(REGION_BITS["horizontal"], 12),
        np.repeat(REGION_BITS["vertical"], 18),
        np.repeat(REGION_BITS["full1"], 9),
        np.repeat(REGION_BITS["full2"], 18),
        # Region bits are counted independently, matching histogram fills.
        np.asarray([REGION_BITS["diagonal"] | REGION_BITS["horizontal"]]),
    )).astype(np.uint32)
    result = legacy_estimate_from_region_masks(masks)

    assert result["raw_mask_prompt"] == 10.0
    assert result["raw_mask_diagonal"] == 7.0
    assert result["raw_mask_horizontal"] == 13.0
    assert result["raw_mask_vertical"] == 18.0
    assert result["raw_mask_full1"] == 9.0
    assert result["raw_mask_full2"] == 18.0

    expected = 7.0 / 6.0 + 0.5 * (18.0 / 6.0 + 13.0 / 6.0) - 0.5 * (
        9.0 / 9.0 + 18.0 / 9.0)
    variance = (7.0 / 36.0 + 0.25 * 18.0 / 36.0 +
                0.25 * 13.0 / 36.0 + 0.25 * 9.0 / 81.0 +
                0.25 * 18.0 / 81.0)
    assert math.isclose(result["raw_mask_box_formula_estimate"], expected,
                        rel_tol=0.0, abs_tol=1.0e-12)
    assert math.isclose(result["raw_mask_box_formula_standard_error"],
                        math.sqrt(variance), rel_tol=0.0, abs_tol=1.0e-12)
    print("legacy timing background comparison test: PASS")


if __name__ == "__main__":
    main()
