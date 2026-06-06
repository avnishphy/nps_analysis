# Two-Stage Smearing Implementation

## Overview

The code has been modified to implement a **two-stage smearing optimization** approach that separates position and energy smearing into distinct stages, each using the most appropriate observable.

## Physical Motivation

### Why Two Stages?

1. **Position smearing** primarily affects the **angular resolution** between photons, which directly impacts the **invariant mass M_γγ** through:
   ```
   M_γγ² = 2E₁E₂(1 - cos θ₁₂)
   ```
   Where θ₁₂ is the opening angle between photons. Position errors change the photon directions, affecting this angle.

2. **Energy smearing** affects the **absolute energy scale**, which is best constrained by **missing mass M_x**:
   ```
   M_x = √[(E_beam + M_p - E_e' - E_γ1 - E_γ2)² - p_miss²]
   ```
   This provides an absolute calibration using the known proton mass.

### Why This is Better

- **Decouples systematics**: Position and energy effects are largely independent
- **Uses appropriate constraints**: Each observable constrains what it's most sensitive to
- **Avoids over-constraining**: Prevents fitting correlated observables simultaneously
- **More physical**: Matches how detectors actually work (position and energy measured separately)

## Implementation Details

### Stage 1: Position Smearing (σ_pos)

**Objective**: Fit `sigma_pos` using **invariant mass M_γγ only**

**Fixed parameters:**
- `mu = 1.0` (no energy scale adjustment)
- `sigma = 0.02` (minimal energy resolution effect)

**Optimization:**
1. Coarse grid scan over `sigma_pos` range (Config::SIGMA_POS_MIN to MAX)
2. Fine refinement around minimum using decreasing step sizes
3. Minimize: `χ²_M_γγ = Σ((sim - data)² / (sim + data))` over invariant mass bins

**Output**: Optimal `sigma_pos` that matches the shape/width of the M_γγ peak

### Stage 2: Energy Smearing (μ, σ)

**Objective**: Fit `mu` and `sigma` using **missing mass M_miss only**

**Fixed parameters:**
- `sigma_pos` = value from Stage 1

**Optimization:**
1. Coarse 2D grid scan over (mu, sigma) ranges
2. Fine refinement around minimum using decreasing step sizes
3. Minimize: `χ²_M_miss = Σ((sim - data)² / (sim + data))` over missing mass bins

**Output**: Optimal energy scale `mu` and resolution `sigma` that matches M_miss peak position

## New Functions Added

### `eval_chi2_mpi0_only()`
Evaluates chi-squared using **only** the invariant mass histogram:
- Used in Stage 1 for position smearing
- Applies both position and energy smearing, but only compares M_γγ
- Returns: `χ²_M_γγ`

### `eval_chi2_mmiss_only()`
Evaluates chi-squared using **only** the missing mass histogram:
- Used in Stage 2 for energy smearing
- Applies both position and energy smearing, but only compares M_miss
- Returns: `χ²_M_miss`

### Modified `fit_section()`
Now implements the two-stage optimization:
1. Stage 1: Position smearing using M_γγ
2. Stage 2: Energy smearing using M_miss
3. Final: Calculate combined χ² for reporting

## Usage

The code works exactly the same from the user's perspective:

```bash
./nps_sim_smearing_new data.root dataTree sim.root simTree out.root 8 8 -30 30 -36 36 0.2 20
```

**No changes to command-line interface!**

## Configuration Options

In the `Config` namespace at the top of the file:

```cpp
// Enable/disable position smearing
const bool ENABLE_POSITION_SMEARING = true;  // Set to false to skip Stage 1

// Observable weights for combined chi2 (used in final reporting only)
const double W_MPI0 = 1.0;   // Weight for M_γγ
const double W_MMISS = 2.0;  // Weight for M_miss
```

If `ENABLE_POSITION_SMEARING = false`:
- Stage 1 is skipped (sigma_pos = 0)
- Only Stage 2 runs (energy smearing)

## Output Interpretation

### Console Output Example:
```
==== TWO-STAGE OPTIMIZATION ====
Stage 1: Position smearing (sigma_pos) using M_γγ
Stage 2: Energy smearing (mu, sigma) using M_miss with fixed sigma_pos

--- STAGE 1: Fitting position smearing (sigma_pos) using M_γγ ---
Stage 1 coarse minimum: sigma_pos=0.45 cm, chi2_M_γγ=125.3
Stage 1: Refining sigma_pos...
Stage 1 final: sigma_pos=0.432 cm, chi2_M_γγ=124.8

--- STAGE 2: Fitting energy smearing (mu, sigma) using M_miss ---
  Using fixed sigma_pos=0.432 cm from Stage 1
Stage 2 coarse minimum: mu=1.025, sigma=0.045, chi2_M_miss=156.7
Stage 2: Refining mu and sigma...
Stage 2 final: mu=1.0237, sigma=0.0443, chi2_M_miss=155.2

--- FINAL RESULTS ---
mu=1.0237, sigma=0.0443, sigma_pos=0.432 cm
chi2_combined=405.2 (w_M_γγ=1.0, w_M_miss=2.0)
```

### Interpretation:
- **sigma_pos = 0.432 cm**: Position resolution that matches M_γγ shape
- **mu = 1.0237**: Photon energies need to be scaled up by 2.37% to match M_miss
- **sigma = 0.0443**: Energy resolution parameter
- **chi2_combined**: Weighted sum for overall quality check

## Visualization

The chi² scan plots show:
- **1D slices for mu**: Uses M_miss chi² (Stage 2 observable)
- **1D slices for sigma**: Uses M_miss chi² (Stage 2 observable)  
- **1D slices for sigma_pos**: Uses M_γγ chi² (Stage 1 observable)
- **2D map (mu vs sigma)**: Uses combined chi² for overview

This clearly shows which observable was used to constrain each parameter.

## Benefits

1. **Physical clarity**: Each stage has a clear physical meaning
2. **Computational efficiency**: 2D + 1D scans instead of 3D scan
3. **Better convergence**: Each stage has well-defined minimum
4. **Interpretability**: Can see how each observable constrains parameters
5. **Debugging**: Can examine each stage independently

## Technical Notes

- The combined chi² is still calculated at the end for quality assessment
- Visualization includes both individual (per-stage) and combined chi² values
- The approach is particularly effective when:
  - Position resolution is important (large detector elements)
  - Missing mass calculation is available (elastic/exclusive reactions)
  - Energy scale needs absolute calibration (known mass constraints)

## Comparison to Old Approach

### Old Method (Simultaneous):
- Fit all parameters (mu, sigma, sigma_pos) together
- Used weighted combination: `χ²_total = w₁·χ²_M_γγ + w₂·χ²_M_miss`
- Required careful weight tuning
- Could have parameter degeneracies

### New Method (Two-Stage):
- Fit position first using M_γγ only
- Then fit energy using M_miss only  
- Natural separation by physics
- No weight tuning needed during optimization
- Clearer convergence properties

## Validation

To validate the new approach:

1. **Check M_γγ shape**: Should match data well after Stage 1
2. **Check M_miss position**: Should match data well after Stage 2
3. **Compare combined chi²**: Should be comparable or better than old method
4. **Check parameter stability**: Run multiple times with different random seeds

## Future Enhancements

Possible improvements:
- Iterate between stages (Stage 1 → Stage 2 → Stage 1 → ...)
- Use different resolution models per stage
- Add Stage 3 for HMS electron momentum calibration
- Implement uncertainty propagation between stages
