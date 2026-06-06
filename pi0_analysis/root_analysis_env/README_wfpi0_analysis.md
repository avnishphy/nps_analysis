# NPS Analysis for Production/Waveform Calibrated Data

## Overview
`nps_analysis_wfpi0.C` is a modified version of `nps_analysis.C` designed to analyze waveform-calibrated production data from the DVCS experiment.

## Key Changes from nps_analysis.C

### 1. Tree Name
- **Original**: `"T"`
- **New**: `"t_prod"`

### 2. NPS Branch Names
- **Original**: `NPS.cal.*`
- **New**: `NPS.prod.*`

Specifically:
- `NPS.cal.nclust` → `NPS.prod.nclust`
- `NPS.cal.clusE` → `NPS.prod.clusE`
- `NPS.cal.clusX` → `NPS.prod.clusX`
- `NPS.cal.clusY` → `NPS.prod.clusY`
- `NPS.cal.clusT` → `NPS.prod.clusT`

### 3. HMS Scaler Branches
- **Original**: `H.BCM2.*`
- **New**: `H.BCM4A.*`

Specifically:
- `H.BCM2.scalerCurrent` → `H.BCM4A.scalerCurrent`
- `H.BCM2.scalerCharge` → `H.BCM4A.scalerCharge`
- `H.1MHz.scalerTime` → `H.1MHz.scaler`

### 4. File Handling
- **Original**: Single file per run (`skim_run*.root`)
- **New**: Multiple split files per run handled via `TChain`

Production files are named: `nps_production_<run>_<split>_wf_calib.root`

Example: `nps_production_3728_0_wf_calib.root`, `nps_production_3728_1_wf_calib.root`, etc.

## Default Paths

### Input Directory
```cpp
/lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS
```

### Output Directory
```cpp
output/plots/production_wfpi0
```

### Runlist File
```cpp
config/runlist_production.txt
```

## Usage

### Basic Usage
```bash
root -l -b nps_analysis_wfpi0.C
```

This will use default parameters and process all runs listed in `config/runlist_production.txt`.

### Custom Parameters
```cpp
root -l
.L nps_analysis_wfpi0.C
nps_analysis_wfpi0("/path/to/production/files", "output/custom_plots", "config/custom_runlist.txt", 10.538)
```

### Environment Variables
You can also set environment variables to override defaults:
```bash
export NPS_SKIM_DIR="/path/to/production/files"
export NPS_OUTPUT_DIR="output/custom_plots"
export NPS_RUNLIST="config/custom_runlist.txt"
export NPS_EBEAM="10.538"
root -l -b nps_analysis_wfpi0.C
```

## Runlist Format

The runlist file should contain one run number per line:
```
3728
3729
3731
3732
...
```

See `config/runlist_production.txt` for an example.

## Production Data Tree Structure

The production data (`t_prod`) contains the following key branches:

### Event Information
- `g.runnum` (Double_t)
- `g.evnum` (Double_t)

### HMS Tracking
- `H.gtr.dp`, `H.gtr.ph`, `H.gtr.th` (Double_t)
- `H.gtr.px`, `H.gtr.py`, `H.gtr.pz` (Double_t)  
- `H.react.x`, `H.react.y`, `H.react.z` (Double_t)
- `H.react.ok`, `H.gtr.ok` (Double_t)

### HMS PID
- `H.hod.beta` (Double_t)
- `H.cal.etottracknorm`, `H.cal.etracknorm`, `H.cal.etotnorm` (Double_t)
- `H.cer.npeSum` (Double_t)

### Scalers
- `H.1MHz.scaler` (Double_t)
- `H.BCM4A.scaler`, `H.BCM4A.scalerCharge`, `H.BCM4A.scalerCurrent` (Double_t)
- `T.helicity.hel` (Double_t)

### NPS Production Data
- `NPS.prod.nclust`, `NPS.prod.nclustCoin`, `NPS.prod.nclustAcc1`, `NPS.prod.nclustAcc2` (Int_t)
- `NPS.prod.clusE`, `NPS.prod.clusX`, `NPS.prod.clusY`, `NPS.prod.clusZ` (vector<double>)
- `NPS.prod.clusXcorr`, `NPS.prod.clusYcorr` (vector<double>)
- `NPS.prod.clusT`, `NPS.prod.clusDepth`, `NPS.prod.clusSize` (vector<double>/vector<int>)
- `NPS.prod.trk.px`, `NPS.prod.trk.py`, `NPS.prod.trk.pz`, `NPS.prod.trk.ene` (vector<double>)
- `NPS.prod.M`, `NPS.prod.Mx2` (vector<double>)

## Output

The analysis produces:

1. **Per-run ROOT files**: `diagnostics_run<run>.root`
   - Contains all histograms and analysis trees for each run
   
2. **Per-run PNG plots**: Various diagnostic plots saved as PNG
   
3. **Summary CSV**: `summary_all_runs.csv`
   - Contains run-by-run summary statistics
   
4. **Per-run text summaries**: `summary_run<run>.txt`

## Notes

- The code uses `TChain` to automatically combine all split files for each run
- All physics calculations and cuts remain identical to `nps_analysis.C`
- The code maintains backward compatibility with helper functions and utilities
- Memory management uses RAII principles with `std::unique_ptr` for automatic cleanup

## Troubleshooting

### "Tree 't_prod' not found"
- Verify the input files are production files with the `t_prod` tree
- Check that files match the expected naming pattern

### "No files found matching pattern"
- Verify the input directory path is correct
- Check that production files exist for the specified run numbers

### Branch reading errors
- Verify branch names match the production tree structure
- Use `root file.root` and `.ls` or `t_prod->Print()` to inspect the tree

## Contact

For questions or issues, contact Avnish Singh.
