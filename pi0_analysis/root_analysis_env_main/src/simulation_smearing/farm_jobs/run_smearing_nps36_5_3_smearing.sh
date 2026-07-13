#!/usr/bin/env bash
set -euo pipefail

# Reload this wrapper once inside the Hall C software environment.
if [[ "${NPS_FARM_ENV_LOADED:-0}" != 1 ]]; then
  export NPS_FARM_ENV_LOADED=1
  exec csh -c 'if ( -f /usr/share/Modules/init/csh ) source /usr/share/Modules/init/csh; source /group/nps/singhav/setup.csh; bash "/work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/simulation_smearing/farm_jobs/run_smearing_nps36_5_3_smearing.sh" 0'
fi

export OMP_NUM_THREADS="1"
cd "/work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main"

if [[ "1" == 1 ]]; then
  printf 'y\n' | /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/simulation_smearing/run_smearing_pipeline.sh --kin KinC_x36_5 --target LH2 --combined-file output/KinC_x36_5_yaopeng_08gev/KinC_x36_5/root/combined_branches_LH2_minus_3728_3729.root --nx 6 --ny 6 --x-min -27 --x-max 27 --y-min -33 --y-max 33 --overlap 0.1 --nsmear 80 --smear-seed 42 --beam-energy-gev 10.538 --z-nps-cm 407 --mode auto --nps-angle-deg 12.2 --root-dir output/KinC_x36_5_yaopeng_08gev/KinC_x36_5/root 
else
  /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/src/simulation_smearing/run_smearing_pipeline.sh --kin KinC_x36_5 --target LH2 --combined-file output/KinC_x36_5_yaopeng_08gev/KinC_x36_5/root/combined_branches_LH2_minus_3728_3729.root --nx 6 --ny 6 --x-min -27 --x-max 27 --y-min -33 --y-max 33 --overlap 0.1 --nsmear 80 --smear-seed 42 --beam-energy-gev 10.538 --z-nps-cm 407 --mode auto --nps-angle-deg 12.2 --root-dir output/KinC_x36_5_yaopeng_08gev/KinC_x36_5/root  </dev/null
fi
