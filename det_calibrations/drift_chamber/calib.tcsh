#!/bin/tcsh

# Check if a file with run numbers is provided
if ( $#argv == 0 ) then
  echo "Usage: $0 <run_numbers_file>"
  exit 1
endif

# Assign the first argument to a variable
set run_numbers_file = $argv[1]

# Loop through each line in the file (each run number)
foreach run_number ( `cat $run_numbers_file` )
  if ( "$run_number" != "" ) then
    echo "Executing hcana for run number: $run_number"
    hcana -l "/u/group/nps/singhav/nps_replay/CALIBRATION/dc_calib/scripts/main_calib.C($run_number)"
  endif
end
