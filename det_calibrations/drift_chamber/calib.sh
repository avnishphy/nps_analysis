#!/bin/bash

# Check if a file with run numbers is provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 <run_numbers_file>"
  exit 1
fi

# Assign the first argument to a variable
run_numbers_file=$1

# Check if the run numbers file exists and is readable
if [[ ! -f "$run_numbers_file" || ! -r "$run_numbers_file" ]]; then
  echo "Error: Cannot read file '$run_numbers_file'"
  exit 1
fi

# Loop through each line in the file (each run number)
while IFS= read -r run_number || [[ -n "$run_number" ]]; do
  # Trim leading/trailing whitespace
  run_number=$(echo "$run_number" | xargs)

  # Skip empty lines or lines starting with a comment character (#)
  if [[ -n "$run_number" && ! "$run_number" =~ ^# ]]; then
    echo "Executing hcana for run number: $run_number"
    hcana -l "/u/group/nps/singhav/nps_replay/CALIBRATION/dc_calib/scripts/main_calib.C($run_number)"
  fi
done < "$run_numbers_file"
