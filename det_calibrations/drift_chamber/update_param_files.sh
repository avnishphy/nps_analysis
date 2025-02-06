#!/bin/bash

# Check if directory argument is provided
if [ $# -eq 0 ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

# Assign input directory to a variable
input_dir=$1

# Loop through all files of the form det_calib_runnumber.param in the input directory
for file in "$input_dir"/det_calib_*.param; do
  # Extract the run number from the file name
  run_number=$(basename "$file" | grep -oP '\d+')

  # Check if the file exists
  if [[ -f "$file" ]]; then
    echo "Processing file: $file"

    # Replace the target lines with the ones containing the specific run number
    sed -i.bak \
    -e "s|#include \"PARAM/HMS/DC/hdc_calib_nps23.param\"|#include \"PARAM/HMS/DC/hdc_calib_$run_number.param\"|g" \
    -e "s|#include \"PARAM/HMS/DC/hdc_tzero_per_wire_nps23.param\"|#include \"PARAM/HMS/DC/hdc_tzero_per_wire_$run_number.param\"|g" \
    "$file"

    echo "Updated file: $file"
  else
    echo "File not found: $file"
  fi
done

echo "Replacement complete."
