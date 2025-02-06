#!/bin/bash

# Specified directory to search for HMS_DC_cardLog_runnumber directories
search_dir="/u/group/nps/singhav/nps_replay/CALIBRATION/dc_calib/scripts/"  # Change this to your actual search directory

# Target directory where the files will be copied
target_dir="/u/group/nps/singhav/nps_replay/PARAM/HMS/DC/"

# Check if the search directory exists
if [ ! -d "$search_dir" ]; then
    echo "The specified search directory does not exist."
    exit 1
fi

# Check if the target directory exists
if [ ! -d "$target_dir" ]; then
    echo "The target directory does not exist. Creating it..."
    mkdir -p "$target_dir"
fi

# Find directories created on October 18th
dirs=$(find "$search_dir" -type d -name "HMS_DC_cardLog_*" -newermt "2024-10-18" ! -newermt "2024-10-19")

# Loop through directories of type HMS_DC_cardLog_runnumber
# for dir in "$search_dir"/HMS_DC_cardLog_*; do
for dir in $dirs; do
    if [ -d "$dir" ]; then
        # Extract the run number from the directory name
        runnumber=$(basename "$dir" | cut -d'_' -f4)

        # Define the filenames to copy
        calib_file="$dir/hdc_calib_${runnumber}.param"
        tzero_file="$dir/hdc_tzero_per_wire_${runnumber}.param"

        # Check if the files exist before copying
        if [ -f "$calib_file" ] && [ -f "$tzero_file" ]; then
            echo "Copying files for run number $runnumber..."
            cp "$calib_file" "$target_dir"
            cp "$tzero_file" "$target_dir"
        else
            echo "Files not found for run number $runnumber. Skipping..."
        fi
    else
        echo "No valid directory found for $dir"
    fi
done

echo "Script execution completed."
