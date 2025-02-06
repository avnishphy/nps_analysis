#!/bin/bash

# List of run numbers
runs=(1692 5244 5428 2057 4015 5413 5992 6475 1593 6428 4845 3098 5255 6873 2990 3017 3010 6794 2195 6991 3942 5286 6965 6861 3804 2411 4668)

# Base directory for the output files
output_dir="/u/group/nps/singhav/nps_replay/DBASE/HMS"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through each run number and create the corresponding file
for run in "${runs[@]}"; do
    # Define the filename
    file_path="$output_dir/det_calib_${run}.param"
    
    # Write the content to the file
    echo "#include \"PARAM/HMS/DC/hdc_calib_${run}.param\"" > "$file_path"
    echo "#include \"PARAM/HMS/DC/hdc_tzero_per_wire_${run}.param\"" >> "$file_path"
    
    echo "Created $file_path"
done

echo "Script execution completed."
