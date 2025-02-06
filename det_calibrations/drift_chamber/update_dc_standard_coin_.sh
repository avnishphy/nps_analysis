#!/bin/bash

# List of run numbers
runs=(1692 5244 5428 2057 4015 5413 5992 6475 1593 6428 4845 3098 5255 6873 2990 3017 3010 6794 2195 6991 3942 5286 6965 6861 3804 2411 4668)

# Path to the configuration file
config_file="/u/group/nps/singhav/nps_replay/DBASE/NPS/standard_coin.database"

# Check if the configuration file exists
if [ ! -f "$config_file" ]; then
    echo "The configuration file does not exist. Please check the path."
    exit 1
fi

# Function to append entries to the config file
add_config_entry() {
    run=$1

    echo "" >> "$config_file"  # Add a newline for separation
    echo "$run-$run" >> "$config_file"
    echo "g_ctp_det_calib_filename   = \"DBASE/HMS/det_calib_${run}.param\"" >> "$config_file"
}

# Loop through each run number to add new entries to the configuration file
for run in "${runs[@]}"; do
    add_config_entry "$run"
done

echo "Script execution completed. Configuration file updated."
