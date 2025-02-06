#!/bin/bash

# File paths
database_file="/u/group/nps/singhav/nps_replay/DBASE/NPS/standard_coin.database"
runlist_file="/u/group/nps/singhav/nps_analysis/det_calibrations/runlists/setreftimes_runlist.txt"

# Check if both files exist
if [ ! -f "$database_file" ]; then
    echo "Database file does not exist."
    exit 1
fi

if [ ! -f "$runlist_file" ]; then
    echo "Runlist file does not exist."
    exit 1
fi

# Read run numbers from the runlist and append the required lines to the database file
while IFS= read -r runnumber; do
    if [ -n "$runnumber" ]; then
        echo "Processing run number: $runnumber"
        
        # Append the required lines to the database file
        echo -e "\n$runnumber - $runnumber" >> "$database_file"
        echo "g_ctp_det_calib_filename   = \"DBASE/HMS/det_calib_${runnumber}.param\"" >> "$database_file"
    fi
done < "$runlist_file"

echo "Script execution completed. Lines added to the database file."
