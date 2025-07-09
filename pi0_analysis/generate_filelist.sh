#!/bin/bash

input_dir="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated"
output_file="filelist.txt"

# Remove old list
rm -f "$output_file"

# Loop through matching files and add to filelist
find "$input_dir" -name "nps_hms_coin_*_*_1_-1.root" | sort > "$output_file"

echo "Created $output_file with $(wc -l < $output_file) entries."
