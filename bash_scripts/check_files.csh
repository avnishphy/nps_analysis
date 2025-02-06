#!/bin/csh

# Set the directory where the files are located
set directory = "/lustre24/expphy/cache/hallc/c-nps/analysis/online/replays/dc_calib_coin"

# Specify the file containing the run numbers
set run_file = "/u/group/nps/singhav/nps_analysis/det_calibrations/runlists/setreftimes_runlist.txt"

# Loop through each run number in the run.txt file
foreach run_number (`cat $run_file`)
    # Define the file pattern(s) for the run number
    set file1 = "${directory}/nps_hms_coin_${run_number}_1_-1.root"
    set file2 = "${directory}/nps_hms_coin_${run_number}_1_-1_1.root"

    # Check if the files exist and display their status
    if ( -f "$file1" ) then
        echo "File found: $file1"
    else
        echo "File not found: $file1"
    endif

    if ( -f "$file2" ) then
        echo "File found: $file2"
    else
        echo "File not found: $file2"
    endif
end
