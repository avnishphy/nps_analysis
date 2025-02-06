#!/bin/tcsh

# Source the setup script
source /group/nps/singhav/setup.csh

# Change to the replay directory
cd $REPLAYDIR

# Define the path to the runlist file
set runlist = "$REPLAYDIR/../nps_analysis/bash_scripts/setreftimes_runlist.txt"

# Loop over each run number in the runlist
# foreach run (`cat $runlist`)
#     # Run the root macro with the run number and 1 million events
#     hcana -b -q "SCRIPTS/NPS/TIMING/no_reference_times_nps_hms.C(${run}, 1000000)"

# end

# Verify that all output files have been created
foreach run (`cat $runlist`)
    set outfile = "$REPLAYDIR/ROOTfiles/NPS/TIMING/nps_hms_noReferenceTime_${run}_1000000.root"
    if (! -e $outfile) then
        echo "Error: Output file $outfile does not exist."
        exit 1
    endif
end

# Change back to the replay directory
cd $REPLAYDIR

# Run the second root script for each run
foreach run (`cat $runlist`)
    set infile = "$REPLAYDIR/ROOTfiles/NPS/TIMING/nps_hms_noReferenceTime_${run}_1000000.root"
    
    # Run the root script with the function and arguments
    # hcana -b -q "CALIBRATION/set_reftimes/reference_time_setup.C(\"run_hms_reference_time_setup('$infile', $run)\")"
    hcana .L CALIBRATION/set_reftimes/reference_time_setup.C << EOF
    run_hms_reference_time_setup("$infile", $run, "$REPLAYDIR/../nps_analysis/det_calibrations/reftimes_14400/move_me_${run}.root")
EOF

end
