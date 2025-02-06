#!/bin/tcsh

# Step 1: Source the setup script
source /group/nps/singhav/setup.csh

# Step 2: Change to the REPLAYDIR directory
cd $REPLAYDIR

# Step 5: Run the replay script with specified parameters
set run_list_file = "$REPLAYDIR/nps_analysis/bash_scripts/temp_runlist.txt"
set run_list = (`cat $run_list_file`)
set num_events = 1000
set replay_script = "$REPLAYDIR/SCRIPTS/NPS/replay_production_coin_NPS_HMS.C"

foreach run_num ($run_list)

    # Step 3: Check that the flag htofusinginvadc = 0
    set flag_file1 = "$REPLAYDIR/PARAM/HMS/HODO/hhodo_cuts_nps23.param"
    set htofusinginvadc = `grep -E '^htofusinginvadc' $flag_file1 | awk -F '=' '{print $2}'`

    if ("$htofusinginvadc" != "0") then
        echo "Error: htofusinginvadc is not set to 0 for run number $run_num. Exiting."
        exit 1
    endif

    # Step 4: Check and modify the flag h_using_tzero_per_wire
    set flag_file2 = "$REPLAYDIR/PARAM/HMS/DC/hdc_cuts_nps23.param"
    set h_using_tzero_per_wire = `grep -i 'h_using_tzero_per_wire' $flag_file2 | awk '{print $3}'`

    if ("$h_using_tzero_per_wire" == "1") then
        sed -i 's/h_using_tzero_per_wire = 1/h_using_tzero_per_wire = 0/' $flag_file2
    endif

   
    # change the rootfilenamepattern dynamically
    set old_string = "nps_hms_coin_%d_%d_%d.root"
    set new_string = "nps_hms_coin_%d_%d_%d_hodo_dc_uncalibrated.root"
    # Replace the old string with the new string in the file
    sed -i "s|$old_string|$new_string|g" "$replay_script"
    
    

    # Run the replay script
    hcana -l -b -q "$REPLAYDIR/SCRIPTS/NPS/replay_production_coin_NPS_HMS.C($run_num, $num_events)"

    # Change directory to hms_hodo_calib
    cd $REPLAYDIR/CALIBRATION/hms_hodo_calib

    
    # Step 6: Run the timeWalkHistos.C script
    set input_rootfile = "$REPLAYDIR/ROOTfiles/COIN/PRODUCTION/nps_hms_coin_${run_num}_${num_events}_1_hodo_dc_uncalibrated.root"
    hcana -l << EOF
    .x timeWalkHistos.C(\"$input_rootfile\", $run_num)
EOF
    
    # hcana -l -b -q "timeWalkHistos.C(\"$input_rootfile\", $run_num)"


    # Step 7: Run the timeWalkCalib.C script
    hcana -l -b -q "timeWalkCalib.C($run_num)"

    # Step 8: Rename the parameter file
    cd $REPLAYDIR/PARAM/HMS/HODO
    if (-e "hhodo_TWcalib_nps23.param") then
        mv hhodo_TWcalib_nps23.param hhodo_TWcalib_nps23_original.param
    endif
    mv hhodo_TWcalib_${run_num}.param hhodo_TWcalib_nps23.param

    # Step 9: Replay the raw data with updated parameters
    
    # change the rootfilenamepattern dynamically
    set old_string = "nps_hms_coin_%d_%d_%d_hodo_dc_uncalibrated.root"
    set new_string = "nps_hms_coin_%d_%d_%d_twcalibed.root"
    # Replace the old string with the new string in the file
    sed -i "s|$old_string|$new_string|g" "$replay_script"
    
    hcana -l -b -q "$REPLAYDIR/SCRIPTS/NPS/replay_production_coin_NPS_HMS.C($run_num, $num_events)"

    # Step 10: Run the fitHodoCalib.C script
    cd $REPLAYDIR/CALIBRATION/hms_hodo_calib
    set input_rootfile = "$REPLAYDIR/ROOTfiles/COIN/PRODUCTION/nps_hms_coin_${run_num}_${num_events}_1_twcalibed.root"
    root -l -b -q "fitHodoCalib.C(\"$input_rootfile\", $run_num)"

    # Step 11: Rename the Vpcalib parameter file
    cd $REPLAYDIR/PARAM/HMS/HODO
    if (-e "hhodo_Vpcalib_nps23.param") then
        mv hhodo_Vpcalib_nps23.param hhodo_Vpcalib_nps23_original.param
    endif
    mv hhodo_Vpcalib_${run_num}.param hhodo_Vpcalib_nps23.param

    # change the rootfilenamepattern dynamically
    set old_string = "nps_hms_coin_%d_%d_%d_twcalibed.root"
    set new_string = "nps_hms_coin_%d_%d_%d_vpcalibed.root"
    # Replace the old string with the new string in the file
    sed -i "s|$old_string|$new_string|g" "$replay_script"

    # Step 12: Replay the raw data one last time
    hcana -l -b -q "$REPLAYDIR/SCRIPTS/NPS/replay_production_coin_NPS_HMS.C($run_num, $num_events)"

    # Step 13: Restore all the name changes
    mv $REPLAYDIR/PARAM/HMS/HODO/hhodo_Vpcalib_nps23.param $REPLAYDIR/PARAM/HMS/HODO/hhodo_Vpcalib_${run_num}.param
    mv $REPLAYDIR/PARAM/HMS/HODO/hhodo_TWcalib_nps23.param $REPLAYDIR/PARAM/HMS/HODO/hhodo_TWcalib_${run_num}.param

end

# Step 14: Restore param *_nps23 param files as they were at the start.
mv $REPLAYDIR/PARAM/HMS/HODO/hhodo_TWcalib_nps23_original.param $REPLAYDIR/PARAM/HMS/HODO/hhodo_TWcalib_nps23.param
mv $REPLAYDIR/PARAM/HMS/HODO/hhodo_Vpcalib_nps23_original.param $REPLAYDIR/PARAM/HMS/HODO/hhodo_Vpcalib_nps23.param

# Step 15: Restore the ROOTFileNamePattern in the replay script.
# change the rootfilenamepattern dynamically
set old_string = "nps_hms_coin_%d_%d_%d_vpcalibed.root"
set new_string = "nps_hms_coin_%d_%d_%d.root"
# Replace the old string with the new string in the file
sed -i "s|$old_string|$new_string|g" "$replay_script"

echo "Script execution completed."
