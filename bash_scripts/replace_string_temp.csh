#!/bin/tcsh

source /group/nps/singhav/setup.csh

set filename = "$REPLAYDIR/SCRIPTS/NPS/replay_production_coin_NPS_HMS.C"

# Escape forward slashes in the strings
# set old_string_escaped = "ROOTFileNamePattern = \"ROOTfiles/COIN/PRODUCTION/nps_hms_coin_%d_%d_%d.root\""
# set new_string_escaped = "ROOTFileNamePattern = \"ROOTfiles/COIN/PRODUCTION/nps_hms_coin_%d_%d_%d_hodo_dc_uncalibrated.root\""

set old_string_escaped = "nps_hms_coin_%d_%d_%d.root"
set new_string_escaped = "nps_hms_coin_%d_%d_%d_hodo_dc_uncalibrated.root"

# Replace the old string with the new string in the file
sed -i "s|$old_string_escaped|$new_string_escaped|g" "$filename"
