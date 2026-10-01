#!/usr/bin/bash

# Check if /etc/profile exists and source it
if [ -f /etc/profile ]; then
    . /etc/profile
fi

# Alternatively, if you only need /etc/profile.d/module.sh
#if [ -f /etc/profile.d/module.sh ]; then
#    . /etc/profile.d/module.sh
#fi

module use /group/halla/modulefiles
module use /group/nps/modulefiles
module load nps_replay
module load panguin
