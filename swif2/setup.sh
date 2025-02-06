#!/usr/bin/bash

#source /apps/modules/5.2.0/init/profile.sh
module use /group/halla/modulefiles
module use /group/nps/modulefiles
module load nps_replay/5.28.24

echo "Modules in use:"
module list
echo "This is the PATH variable:"
echo $PATH
echo "This is LD_LIBRARY_PATH"
echo $LD_LIBRARY_PATH
# Source setup scripts
curdir=`pwd`