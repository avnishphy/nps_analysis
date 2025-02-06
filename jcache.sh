#!/usr/bin/bash


#arguments pass to this bash script
#Use: on terminal type (for example) >> ./jcacheGet.sh shms 1791

#User input
# spec=$1    #hms, shms or coin
# run=$2


#for run in {1149..1171}
#do
    

#Use: on terminal type (for example) >> jcacheGet.sh shms 1791
#spec=$1 
#runNUM=$2

#for run in {1149..1171}
#do

#filename='DEUTERON_ANALYSIS/hcswif/runlists/d2_full.dat'


#for run in $(cat $filename) ; do    
# for run in {3368..3379} ; do    
# # echo ${run}
# #mss="/mss/hallc/spring17/raw/${spec}_all_0${run}.dat"
# #mss="/mss/hallc/c-polhe3/raw/${spec}_all_${run}.dat"                                
# #mss="/mss/hallc/c-cafe-2022/raw/${spec}_all_${run}.dat"
#     mss="/mss/hallc/c-nps/raw/${spec}_all_${run}.dat"
# #mss="/mss/hallc/c-pionlt/raw/${spec}_all_${run}.dat" 
# #mss="/mss/hallc/xem2/raw/${spec}_all_${run}.dat"

#     jcacheCMD="jcache get ${mss} -e singhav@cua.edu -x"

#     echo "Executing command: ${jcacheCMD}" 
#     eval ${jcacheCMD}    

for run in 3013 3017 2023 2057 6991 2189 2195 6932 3991 4015 4668 4845 5992 6794 6428 6475 2976 2990 3010 3063 3098 1598 1593 1692 6965 3796 3804 3823 3828 5428 6861 5244 5255 5413 5286 2340 2411 6921 3913 3942 3162 3302 6873  ; do
    for file in /mss/hallc/c-nps/raw/nps_coin_${run}.dat.* ; do
        if [ -e "$file" ]; then
            jcacheCMD="jcache get ${file} -e as8oct@gmail.com -x"

            echo "Executing command: ${jcacheCMD}" 
            eval ${jcacheCMD}
        else
            echo "No files found for run ${run}"
        fi
    done
done

#done

  
