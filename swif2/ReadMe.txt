!!! Memory Required: 2442	Memory used: 4.02	Perentage: 0.16%
du -h --max-depth=1

Instructions:
** setup.sh : Set up the environment for running the script
** run_myscript.sh : Execute the main analysis macro
** make_script_runlist.sh: Generate a run number list which is going to be analyzed
** myswif.py: Generate the json file for Slurm (Usage: python3 myswif.py.py job_name job_nunber)

change: run-list

Run list for reference time check:

1753    KinC_x60_3      2566591     ps6=0
1823    KinC_x60_3      1043734     ps3=1

2770    KinC_x60_3'     4490075     ps6=0
2776    KinC_x60_3'     1436714     ps3=2
2834    KinC_x60_3'     1515807     ps4=0

4064    KinC_x60_3a     4042037     ps6=0
4122    KinC_x60_3a     697094      ps4=0

4088    KinC_x60_3b     4332561     ps6=0
4180    KinC_x60_3b     421882      ps4=0

4226    KinC_x60_4a     1480273     ps6=0
4458    KinC_x60_4a     2758804     ps4=0

4254    KinC_x60_4b     1132375     ps6=0
4551    KinC_x60_4b     1025818     ps4=0

4972    KinC_x25_4      1852061     ps6=0
5098    KinC_x25_4      2946995     ps4=3

5540    KinC_x60_1      381877      ps4=0
5894    KinC_x60_1      401478      ps4=0

5904    KinC_x25_1      1446546     ps4=0
5974    KinC_x25_1      1405862     ps4=0

6190    KinC_x60_2b     902166      ps4=0
6402    KinC_x60_2b     1009556     ps4=0

6571    KinC_x25_3      2472220     ps4=0
6687    KinC_x25_3      2420730     ps4=0

3389    KinC_x60_2      1310704     ps6=0
3635    KinC_x60_2      315840      ps4=0

3430    KinC_x60_2'     1670411     ps6=0
3654    KinC_x60_2'     522809      ps4=0




KinC_x25_2

KinC_x25_4


Steps:
1.
cd nps_replay/ && tar -czf ../nps_replay.tar.gz . && cd -
cd swif2
2.
./make_script_runlist.sh 6810 6810 run-lists/runlist_1.dat
./make_script_runlist.sh run-lists/runlist_2.dat

3.
python3 myswif.py RefTime 1
python3 myswif.py RefTime 2

4.
swif2 import -file jsons/job1.json
swif2 import -file jsons/job2.json
swif2 list

5.
swif2 run Yaopeng_RefTime_1
swif2 status Yaopeng_RefTime_1
swif2 cancel Yaopeng_RefTime_1 --delete

swif2 run Yaopeng_RefTime_2
swif2 status Yaopeng_RefTime_2
swif2 cancel Yaopeng_RefTime_2 --delete




swif2 list
swif2 status Yaopeng_DeadTime_1
swif2 status Yaopeng_DeadTime_1 -problems
swif2 retry-jobs Yaopeng_DeadTime_1 -problems SLURM_NODE_FAIL 
swif2 cancel Yaopeng_DeadTime_1 --delete

1   problems = 7, 69
SLURM_NODE_FAIL ,SLURM_FAILED

2   problems = 1, 1, 65
SITE_LAUNCH_FAIL,SLURM_CANCELLED,SLURM_FAILED

3   problems = 3, 3, 109
SITE_LAUNCH_FAIL,SLURM_NODE_FAIL,SLURM_FAILED

4   problems = 4, 58
SLURM_NODE_FAIL,SLURM_FAILED

5   problems = 4, 14
SLURM_NODE_FAIL,SLURM_FAILED

6   problems = 1, 3, 30
SITE_LAUNCH_FAIL,SLURM_NODE_FAIL,SLURM_FAILED

7   problems = 9, 48
SLURM_NODE_FAIL,SLURM_FAILED


403 bad runs in 3000 runs