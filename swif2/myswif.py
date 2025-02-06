import json
import sys

def read_txt_file(file_path, job_name):
    data = []
    with open(file_path, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                disk_bytes = int(parts[1]) * 3
                job_info = {
                    "account": "hallc",
                    "command": [
                        f"/u/group/nps/singhav/nps_analysis/swif2/run_myscript.sh SCRIPTS/HMS/CALIBRATION/replay_hodo_coin_NPS_HMS_beta.C {parts[0]} -1"
                    ],
                    # "constraint": "centos79",
                    "cpu_cores": 1,
                    "disk_bytes": disk_bytes,
                    "inputs": [
                        {
                            "local": "nps_replay.tar.gz",
                            "remote": "/group/nps/singhav/nps_replay.tar.gz"
                        },
                        {
                            "local": f"nps_coin_{parts[0]}.dat.0",
                            "remote": f"/mss/hallc/c-nps/raw/nps_coin_{parts[0]}.dat.0"
                        }
                    ],
                    "name": f"{job_name}_{parts[0]}",
                    "partition": "production",
                    "ram_bytes": 2500000000,
                    "stderr": f"/farm_out/singhav/{job_name}_stderr/test_{parts[0]}.err",
                    "stdout": f"/farm_out/singhav/{job_name}_stdout/test_{parts[0]}.out"
                }
                data.append(job_info)
    return data

def generate_json_file(data, output_file, job_name, job_nunber):
    json_data = {
        "jobs": data,
        "name": "singhav_" + job_name + "_" + job_nunber
    }

    with open(output_file, 'w') as json_file:
        json.dump(json_data, json_file, indent=2)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 myswif.py.py job_name job_nunber")
        sys.exit(1)

    job_name = sys.argv[1]
    job_nunber = sys.argv[2]
    runlist_path = "run-lists/runlist_" + sys.argv[2] + ".dat"
    output_json_file = "jsons/job" + sys.argv[2] + ".json"

    job_data = read_txt_file(runlist_path, job_name)
    generate_json_file(job_data, output_json_file, job_name, job_nunber)
