#!/bin/bash

number=$1

input_file="./run-lists/runlist_$number.dat"

output_file="./run-lists/runlist_failed_$number.dat"

folder_path="../rootfiles/"

while IFS= read -r line; do

    number=$(echo "$line" | awk '{print $1}')

    if [[ ! -e "$folder_path/DeadTime_$number.root" ]]; then
        echo "$line" >> "$output_file"
    fi
done < "$input_file"
