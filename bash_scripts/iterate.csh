#!/bin/tcsh

#Check if the correct number of arguments is provided
if ($#argv != 1) then
	echo "Usage: $0 <input_dat_file>"
	exit 1
endif

set input_dat_file = $1

# Check if the file exists
if (! -e "$input_dat_file") then
     echo "Error: File not found: $input_dat_file"
         exit 1
endif

set entries = "awk '{print $1}' $input_dat_file"

# Loop through each entry in the first column and run setreftimes.csh
foreach entry ($entries )
	./setreftimes.csh $entry
end
