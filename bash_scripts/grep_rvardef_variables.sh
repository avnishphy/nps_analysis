#!/bin/bash

# Function to search for RVarDef patterns in all files within a directory
search_rvardef_in_all_files() {
  local search_directory="$1"
  local dir_name=$(basename "$search_directory")
  echo "===== Searching in directory: $dir_name =====" >> "$output_file"
  find "$search_directory" -type f -name '*.cxx' -print0 | while IFS= read -r -d '' file; do
    local filename=$(basename "$file")
    echo "Searching in $file..."
    echo -e "\n--- File: $filename ---" >> "$output_file"
    grep -Pzo "$pattern" "$file" | sed 's/^/    /' >> "$output_file"
    echo >> "$output_file"  # Add a newline for separation
  done
}

# Define the directories containing the files
primary_directory="/group/nps/apps/NPSlib/fddb8e9/src/src"
secondary_directory="/group/nps/apps/hcana/1.0.2"

# Define the output file
output_file="grep_rvardef_variables_output.txt"

# Clear the output file if it exists
> "$output_file"

# Grep pattern to search for RVarDef structure
pattern='RVarDef\s+\w+\[\]\s*=\s*\{([\s\S]*?)\{\s*0\s*\}\s*\};'

# Search in all files in the primary directory
echo "Searching in the primary directory: $primary_directory"
search_rvardef_in_all_files "$primary_directory"

# Prompt the user whether to search in the secondary directory
read -p "Do you also want to search in the $secondary_directory directory? (y/n): " user_response

if [[ "$user_response" =~ ^[Yy]$ ]]; then
  echo "Searching in the $secondary_directory directory: $secondary_directory"
  search_rvardef_in_all_files "$secondary_directory"
else
  echo "Skipping search in the $secondary_directory directory."
fi

echo "Results have been saved to $output_file"
