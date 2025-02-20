#!/bin/bash

# Directory containing the input files
input_dir="."

# Path to the executable
executable="/home/ui330360/software/SusHi-1.6.1/bin/sushi"

# Number of generated files
num_files=20  # Adjust this to match the number of your generated input files

# Loop over each file
for i in $(seq 1 $num_files); do
    # Define input and output filenames based on loop index
    input_file="${input_dir}/scale_scan_${i}.in"
    output_file="${input_dir}/scale_scan_${i}.out"

    # Call the executable with the current input and output files
    echo "Processing: $input_file -> $output_file"
    $executable "$input_file" "$output_file"
done

echo "All files processed."