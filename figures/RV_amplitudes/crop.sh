#!/bin/bash

# Define the output directory
output_dir="../../LaTeX/Images/RV_amplitudes"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop over all PDF files in the current directory
for file in *.pdf; do
    # Check if the current item is a file
    if [ -f "$file" ]; then
        echo "Processing '$file'..."

        # Construct the output file path
        output_file="$output_dir/$file"

        # Execute the crop command
        python crop.py "$file" "$output_file"

        echo "Saved cropped PDF to '$output_file'"
    fi
done

echo "All PDF files have been processed."