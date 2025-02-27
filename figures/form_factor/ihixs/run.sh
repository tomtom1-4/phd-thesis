#!/bin/bash

# Loop through all .card files in the current directory
for file in *.card; do
  # Check if there are no .card files (prevents running if none exist)
  [ -e "$file" ] || continue

  # Extract the filename without the extension
  base_name="${file%.card}"

  # Run the ihixs command with input and output
  /home/tom/Documents/software/software/ihixs/ihixs/build/ihixs -i "$file" -o "${base_name}.out"
done

echo "Processing complete!"
