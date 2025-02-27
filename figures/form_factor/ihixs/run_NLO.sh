#!/bin/bash

# Define the combinations for fi and fj
fi_fj_list=("t t" "t b" "t c" "b b" "b c" "c c")

# Define the combinations for mur and muf
mur_muf_list=("1 1" "0.5 1" "0.5 0.5" "1 0.5" "2 1" "2 2" "1 2")

# Loop over each fi and fj combination
for fi_fj in "${fi_fj_list[@]}"; do
    # Extract fi and fj from the combination
    read fi fj <<< "$fi_fj"

    # Set y_top to 1.0 if either fi or fj is 't', else 0.0
    if [ "$fi" == "t" ] || [ "$fj" == "t" ]; then
        y_top=1.0
    else
        y_top=0.0
    fi

    # Set y_bottom to 1.0 if either fi or fj is 'b', else 0.0
    if [ "$fi" == "b" ] || [ "$fj" == "b" ]; then
        y_bottom=1.0
    else
        y_bottom=0.0
    fi

    # Set y_charm to 1.0 if either fi or fj is 'c', else 0.0
    if [ "$fi" == "c" ] || [ "$fj" == "c" ]; then
        y_charm=1.0
    else
        y_charm=0.0
    fi

    # Loop over each mur and muf combination
    for mur_muf in "${mur_muf_list[@]}"; do
        # Extract mur and muf from the combination
        read mur muf <<< "$mur_muf"


        # Calculate mur_value and muf_value by multiplying by 62.5
        mur_value=$(echo "$mur * 62.5" | bc -l)
        muf_value=$(echo "$muf * 62.5" | bc -l)

        # Replace dots with underscores for file naming
        mur_safe=$(echo "$mur" | sed 's/\.//g')
        muf_safe=$(echo "$muf" | sed 's/\.//g')

        # Construct the output file name
        output_file="${fi}${fj}_${mur_safe}_${muf_safe}.out"

        # Execute the command with the specified options
        /home/tom/Documents/software/software/ihixs/ihixs/build/ihixs \
            -i NLO.card \
            --y_top $y_top \
            --y_bot $y_bottom \
            --y_charm $y_charm \
            --mur $mur_value \
            --muf $muf_value \
            -o $output_file
    done
done