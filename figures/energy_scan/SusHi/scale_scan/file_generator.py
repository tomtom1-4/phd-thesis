import os

def modify_and_save_file(original_file, output_dir, num_files):
    # Ensure output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Calculate step size for values between 0.05 and 3
    start_value = 0.05
    end_value = 3.0
    step_size = (end_value - start_value) / (num_files - 1)

    with open(original_file, 'r') as file:
        lines = file.readlines()

    for i in range(num_files):
        value = start_value + i * step_size
        # Modify the 39th line (index 38)
        modified_line = f"  2   {value:.2f} \t# factorization scale muF/mh\n"
        lines[38] = modified_line

        # Create new filename
        new_filename = f"{output_dir}/scale_scan_{i+1}.in"

        # Write modified content to a new file
        with open(new_filename, 'w') as new_file:
            new_file.writelines(lines)

        print(f"Created: {new_filename}")

# Usage example:
original_file_path = 'scale_scan.in'
output_directory = '.'
number_of_files_to_generate = 20

modify_and_save_file(original_file_path, output_directory, number_of_files_to_generate)