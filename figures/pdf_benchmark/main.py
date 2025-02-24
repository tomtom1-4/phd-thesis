import numpy as np
import matplotlib.pyplot as plt
import glob

import csv
import os
import random as rdm
import sys

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

def extract_value(line):
  """
  Extracts the numerical value before the square bracket '[' in a line.
  """
  # Split the line at '='
  parts = line.split('=')
  if len(parts) > 1:
    rhs = parts[1]
    # Split the right-hand side at '[' to exclude uncertainties
    value_part = rhs.split('[')[0]
    # Remove any whitespace
    value_str = value_part.strip()
    # Convert to float
    value = float(value_str)
    return value
  else:
    # Return None if the line doesn't contain '='
    return None

def create_latex_table(filenames, central_values, up_errors, down_errors, normalized=True, output_file=None):
  """
  Generate a LaTeX table from the given data.

  Parameters:
  - filenames: List of filenames or labels.
  - central_values: Array of central values.
  - up_errors: Array of upward uncertainties.
  - down_errors: Array of downward uncertainties.
  - normalized: Boolean indicating if the values are normalized.
  - output_file: Optional path to save the LaTeX table to a file.

  Returns:
  - tex_table: String containing the LaTeX table.
  """
  # Begin the LaTeX table
  tex_table = "\\begin{table}[h!]\n"
  tex_table += "\\centering\n"
  if normalized:
    tex_table += "\\caption{Normalized Central Values with Uncertainties}\n"
    tex_table += "\\label{tab:normalized_values}\n"
  else:
    tex_table += "\\caption{Central Values with Uncertainties}\n"
    tex_table += "\\label{tab:central_values}\n"
  tex_table += "\\begin{tabular}{l c}\n"
  tex_table += "\\hline\n"
  tex_table += "File Name &; Central Value \\\\\n"
  tex_table += "\\hline\n"

  # Format each row
  for fname, value, err_up, err_down in zip(filenames, central_values, up_errors, down_errors):
    # Format the uncertainties
    # Use asymmetric uncertainties if they are different
    if np.isclose(err_up, err_down):
      uncertainty = f"\\pm {err_up:.2g}"
    else:
      uncertainty = f"^{{{{+{err_up:.2g}}}}}_{{{{-{err_down:.2g}}}}}"
    # Format the central value with uncertainties
    # For LaTeX, we might want to use \num{} or \SI{}{} from siunitx if desired
    formatted_value = f"${value:.4g} {uncertainty}$"
    # Extract base filename or use a label
    label = fname  # Use fname or modify as needed
    tex_table += f"{label} &; {formatted_value} \\\\\n"
  tex_table += "\\hline\n"
  tex_table += "\\end{tabular}\n"
  tex_table += "\\end{table}\n"

  # Save to file if an output path is provided
  if output_file:
    with open(output_file, 'w') as f:
      f.write(tex_table)
    print(f"LaTeX table saved to {output_file}")

  return tex_table

def main():
  with open('../style.py', 'r') as f:
    s = f.read()
    plot = eval(s)
  # Collect all input files (adjust the pattern according to your files)
  #files = glob.glob('ihixs/*.out')  # Replace 'input_files/*.txt' with your file path pattern
  files = ['ihixs/NNPDF31_nnlo_as_0118.out', 'ihixs/NNPDF40_an3lo_as_01180.out', 'ihixs/NNPDF40_an3lo_as_01180_mhou.out', 'ihixs/NNPDF40_an3lo_as_01180_pdfas.out',
           'ihixs/NNPDF40_an3lo_as_01180_qed.out', 'ihixs/MSHT20an3lo_as118.out', 'ihixs/MSHT20xNNPDF40_aN3LO.out', 'ihixs/MSHT20xNNPDF40_aN3LO_qed.out']
  filenames = []
  central_values = []
  up_errors = []
  down_errors = []

  for filename in files:
    filenames.append("$\mathtt{" + filename[6:-4].replace("_", "\_") + "}$")
    with open(filename, 'r') as f:
      lines = f.readlines()
      # Ensure the file has enough lines
      if len(lines) >= 65:
        # Extract central value from line 54 (index 53)
        central_line = lines[53]
        central_value = extract_value(central_line)
        central_values.append(central_value)
        # Extract upward error from line 64 (index 63)
        up_error_line = lines[63]
        up_error = extract_value(up_error_line)
        up_errors.append(up_error)
        # Extract downward error from line 65 (index 64)
        down_error_line = lines[64]
        down_error = extract_value(down_error_line)
        down_errors.append(down_error)
      else:
        print(f"File {filename} does not have enough lines.")
        central_values.append(np.nan)
        up_errors.append(0)
        down_errors.append(0)

  # Convert lists to numpy arrays for plotting
  central_values = np.array(central_values)
  up_errors = np.array(up_errors)
  down_errors = np.abs(np.array(down_errors))


  # Extract reference values
  print(filenames)
  ref_index = filenames.index("$\\mathtt{NNPDF31\\_nnlo\\_as\\_0118}$")
  ref_central = central_values[ref_index]
  ref_up_error = up_errors[ref_index]
  ref_down_error = down_errors[ref_index]

  # Create x positions for categorical x-axis
  x = np.arange(len(filenames))

  # Prepare asymmetric error bars
  yerr = np.array([down_errors/ref_central, up_errors/ref_central])

  # Plotting
  fig, ax = plt.subplots(figsize=(plot['singleplot']['width'], plot['singleplot']['height']))
  colors = [plot['color'][0], plot['color'][1], plot['color'][1], plot['color'][1], plot['color'][1], plot['color'][2], plot['color'][3], plot['color'][3]]

  # Plot central values with error bars
  for i in range(0, len(central_values)):
    ax.errorbar([x[i]], [central_values[i]/ref_central], yerr=[[yerr[0][i]], [yerr[1][i]]], fmt='o', capsize=5, color=colors[i])

  # Set x-axis labels to filenames
  ax.set_xticks(x)
  ax.set_xticklabels(filenames, rotation=45, ha='right', rotation_mode='anchor')

  # Set axis labels and title
  ax.set_ylabel("$\sigma/\sigma (\mathtt{NNPDF31\_nnlo\_as\_0118})$")
  ax.text(5, 0.985, '$\mathrm{LHC@13TeV}$\n$\mu_R=\mu_F=m_H/2$')

  # Adjust layout to prevent label cutoff
  plt.tight_layout()

  plt.grid()

  # Show the plot
  plt.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/pdf_benchmark.pdf")

  # Call the function to create the LaTeX table
  latex_table = create_latex_table(
      filenames=filenames,  # Use the labels or filenames as desired
      central_values=central_values,
      up_errors=up_errors,
      down_errors=down_errors,
      normalized=False,  # Set to False if not normalized
      output_file='pdf_benchmark.tex'  # Optional: specify output file path
  )

  # Print the LaTeX table (optional)
  print(latex_table)


if __name__=="__main__":
  main()