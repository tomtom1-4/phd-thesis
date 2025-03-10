import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from scipy.interpolate import griddata

import csv
import os
import random as rdm
import sys

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

# Define the symbolic variables
nl_sym, nt_sym, nb_sym = sp.symbols('nl nt nb')

# Function to get the total number of lines in a file
def get_total_lines(filename):
  with open(filename, 'r') as f:
    for i, _ in enumerate(f):
      pass
  return i + 1  # Since 'i' is zero-based

# Function to read values and their original indices
def read_values_with_indices(filename, skip=1, limits=None):
  values = []
  indices = []
  with open(filename, 'r') as f:
    for idx, line in enumerate(f):
      if idx % skip == 0:
        try:
          value = float(line.strip())
          if limits is None:
            values.append(value)
            indices.append(idx)
          elif( value > limits[0] and value < limits[1]):
            values.append(value)
            indices.append(idx)
        except ValueError:
          print(f"Could not convert line {idx} to float: {line.strip()}")
  return values, indices

# Updated function to read and evaluate z-values
def read_and_evaluate_z(filename, nl_value, nt_value, nb_value, N_x, N_y, x_indices_set, y_indices_set):
  z_values = []
  with open(filename, 'r') as f:
    x_index = 0
    line_num = 0  # Line number in z file
    while x_index < N_x:
      if x_index in x_indices_set:
        y_index = 0
        while y_index < N_y:
          line = f.readline()
          if not line:
            break  # EOF
          line_num += 1
          if y_index in y_indices_set:
            expr_str = line.strip()
            # Check if 'NAN' is in expr_str
            if 'NAN' in expr_str.upper():
              print(f"'NAN' detected at line {line_num} (x_index={x_index}, y_index={y_index})")
              z_values.append(np.nan)
            else:
              try:
                expr = sp.sympify(expr_str)
                expr_sub = expr.subs({nl_sym: nl_value, nt_sym: nt_value, nb_sym: nb_value})
                z = float(expr_sub.evalf())
                z_values.append(z)
              except Exception as e:
                print(f"Error evaluating line {line_num} at x_index={x_index}, y_index={y_index}: {e}")
                z_values.append(np.nan)
          else:
            # We've already read the line; no need to read it again
            pass  # Do nothing
          y_index += 1
        x_index += 1
      else:
        # Skip over N_y lines since we're skipping this x_index
        for _ in range(N_y):
          line = f.readline()
          if not line:
            break  # EOF
          line_num += 1
        x_index += 1
  return z_values

# Interpolate to fill NaN values
def fill_nan_nearest(Z, X, Y):
  # Flatten the arrays
  X_flat = X.flatten()
  Y_flat = Y.flatten()
  Z_flat = Z.flatten()

  # Create mask of valid values
  valid_mask = ~np.isnan(Z_flat)

  # Points with valid data
  points_valid = np.column_stack((X_flat[valid_mask], Y_flat[valid_mask]))
  values_valid = Z_flat[valid_mask]

  # Points to interpolate
  invalid_mask = np.isnan(Z_flat)
  points_invalid = np.column_stack((X_flat[invalid_mask], Y_flat[invalid_mask]))

  # Interpolate using 'linear' method
  Z_interp_values = griddata(points_valid, values_valid, points_invalid, method='linear')

  # Where interpolation result is still NaN (e.g., outside convex hull), use 'nearest'
  nan_interp_mask = np.isnan(Z_interp_values)
  if np.any(nan_interp_mask):
    Z_interp_values_nearest = griddata(points_valid, values_valid, points_invalid[nan_interp_mask], method='nearest')
    Z_interp_values[nan_interp_mask] = Z_interp_values_nearest

  # Replace NaNs with interpolated values
  Z_flat[invalid_mask] = Z_interp_values

  # Reshape
  Z_filled = Z_flat.reshape(Z.shape)

  return Z_filled

def plot_amp(x_filename, y_filename, z_filename, out_filename, ct_filename=None,
             nl_value=5, nt_value=1, nb_value=0, subtract_ct=False,
             x_lim=[0.02, 1], y_lim=[0.5, 0.98], z_label="$\mathrm{Amp}$", skip_x=10, skip_y=10, elev=10):
  # Get the total number of x and y values before skipping
  N_x = get_total_lines(x_filename)
  N_y = get_total_lines(y_filename)

  # Read the x and y values along with their original indices
  x_values, x_indices = read_values_with_indices(x_filename, skip=skip_x, limits=x_lim)
  y_values, y_indices = read_values_with_indices(y_filename, skip=skip_y, limits=y_lim)

  # Convert indices to sets for faster lookup
  x_indices_set = set(x_indices)
  y_indices_set = set(y_indices)

  # Generate the meshgrid
  X, Y = np.meshgrid(x_values, y_values, indexing='ij')

  # Read and evaluate z-values
  z_values = read_and_evaluate_z(
    z_filename,
    nl_value=nl_value,
    nt_value=nt_value,
    nb_value=nb_value,
    N_x=N_x,
    N_y=N_y,
    x_indices_set=x_indices_set,
    y_indices_set=y_indices_set
  )
  # Convert z_values to a NumPy array and reshape
  Z = np.array(z_values).reshape(len(x_values), len(y_values))

  if subtract_ct:
    ct_values = read_and_evaluate_z(
      ct_filename,
      nl_value=nl_value,
      nt_value=nt_value,
      nb_value=nb_value,
      N_x=N_x,
      N_y=N_y,
      x_indices_set=x_indices_set,
      y_indices_set=y_indices_set
    )
    CT = np.array(ct_values).reshape(len(x_values), len(y_values))
    Z = Z - CT


  #Z[Z>1.e5]= np.nan
  #Z[Z<-1.e5]= np.nan
  # Create a masked array where NaNs are masked
  #Z_masked = np.ma.masked_invalid(Z)

  # Apply the function to fill NaNs
  Z_filled = fill_nan_nearest(Z, X, Y)

  # Now, plot using Z_filled
  fig = plt.figure(figsize=(7,6))
  ax = fig.add_subplot(111, projection='3d')

  # Plot the surface
  surf = ax.plot_surface(X, Y, Z_filled, cmap='viridis', alpha=0.8, antialiased=True, edgecolor='gray', linewidth=0.25, rcount=100, ccount=50)

  # Add a color bar
  #fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
  ax.set_xlim3d(0, 1)
  ax.set_ylim3d(0.5, 1)
  if(np.min(Z_filled) > 0):
    ax.set_zlim3d(0, np.max(Z_filled)*1.05)

  ax.set_xlabel("$z$", fontsize=18)
  ax.yaxis.set_rotate_label(False)  # disable automatic rotation
  ax.set_ylabel("$\lambda$", fontsize=18)
  ax.zaxis.set_rotate_label(False)  # disable automatic rotation
  ax.set_zlabel(z_label, fontsize=18, labelpad=25, rotation=0)

  # Adjust label rotations
  ax.xaxis.label.set_rotation(0)
  ax.yaxis.label.set_rotation(0)
  ax.zaxis.label.set_rotation(90)

  # Increase the font size of the tick labels for all axes
  ax.tick_params(axis='x', which='major', labelsize=14)
  ax.tick_params(axis='y', which='major', labelsize=14)
  ax.tick_params(axis='z', which='major', labelsize=14, pad=10)

  ax.view_init(elev, -112.8044)

  plt.tight_layout()
  plt.savefig(out_filename)
  plt.show()
  print('ax.azim {}'.format(ax.azim))
  print('ax.elev {}'.format(ax.elev))

def main():
  # File paths
  base_dir = '/home/tom/Documents/software/software/Stripper/data4hardcoded/ggh/Grids'
  # Known values
  nl_value = 5
  nt_value = 1
  nb_value = 0

  mH2omt2 = 12./23.
  z_threshold = 1. - 1./4.*mH2omt2

  #x_filename = base_dir + '/mb2=mH2_684/z_values.grid'
  #y_filename = base_dir + '/mb2=mH2_684/lam_values.grid'
  #z_filename = base_dir + '/mb2=mH2_684/amp-ct_gg_6556.grid'
  #plot_amp(x_filename, y_filename, z_filename,
  #         out_filename="tOSbOS_gg.pdf", x_lim=[0.02, 0.98], y_lim=[0.5, 0.98],
  #         z_label="$2 \mathrm{Re} \overline{\langle \mathcal{M}_{gg \\rightarrow gH}^{(0)}|\mathcal{M}_{gg \\rightarrow gH}^{(1)} \\rangle}\\bigg \\vert_{\mathrm{regulated}}$",
  #         elev=20)

  #x_filename = base_dir + '/mb2=mH2_684/z_values.grid'
  #y_filename = base_dir + '/mb2=mH2_684/lam_values.grid'
  #z_filename = base_dir + '/mb2=mH2_684/amp1-ct_qq_6556.grid'
  #plot_amp(x_filename, y_filename, z_filename,
  #         out_filename="tOSbOS_qBq.pdf", x_lim=[0.02, 0.98], y_lim=[0.5, 0.98],
  #         z_label="$2 \mathrm{Re} \overline{\langle \mathcal{M}_{\\bar{q}q \\rightarrow gH}^{(0)}|\mathcal{M}_{\\bar{q}q \\rightarrow gH}^{(1)} \\rangle}\\bigg \\vert_{\mathrm{regulated}}$",
  #         elev=20)

  #x_filename = base_dir + '/mb2=mH2_684/z_values.grid'
  #y_filename = base_dir + '/mb2=mH2_684/lam_values.grid'
  #z_filename = base_dir + '/mb2=mH2_684/amp2-ct_qg_6556.grid'
  #plot_amp(x_filename, y_filename, z_filename,
  #         out_filename="tOSbOS_gqB.pdf", x_lim=[0.02, 0.98], y_lim=[0.5, 0.98],
  #         z_label="$2 \mathrm{Re} \overline{\langle \mathcal{M}_{g\\bar{q} \\rightarrow \\bar{q}H}^{(0)}|\mathcal{M}_{g\\bar{q} \\rightarrow \\bar{q}H}^{(1)} \\rangle}\\bigg \\vert_{\mathrm{regulated}}$",
  #         elev=20)

  #x_filename = base_dir + '/mb2=mH2_684/z_values.grid'
  #y_filename = base_dir + '/mb2=mH2_684/lam_values.grid'
  #z_filename = base_dir + '/mb2=mH2_684/amp3-ct_qg_6556.grid'
  #plot_amp(x_filename, y_filename, z_filename,
  #         out_filename="tOSbOS_qg.pdf", x_lim=[0.02, 0.98], y_lim=[0.5, 0.98],
  #         z_label="$2 \mathrm{Re} \overline{\langle \mathcal{M}_{qg \\rightarrow qH}^{(0)}|\mathcal{M}_{qg \\rightarrow qH}^{(1)} \\rangle}\\bigg \\vert_{\mathrm{regulated}}$",
  #         elev=20.)

  #x_filename = base_dir + '/mtOSmtOS/z_values_tt.grid'
  #y_filename = base_dir + '/mtOSmtOS/lam_values.grid'
  #z_filename = base_dir + '/mtOSmtOS/amp_gg_tt.grid'
  #ct_filename = base_dir + '/mtOSmtOS/ct_gg_tt.grid'
  #plot_amp(x_filename, y_filename, z_filename, out_filename="tOStOS_gg.pdf",
  #         ct_filename=ct_filename, subtract_ct=True, x_lim=[0.002, 0.965], y_lim=[0.5,0.98],
  #         z_label="$2 \mathrm{Re} \overline{\langle \mathcal{M}_{gg \\rightarrow gH}^{(0)}|\mathcal{M}_{gg \\rightarrow gH}^{(1)} \\rangle}\\bigg \\vert_{\mathrm{regulated}}$",
  #         skip_x=3, skip_y=4, elev=10)

  #x_filename = base_dir + '/mtOSmtOS/z_values_tt.grid'
  #y_filename = base_dir + '/mtOSmtOS/lam_values.grid'
  #z_filename = base_dir + '/mtOSmtOS/amp1_qq_tt.grid'
  #ct_filename = base_dir + '/mtOSmtOS/ct1_qq_tt.grid'
  #plot_amp(x_filename, y_filename, z_filename, out_filename="tOStOS_qBq.pdf",
  #         ct_filename=ct_filename, subtract_ct=True, x_lim=[0.02, 0.98], y_lim=[0.5,0.98],
  #         z_label="$2 \mathrm{Re} \overline{\langle \mathcal{M}_{\\bar{q}q \\rightarrow gH}^{(0)}|\mathcal{M}_{\\bar{q}q \\rightarrow gH}^{(1)} \\rangle}\\bigg \\vert_{\mathrm{regulated}}$",
  #         skip_x=11, skip_y=10, elev=20)

  x_filename = base_dir + '/mtOSmtOS/z_values_tt.grid'
  y_filename = base_dir + '/mtOSmtOS/lam_values.grid'
  z_filename = base_dir + '/mtOSmtOS/amp2_qq_tt.grid'
  ct_filename = base_dir + '/mtOSmtOS/ct2_qq_tt.grid'
  plot_amp(x_filename, y_filename, z_filename, out_filename="tOStOS_gqB.pdf",
           ct_filename=ct_filename, subtract_ct=True, x_lim=[0.02, 0.98], y_lim=[0.5,0.98],
           z_label="$2 \mathrm{Re} \overline{\langle \mathcal{M}_{g \\bar{q} \\rightarrow \\bar{q}H}^{(0)}|\mathcal{M}_{g\\bar{q} \\rightarrow \\bar{q}H}^{(1)} \\rangle}\\bigg \\vert_{\mathrm{regulated}}$",
           skip_x=5, skip_y=10, elev=20)

  #x_filename = base_dir + '/mtOSmtOS/z_values_tt.grid'
  #y_filename = base_dir + '/mtOSmtOS/lam_values.grid'
  #z_filename = base_dir + '/mtOSmtOS/amp3_qq_tt.grid'
  #ct_filename = base_dir + '/mtOSmtOS/ct3_qq_tt.grid'
  #plot_amp(x_filename, y_filename, z_filename, out_filename="tOStOS_qg.pdf",
  #         ct_filename=ct_filename, subtract_ct=True, x_lim=[0.02, 0.98], y_lim=[0.5,0.98],
  #         z_label="$2 \mathrm{Re} \overline{\langle \mathcal{M}_{qg \\rightarrow qH}^{(0)}|\mathcal{M}_{qg \\rightarrow qH}^{(1)} \\rangle}\\bigg \\vert_{\mathrm{regulated}}$",
  #         skip_x=10, skip_y=10, elev=20)

if __name__=="__main__":
  main()