import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.transforms as mtransforms
import numpy as np

import csv
import os
import random as rdm
import sys

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

def parse_line(line):
  line = line.strip()
  main_str, rest = line.split(',', 1)
  main_val = float(main_str.strip())

  # Extract contents inside parentheses
  pairs = re.findall(r'-?\d+(?:\.\d+)?(?:e-?\d+)?', rest)

  complex_numbers = []
  i = 0
  while i < len(pairs) - 1:
    real_str = pairs[i]
    imag_str = pairs[i + 1]
    real_part = float(real_str.strip())
    imag_part = float(imag_str.strip())
    complex_numbers.append(complex(real_part, imag_part))
    i += 2

  return main_val, complex_numbers

def parse_file(filepath):
  data = []
  with open(filepath, 'r') as f:
    for line in f:
      # Skip any empty or whitespace-only lines
      if not line.strip():
        continue
      main_val, comps = parse_line(line)
      data.append((main_val, comps))
  return data

def main():
  with open('../style.py', 'r') as f:
    s = f.read()
    plot = eval(s)

  results_5FS = parse_file('5FS_nl0.csv') # Cb
  results_5FS_nl1 = parse_file('5FS_nl1.csv') # Cb + Cbl
  results_4FS = parse_file('4FS_nl0.csv') # Cb + Cbb
  results_bt = parse_file('CMixed_NNLO.csv' )
  print(results_bt)
  z = []
  Cb_NNLO = []
  Cbb_NNLO = []
  Cbl_NNLO = []
  Ctptt_NNLO = []
  Ctl_NNLO = []
  CMixed_NNLO = []
  for i in range(0, len(results_5FS)):
    z.append(results_5FS[i][0])
    Cbl_NNLO.append(results_5FS_nl1[i][1][2] - results_5FS[i][1][2])
    Cb_NNLO.append(results_5FS[i][1][2])
    Cbb_NNLO.append(results_4FS[i][1][2] - Cb_NNLO[-1])
    Ctptt_NNLO.append(results_5FS[i][1][2])
    Ctl_NNLO.append(results_5FS_nl1[i][1][2] - results_5FS[i][1][2])
    CMixed_NNLO.append(results_bt[i][1][0])

  z = np.array(z)
  Cb_NNLO = np.array(Cb_NNLO)
  Cbb_NNLO = np.array(Cbb_NNLO)
  Cbl_NNLO = np.array(Cbl_NNLO)
  Ctptt_NNLO = np.array(Ctptt_NNLO)
  Ctl_NNLO = np.array(Ctl_NNLO)
  CMixed_NNLO = np.array(CMixed_NNLO)
  print(z)
  print(Cb_NNLO)

  z_thresh = 1
  i_thresh = 0
  for i in range(0, len(z)):
    if z[i] > z_thresh:
      i_thresh = i
      break

  fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(plot['subplot']['width'], plot['subplot']['height']))

  #ax1.plot(z[:i_thresh], Cb_NNLO.real[:i_thresh], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'])
  ax1.plot(z[i_thresh:], Cb_NNLO.real[i_thresh:], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{b}$")
  ax1.plot(z[i_thresh:], Cbb_NNLO.real[i_thresh:], color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{bb}$")
  ax1.plot(z[i_thresh:], Cbl_NNLO.real[i_thresh:], color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{bl}$")
  ax1.plot(z[i_thresh:], CMixed_NNLO.real[i_thresh:], color=plot['color'][3], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{bt}$")
  ax1.set_ylabel(r'$\mathrm{Re}(\mathcal{C})$')
  ax1.legend()

  # Show grid for better readability
  ax1.grid()
  ax1.set_xlim(1, 500)

  #ax2.plot(z[:i_thresh], Cb_NNLO.imag[:i_thresh], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'])
  ax2.plot(z[i_thresh:], Cb_NNLO.imag[i_thresh:], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'])
  ax2.plot(z[i_thresh:], Cbb_NNLO.imag[i_thresh:], color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'])
  ax2.plot(z[i_thresh:], Cbl_NNLO.imag[i_thresh:], color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'])
  ax2.plot(z[i_thresh:], CMixed_NNLO.imag[i_thresh:], color=plot['color'][3], linewidth=plot['linewidth'], linestyle=plot['linestyle'])
  ax2.set_xlabel(r'$z$')
  ax2.set_ylabel(r'$\mathrm{Im}(\mathcal{C})$')
  ax2.grid()
  fig.subplots_adjust(top=0.98, bottom=0.12, left=0.11, right=0.98, hspace=plot['subplot']['hspace'])

  fig.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/form_factor_coefficients_b.pdf")

  fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(plot['subplot']['width'], plot['subplot']['height']))

  ax1.plot(z[:i_thresh], Ctptt_NNLO.real[:i_thresh], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{t} + \mathcal{C}^{(2)}_{tt}$")
  ax1.plot(z[:i_thresh], Ctl_NNLO.real[:i_thresh], color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{tl}$")
  ax1.plot(z[:i_thresh], CMixed_NNLO.real[:i_thresh], color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{tb}$")
  ax1.set_ylabel(r'$\mathrm{Re}(\mathcal{C})$')
  ax1.legend()

  # Show grid for better readability
  ax1.grid()
  ax1.set_xlim(0, 1)

  ax2.plot(z[:i_thresh], Ctptt_NNLO.imag[:i_thresh], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{t} + \mathcal{C}^{(2)}_{tt}$")
  ax2.plot(z[:i_thresh], Ctl_NNLO.imag[:i_thresh], color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{tl}$")
  ax2.plot(z[:i_thresh], CMixed_NNLO.imag[:i_thresh], color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}_{tb}$")
  ax2.set_xlabel(r'$z$')
  ax2.set_ylabel(r'$\mathrm{Im}(\mathcal{C})$')
  ax2.grid()
  fig.subplots_adjust(top=0.98, bottom=0.12, left=0.11, right=0.98, hspace=plot['subplot']['hspace'])

  fig.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/form_factor_coefficients_t.pdf")


if __name__ == "__main__":
  main()