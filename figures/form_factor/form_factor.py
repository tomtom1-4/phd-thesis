import matplotlib.pyplot as plt
import numpy as np

import csv
import os
import random as rdm
import sys

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

def read_data_pairs(file_name):
  data_pairs = []

  # Read the CSV file
  with open(file_name, mode='r') as file:
    reader = csv.reader(file)
    for row in reader:
      # Convert strings to float and append to list
      data_pairs.append((float(row[0]), float(row[1])))
  data_pairs = np.transpose(np.array(data_pairs))
  return data_pairs

def main():
  with open('../style.py', 'r') as f:
    s = f.read()
    plot = eval(s)

  LO_Re = read_data_pairs('LO_Re.csv')
  NLO_Re = read_data_pairs('NLO_Re.csv')
  NNLO_Re = read_data_pairs('NNLO_Re.csv')

  LO_Im = read_data_pairs('LO_Im.csv')
  NLO_Im = read_data_pairs('NLO_Im.csv')
  NNLO_Im = read_data_pairs('NNLO_Im.csv')

  # Create the plot with corrected function name 'subplots'
  fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(plot['singleplot']['width'], 3))

  ax1.plot(LO_Re[0], LO_Re[1], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{LO}$")
  ax1.plot(NLO_Re[0], NLO_Re[1], color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{NLO}$")
  ax1.plot(NNLO_Re[0], NNLO_Re[1]/10, color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{NNLO}/10$")
  ax1.set_ylabel(r'$\mathrm{Re}(\mathcal{C})$')
  ax1.legend()

  # Show grid for better readability
  ax1.grid()
  ax1.set_xlim(LO_Re[0][0], LO_Re[0][-1])
  #ax1.set_ylim(-1, 3)
  # plt.subplots_adjust(left=0.2, right=0.97, top=0.97, bottom=0.1)

  ax2.plot(LO_Im[0], LO_Im[1], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{LO}$")
  ax2.plot(NLO_Im[0], NLO_Im[1], color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{NLO}$")
  ax2.plot(NNLO_Im[0], NNLO_Im[1]/10, color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{NNLO}/10$")
  ax2.set_xlabel(r'$z$')
  ax2.set_ylabel(r'$\mathrm{Im}(\mathcal{C})$')
  ax2.grid()

  fig.subplots_adjust(top=0.98, bottom=0.12, left=0.12, right=0.98, hspace=0)

  # Optional: Adjust the y-axis labels to prevent overlapping
  ax1.yaxis.set_label_coords(-0.1, 0.5)
  ax2.yaxis.set_label_coords(-0.1, 0.5)

  #fig.tight_layout()
  fig.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/form_factor.pdf")


if __name__=="__main__":
  main()