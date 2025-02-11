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

zTop = 12./23./4.
zBottom = 125**2/4.78**2/4.

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
  fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(5, 4))
  #ax1.tick_params(axis='both', direction='in')
  #ax2.tick_params(axis='both', direction='in')

  ax1.plot(LO_Re[0], LO_Re[1], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(0)}$")
  ax1.plot(NLO_Re[0], NLO_Re[1], color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(1)}$")
  ax1.plot(NNLO_Re[0], NNLO_Re[1]/10, color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathcal{C}^{(2)}/10$")
  ax1.vlines([zTop, zBottom], ymin=-1.5, ymax=3.5, colors='k', linestyles='dashed')
  ax1.set_ylim(-1.5, 3.5)
  ax1.set_ylabel(r'$\mathrm{Re}(\mathcal{C})$')
  ax1.legend()
  ax1.text(0.2, -1., "$\mathrm{top}$\n$\mathrm{quark}$")
  ax1.text(70, -1., "$\mathrm{bottom}$\n$\mathrm{quark}$")

  # Show grid for better readability
  ax1.grid()
  ax1.set_xlim(LO_Re[0][0], LO_Re[0][-1])
  ax1.set_xscale('symlog', linthresh=1)
  #ax1.set_ylim(-1, 3)
  # plt.subplots_adjust(left=0.2, right=0.97, top=0.97, bottom=0.1)

  ax2.plot(LO_Im[0], LO_Im[1], color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{LO}$")
  ax2.plot(NLO_Im[0], NLO_Im[1], color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{NLO}$")
  ax2.plot(NNLO_Im[0], NNLO_Im[1]/10, color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=r"$\mathrm{NNLO}/10$")
  ax2.set_xlabel(r'$z$')
  ax2.set_ylabel(r'$\mathrm{Im}(\mathcal{C})$')
  ax2.vlines([zTop, zBottom], ymin=-1, ymax=2.5, colors='k', linestyles='dashed')
  ax2.set_ylim(-1, 2.5)
  ax2.grid()

  # Get current ticks and add a new one at position 5
  current_ticks = ax1.get_xticks()
  new_tick = 0.5

  # Combine current ticks with new tick and set them
  all_ticks = np.append(current_ticks, new_tick)
  ax1.set_xticks(all_ticks)
  # Optionally set custom labels (e.g., if needed)
  ticklabels = ["$0$", "$1$", "$10$", "$100$", "$1/2$"]
  ax1.set_xticklabels([tick for tick in ticklabels])

  ax1.yaxis.set_label_coords(-0.04, 0.5)
  ax2.yaxis.set_label_coords(-0.04, 0.5)


  fig.subplots_adjust(top=0.98, bottom=0.12, left=0.08, right=0.98, hspace=0.05)

  #fig.tight_layout()
  fig.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/form_factor.pdf")

  # Create a 3x3 array of numbers
  data = np.array([[16.3038, 1.70219, 0.344081],
                  [1.70219, 0.122182, 0.0374345],
                  [0.344081, 0.0374345, 0.0030346]])
  data_labels = np.array([["$16.30\ \mathrm{pb}$", "", ""],
                 ["$-1.70\ \mathrm{pb}$", "$0.12\ \mathrm{pb}$", ""],
                 ["$-0.34\ \mathrm{pb}$", "$0.04\ \mathrm{pb}$", "$0.003\ \mathrm{pb}$"]])
  # Create the plot
  fig, ax = plt.subplots(figsize=(2.5,2.5))

  # Display the data as an image with a heatmap color scheme
  cax = ax.matshow(data, cmap='viridis')

  # Annotate each cell with the numeric value
  for (i, j), val in np.ndenumerate(data_labels):
    if(i==0 and j==0):
      col = "black"
    else:
      col = "white"
    ax.text(j, i, val, ha='center', va='center', color=col)

  # Set ticks and labels for clarity (optional)
  ax.set_xticks(np.arange(0., len(data), step=1))
  ax.set_yticks(np.arange(0., len(data), step=1))
  ax.set_xticklabels(["$\mathrm{top}$", "$\mathrm{bottom}$", "$\mathrm{charm}$"])
  ax.set_yticklabels(["$\mathrm{top}$", "$\mathrm{bottom}$", "$\mathrm{charm}$"])

  ax.xaxis.set_ticks_position('bottom')
  ax.tick_params(axis='both', which='both', length=0)

  # Rotate tick labels for better readability
  plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
  plt.setp(ax.get_yticklabels(), rotation=45)

 # Grey out specific cells by adding rectangle patches over them
  grey_cells = [(0, 1), (0, 2), (1, 2)]
  for cell in grey_cells:
    rect = patches.Rectangle((cell[1] - 0.5, cell[0] - 0.5), width=1, height=1,
                            transform=ax.transData,
                            color='black', alpha=1)
    ax.add_patch(rect)

  # Adjust colorbar to match the height of the plot box
  #divider = make_axes_locatable(ax)
  #cax_cb = divider.append_axes("right", size="5%", pad=0.05)
  #color_bar = plt.colorbar(cax, cax=cax_cb)
  #color_bar.ax.set_ylabel("$\sigma [\mathrm{pb}]$", rotation=-90, va="bottom")
  #def fmt(x):
  #  return f'{x:.0f} pb'
  #color_bar.ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda x,_: fmt(x)))


  # Create a secondary y-axis on the right-hand side with additional tick labels
  sec_ax = ax.twinx()

  sec_ax.set_ylim(ax.get_ylim())
  sec_ax.set_yticks(np.arange(3))
  sec_ax.set_yticklabels(["$+16.30\ \mathrm{pb}$", "$-\ 1.58\ \mathrm{pb}$", "$-\ 0.30\ \mathrm{pb}$"])
  sec_ax.tick_params(length=0, pad=12)

  # Force the figure to draw to get the updated positions of tick labels
  fig.canvas.draw()

  # Get the renderer
  renderer = fig.canvas.get_renderer()

  # Get the bounding boxes of the tick labels in display coordinates
  tick_labels = sec_ax.get_yticklabels()
  bboxes = [label.get_window_extent(renderer=renderer) for label in tick_labels]

  # Combine the bounding boxes to get a single bounding box that encompasses all tick labels
  bbox = mtransforms.Bbox.union(bboxes)

  # Adjust padding (increase 'pad_x' to add space between the box and the plot)
  pad_x = -4  # horizontal padding in points (distance from the plot)
  pad_y = 0   # vertical padding in points around the labels

  # Expand the bounding box for padding and translate it to add space between box and plot
  bbox_expanded = bbox.expanded(1.0, 1.0).translated(pad_x, 0)

  # Convert the tick labels' bounding box to figure coordinates
  bbox_fig = bbox_expanded.transformed(fig.transFigure.inverted())

  # Get the axes bounding box in figure coordinates
  ax_bbox = sec_ax.get_position()

  # Adjust the vertical position to match the plot, with optional padding
  # Add vertical padding by increasing 'pad_y' above
  bbox_fig.y0 = ax_bbox.y0 - pad_y / fig.bbox.height + 0.02 - 0.15
  bbox_fig.y1 = ax_bbox.y1 + pad_y / fig.bbox.height - 0.02
  # 'bbox_fig.height' is automatically updated from 'y0' and 'y1'

  # Create a FancyBboxPatch with rounded edges
  fancy_box = patches.FancyBboxPatch(
      (bbox_fig.x0, bbox_fig.y0),
      width=bbox_fig.width*1.1,
      height=bbox_fig.height,
      transform=fig.transFigure,
      boxstyle="round,pad=0.02",
      facecolor='none',
      edgecolor='red',
      linewidth=1
  )
  fig.text(0.97, 0.04, "$+14.42\ \mathrm{pb}$")
  fig.text(0.95, 0.10, "$\overline{\hspace{1.7cm}}$")
  # Add the FancyBboxPatch to the figure
  fig.patches.append(fancy_box)

  #plt.tight_layout()
  #fig.subplots_adjust(top=0.99, bottom=0.13, left=0.13, right=0.8)
  fig.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/quark_effects_LO.pdf", bbox_inches='tight')


if __name__=="__main__":
  main()