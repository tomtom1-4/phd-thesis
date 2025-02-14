import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import os
import csv
import cmath

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

def extract_line(file_path, line_number):
  try:
    with open(file_path, 'r') as file:
      lines = file.readlines()
      if line_number > 0 and line_number <= len(lines):
        return lines[line_number - 1].strip()  # Returning the specified line, line_number is 1-based
      else:
        return f"Line number {line_number} is out of range."
  except FileNotFoundError:
    return "File not found."
  except Exception as e:
    return f"An error occurred: {e}"

def C0(z):
  t_1 = 1/z
  if (t_1 > 1.):
    ft_1 = np.arcsin(1./np.sqrt(t_1))**2
  else:
    ft_1 = -.25*(np.log((1.+np.sqrt(1.-t_1))/(1.-np.sqrt(1.-t_1))) -np.pi*(0+1j))**2;
  C0_1 = 0.5*t_1*(1.+(1.-t_1)*ft_1)
  return abs(C0_1)

def main():
  with open('../style.py', 'r') as f:
    s = f.read()
    plot = eval(s)

  mH = 125
  mt = 173.06
  rescaling = C0(mH**2/4/mt**2)**2*9
  print("rescaling = ", rescaling)

  E = np.arange(1, 21, 1)
  LO = []
  NLO = []
  HTL_LO = []
  HTL_NLO = []
  HTL_NNLO = []
  HTL_N3LO = []
  HTL_gg_NLO = []
  HTL_qg_NLO = []
  HTL_qq_NLO = []
  HTL_soft_virtual_NLO = []

  for e in E:
    LO_central = float(extract_line("SusHi/HEFT_LO_" + str(e) + "000.out", 18)[6:22])
    LO_down = float(extract_line("SusHi/HEFT_LO_" + str(e) + "000.out", 20)[6:22])
    LO_up = float(extract_line("SusHi/HEFT_LO_" + str(e) + "000.out", 21)[6:22])
    LO.append([LO_central, LO_down, LO_up])

    NLO_central = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 21)[6:22])
    NLO_down = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 23)[6:22])
    NLO_up = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 24)[6:22])
    NLO.append([NLO_central, NLO_down, NLO_up])

    HTL_NLO_central = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 38)[6:22])
    HTL_NLO_down = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 23)[6:22])/rescaling
    HTL_NLO_up = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 24)[6:22])/rescaling
    HTL_NLO.append([HTL_NLO_central, HTL_NLO_down, HTL_NLO_up])

    HTL_gg_NLO_central = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 43)[6:22])
    HTL_gg_NLO.append([HTL_gg_NLO_central, 0, 0])

    HTL_qg_NLO_central = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 44)[6:22])
    HTL_qg_NLO.append([HTL_qg_NLO_central, 0, 0])

    HTL_qq_NLO_central = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 45)[6:22])
    HTL_qq_NLO.append([HTL_qq_NLO_central, 0, 0])

    HTL_soft_virtual_NLO_central = float(extract_line("SusHi/HEFT_NLO_" + str(e) + "000.out", 42)[6:22])
    HTL_soft_virtual_NLO.append([HTL_soft_virtual_NLO_central, 0, 0])

    HTL_NNLO_central = float(extract_line("SusHi/HEFT_NNLO_" + str(e) + "000.out", 50)[6:22])
    HTL_NNLO_down = float(extract_line("SusHi/HEFT_NNLO_" + str(e) + "000.out", 23)[6:22])/rescaling
    HTL_NNLO_up = float(extract_line("SusHi/HEFT_NNLO_" + str(e) + "000.out", 24)[6:22])/rescaling
    HTL_NNLO.append([HTL_NNLO_central, HTL_NNLO_down, HTL_NNLO_up])

    HTL_N3LO_central = float(extract_line("SusHi/HEFT_N3LO_" + str(e) + "000.out", 67)[6:22])
    HTL_N3LO_down = float(extract_line("SusHi/HEFT_N3LO_" + str(e) + "000.out", 26)[6:22])
    HTL_N3LO_up = float(extract_line("SusHi/HEFT_N3LO_" + str(e) + "000.out", 27)[6:22])
    HTL_N3LO.append([HTL_N3LO_central, HTL_N3LO_down, HTL_N3LO_up])

  print(LO)

  LO = np.array(LO)
  HTL_LO = LO/rescaling
  NLO = np.array(NLO)
  HTL_NLO = np.array(HTL_NLO)
  HTL_gg_NLO = np.array(HTL_gg_NLO)
  HTL_qg_NLO = np.array(HTL_qg_NLO)
  HTL_qq_NLO = np.array(HTL_qq_NLO)
  HTL_soft_virtual_NLO = np.array(HTL_soft_virtual_NLO)
  HTL_NNLO = np.array(HTL_NNLO)
  HTL_N3LO = np.array(HTL_N3LO)

  HTL_NLO_control = HTL_gg_NLO + HTL_qg_NLO + HTL_qq_NLO + HTL_soft_virtual_NLO
  print(HTL_NLO_control)
  print(HTL_NLO)

  E_fine = np.linspace(min(E), max(E), len(E)*10)
  spline_LO = CubicSpline(E, HTL_LO[:, 0])
  spline_soft_virtual_NLO = CubicSpline(E, HTL_soft_virtual_NLO[:, 0])
  spline_gg_NLO = CubicSpline(E, HTL_gg_NLO[:, 0])
  spline_qg_NLO = CubicSpline(E, HTL_qg_NLO[:, 0])
  spline_qq_NLO = CubicSpline(E, HTL_qq_NLO[:, 0])

  HTL_LO_spline = spline_LO(E_fine)
  HTL_soft_virtual_NLO_spline = spline_soft_virtual_NLO(E_fine)
  HTL_gg_NLO_spline = spline_gg_NLO(E_fine)
  HTL_qg_NLO_spline = spline_qg_NLO(E_fine)
  HTL_qq_NLO_spline = spline_qq_NLO(E_fine)

  # Compute contributions (not yet stacked)
  soft = HTL_soft_virtual_NLO_spline
  gg = HTL_gg_NLO_spline
  qg = HTL_qg_NLO_spline
  qq = HTL_qq_NLO_spline

  print(qg)
  print(qq)

  # Compute total sum
  HTL_total = soft + gg + qg + qq

  # Sort contributions dynamically at each energy value
  contributions = np.array([soft, gg, qg, qq])  # Shape (4, len(E))
  channels = ["$\mathrm{soft}+\mathrm{virtual}$", "$gg$", "$qg$", "$q \\bar{q}$"]
  sorted_indices = np.argsort(contributions, axis=0)  # Indices that would sort each column
  sorted_contributions = np.take_along_axis(contributions, sorted_indices, axis=0)
  print(sorted_indices)

  # Compute cumulative sum for stacking
  cumulative = np.cumsum(sorted_contributions, axis=0)

  # Create figure and subplots
  fig, axs = plt.subplots(2, 1, figsize=(plot['subplot']['width'], plot['subplot']['height']), sharex=True)

  # Upper plot (stacked contributions)
  #axs[0].plot(E_fine, HTL_total, color='black', linewidth=plot['linewidth'] + 1.5, linestyle='-', label="Total")

  # Fill contributions dynamically in sorted order
  colors = [plot['color'][i] for i in range(4)]
  # WARNING: We do not assign the color dynamically. But rely on a fixed hierarchy
  for i in range(4):
    sorted_color = [colors[idx] for idx in sorted_indices[i]]  # Dynamic color assignment
    axs[0].plot(E_fine, cumulative[i], color=sorted_color[-1], label=channels[sorted_indices[i, -1]])
    if i == 0:
      axs[0].fill_between(E_fine, 0, cumulative[i], color=sorted_color[-1], alpha=0.5)
    else:
      axs[0].fill_between(E_fine, cumulative[i-1], cumulative[i], color=sorted_color[-1], alpha=0.5)

  axs[0].set_ylabel("$\sigma\ [\mathrm{pb}]$")
  axs[0].yaxis.set_label_coords(-0.1, 0.5)
  handles, labels = axs[0].get_legend_handles_labels()  # Get current legend items
  axs[0].legend(handles[::-1], labels[::-1], loc='upper left')
  axs[0].grid()
  axs[0].set_xlim(0, max(E_fine))
  axs[0].set_ylim(0, max(HTL_total)*1.05)
  axs[0].set_yticks([25, 50])
  axs[0].text(6.7, 55, "$\mathrm{NNPDF31}$\n$\mu_R=\mu_F=m_H/2$")


  '''
  # First subplot (linear y-scale)
  axs[0].plot(E_fine, HTL_soft_virtual_NLO_spline, color=plot['color'][0], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label="$\mathrm{soft}+\mathrm{virtual}$")
  axs[0].plot(E_fine, HTL_gg_NLO_spline, color=plot['color'][1], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label="$gg$")
  axs[0].plot(E_fine, HTL_qg_NLO_spline, color=plot['color'][2], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label="$qg$")
  axs[0].plot(E_fine, HTL_qq_NLO_spline, color=plot['color'][3], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label="$q\\bar{q}$")
  axs[0].set_ylabel("$\sigma\ [\mathrm{pb}]$")
  axs[0].yaxis.set_label_coords(-0.1, 0.5)
  axs[0].legend()
  axs[0].grid()
  axs[0].set_xlim(0, max(E))
  '''

  # Second subplot (logarithmic y-scale)
  axs[1].plot(E_fine, contributions[sorted_indices[0, 0]]/HTL_total, color=plot['color'][sorted_indices[0, 0]], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=channels[sorted_indices[i, 0]])
  axs[1].plot(E_fine, contributions[sorted_indices[1, 0]]/HTL_total, color=plot['color'][sorted_indices[1, 0]], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=channels[sorted_indices[i, 0]])
  axs[1].plot(E_fine, contributions[sorted_indices[2, 0]]/HTL_total, color=plot['color'][sorted_indices[2, 0]], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=channels[sorted_indices[i, 0]])
  axs[1].plot(E_fine, contributions[sorted_indices[3, 0]]/HTL_total, color=plot['color'][sorted_indices[3, 0]], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=channels[sorted_indices[i, 0]])
  axs[1].set_xlabel("$\sqrt{S}\ [\mathrm{TeV]}$")
  axs[1].set_ylabel("$\mathrm{Ratio\ to\ NLO}$")
  axs[1].yaxis.set_label_coords(-0.1, 0.5)
  axs[1].set_yscale("log")
  axs[1].grid()
  axs[1].set_xlim(0, max(E))
  axs[1].set_ylim(5.e-4, 1.)
  axs[1].set_xticks([0, 5, 10, 15, 20])
  # Add labels, title, and legend
  #plt.ylabel('$x f_{P,i}(x, m_H/2)$')

  # Show grid for better readability
  plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.11, hspace=plot['subplot']['hspace'])
  plt.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/channel_comparison_HTL_NLO.pdf")

if __name__=="__main__":
  main()