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

  #constant muF
  # Load data from file
  data = np.loadtxt('SusHi/energy_scan/HEFT_N3LO_13000.out_murdep', comments='#')

  # Extract columns
  muR = np.array(data[:, 0])
  muR_LO = np.array(data[:, 1])
  muR_NLO = np.array(data[:, 2])
  muR_NNLO = np.array(data[:, 3])
  muR_N3LO = np.array(data[:, 4])

  #constant muR
  muF = []
  muF_LO = []
  muF_NLO = []
  muF_NNLO = []
  muF_N3LO = []
  for i in range(1, 21):
    muF.append(float(extract_line("SusHi/scale_scan/scale_scan_" + str(i) + ".out", 152)[6:22])*mH)
    muF_LO.append(float(extract_line("SusHi/scale_scan/scale_scan_" + str(i) + ".out", 30)[6:22]))
    muF_NLO.append(float(extract_line("SusHi/scale_scan/scale_scan_" + str(i) + ".out", 41)[6:22])*rescaling)
    muF_NNLO.append(float(extract_line("SusHi/scale_scan/scale_scan_" + str(i) + ".out", 53)[6:22])*rescaling)
    muF_N3LO.append(float(extract_line("SusHi/scale_scan/scale_scan_" + str(i) + ".out", 67)[6:22])*rescaling)

  muF = np.array(muF)
  muF_LO = np.array(muF_LO)
  muF_NLO = np.array(muF_NLO)
  muF_NNLO = np.array(muF_NNLO)
  muF_N3LO = np.array(muF_N3LO)

  spline_muF_LO = CubicSpline(muF, muF_LO)
  spline_muF_NLO = CubicSpline(muF, muF_NLO)
  spline_muF_NNLO = CubicSpline(muF, muF_NNLO)
  spline_muF_N3LO = CubicSpline(muF, muF_N3LO)

  # Create plot
  fig, ax = plt.subplots(2, 1, figsize=(plot['subplot']['width'], plot['subplot']['height']), sharex=True)
  ax[0].plot(muR/mH, muR_LO, label="$\mathrm{rHTL} (\mathrm{LO})$", color=plot['color'][0], linewidth=plot['linewidth'])
  ax[0].plot(muR/mH, muR_NLO, label="$\mathrm{rHTL} (\mathrm{NLO})$", color=plot['color'][1], linewidth=plot['linewidth'])
  ax[0].plot(muR/mH, muR_NNLO, label="$\mathrm{rHTL} (\mathrm{NNLO})$", color=plot['color'][2], linewidth=plot['linewidth'])
  ax[0].plot(muR/mH, muR_N3LO, label="$\mathrm{rHTL} (\mathrm{N}^3\mathrm{LO})$", color=plot['color'][3], linewidth=plot['linewidth'])

  muF = np.linspace(muF[0], muF[-1], 100)
  ax[1].plot(muF/mH, spline_muF_LO(muF), label="$\mathrm{rHTL} (\mathrm{LO})$", color=plot['color'][0], linewidth=plot['linewidth'])
  ax[1].plot(muF/mH, spline_muF_NLO(muF), label="$\mathrm{rHTL} (\mathrm{NLO})$", color=plot['color'][1], linewidth=plot['linewidth'])
  ax[1].plot(muF/mH, spline_muF_NNLO(muF), label="$\mathrm{rHTL} (\mathrm{NNLO})$", color=plot['color'][2], linewidth=plot['linewidth'])
  ax[1].plot(muF/mH, spline_muF_N3LO(muF), label="$\mathrm{rHTL} (\mathrm{N}^3\mathrm{LO})$", color=plot['color'][3], linewidth=plot['linewidth'])


  # Add labels and legend
  ax[0].set_ylabel("$\sigma_{pp \\rightarrow HX}$")
  ax[0].set_yticks([0,20, 40, 60])
  ax[0].set_ylim(0, 60)
  ax[1].set_ylim(0, 60)
  ax[0].legend()
  ax[0].grid()
  ax[0].text(2.4, 47, "$\mathrm{NNPDF31}$\n$\mu_F=m_H/2$")


  ax[1].set_yticks([0,20, 40, 60])
  ax[1].set_ylabel("$\sigma_{pp \\rightarrow HX}$")
  ax[1].grid()
  ax[1].set_xlim(0, 3)
  ax[1].set_xticks([0,1,2,3])
  ax[1].set_xlabel("$\mu/m_H$")
  ax[1].text(2.4, 5, "$\mathrm{NNPDF31}$\n $\mu_R=m_H/2$")


  # Show grid for better readability
  plt.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.11, hspace=plot['subplot']['hspace']*1.1)
  plt.savefig("../../LaTeX/Images/scale_scan.pdf")
  plt.close(fig)

if __name__=="__main__":
  main()