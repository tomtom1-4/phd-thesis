import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, interp1d
from scipy.integrate import quad
import os
import csv
import cmath

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

alphas = 0.12515670484984653 # alphas @ mH/2 extracted from NNPDF31 using LHAPDF
GF = 1.16637e-5 # Fermi constant in GeV^-2
v = 1./np.sqrt(np.sqrt(2.)*GF) # vacuum expectation value in #GeV
mH = 125.
mu = mH/2.
GeV2pb = 3.89379e8 # conversion factor from 1/GeV^2 to pb calculated in figures/form_factor/form_factor.nb

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

def integrate(x, y):
  """
  Computes the area under the curve defined by arrays of x and y values
  using the trapezoidal rule.

  Parameters:
    x (list or array-like): The x-coordinates (must be sorted in ascending order).
    y (list or array-like): The corresponding y-coordinates.

  Returns:
    float: The computed area under the curve.
  """
  if len(x) != len(y):
    raise ValueError("The arrays x and y must have the same length.")
  if any(x[i] >= x[i + 1] for i in range(len(x) - 1)):
    raise ValueError("The array x must be sorted in ascending order.")

  # Use the trapezoidal rule for integration
  area = 0.0
  for i in range(1, len(x)):
    area += 0.5 * (x[i] - x[i - 1]) * (y[i] + y[i - 1])

  return area
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
  LO = HTL_LO_spline
  soft = HTL_soft_virtual_NLO_spline - LO
  gg = HTL_gg_NLO_spline
  qg = HTL_qg_NLO_spline
  qq = HTL_qq_NLO_spline


  # Compute total sum
  HTL_total = LO + soft + gg + qg + qq
  print("qg:")
  print(qg/HTL_total)
  print("gg:")
  print(gg/HTL_total)
  print("LO:")
  print(LO/HTL_total)
  print("qq:")
  print(qq/HTL_total)
  print("soft + virtual:")
  print(soft/HTL_total)

  # Sort contributions dynamically at each energy value
  contributions = np.array([LO, qq, qg, gg, soft])  # Shape (4, len(E))
  channels = ["$\mathrm{LO}$", "$q \\bar{q}$", "$qg$", "$gg$", "$\mathrm{soft}+\mathrm{virtual}$"]

  # Compute cumulative sum for stacking
  cumulative = np.cumsum(contributions, axis=0)

  # Create figure and subplots
  fig, axs = plt.subplots(2, 1, figsize=(plot['subplot']['width'], plot['subplot']['height']), sharex=True)

  # Upper plot (stacked contributions)
  #axs[0].plot(E_fine, HTL_total, color='black', linewidth=plot['linewidth'] + 1.5, linestyle='-', label="Total")

  # Fill contributions dynamically in sorted order
  colors = [plot['color'][i] for i in range(len(contributions))]
  # WARNING: We do not assign the color dynamically. But rely on a fixed hierarchy
  for i in range(len(contributions)):
    axs[0].plot(E_fine, cumulative[i], color=colors[i], label=channels[i])
    if i == 0:
      axs[0].fill_between(E_fine, 0, cumulative[i], color=colors[i], alpha=0.5)
    else:
      axs[0].fill_between(E_fine, cumulative[i-1], cumulative[i], color=colors[i], alpha=0.5)

  axs[0].set_ylabel("$\sigma\ [\mathrm{pb}]$")
  axs[0].yaxis.set_label_coords(-0.1, 0.5)
  handles, labels = axs[0].get_legend_handles_labels()  # Get current legend items
  axs[0].legend(handles[::-1], labels[::-1], loc='upper left')
  axs[0].grid()
  axs[0].set_xlim(0, max(E_fine))
  axs[0].set_ylim(0, max(HTL_total)*1.05)
  axs[0].set_yticks([25, 50])
  axs[0].text(6.7, 55, "$\mathrm{NNPDF31}$\n$\mu_R=\mu_F=m_H/2$")

  # Second subplot (logarithmic y-scale)
  for i in range(0, len(contributions)):
    axs[1].plot(E_fine, contributions[i]/HTL_total, color=plot['color'][i], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label=channels[i])

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

  # By hand approach
  # QQ Channel
  # Import partonic luminosity
  with open('../luminosity/luminosity_qqB.csv', 'r') as f:
    reader = csv.reader(f)
    lum_data = np.array([tuple(row) for row in reader])
    luminosity_qq = CubicSpline(lum_data[:,0], lum_data[:,1])

  with open('../luminosity/luminosity_qg.csv', 'r') as f:
    reader = csv.reader(f)
    lum_data = np.array([tuple(row) for row in reader])
    luminosity_qg = CubicSpline(lum_data[:,0], lum_data[:,1])

  with open('../luminosity/luminosity_gg.csv', 'r') as f:
    reader = csv.reader(f)
    lum_data = []
    for row in reader:
      lum_data.append([float(row[0]), float(row[1])])
    lum_data = np.array(lum_data)
    luminosity_gg = CubicSpline(lum_data[:,0], lum_data[:,1])

  #print(lum_data[:,1])

  def qgHq_integrand(S, tau, mu):
    s = S*tau
    xi = mH**2/s
    #if(1. - xi < 0):
    #  print("below threshold")
    #  return 0.
    partonic_xSec = 1./576./(np.pi**2)*alphas**3/(v**2)\
      *((1. - xi)*(3.*xi - 7.)/3. + 1./2.*xi*4./3.*((1. + (1. - xi)**2)/xi)*(1. + np.log(mH**2/(mu**2)) + np.log((1. - xi)**2/xi)))

    return partonic_xSec*luminosity_qg(tau)/tau

  def qqHg_integrand(S, tau):
    s = S*tau
    xi = mH**2/s
    if(1. - xi < 0):
      print("below threshold")
      return 0.
    partonic_xSec = 1./486./(np.pi**2)*alphas**3/(v**2)*(1. - xi)**3

    return partonic_xSec*luminosity_qq(tau)/tau

  def ggHg_integrand(S, tau, mu):
    s = S*tau
    xi = mH**2/s
    #if(1. - xi < 0):
    #  print("below threshold")
    #  return 0.

    partonic_xSec_woDist = 1./576./(np.pi**2)*alphas**3/(v**2) \
      *(-11./2.*(1. - xi)**3 + 6.*(1. - xi)*(1. + xi**2)*np.log(s/mu**2))

    result = partonic_xSec_woDist*luminosity_gg(tau)/tau

    return result

  def ggH_integrand_delta(S, tau):
    s = S*tau
    xi = mH**2/s
    partonic_xSec = alphas**3/576./np.pi**2/(v**2)*(11./2. + np.pi**2)*tau
    return partonic_xSec*luminosity_gg(tau)/tau

  def ggH_integrand_PlusDist(S, tau, mu):
    s = S*tau
    xi = mH**2/s
    #if(1. - xi < 0):
    #  print("below threshold")
    #  return 0.

    partonic_xSec = 1./576./(np.pi**2)*alphas**3/(v**2) \
      *(6.*(1. + xi**4 + (1. - xi)**4)*np.log(1. - xi)/(1. - xi) \
        + 6.*(xi**2/(1. - xi))*np.log(s/mu**2))

    result = partonic_xSec*luminosity_gg(tau)/tau

    # plus distribution contribution
    tau0 = mH**2/S
    result = result - alphas**3/576./np.pi**2/v**2*(6.*2.*np.log(1. - xi)/(1. - xi) + 6.*1./(1. - xi)*np.log(mH**2/mu**2))*luminosity_gg(tau)/tau

    return result

  def LO_integrand(S, tau):
    s = S*tau
    xi = mH**2/s
    partonic_xSec = alphas**2/576./np.pi/(v**2)*tau
    return partonic_xSec*luminosity_gg(tau)/tau

  COM = 13000.
  xSec_qq, error_qq = quad(lambda tau: qqHg_integrand(COM**2, tau)*GeV2pb, mH**2/COM**2, 1)
  xSec_qg, error_qg = quad(lambda tau: qgHq_integrand(COM**2, tau, mu)*GeV2pb, mH**2/COM**2, 1)
  xSec_gg, error_gg = quad(lambda tau: ggHg_integrand(COM**2, tau, mu)*GeV2pb, mH**2/COM**2, 1)
  xSec_soft, error_soft = (ggH_integrand_delta(COM**2, mH**2/COM**2)*GeV2pb, 0)
  xSec_gg = xSec_gg + quad(lambda tau: ggH_integrand_PlusDist(COM**2, tau, mu)*GeV2pb, mH**2/COM**2 + 1.e-9, 1)[0]
  xSec_LO = LO_integrand(COM**2, mH**2/COM**2)*GeV2pb

  print("xSec_qq:", xSec_qq, " +- ", error_qq)
  print("Control with SusHi = ", spline_qq_NLO(COM/1000.))
  print("ratio = ", spline_qq_NLO(COM/1000.)/xSec_qq)

  print("\nxSec_qg:", xSec_qg, " +- ", error_qg)
  print("Control with SusHi = ", spline_qg_NLO(COM/1000.))
  print("ratio = ", spline_qg_NLO(COM/1000.)/xSec_qg)

  print("\nxSec_soft:", xSec_soft, " +- ", error_soft)
  print("Control with SusHi = ", spline_soft_virtual_NLO(COM/1000.) - spline_LO(COM/1000.))
  print("ratio = ", (spline_soft_virtual_NLO(COM/1000.) - spline_LO(COM/1000.))/xSec_soft)

  print("\nxSec_LO:", xSec_LO, " +- ", 0)
  print("Control with SusHi = ", spline_LO(COM/1000.))
  print("ratio = ", spline_LO(COM/1000.)/xSec_LO)

  print("\nxSec_gg:", xSec_gg, " +- ", error_gg)
  print("Control with SusHi = ", spline_gg_NLO(COM/1000.))
  print("ratio = ", xSec_gg/spline_gg_NLO(COM/1000.))

  print("\nxSec_gg (combined):", xSec_gg + xSec_soft, " +- ", error_gg)
  print("Control with SusHi = ", spline_gg_NLO(COM/1000.) + spline_soft_virtual_NLO(COM/1000.) - spline_LO(COM/1000.))
  print("ratio = ", (xSec_gg + xSec_soft)/(spline_gg_NLO(COM/1000.) + spline_soft_virtual_NLO(COM/1000.) - spline_LO(COM/1000.)))

  #print(lum_data[:,0])

  OneMinusXi = np.logspace(-7, 0, 10000)[:-1]
  #xi = 1 - OneMinusXi
  tau = mH**2/COM**2 + OneMinusXi
  #print("min = ", min(tau), ", max = ", max(tau))
  #check = integrate(tau, ggHg_integrand(COM**2, tau, mu)*GeV2pb)
  #print(xi)
  #print(ggHg_integrand(COM**2, tau, mu)*GeV2pb)
  #print("trapezoid integration = ", check)

  #plt.figure(figsize=(5,6))
  #plt.plot(tau, abs(ggHg_integrand(COM**2, tau, mu) + ggH_integrand_PlusDist(COM**2, tau, mu))*GeV2pb, label="gg")
  #plt.plot(tau, abs(qgHq_integrand(COM**2, tau, mu))*GeV2pb, label="qg")
  #plt.yscale('log')
  #plt.xscale('log')
  #plt.legend()
  #plt.savefig("ggHg_integrand.pdf")

if __name__=="__main__":
  main()