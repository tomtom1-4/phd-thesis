import lhapdf
import numpy as np
import matplotlib.pyplot as plt

import csv
import os
import random as rdm
import sys

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

mH = 125
S = 13000**2
#labels = ['$g$', '$d$', '$u$','$s$', '$c$', '$b$', '$\overline{d}$', '$\overline{u}$', '$\overline{s}$', '$\overline{c}$', '$\overline{b}$']
#pdg_id = [21, 1, 2, 3, 4, 5, -1, -2, -3, -4, -5]
pdg_ids = [[21,21], [21, 2], [21, 3], [2,-2], [2,2], [3, 3]]

labels = {
  21: "g",
  1: "d",
  2: "u",
  3: "s",
  4: "c",
  5: "b",
  -1: "\overline{d}",
  -2: "\overline{u}",
  -3: "\overline{s}",
  -4: "\overline{c}",
  -5: "\overline{b}",
}

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
  p = lhapdf.mkPDF("NNPDF31_nnlo_as_0118", 0)

  pdf = []
  mu = 125./2.
  xs = np.logspace(-3, 0, 200)
  for pid in p.flavors():
    pdf_i = []
    for x in xs:
      pdf_i.append(p.xfxQ(pid, x, mu)/x)
    pdf.append(pdf_i)
  pdf = np.array(pdf)
  # Create the plot
  plt.figure(figsize=(plot['singleplot']['width']/2, 3))
  for i in range(0, len(p.flavors())):
    plt.plot(xs, pdf[i-1]*xs, color=plot['color'][i], linewidth=plot['linewidth'], linestyle=plot['linestyle'], label="$x f_{" + labels[p.flavors()[i-1]] + "}(x)$")
  # Set logarithmic scaling for x-axis
  #plt.xscale('log')

  # Add labels, title, and legend
  plt.xlabel('$x$')
  #plt.ylabel('$x f_{P,i}(x, m_H/2)$')
  plt.legend()

  # Show grid for better readability
  plt.grid()
  plt.xlim(xs[0], 1)
  plt.ylim(0, 1)
  plt.text(0.2, 0.82, '$\mathrm{NNPDF31}$\n$\mu_F=m_H/2$')
  plt.subplots_adjust(left=0.2, right=0.97, top=0.97, bottom=0.1)
  plt.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/PDF.pdf")

  # Luminosities
  plt.figure(figsize=(plot['singleplot']['width']/2, 3))
  plt.ylim(0.01, 500)

  ts = np.logspace(-5, -0.000001, 300)
  plt.xlim(ts[0], 1)

  for i in range(0, len(pdg_ids)):
    luminosity = []
    combinatorics = 1.
    if pdg_ids[i][0] != pdg_ids[i][1]:
      combinatorics = 2.
    for t in ts:
      integrand = []
      xrange = np.linspace(t, 1, 1000)
      for x in xrange:
        integrand.append(combinatorics*p.xfxQ(pdg_ids[i][0], x, mu)/x**2 * p.xfxQ(pdg_ids[i][1], t/x, mu)*x/t*t)

      L = integrate(xrange, integrand)
      luminosity.append(L)

    plot_label = "$\mathcal{L}_{" + labels[pdg_ids[i][0]] + labels[pdg_ids[i][1]] + "}"
    if combinatorics != 1:
      plot_label += r"\times" + str(int(combinatorics))
    plot_label += "$"
    plt.plot(ts, luminosity, label=plot_label, \
      color=plot['color'][i], linewidth=plot['linewidth'], linestyle=plot['linestyle'])

  # Add labels, title, and legend
  plt.xlabel('$\\tau$')
  plt.ylabel('$\mathcal{L}$')
  plt.xscale('log')
  plt.yscale('log')
  plt.legend()
  plt.grid()
  #plt.grid(which='both', linestyle='--', linewidth=0.5)
  #plt.tight_layout()
  plt.subplots_adjust(left=0.2, right=0.97, top=0.97, bottom=0.1)
  plt.savefig("/home/tom/Uni/phd/PhD_thesis/thesis/Images/luminosity.pdf")




if __name__=="__main__":
  main()

