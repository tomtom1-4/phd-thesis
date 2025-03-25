import numpy as np
import matplotlib.pyplot as plt

import csv
import os
import random as rdm
import sys

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

# Constants
Gamma_SM = 0.00407  # Standard Model Higgs total width in GeV
GF = 1.16637e-5
mH = 125
vev = 1./np.sqrt(np.sqrt(2.)*GF)
print(vev)

# Function to calculate lambda_HSS as a function of m_S
def lambda_HSS_max(m_S, Gamma_inv_max):
    if m_S >= mH / 2:
        return np.nan  # Decay is kinematically forbidden
    sqrt_term = np.sqrt(1 - 4 * m_S**2 / mH**2)
    numerator = 128 * np.pi * mH * Gamma_inv_max
    denominator = vev**2 * sqrt_term
    lambda_squared = numerator / denominator
    return np.sqrt(lambda_squared)

def Gamma_HSS(M_S, Lambda_S):
  m_S = np.sqrt(M_S**2 + 1./4.*Lambda_S*vev**2)
  if m_S >= mH / 2:
    return np.nan  # Decay is kinematically forbidden

  output = Lambda_S**2*vev**2/128./np.pi/mH*np.sqrt(1. - 4.*m_S**2/mH**2)
  return output

# Function to calculate lambda_HVV as a function of m_V
def lambda_HVV_max(m_V, Gamma_inv_max):
    if m_V >= mH / 2:
        return np.nan  # Decay is kinematically forbidden
    sqrt_term = np.sqrt(1 - 4 * m_V**2 / mH**2)
    bracket_term = 1 - 4 * m_V**2 / mH**2 + 12 * m_V**4 / mH**4
    numerator = 512 * np.pi * m_V**4 * Gamma_inv_max
    denominator = vev**2 * mH**3 * sqrt_term * bracket_term
    lambda_squared = numerator / denominator
    return np.sqrt(lambda_squared)

def Gamma_HVV(M_V, Lambda_V):
  if m_V >= mH / 2:
    return np.nan  # Decay is kinematically forbidden

# Function to calculate lambda_Hchi as a function of m_chi and Lambda
def lambda_Hchi_max(m_chi, Lambda, Gamma_inv_max):
    if m_chi >= mH / 2.:
        return np.nan  # Decay is kinematically forbidden
    sqrt_term = np.sqrt(1. - 4. * m_chi**2 / mH**2)**3
    numerator = 32. * np.pi * Lambda**2 * Gamma_inv_max
    denominator = vev**2 * mH * sqrt_term
    lambda_squared = numerator / denominator
    return np.sqrt(lambda_squared)

def main():
  with open('../style.py', 'r') as f:
    s = f.read()
    plot = eval(s)

  # Mass ranges
  data_points = 1000
  mass_values = np.linspace(0, mH / 2, data_points)
  Lambda = 1000.0  # Example value for Lambda in GeV

  lambda_HSS_vals = []
  lambda_HVV_vals = []
  lambda_Hchi_vals = []
  # Invisible branching ratio limit
  for BR_inv_limit in [0.20, 0.10, 0.05]:  # 20%
    # Calculate the maximum partial width allowed
    Gamma_inv_max = BR_inv_limit*Gamma_SM/(1. - BR_inv_limit)
    # Calculate lambda_max values
    lambda_HSS_vals.append(np.array([lambda_HSS_max(m_S, Gamma_inv_max) for m_S in mass_values]))
    lambda_HVV_vals.append(np.array([lambda_HVV_max(m_V, Gamma_inv_max) for m_V in mass_values]))
    lambda_Hchi_vals.append( np.array([lambda_Hchi_max(m_chi, Lambda, Gamma_inv_max) for m_chi in mass_values]))

  data = [lambda_HSS_vals, lambda_HVV_vals, lambda_Hchi_vals]
  xlabels = ["$m_S\ \mathrm{[GeV]}$", "$m_V\ \mathrm{[GeV]}$", "$m_\chi\ \mathrm{[GeV]}$"]
  ylabels = ["$\lambda_{HSS}$", "$\lambda_{HVV}$", "$\lambda_{H\chi \chi}$"]
  # Plotting
  fig, axs = plt.subplots(3, 1, figsize=(5, 6), sharex=True)

  for i in range(0, len(data)):
    axs[i].plot(mass_values, data[i][0], color=plot['color'][0])
    axs[i].fill_between(mass_values,data[i][0],y2=1e2,where=~np.isnan(data[i][0]),color=plot['color'][0],alpha=0.3)
    axs[i].plot(mass_values, data[i][1], color=plot['color'][1])
    axs[i].fill_between(mass_values,data[i][1],y2=data[i][0],where=~np.isnan(data[i][1]),color=plot['color'][1],alpha=0.3)
    axs[i].plot(mass_values, data[i][2], color=plot['color'][2])
    axs[i].fill_between(mass_values,data[i][2],y2=data[i][1],where=~np.isnan(data[i][2]),color=plot['color'][2],alpha=0.3)

    #axs[i].plot(m_V_vals, lambda_HVV_vals, label=r'$\lambda_{HVV}^{\mathrm{max}}$')
    #axs[i].plot(m_chi_vals, lambda_Hchi_vals, label=r'$\lambda_{H\chi\chi}^{\mathrm{max}}$ (with $\Lambda=1000$ GeV)')

    axs[i].set_xlabel(xlabels[i])
    axs[i].set_ylabel(ylabels[i])
    axs[i].grid(True)
    axs[i].set_xlim(0, mH / 2)
    axs[i].set_ylim(1e-4, 1)
    axs[i].set_yscale('log')

  plt.tight_layout()
  plt.subplots_adjust(hspace=0.2)
  plt.savefig("../../Images/Higgs_portal_exclusion_bounds.pdf")
  #plt.show()

if __name__=="__main__":
  main()