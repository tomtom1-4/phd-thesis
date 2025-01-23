import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')


def read_csv_and_return_arrays(file_path):
    # Read the CSV file into a pandas DataFrame
    df = pd.read_csv(file_path)

    # Extract each column as an array
    column1 = df.iloc[:, 0].values
    column2 = df.iloc[:, 1].values
    column3 = df.iloc[:, 2].values
    column4 = df.iloc[:, 3].values
    column5 = df.iloc[:, 4].values
    column6 = df.iloc[:, 5].values
    column7 = df.iloc[:, 6].values
    column8 = df.iloc[:, 7].values


    return column1, column2, column3, column4, column5, column6, column7, column8

def alphas_running_nf(mu, mu0, alphas0, looporder, nf, mt, decoupling_scale):
  beta0 = 1./4./np.pi*(11. - 2./3.*nf)
  beta1 = (1./4./np.pi)**2*(102. - 38./3.*nf)

  decoupling = 0
  decoupling = 4./3.*0.5*np.log(mt**2/mu0**2)
  alphas = 0
  if(nf == 6):
    match looporder:
      case 1:
        alphas =  alphas0/(1. + alphas0*beta0*np.log(mu**2/mu0**2))
      case _:
        print("Running not implemented at loop order: ", looporder)
        alphas = 0
  elif(nf == 5):
    match looporder:
      case 1:
        alphas_nfp1 = alphas_running_nf(decoupling_scale, mu0, alphas0, 1, 6, mt, decoupling_scale)
        #alphas_nfp1 = alphas0
        alphas_d = alphas_nfp1
        for i in range(0, 20):
          decoupling = 1. - alphas_d/4./np.pi*(4./3.*0.5*np.log(mt**2/decoupling_scale**2))
          alphas_d = alphas_nfp1/decoupling
        alphas = alphas_d/(1. + alphas_d*beta0*np.log(mu**2/decoupling_scale**2))
      case _:
        print("Running not implemented at loop order: ", looporder)
        alphas = 0

  return alphas


def alphas_running(mu, mu0, alphas0, looporder, mt, decoupling_scale):
  if(mu < decoupling_scale):
    return alphas_running_nf(mu, mu0, alphas0, looporder, 5, mt, decoupling_scale)
  else:
    return alphas_running_nf(mu, mu0, alphas0, looporder, 6, mt, decoupling_scale)




if __name__=='__main__':
  alphas_MZ = 0.118002
  MZ = 91.1876
  mt0 = 162.7

  # Example usage
  file_path = 'running4l.txt'  # Replace with the path to your CSV file

  muR, alphas, mt, mb, alphas_LHAPDF, mbRunDec, mtRunDec, alphasRunDec = read_csv_and_return_arrays(file_path)
  muR1l, alphas1l, mt1l, mb1l, alphas_LHAPDF3, mbRunDec, mtRunDec, alphasRunDec = read_csv_and_return_arrays('running1l.txt')

  # alphas runnging at 1 loop
  fig1 = plt.figure(figsize=(9, 6))
  plt.plot(muR1l, 1./np.array(alphas1l), label='REvolver')
  #plt.plot(muR1l, alphas_running(muR1l, MZ, alphas_MZ, 6, 1))
  alphas_mu = []
  for mu in muR1l:
    alphas_mu.append(alphas_running(mu, MZ, alphas_MZ, 1, mt0, mt0*0.5))
  plt.plot(muR1l, 1./np.array(alphas_mu), label='Analytic runnging')
  #plt.plot(muR1l, 1./np.array(alphas_running_nf(muR1l, MZ, alphas_MZ, 1, 6, mt0)), label='nf = 6')
  #plt.plot(muR1l, 1./np.array(alphas_running_nf(muR1l, MZ, alphas_MZ, 1, 5, mt0)), label='nf = 5')


  plt.xlabel(r'$\mu_R$/GeV')
  plt.ylabel(r'$\alpha_s$')
  plt.yscale('log')
  plt.xscale('log')
  plt.axvline(mt0*0.5, linestyle='dashed', color='black', linewidth = 0.5)
  plt.legend()
  fig1.savefig("../../../figures/running/alphas_running_1l.pdf")

  # alphas running
  fig1, ax = plt.subplots(2, 1, figsize=(9, 6))
  ax[0].plot(muR, alphas, label="REvolver")
  ax[0].plot(muR, alphasRunDec, label="RunDec")
  ax[0].grid()
  ax[1].set_xlabel(r'$\mu_R$/GeV')
  ax[0].set_ylabel(r'$\alpha_s$')
  ax[0].legend()
  ax[1].plot(muR, (alphasRunDec-alphas_LHAPDF)/alphasRunDec)
  ax[0].sharex(ax[1])
  ax[1].vlines([62.5/2.], ymin = -1, ymax = 1, linestyle='dashed', color='gray')
  ax[1].set_ylim(-0.0005, 0.001)
  ax[1].grid()
  ax[1].set_ylabel(r'$1 - \frac{\alpha_s(LHAPDF)}{\alpha_s(RunDec)}$')
  fig1.savefig("../../../figures/running/alphas_running.pdf")

  # mt running
  fig2, ax = plt.subplots(2, 1, figsize=(9, 6))
  ax[0].plot(muR, mt, label="REvolver")
  ax[0].plot(muR, mtRunDec, label="RunDec")
  ax[0].grid()
  ax[1].set_xlabel(r'$\mu_R$/GeV')
  ax[0].set_ylabel(r'$m_t$/GeV')
  ax[0].legend()
  ax[1].plot(muR, (mt-mtRunDec)/mtRunDec)
  ax[0].sharex(ax[1])
  ax[1].vlines([62.5/2., 162.7], ymin = -1, ymax = 1, linestyle='dashed', color='gray')
  ax[1].set_ylim(-0.009, 0.001)
  ax[1].set_ylabel(r'$\frac{m_t(REvolver)}{m_t(RunDec)} - 1$')
  ax[1].grid()
  plt.subplots_adjust(hspace=0.05)
  fig2.savefig("../../../figures/running/mt_running.pdf")

  # mb running
  fig3, ax = plt.subplots(2, 1, figsize=(9, 6))
  ax[0].plot(muR, mb, label="REvolver")
  ax[0].plot(muR, mbRunDec, label="RunDec")
  ax[0].grid()
  ax[1].set_xlabel(r'$\mu_R$/GeV')
  ax[0].set_ylabel(r'$m_b$/GeV')
  ax[0].legend()
  ax[1].plot(muR, (mb-mbRunDec)/mbRunDec)
  ax[0].sharex(ax[1])
  ax[1].vlines([62.5/2.], ymin = -1, ymax = 1, linestyle='dashed', color='gray')
  ax[1].set_ylim(-0.009, 0.001)
  ax[1].set_ylabel(r'$\frac{m_b(REvolver)}{m_b(RunDec)} - 1$')
  ax[1].grid()
  plt.subplots_adjust(hspace=0.05)
  #ax[0].tight_layout()
  fig3.savefig("../../../figures/running/mb_running.pdf")

  #plt.xkcd()
  fig3, ax = plt.subplots(1, 1, figsize=(3, 1.5))
  #ax.plot(muR, mb, label="REvolver")
  ax.plot(muR, mbRunDec, label="RunDec")
  ax.grid()
  ax.set_xlabel(r'$\mu_R$/GeV')
  ax.set_ylabel(r'$m_b$/GeV')
  plt.tight_layout()
  fig3.savefig("mb_running.pdf")