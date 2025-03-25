import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad

CA = 3.
nf = 6.
CF = 4./3.

mt = 90 #152.5
alphas_mZ = 0.118
gs_mZ = np.sqrt(4.*np.pi*alphas_mZ)
betas = (-11./3.*CA + 2./3.*nf)
beta2 = -19./6.
betaY = 41./6.
mH = 0.1 #100
mZ = 91.
sinThetaWsq_mZ = 0.23120
mW = mZ*np.sqrt(1. - sinThetaWsq_mZ)
alpha_mZ = 1./127.937
e_mZ = np.sqrt(alpha_mZ*4.*np.pi)
g2_mZ = e_mZ/np.sqrt(sinThetaWsq_mZ)
gY_mZ = e_mZ/np.sqrt(1. - sinThetaWsq_mZ)
alpha2_mZ = g2_mZ**2/(4.*np.pi)
alphaY_mZ = gY_mZ**2/(4.*np.pi)
vev = mZ/g2_mZ*2.*np.sqrt(1. - sinThetaWsq_mZ) # GeV
lam_mZ = mH**2/vev**2/2.
mu2_mZ = -mH**2./2.
Y_mZ = mt/vev*np.sqrt(2)
GF = 1./vev**2/np.sqrt(2.)

def B0(ma, mb, mc, mu):
  output, err = quad(lambda x: -np.log(((1. - x)*ma**2 + x*mb**2 - x*(1 - x)*mc**2)/mu**2), 0, 1)
  return output

def delta_lam(mu):
  LH = np.log(mH**2/mu**2)
  LW = np.log(mW**2/mu**2)
  LT = np.log(mt**2/mu**2)
  LZ = np.log(mZ**2/mu**2)
  del_lam = 1./2.*GF**2/(4.*np.pi)**2*(6*(LH - LW)*mH**6/(mH**2 - mW**2) - 8.*(2.*mW**4 + mZ**4) - 2.*(-3 + 6.*LT)*mH**2*mt**2
                                     + mH**4*(19 - 15*LH + 6*LW - 3.*np.sqrt(3)*np.pi) + 12*(mH**2 - 4*mt**2)*mt**2*B0(mt, mt, mH, mu)
                                     + 2*(mH**4 - 4*mH**2*mW**2 + 12*mW**4)*B0(mW, mW, mH, mu)
                                     + (mH**4 - 4*mH**2*mZ**2 + 12*mZ**4)*B0(mZ, mZ, mH, mu)
                                     + mH**2*(2.*(8*LW - 7)*mW**2 + (8*LZ - 7)*mZ**2 - 6*mZ**2*mW**2/(mZ**2 - mW**2)*(LZ - LW)))
  print("del_lam = ", del_lam)
  return del_lam

lam_mZ += delta_lam(mZ)
#gY_mZ = gY_mZ*np.sqrt(5./3.)

print("g2_mZ = ", g2_mZ)
print("gY_mZ = ", gY_mZ)
print("lam_mZ = ", lam_mZ)
print("Y_mZ = ", Y_mZ)

def Veff(phi, t):
  mu = mZ*np.exp(t)
  print("mu = ", mu)
  V0 = 1./2.*mu2_mZ*phi**2 + 1./4.*lam_mZ*phi**4
  V1 = 1./64./np.pi**2*((3./16.*(3.*g2_mZ**4 + 2.*g2_mZ**2*gY_mZ**2 + gY_mZ**4) - 3.*Y_mZ**4)*phi**4*np.log(phi**2/mu**2) + (mu2_mZ + 3.*lam_mZ*phi**2)**2*np.log(abs(mu2_mZ + 3.*lam_mZ*phi**2)/(mu**2)) + 3.*(mu2_mZ + lam_mZ*phi**2)**2*np.log(abs(mu2_mZ + lam_mZ*phi**2)/mu**2))
  return V0+V1

#def DVeff(phi, t):


def gs_running(t, gs): # Derivative for numerical integration
  return -7*gs**3/16./np.pi**2

def gs_running_analytic(t):
  return np.sqrt(4.*np.pi*alphas_mZ/(1. - 2.*alphas_mZ/(4.*np.pi)*betas*t))

def g2_running_analytic(t):
  return np.sqrt(4.*np.pi*alpha2_mZ/(1. - 2.*alpha2_mZ/(4.*np.pi)*beta2*t))

def gY_running_analytic(t):
  return np.sqrt(4.*np.pi*alphaY_mZ/(1. - 2.*alphaY_mZ/(4.*np.pi)*betaY*t))

def Y_running(t, Y):
  beta0 = (9./2.*Y**2 - 8.*gs_running_analytic(t)**2 - 9./4.*g2_running_analytic(t)**2 - 17./12.*gY_running_analytic(t)**2)
  return beta0*Y/16./np.pi**2

def delta_Y(mu):
  t = np.log(mu/mZ)
  LT = np.log(mt**2/mu**2)
  output = np.sqrt(2)*GF*mt**2*(8./3./(4.*np.pi)**2*gs_running_analytic(t)**2*(3*LT - 4) + 1./(4*np.pi)**2*np.sqrt(2)*GF*mt**2*(-9.*LT + 11.))
  #full =  2.*np.sqrt(2)*GF*mt**2*(1 + 8./3./(4.*np.pi)**2*gs_running_analytic(t)**2*(3*LT - 4) + 1./(4*np.pi)**2*np.sqrt(2)*GF*mt**2*(-9.*LT + 11.))
  return output/np.sqrt(2*np.sqrt(2)*GF*mt**2)

Y_mZ += delta_Y(mZ)
print("Y_mZ = ", Y_mZ)

def main():
  muRange = (mZ, 1e19)
  tRange = (0, np.log(muRange[1]/muRange[0]))

  '''
  gs_running_numerical = solve_ivp(gs_running, tRange, y0=[gs_mZ], dense_output=True, rtol=1.e-8)
  ts = np.linspace(0, tRange[-1], 100)
  gss = []
  for t in ts:
    gss.append(gs_running_numerical.sol(t)[0])

  plt.plot(np.exp(ts), gss)
  plt.plot(np.exp(ts), gs_running_analytic(ts), '--')
  plt.xscale('log')
  plt.show()
  return
  '''

  Y_running_sol = solve_ivp(Y_running, tRange, y0=[Y_mZ], dense_output=True, rtol=1.e-8)

  def delta_mu2(mu):
    lh = np.log(mH**2/mu**2)
    lW = np.log(mW**2/mu**2)
    lZ = np.log(mZ**2/mu**2)
    t = np.log(mu/mZ)
    G2 = g2_running_analytic(t)**2 + gY_running_analytic(t)**2
    del_mu2 = 1./(4.*np.pi)**2*(3.*Y_running_sol.sol(t)[0]**2*(4*mt**2 - mH**2)*B0(mt, mt, mH, mu) + 6.*lam_mZ**2*vev**2*(3*lh - 6 + np.pi*np.sqrt(3))
                                -vev**2/4.*(3*g2_running_analytic(t)**4 - 8*lam_mZ*g2_running_analytic(t)**2 + 16*lam_mZ**2)*B0(mW, mW, mH,mu)
                                -vev**2/8.*(3.*G2**2 - 8*lam_mZ*G2 + 16*lam_mZ**2)*B0(mZ, mZ, mH, mu)
                                +2*mW**2*(g2_running_analytic(t)**2 - 2*lam_mZ*(lW - 1.)) + mZ**2*(G2 - 2*lam_mZ*(lZ - 1)))

    return del_mu2

  mu2_mZ = -mH**2./2. + delta_mu2(mZ)
  print("mu2_mZ = ", mu2_mZ)
  print("delta_mu2 = " , delta_mu2(mZ))

  def gamma(t):
    return (-9*g2_running_analytic(t)**2 - 3*gY_running_analytic(t)**2 + 12.*Y_running_sol.sol(t)[0]**2)/64./np.pi**2

  def G(t):
    exponent, err = quad(gamma, 0, t, epsrel=1e-6)
    return np.exp(-exponent)

  #def B(mu):
  #  return 3./64./np.pi**2*((3.*g2_running_analytic(mu)**4 + 2.*g2_running_analytic(mu)**2*gY_running_analytic(mu)**2 + gY_running_analytic(mu)**4)/16. - gY_running_analytic(mu)**4)

  def lam_running(t, lam):
    #beta_lam = 4.*lam*gamma(mu) + (12.*lam**2 + B(mu))/8./np.pi**2
    beta_lam = 4*lam*gamma(t) + (24.*lam**2 + 3./4.*g2_running_analytic(t)**4 + 3./8.*(gY_running_analytic(t)**2 + g2_running_analytic(t)**2)**2 - 6.*Y_running_sol.sol(t)[0]**4)/16./np.pi**2
    return beta_lam

  lam_running_sol = solve_ivp(lam_running, tRange, y0=[lam_mZ], dense_output=True, rtol=1.e-8)

  #def beta_lam(mu):
  #  return 4.*lam_running_sol.sol(mu)[0]*gamma(mu) + (12.*lam_running_sol.sol(mu)[0]**2 + B(mu))/8./np.pi**2

  def mu2_running(t, mu2):
    #beta_mu2 = 2.*gamma(mu) + 3.*lam_running_sol.sol(mu)[0]/4./np.pi**2
    beta_mu2 = 6.*lam_running_sol.sol(t)[0]/16./np.pi**2 + gamma(t)
    return beta_mu2*2.*mu2

  mu2_running_sol = solve_ivp(mu2_running, tRange, y0=[mu2_mZ], dense_output=True, rtol = 0.0001)

  #def beta_m2(mu):
  #  return 2.*gamma(mu) + 3.*lam_running_sol.sol(mu)[0]/4./np.pi

  def Veff_resummed(phi, t):
    output = 1./2.*mu2_running_sol.sol(t)[0]*G(t)**2*phi**2 + 1./4.*lam_running_sol.sol(t)[0]*G(t)**4*phi**4
    return output

  def LogVeff_resummed(phi, t):
    output = 2*np.log(G(t)) + np.log(abs(1./2.*mu2_running_sol.sol(t)[0]*phi**2 + 1./4.*lam_running_sol.sol(t)[0]*G(t)**2*phi**4))
    return output

  print("g2(mT) = ", g2_running_analytic(np.log(173./mZ)))
  print("gY(mT) = ", gY_running_analytic(np.log(173./mZ)))
  print("Y(mT) = ", Y_running_sol.sol(np.log(173./mZ))[0])
  print("lam(mT) = ", lam_running_sol.sol(np.log(173./mZ))[0])
  print("mH(mT) = ", np.sqrt(-2.*mu2_running_sol.sol(np.log(173./mZ))[0]))

  lams = []
  mu2s = []
  gYs = []
  g2s = []
  gss = []
  Ys = []
  ts = np.linspace(tRange[0], tRange[1], 100)
  mus = []
  for t in ts:
    lams.append(lam_running_sol.sol(t)[0])
    mu2s.append(mu2_running_sol.sol(t)[0])
    gss.append(gs_running_analytic(t))
    g2s.append(g2_running_analytic(t))
    gYs.append(np.sqrt(5./3.)*gY_running_analytic(t))
    Ys.append(Y_running_sol.sol(t)[0])
    mus.append(mZ*np.exp(t))
  plt.figure()
  plt.plot(mus, lams, label="lam")
  #plt.plot(mus, np.sqrt(np.array(mu2s)*(-2.)), label="mu2")
  plt.plot(mus, gss, label="gs")
  plt.plot(mus, g2s, label="g2")
  plt.plot(mus, gYs, label="g1")
  plt.plot(mus, Ys, label="Yt")
  plt.xscale('log')
  plt.xlim(1e2, muRange[-1])
  plt.grid()
  plt.legend()
  plt.show()

  plt.figure()
  phis = np.logspace(0, 18, 100)
  Vs = []
  Vs_resummed = []
  mus = [vev]
  for mu in mus:
    Vs_resummed_mu = []
    Vs_mu = []
    for phi in phis:
      Vs_resummed_mu.append(Veff_resummed(phi, np.log(phi/mu)))
      Vs_mu.append(Veff(phi, np.log(phi/mu)))
    Vs.append(Vs_mu)
    Vs_resummed.append(Vs_resummed_mu)

  print(Vs_resummed)
  color = ["red", "blue", "green", "purple"]
  for i in range(0, len(mus)):
    plt.plot(phis, Vs[i], label=str(mus[i]), color=color[i])
    plt.plot(phis, np.array(Vs_resummed[i]), label=str(mus[i]) + " resummed", color=color[i], linestyle='dashed')
  plt.legend()
  plt.xlabel("phi")
  plt.ylabel("V")
  plt.yscale('symlog', linthresh=1e2, linscale=5)
  plt.xscale('log')
  plt.xlim(1e0, 1e13)

  plt.show()

if __name__=="__main__":
  main()