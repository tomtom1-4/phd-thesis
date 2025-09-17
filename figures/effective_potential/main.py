#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad
import math


import csv
import os
import random as rdm
import sys

plt.rcParams['text.usetex'] = True
os.environ['PATH'] = os.environ['PATH'] + ':/Library/TeX/texbin'
plt.close('all')

###############################################################################
# 1) Global constants and utility functions
###############################################################################

CA = 3.0
nf = 6.0
CF = 4.0/3.0

# Measured or fixed parameters at mZ scale
alphas_mZ = 0.118
gs_mZ = np.sqrt(4.0 * np.pi * alphas_mZ)
betas = (-11.0 / 3.0 * CA + 2.0 / 3.0 * nf)
beta2 = -19.0 / 6.0
betaY = 41.0 / 6.0

mZ = 91.0
sinThetaWsq_mZ = 0.23120
mW = mZ * np.sqrt(1.0 - sinThetaWsq_mZ)
alpha_mZ = 1.0 / 127.937
e_mZ = np.sqrt(alpha_mZ * 4.0 * np.pi)
g2_mZ_default = e_mZ / np.sqrt(sinThetaWsq_mZ)
gY_mZ_default = e_mZ / np.sqrt(1.0 - sinThetaWsq_mZ)
alpha2_mZ = g2_mZ_default**2 / (4.0 * np.pi)
alphaY_mZ = gY_mZ_default**2 / (4.0 * np.pi)

# Default vacuum expectation value
vev = (
    mZ / g2_mZ_default * 2.0 * np.sqrt(1.0 - sinThetaWsq_mZ)
)  # ~246 GeV

def B0(ma, mb, mc, mu):
    """
    One-loop Passarino–Veltman function B0, adapted from the original script.
    """
    val, err = quad(
        lambda x: -np.log(((1.0 - x) * ma**2 + x * mb**2 - x * (1.0 - x) * mc**2) / mu**2),
        0,
        1
    )
    return val

###############################################################################
# 2) RGE setup and running definitions
###############################################################################

def gs_running_analytic(t):
    """
    Leading-order RGE for g_s(μ).
    μ is related to t by μ = mZ * e^t.
    """
    numerator = 4.0 * np.pi * alphas_mZ
    denominator = 1.0 - 2.0 * alphas_mZ / (4.0 * np.pi) * betas * t
    return np.sqrt(numerator / denominator)

def g2_running_analytic(t):
    """
    Leading-order RGE for g2(μ).
    """
    numerator = 4.0 * np.pi * alpha2_mZ
    denominator = 1.0 - 2.0 * alpha2_mZ / (4.0 * np.pi) * beta2 * t
    return np.sqrt(numerator / denominator)

def gY_running_analytic(t):
    """
    Leading-order RGE for gY(μ).
    """
    numerator = 4.0 * np.pi * alphaY_mZ
    denominator = 1.0 - 2.0 * alphaY_mZ / (4.0 * np.pi) * betaY * t
    return np.sqrt(numerator / denominator)

def Y_running(t, Y):
    """
    1-loop Beta function for the top Yukawa coupling:
      dY/dt = beta0 * Y / (16π^2).
    """
    term = (
        9.0 / 2.0 * Y**2
        - 8.0 * gs_running_analytic(t)**2
        - 9.0 / 4.0 * g2_running_analytic(t)**2
        - 17.0 / 12.0 * gY_running_analytic(t)**2
    )
    return term * Y / (16.0 * np.pi**2)

def gamma_function(t, Y_sol):
    """
    Wavefunction renormalization factor:
      gamma(μ) = [ -9*g2^2(μ) - 3*gY^2(μ) + 12*Y^2(μ) ] / (64π^2).
    """
    g2_val = g2_running_analytic(t)
    gY_val = gY_running_analytic(t)
    Y_val = Y_sol(t)[0]
    return (
        -9.0 * g2_val**2
        - 3.0 * gY_val**2
        + 12.0 * Y_val**2
    ) / (64.0 * np.pi**2)

def lam_running(t, lam, Y_sol):
    """
    1-loop Beta function for λ(μ), including wavefunction factor:
      dλ/dt = 4*γ*λ + poly(...) / (16π^2).
    """
    g = gamma_function(t, Y_sol)
    g2_val = g2_running_analytic(t)
    gY_val = gY_running_analytic(t)

    lam_part = (
        24.0 * lam**2
        + 3.0 / 4.0 * g2_val**4
        + 3.0 / 8.0 * (g2_val**2 + gY_val**2)**2
        - 6.0 * (Y_sol(t)[0]**4)
    ) / (16.0 * np.pi**2)

    return 4.0 * lam * g + lam_part

def mu2_running(t, mu2, lam_sol, Y_sol):
    """
    Very approximate 1-loop RGE for mu^2(μ).
    """
    g = gamma_function(t, Y_sol)
    lam_val = lam_sol(t)[0]
    return (6.0 * lam_val / (16.0 * np.pi**2) + g) * 2.0 * mu2

###############################################################################
# 3) Setup and run the RGEs for a given (mt, mH)
###############################################################################

def setup_and_run_RGE(mt, mH):
    """
    1) Initialize parameters for a given top mass and Higgs mass.
    2) Solve RGEs for Y(t), λ(t), μ²(t).
    3) Return the solutions for further analysis.
    """
    mu2_mZ = - (mH**2) / 2.0
    lam_mZ = (mH**2) / (2.0 * vev**2)
    Y_mZ = (mt / vev) * np.sqrt(2.0)

    mu2_found, lam_found = solve_mu2_lam(
        phi_target=vev, mH_target=mH,
        g2=g2_mZ_default, gY=gY_mZ_default, Yt=Y_mZ, muR=mZ
    )
    #print(f"Solved (mu2, lam) = ({mu2_found:.5g}, {lam_found:.5g}) at one-loop")

    # Check the conditions
    #eqs = residual_equations([mu2_found, lam_found], vev, mH, g2_mZ_default, gY_mZ_default, np.sqrt(2)*mt/vev, mZ, order=1)
    #dV_phi = eqs[0]
    #d2V_phi = eqs[1] + 125**2  # eq2 = d2V - mH^2
    mu2_mZ = mu2_found
    lam_mZ = lam_found
    #print(f"Residuals: dV/dphi(v)={dV_phi:.5g}, [d²V/dphi²(v) - mH²]={eqs[1]:.5g}")
    #print(f"Thus, d²V/dphi²(v) = {d2V_phi:.5g} GeV^2 vs. mH^2={mH**2:.5g} GeV^2")

    # Integration range in t (μ from mZ to ~1e19):
    tRange = (0.0, np.log(1.0e19 / mZ))

    # Solve for Y(t)
    Y_sol_ivp = solve_ivp(
        lambda t, Y: Y_running(t, Y),
        tRange,
        y0=[Y_mZ],
        dense_output=True,
        rtol=1.e-7,
        atol=1.e-9
    )

    # Solve for λ(t)
    def lam_ode(t, lam_val):
        return lam_running(t, lam_val, Y_sol_ivp.sol)

    lam_sol_ivp = solve_ivp(
        lam_ode,
        tRange,
        y0=[lam_mZ],
        dense_output=True,
        rtol=1.e-7,
        atol=1.e-9
    )

    # Solve for μ²(t)
    def mu2_ode(t, mu2_val):
        return mu2_running(t, mu2_val, lam_sol_ivp.sol, Y_sol_ivp.sol)

    mu2_sol_ivp = solve_ivp(
        mu2_ode,
        tRange,
        y0=[mu2_mZ],
        dense_output=True,
        rtol=1.e-7,
        atol=1.e-9
    )

    return {
        "Y_sol":   Y_sol_ivp,
        "lam_sol": lam_sol_ivp,
        "mu2_sol": mu2_sol_ivp
    }

###############################################################################
# 4) Simple stability check
###############################################################################

def check_stability(rge_solutions):
    """
    Returns 'unstable', 'metastable', or 'stable' by checking the sign of λ(μ).
    """
    lam_sol = rge_solutions["lam_sol"].sol
    t_vals = rge_solutions["lam_sol"].t
    t_end = t_vals[-1]

    # Evaluate λ at largest scale
    lam_high = lam_sol(t_end)[0]
    if lam_high < 0:
        if True: #lam_high < -0.05:
            return "unstable"
        else:
            return "metastable"

    # Check for negative crossing anywhere else
    sample_ts = np.linspace(t_vals[0], t_end, 200)
    lam_vals = [lam_sol(tt)[0] for tt in sample_ts]
    if any(lv < 0.0 for lv in lam_vals):
        return "metastable"

    return "stable"

###############################################################################
# 5) Extended scanning function for domain plot
###############################################################################

def scan_mt_mH_for_plot(mt_values, mH_values):
    """
    Builds a 2D grid over the arrays mH_values and mt_values and classifies each point.

    Returns:
        class_array: 2D numpy array with shape (len(mH_values), len(mt_values)).
                     The codes are 0=unstable, 1=metastable, 2=stable.
    """
    class_array = np.zeros((len(mH_values), len(mt_values)))

    for i, mH in enumerate(mH_values):
        for j, mt in enumerate(mt_values):
            # Improve speed by ignoring some simple regions and focus on the complicated boundary
            if mt < 50 + 0.6*mH:
                status = "stable"
                class_array[i, j] = 2
                continue
            elif mt > 150 + 0.5*mH:
                status = "unstable"
                class_array[i, j] = 0
                continue

            rge_sol = setup_and_run_RGE(mt, mH)
            status = check_stability(rge_sol)
            if status == "unstable":
                class_array[i, j] = 0
            elif status == "metastable":
                class_array[i, j] = 1
            else:  # "stable"
                class_array[i, j] = 2

    return class_array

###############################################################################
# 6) Plotting function for couplings (optional)
###############################################################################

def plot_couplings(rge_solutions, plot):
    """
    Plots the running of g_s, g2, gY, top Yukawa coupling Y_top, and λ as functions of μ.
    """
    Y_sol = rge_solutions["Y_sol"].sol
    lam_sol = rge_solutions["lam_sol"].sol
    mu2_sol = rge_solutions["mu2_sol"].sol

    t_vals = rge_solutions["lam_sol"].t
    t_start = t_vals[0]
    t_end = t_vals[-1]
    t_arr = np.linspace(t_start, t_end, 300)
    mu_arr = mZ * np.exp(t_arr)

    gs_vals = [gs_running_analytic(t) for t in t_arr]
    g2_vals = [g2_running_analytic(t) for t in t_arr]
    gY_vals = [gY_running_analytic(t) for t in t_arr]
    Y_vals = [Y_sol(tt)[0]/np.sqrt(2) for tt in t_arr]
    lam_vals = [lam_sol(tt)[0] for tt in t_arr]
    mu2_vals = [-mu2_sol(tt)[0]/125.**2*2. for tt in t_arr]

    plt.figure(figsize=(plot['singleplot']['height'], plot['singleplot']['height']))
    plt.plot(mu_arr, gs_vals, label='$g_s$', color=plot['color'][0])
    plt.plot(mu_arr, g2_vals, label='$g_2$', color=plot['color'][1])
    plt.plot(mu_arr, gY_vals, label='$g_Y$', color=plot['color'][2])
    plt.plot(mu_arr, Y_vals, label='$Y_{t}$', color=plot['color'][3])
    plt.plot(mu_arr, lam_vals, label='$\\lambda$', color=plot['color'][4])
    plt.plot(mu_arr, mu2_vals, label='$2 \\mu^2/(125\ \mathrm{GeV})^2$', color=plot['color'][5])

    plt.xscale('log')
    plt.xlabel(r'$\mu \ [\mathrm{GeV}]$')
    plt.xlim(1e2, 1e19)
    #plt.ylabel(r'$\mathrm{Coupling\ value}$')
    plt.legend(loc='best')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../../LaTeX/Images/running_couplings.pdf')
    plt.show()

###############################################################################
# 7) Effective Potential
###############################################################################
'''
def Veff(phi, mu2, lam, g2, gY, Yt, muR, order=1):
    #print("mu2 = ", mu2)
    #print("lam = ", lam)
    #print("sqrt(mu2/lam) = ", np.sqrt(mu2/lam))
    #print("mH = ", np.sqrt(2*mu2))
    CW = CZ = 5./6.
    Ct = CH = 3./2.
    nW = 6.
    nZ = 3.
    nt = -12.
    nH = 1.
    mWeff2 = 1./4.*g2**2*phi**2
    mZeff2 = 1./4.*(g2**2 + gY**2)*phi**2
    mteff2 = Yt**2*phi**2/np.sqrt(2)
    mHeff2 = 3*lam*phi**2 + mu2
    Ci = [Ct, CW, CZ, CH]
    ni = [nt, nW, nZ, nH]
    mieff2 = [mteff2, mWeff2, mZeff2, mHeff2]
    print("sqrt mteff2 = ", np.sqrt(mteff2))
    V0 = mu2*phi**2/2 + lam*phi**4/4
    V1 = 0
    for i in range(0, 4):
        V1 += ni[i]/64./np.pi**2*mieff2[i]**2*(np.log(mieff2[i]/muR**2) - 0.*Ci[i])
    if order==0:
        return V0
    else:
        return V0 + V1
'''
def Veff(
    phi, mu2, lam, g2, gY, Yt, muR
):
    """
    One-loop Higgs potential with the omitted C terms included:
      V1 = sum_i [ n_i/(64π²) * (m^2_i)^2 * ( ln(m^2_i / muR^2 ) - C_i ) ].

    The default C_i are:
      Top (Ct) = 3/2,
      W   (CW) = 5/6,
      Z   (CZ) = 5/6,
      Higgs (CH)=3/2.
    """

    # Tree-level piece:
    #   V0 = (mu2/2)*phi^2 + (lam/4)*phi^4
    V0 = 0.5 * mu2 * phi**2 + 0.25 * lam * phi**4

    # Field-dependent masses
    # same definitions from your snippet
    mW2 = 0.25 * g2**2         * phi**2
    mZ2 = 0.25 * (g2**2+gY**2) * phi**2
    mt2 = (Yt**2 / 2.) * phi**2
    mH2 = 3.0 * lam * phi**2 + mu2

    # Multiplicities
    nW, nZ, nt, nH = 6.0, 3.0, -12.0, 1.0
    # Corresponding C factors:
    #   i=0 => top, i=1 => W, i=2 => Z, i=3 => Higgs
    #   so let's keep that order: (top, W, Z, Higgs)
    CW=0 #5.0/6.0
    CZ=0 #5.0/6.0
    Ct=0 #3.0/2.0
    CH=0 #3.0/2.0
    C_array = [Ct, CW, CZ, CH]
    n_array = [nt, nW, nZ, nH]
    m2_array= [mt2, mW2, mZ2, mH2]

    # One-loop correction
    V1 = 0.0
    for (m2_eff, n_fac, c_fac) in zip(m2_array, n_array, C_array):
        if m2_eff <= 0.0:
            # skip or handle negative mass^2 carefully if needed
            continue
        # The piece: n_fac/(64 π²) * m2_eff^2 [ ln(m2_eff / muR^2 ) - c_fac ]
        V1 += n_fac/(64.0*math.pi**2) * (m2_eff**2) * (
            math.log(m2_eff / muR**2) - c_fac
        )

    return V0 + V1

def dVdphi(
    phi, mu2, lam, g2, gY, Yt, muR,
):
    """
    First derivative dV/dphi for the one-loop potential including C terms.
    We'll do an explicit chain rule for each mass sector.
    """

    # We'll do:  d2V/dphi2 ≈ [dV/dphi(phi+eps) - dV/dphi(phi-eps)] / (2eps)
    eps = phi * 1e-4

    dV_minus = Veff(phi - eps, mu2, lam, g2, gY, Yt, muR)
    dV_plus  = Veff(phi + eps, mu2, lam, g2, gY, Yt, muR)
    return (dV_plus - dV_minus)/(2.0*eps)

def d2Vdphi2(
    phi, mu2, lam, g2, gY, Yt, muR
):
    """
    Second derivative d²V/dphi² for the one-loop potential with C terms.
    For brevity, we do a small finite-difference of dV/dphi rather
    than writing out the full symbolic expression (which is lengthy).
    """

    # We'll do:  d2V/dphi2 ≈ [dV/dphi(phi+eps) - dV/dphi(phi-eps)] / (2eps)
    eps = phi * 1e-4

    dV_minus = dVdphi(phi - eps, mu2, lam, g2, gY, Yt, muR)
    dV_plus  = dVdphi(phi + eps, mu2, lam, g2, gY, Yt, muR)
    return (dV_plus - dV_minus)/(2.0*eps)

def residual_equations(vars, phi_target, mH_target, g2, gY, Yt, muR):
    """
    We want:
      eq1 = dV/dphi(phi_target) = 0
      eq2 = d2V/dphi^2(phi_target) - mH_target^2 = 0

    Inputs:
      vars = [mu2, lam]
      phi_target : float, typically v ~ 246 GeV
      mH_target  : float, typically Higgs mass ~ 125 GeV
      g2, gY, Yt, muR, order: see above.
    """
    mu2_guess, lam_guess = vars
    eq1 = dVdphi(phi_target, mu2_guess, lam_guess, g2, gY, Yt, muR)
    eq2 = d2Vdphi2(phi_target, mu2_guess, lam_guess, g2, gY, Yt, muR) - mH_target**2
    return [eq1, eq2]

def solve_mu2_lam(phi_target=246.0, mH_target=125.0,
                  g2=0.65, gY=0.36, Yt=0.95, muR=91):
    """
    Numerically solves for (mu2, lam) so that:
      1) dV/dphi(phi=v) = 0
      2) d²V/dphi²(phi=v) = mH^2
    for the one-loop potential described above.

    Returns (mu2_sol, lam_sol).
    """
    from scipy.optimize import root

    # A decent initial guess from tree-level relations:
    # at tree level, mu^2 ~ -mH^2/2, lam ~ mH^2/(2v^2).
    mu2_init = -0.5 * mH_target**2
    lam_init = mH_target**2 / (2.0 * phi_target**2)
    guess = [mu2_init, lam_init]

    # solve
    sol = root(
        lambda x: residual_equations(x, phi_target, mH_target, g2, gY, Yt, muR),
        guess,
        method='hybr',
        options={'xtol': 1e-7, 'factor':  0.2}
    )
    if not sol.success:
        print("mH = ", mH_target, "mt = ", vev*Yt/np.sqrt(2.))
        return mu2_init, lam_init
        #raise RuntimeError("Could not solve for (mu2, lam). Message: " + sol.message)

    mu2_sol, lam_sol = sol.x
    return mu2_sol, lam_sol

def Veff_resummed(phi, t,
                  mu2_running_sol, lam_running_sol,
                  Y_sol,
                  g2_running_analytic, gY_running_analytic):
    """
    Compute the RG-resummed effective potential at field value phi and scale parameter t,
    given:
      • mu2_running_sol: solution object for mu²(t) from solve_ivp
      • lam_running_sol: solution object for λ(t)
      • Y_sol: solution object for the top Yukawa coupling Y_top(t)
      • g2_running_analytic(t), gY_running_analytic(t): gauge coupling runnings
      • t = ln(μ / μ₀) (for example, μ₀ = mZ)

    The resummed potential takes the form:
        V_resummed(phi, t) = ½ * mu²(t)*G(t)² * φ²  +  ¼ * λ(t)*G(t)⁴ * φ⁴,
    where G(t) = exp[- ∫₀ᵗ γ(x) dx], and
      γ(t) = [ -9*g2(t)² - 3*gY(t)² + 12*Y_top(t)² ] / (64 π²).

    Parameters
    ----------
    phi : float
        The field value at which to evaluate the resummed potential.
    t : float
        The RG time, t = ln(μ / μ₀).
    mu2_running_sol : OdeSolution
        From solve_ivp for d(mu²)/dt. Access mu²(t) via mu2_running_sol.sol(t)[0].
    lam_running_sol : OdeSolution
        From solve_ivp for d(λ)/dt. Access λ(t) via lam_running_sol.sol(t)[0].
    Y_sol : OdeSolution
        From solve_ivp for the top Yukawa coupling. Y_sol.sol(t)[0] gives Y_top(t).
    g2_running_analytic : callable
        A function g2_running_analytic(t) -> g₂(t).
    gY_running_analytic : callable
        A function gY_running_analytic(t) -> gʸ(t).

    Returns
    -------
    float
        The value of the RG-improved effective potential V_resummed(phi, t).
    """

    # 1) Extract mu²(t) and λ(t) from the solutions
    mu2_t = mu2_running_sol.sol(t)[0]
    lam_t = lam_running_sol.sol(t)[0]

    # 2) Define a local gamma function that depends on Y_sol, g2, gY
    def gamma_function_local(tval):
        """
        gamma(t) = [ -9*g2(t)² - 3*gY(t)² + 12*Y_top(t)² ] / (64 π²).
        """
        g2_val = g2_running_analytic(tval)
        gY_val = gY_running_analytic(tval)
        Y_top_val = Y_sol.sol(tval)[0]
        return ( -9.0*g2_val**2 - 3.0*gY_val**2 + 12.0*Y_top_val**2 ) / (64.0*np.pi**2)

    # 3) Compute the wavefunction factor G(t) = exp(- ∫₀ᵗ gamma_function_local(x) dx)
    def G_of_t(tval):
        integral_val, _ = quad(gamma_function_local, 0.0, tval, epsrel=1e-8)
        return np.exp(-integral_val)

    Gval = G_of_t(t)
    print("t = ", t, ", G = ", Gval)
    # 4) Construct the resummed potential:
    #    V_resummed(phi,t) = ½ mu²(t)*G(t)² * φ² + ¼ λ(t)*G(t)⁴ * φ⁴
    part1 = 0.5 * mu2_t * (Gval**2) * (phi**2)
    part2 = 0.25 * lam_t * (Gval**4) * (phi**4)

    return part1 + part2

###############################################################################
# 7) Example main: domain plot + couplings
###############################################################################

def main():
    with open('../style.py', 'r') as f:
        s = f.read()
        plot = eval(s)

    mH = 125
    mt = 173
    mu2_tree = -mH**2/2.
    lam_tree = mH**2 / (2.0 * vev**2)
    muR = mZ
    mu2_found, lam_found = solve_mu2_lam(
        phi_target=vev, mH_target=mH,
        g2=g2_mZ_default, gY=gY_mZ_default, Yt=np.sqrt(2)*mt/vev, muR=mZ
    )
    print(f"Solved (mu2, lam) = ({mu2_found:.5g}, {lam_found:.5g}) at one-loop")
    print(f"mu2_initial , lam_initial) = ({(mH**2/2.):.5g}, {(mH**2/2./vev**2):.5g}) at one-loop")

    # Check the conditions
    eqs = residual_equations([mu2_found, lam_found], vev, mH, g2_mZ_default, gY_mZ_default, np.sqrt(2)*mt/vev, mZ)
    dV_phi = eqs[0]
    d2V_phi = eqs[1] + 125**2  # eq2 = d2V - mH^2
    print(f"Residuals: dV/dphi(v)={dV_phi:.5g}, [d²V/dphi²(v) - mH²]={eqs[1]:.5g}")
    print(f"Thus, d²V/dphi²(v) = {d2V_phi:.5g} GeV^2 vs. mH^2={mH**2:.5g} GeV^2")

    phi_min, phi_max = 1.0, 1.0e20
    n_pts = 200
    phi_vals = np.logspace(np.log10(phi_min), np.log10(phi_max), n_pts)
    phi_vals = np.linspace(0, 400, n_pts)

    V_tree_vals      = []
    V_notres_vals    = []
    V_resum_vals     = []
    single_sol = setup_and_run_RGE(mt, mH)
    for phi in phi_vals:
        # Tree-level
        V_t = mu2_tree/2*phi**2 + lam_tree/4*phi**4
        V_tree_vals.append(V_t)

        # 1-loop not-resummed
        V_nr = Veff(phi, mu2_found, lam_found, g2_mZ_default, gY_mZ_default, np.sqrt(2)*mt/vev, muR)
        V_notres_vals.append(V_nr)

        # Resummed: define t = ln(phi / mZ)
        if phi>0.0:
            t_val = np.log(phi/mZ)
        else:
            t_val = -10.0  # some fallback

        V_r  = Veff_resummed(
            phi, np.log(phi/muR),
            single_sol["mu2_sol"], single_sol["lam_sol"],
            single_sol["Y_sol"],
            g2_running_analytic, gY_running_analytic
        )
        V_resum_vals.append(V_r)

    #------------------------------------------------------------------
    # (D) Plot
    #------------------------------------------------------------------
    plt.figure(figsize=(8,6))

    plt.plot(phi_vals, V_tree_vals,      label='Tree-level')
    plt.plot(phi_vals, V_notres_vals,    label='1-loop not resummed')
    plt.plot(phi_vals, V_resum_vals,     label='RG-resummed')

    # Log scales help if the potential spans many magnitudes
    #plt.xscale('log')
    #plt.yscale('symlog')
    # Potential can sometimes cross zero, so you could do plt.yscale('symlog') if it crosses zero
    # or just do a normal linear y-scale:
    # For demonstration, let's do a normal scale.
    # If your potential gets huge, you might choose 'symlog'.
    # plt.yscale('symlog', linthresh=1e2)

    plt.xlabel(r'$\phi$ (GeV)')
    plt.ylabel(r'$V(\phi)$ (GeV$^4$)')
    plt.title("Comparison of Tree-level, 1-loop, and RG-resummed Potentials")
    plt.legend(loc='best')
    plt.grid(True)

    plt.tight_layout()
    plt.show()

    '''
    mH = 125.
    mt = 0.01
    vev = vev
    mu2_mZ = 0.#(mH**2) / 2.0
    lam_mZ = (mH**2) / (2.0 * vev**2)
    Y_mZ = (mt / vev)
    print("lam_mZ = ", lam_mZ)
    phis = np.linspace(0, 600, 100)
    Veffs1 = []
    Veffs0 = []
    for phi in phis:
        Veffs1.append(Veff(phi, mu2_mZ, lam_mZ, g2_mZ_default, gY_mZ_default, Y_mZ, mZ, order=1))
        Veffs0.append(Veff(phi, mu2_mZ, lam_mZ, g2_mZ_default, gY_mZ_default, Y_mZ, mZ, order=0))

    plt.figure(figsize=(plot['singleplot']['width'], plot['singleplot']['height']))
    plt.plot(phis, Veffs0, color=plot['color'][0], label="V0")
    plt.plot(phis, Veffs1, color=plot['color'][1], label="V0 + V1")
    plt.legend()
    plt.show()

    return
    '''


    # Choose mt and mH ranges
    mt_vals = np.linspace(50, 250, 40)   # 7 points from 150 to 180
    mH_vals = np.linspace(0, 250, 50)   # 7 points from 110 to 140

    # Build clasification grid
    class_array = scan_mt_mH_for_plot(mt_vals, mH_vals)

    # Create a custom colormap for imshow
    cmap_custom = plt.matplotlib.colors.ListedColormap(['red', 'yellow', 'limegreen'])
    bounds = [0, 1, 2, 3]
    norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap_custom.N)

    # Plot the domain
    plt.figure(figsize=(plot['singleplot']['height'], plot['singleplot']['height']))
    im = plt.imshow(
        class_array.transpose(),
        origin='lower',
        extent=[mH_vals[0], mH_vals[-1], mt_vals[0], mt_vals[-1]],
        interpolation='nearest',
        cmap=cmap_custom,
        norm=norm,
        aspect='auto'
    )
    #plt.colorbar(
    #    im,
    #    boundaries=bounds,
    #    ticks=[0.5, 1.5, 2.5],
    #    label='Stability (0=unstable, 1=metastable, 2=stable)'
    #)
    plt.ylabel(r'$m_t\ [\mathrm{GeV}]$')
    plt.xlabel(r'$m_H\ [\mathrm{GeV}]$')
    plt.tight_layout()
    plt.savefig('../../LaTeX/Images/stability_scan.pdf')
    plt.show()



    # Demonstrate couplings for a reference point
    mt_ref, mH_ref = 173.0, 125.0
    single_sol = setup_and_run_RGE(mt_ref, mH_ref)
    classification = check_stability(single_sol)
    print(f"For (mt={mt_ref:.1f}, mH={mH_ref:.1f}), classification = {classification}")
    plot_couplings(single_sol, plot)


if __name__ == "__main__":
    main()