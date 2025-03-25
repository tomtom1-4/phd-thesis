#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, quad


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
vev_default = (
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
    vev = vev_default
    mu2_mZ = - (mH**2) / 2.0
    lam_mZ = (mH**2) / (2.0 * vev**2)
    Y_mZ = (mt / vev) * np.sqrt(2.0)

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
    if lam_high < 0.0:
        return "unstable"

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
            if mt < 75 + 0.6*mH:
                status = "stable"
                class_array[i, j] = 2
                continue
            elif mt > 125 + 0.45*mH:
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

    t_vals = rge_solutions["lam_sol"].t
    t_start = t_vals[0]
    t_end = t_vals[-1]
    t_arr = np.linspace(t_start, t_end, 300)
    mu_arr = mZ * np.exp(t_arr)

    gs_vals = [gs_running_analytic(t) for t in t_arr]
    g2_vals = [g2_running_analytic(t) for t in t_arr]
    gY_vals = [gY_running_analytic(t) for t in t_arr]
    Y_vals = [Y_sol(tt)[0] for tt in t_arr]
    lam_vals = [lam_sol(tt)[0] for tt in t_arr]

    plt.figure(figsize=(plot['singleplot']['width'], plot['singleplot']['height']))
    plt.plot(mu_arr, gs_vals, label='$g_s$', color=plot['color'][0])
    plt.plot(mu_arr, g2_vals, label='$g_2$', color=plot['color'][1])
    plt.plot(mu_arr, gY_vals, label='$g_Y$', color=plot['color'][2])
    plt.plot(mu_arr, Y_vals, label='$Y_{t}$', color=plot['color'][3])
    plt.plot(mu_arr, lam_vals, label='$\\lambda$', color=plot['color'][4])

    plt.xscale('log')
    plt.xlabel(r'$\mu \ [\mathrm{GeV}]$')
    plt.xlim(1e2, 1e19)
    #plt.ylabel(r'$\mathrm{Coupling\ value}$')
    plt.legend(loc='best')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig('../../Images/running_couplings.pdf')
    plt.show()

###############################################################################
# 7) Example main: domain plot + couplings
###############################################################################

def main():
    with open('../style.py', 'r') as f:
        s = f.read()
        plot = eval(s)
    # Choose mt and mH ranges
    mt_vals = np.linspace(50, 250, 200)   # 7 points from 150 to 180
    mH_vals = np.linspace(0, 250, 200)   # 7 points from 110 to 140

    # Build clasification grid
    class_array = scan_mt_mH_for_plot(mt_vals, mH_vals)

    # Create a custom colormap for imshow
    cmap_custom = plt.matplotlib.colors.ListedColormap(['red', 'yellow', 'limegreen'])
    bounds = [0, 1, 2, 3]
    norm = plt.matplotlib.colors.BoundaryNorm(bounds, cmap_custom.N)

    # Plot the domain
    plt.figure(figsize=(plot['singleplot']['width'], plot['singleplot']['height']))
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
    plt.savefig('../../Images/stability_scan.pdf')
    plt.show()

    # Demonstrate couplings for a reference point
    mt_ref, mH_ref = 173.0, 125.0
    single_sol = setup_and_run_RGE(mt_ref, mH_ref)
    classification = check_stability(single_sol)
    print(f"For (mt={mt_ref:.1f}, mH={mH_ref:.1f}), classification = {classification}")
    plot_couplings(single_sol, plot)

if __name__ == "__main__":
    main()