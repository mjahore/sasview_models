import numpy as np  # type: ignore
from numpy import cos, pi, inf, power, errstate
from sasmodels.special import sas_gamma, sas_gammainc, sas_sinx_x, sas_3j1x_x

name = "csc"
title = "Core Shell Chain (CSC)"
description = """
"""
category = "shape:sphere"

#             [ "name",       "units",         default, [lower, upper], "type",   "description"],
parameters = [["sld",         "1e-6/Ang^2",    3.5,     [-inf, inf],    "sld",    "Core scattering length density"],
              ["sld_shell",   "1e-6/Ang^2",    -0.022,  [-inf, inf],    "sld",    "Shell scattering length density"],
	      ["sld_poly",    "1e-6/Ang^2",    1.269,   [-inf, inf],    "sld",    "Grafted polymer scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",    6.37,    [-inf, inf],    "sld",    "Solvent scattering length density"],
              ["radius",      "Ang",           70,      [0, inf],       "volume", "Core radius"],
              ["t_shell",     "Ang",           20,      [0, inf],       "volume", "Shell thickness"],
              ["poly_sig",    "chains/nm^2",   0.33,    [0, inf],       "",       "Polymer grafting density"],
              ["C_infty",     "None",          12,      [1, inf],       "",       "Characteristic ratio"],
              ["M0",          "g/mol",         113,     [1, inf],       "",       "Monomer molar mass"],
              ["Mn",          "g/mol",         11.18,   [0, inf],       "",       "Polymer molar mass"],
              ["nu",          "None",          0.50,    [0.25, 1.0],    "",       "Flory exponent"],
              ["v",           "Ang^3",         162,     [0, inf]   ,    "volume", "Kuhn monomer volume"],
             ]

def Iq(q,
       sld,
       sld_shell,
       sld_poly,
       sld_solvent,
       radius,
       t_shell,
       poly_sig,
       C_infty,
       M0,
       Mn,
       nu,
       v):

    # Bond angles
    theta0 = 68.0 * pi/180.0

    # Kuhn length
    b = C_infty * 1.54 / cos(theta0/2.0)

    # Deg. of polymerization
    N = (Mn/M0) * cos(theta0/2.0)/C_infty

    # Volume of core regions:
    Rcoreshell = radius + t_shell
    Vcore      = 4.0/3.0 * pi * radius**3
    Vcoreshell = 4.0/3.0 * pi * (radius+t_shell)**3

    # Number of grafted chains per core:
    Ng = poly_sig * 4.00 * pi * (0.1 * Rcoreshell) * (0.1 * Rcoreshell)

    
    Vtotal = Vcoreshell + Ng*N*v;

    # One over excl. vol. parm.:
    onu  = 1.0 / nu
    o2nu = 1.0 / 2.0 / nu

    # Propagator function:
    Ea = sas_sinx_x(q*(Rcoreshell))

    # Polymer size variable
    Usub = (q*b)**2 * N**(2*nu) / 6.0

    # Form factor amplitude of core-shell sphere:
    with errstate(divide='ignore'):
        Fs = (sld - sld_shell)*Vcore*sas_3j1x_x(q*radius) + (sld_shell - sld_solvent)*Vcoreshell*sas_3j1x_x(q*Rcoreshell)


    # Form factor amplitude of the polymer:
    with errstate(divide='ignore'):
        Fp = N*v*o2nu*power(Usub, -o2nu) * sas_gamma(o2nu) * sas_gammainc(o2nu, Usub)

    # Form factor of the polymer (Pp(q) is not simply Fp(q)^2!!):
    with errstate(divide='ignore'):
        Pp = (N*v)**2 * (onu * power(Usub, -o2nu)*sas_gamma(o2nu)*sas_gammainc(o2nu, Usub) - onu * power(Usub, -onu)*sas_gamma(onu)*sas_gammainc(onu,Usub))

    # Combine all terms to form intensity:
    #
    # Term 1: Core-shell particle:
    inten = Fs * Fs

    # Term 2: Polymer
    inten = inten + Ng * (sld_poly - sld_solvent) * (sld_poly - sld_solvent) * Pp

    # Term 3: Particle/polymer crossterm:
    inten = inten + 2.0 * Ng * (sld_poly - sld_solvent) * Fs * Ea * Fp

    # Term 4: Polymer/polymer crossterm:
    inten = inten + Ng * (Ng - 1) * (sld_poly - sld_solvent) * (sld_poly - sld_solvent) * Fp * Ea * Ea * Fp
    with errstate(divide='ignore'):
        inten = inten * 1.0e-4 / Vtotal

    return inten

Iq.vectorized = True  # Iq accepts an array of q values

# Effective radius for S(q) calculations:
def ER(radius, t_shell):
    return radius + t_shell

# Volume ratio
def VR(radius, t_shell):

    whole = 4.0/3.0*pi*(radius+t_shell)**3
    core  = 4.0/3.0*pi*radius**3
    return whole, whole-core

def random():
    pars = dict(
	radius   = np.random.uniform(20,200),
        t_shell  = np.random.uniform(10,100),
        poly_sig = np.random.uniform(0,2),
        rg       = np.random.unifomr(20,150),
        nu       = np.random.uniform(0.3,0.6),
        v_poly   = np.random.uniform(10,50),
    )
    return pars

#demo = dict(scale=1, background,0,
#            sld=3.0, sld_shell=1.0, sld_poly = 1.0, sld_solvent=4.3,
#            radius=50, t_shell=20, poly_sig=0.50, rg=70, nu=0.5, v_poly=30)
