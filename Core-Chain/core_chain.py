r"""


Definition
----------

(Coming soon)

References
----------

J. S. Pedersen, Adv. Colloid Interface Sci. 70, 171-210 (1997).
M.J.A. Hore, J. Ford, K. Ohno, R. J. Composto, B. Hammouda, Macromolecules 46, 9341-9348 (2013).
"""

import numpy as np  # type: ignore
from numpy import pi, inf, power, errstate
from sasmodels.special import sas_gammainc, sas_sinx_x, sas_3j1x_x

name = "core_chain"
title = "Spherically symmetric core with grafted polymer chains."
description = """
"""
category = "shape:sphere"

#             [ "name",       "units",         default, [lower, upper], "type",   "description"],
parameters = [["sld",         "1e-6/Ang^2",    3.5,     [-inf, inf],    "sld",    "Core scattering length density"],
	      ["sld_poly",    "1e-6/Ang^2",    1.0,     [-inf, inf],    "sld",    "Grafted polymer scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",    4.4,     [-inf, inf],    "sld",    "Solvent scattering length density"],
              ["radius",      "Ang",           60,      [0, inf],       "volume", "Core radius"],
              ["poly_sig",    "chains/nm^2",   0.50,    [0, inf],       "",       "Polymer grafting density"],
              ["rg",          "Ang",           40,      [0, inf],       "volume", "Grafted polymer radius of gyration"],
              ["nu",          "None",          0.50,    [0.25, 1.0],    "",       "Grafted polymer excluded volume parameter"],
              ["v_poly",      "1/Ang^3",       30,      [0, inf]   ,    "volume", "Volume of one polymer"],
             ]

def Iq(q,
       sld,
       sld_poly,
       sld_solvent,
       radius=60,
       poly_sig=0.50,
       rg=40,
       nu=0.5,
       v_poly=30):
    """
    :param q:              Input q-value
    :param sld:		   Core scattering length density
    :param sld_poly:       Polymer scattering length density
    :param sld_solvent:    Solvent scattering length density
    :param radius:         Core radius
    :param poly_sig:       Polymer grafting density
    :param rg:             Grafted polymer radius of gyration
    :param nu:             Grafted polymer excluded volume parameter
    :param v_poly:         Volume of one polymer 
    :return:               Calculated intensity
    """

    # Volume of core regions:
    Vcore      = 4.0/3.0 * pi * radius**3
    Vtotal     = Vcore + v_poly

    # One over excl. vol. parm.:
    onu  = 1.0 / nu
    o2nu = 1.0 / 2.0 / nu

    # Propagator function:
    Ea = sas_sinx_x(q*radius)

    # Number of grafted chains per core:
    Ng = poly_sig * 4.00 * pi * (0.1 * radius) * (0.1 * radius)

    # Polymer size variable
    Usub = (q*rg)**2 * (2.0*nu + 1.0) * (2.0*nu + 2.0) / 6.0

    # Form factor amplitude of core-shell sphere:
    with errstate(divide='ignore'):
        Fs = 3.0*(sld - sld_solvent)*Vcore*sas_3j1x_x(q*radius)


    # Form factor amplitude of the polymer:
    with errstate(divide='ignore'):
        Fp = o2nu*power(Usub, -o2nu) * sas_gammainc(o2nu, Usub) 


    # Form factor of the polymer (Pp(q) is not simply Fp(q)^2!!):
    with errstate(divide='ignore'):
        Pp = onu * power(Usub, -o2nu)*sas_gammainc(o2nu, Usub) - onu * power(Usub, -onu)*sas_gammainc(onu,Usub)

    # Combine all terms to form intensity:
    #
    # Term 1: Core-shell particle:
    inten = Fs * Fs

    # Term 2: Polymer
    inten = inten + Ng * v_poly * v_poly * (sld_poly - sld_solvent) * (sld_poly - sld_solvent) * Pp

    # Term 3: Particle/polymer crossterm:
    inten = inten + 2.0 * Ng * v_poly * (sld_poly - sld_solvent) * Fs * Ea * Fp

    # Term 4: Polymer/polymer crossterm:
    inten = inten + Ng * (Ng - 1) * v_poly * v_poly * (sld_poly - sld_solvent) * (sld_poly - sld_solvent) * Fp * Ea * Ea * Fp
    with errstate(divide='ignore'):
        inten = inten * 1.0e-6 * 1.0e-6 * 1.0e8 / Vtotal

    return inten

Iq.vectorized = True  # Iq accepts an array of q values

# Effective radius for S(q) calculations:
def ER(radius, rg, v_poly):
    return radius + rg

# VR defaults to 1.0

def random():
    pars = dict(
	radius   = np.random.uniform(20,200),
        poly_sig = np.random.uniform(0,2),
        rg       = np.random.unifomr(20,150),
        nu       = np.random.uniform(0.3,0.6),
        v_poly   = np.random.uniform(10,50),
    )
    return pars

#demo = dict(scale=1, background,0,
#            sld=3.0, sld_shell=1.0, sld_poly = 1.0, sld_solvent=4.3,
#            radius=50, t_shell=20, poly_sig=0.50, rg=70, nu=0.5, v_poly=30)
