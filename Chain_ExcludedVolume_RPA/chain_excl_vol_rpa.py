r"""
Definition
----------

This model combines the random phase approximation (RPA) with the form factor for a (linear) polymer with excluded volume
to analyze a dilute polymer solution.

References
----------

B. Hammouda, "Form Factors for Branched Polymers with Excluded Volume", J. of Research of NIST, 121, 139-164 (2016).

Authorship and Verification
----------------------------

* **Author:** Michael J. A. Hore **Date:** 5 Mar 2018
"""

import numpy as np
from numpy import inf, errstate, power,exp, sqrt
from sasmodels.special import sas_gammainc, sas_gamma

name = "poly_excl_vol_rpa"
title = "Polymer with excluded volume, RPA"
description = """\
      List of default parameters:
      scale = Scaling factor
      nu = Excluded volume parameter
      b = Kuhn length
      n = Degree of polymerization
      background = Incoherent background"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default,       [lower, upper], "type", "description"],
parameters = [
              ["phi_p",      "",                 0.01,      [0, 1],         "",     "Polymer volume fraction"],
              ["nu",       "",                 0.5,       [0, 1],         "",     "Excluded volume parameter"],
              ["b",        "Ang",              7.0,       [1, inf],       "",     "Kuhn length"],
              ["n",        "",                30.0,       [1, inf],       "",     "Degree of polymerization"],
              ["sldp",     "1e-6/Ang^2",       1.4,       [-inf,inf],     "sld",  "Polymer SLD"],
              ["slds",     "1e-6/Ang^2",       6.7,       [-inf,inf],     "sld",  "Solvent SLD"],
              ["vm",       "Ang^3",            178,       [1, inf],       "",      "Monomer volume"],
              ["vs",       "Ang^3",            179,       [1, inf],       "",     "Solvent volume"],
              ["chi",      "",                 0.5,       [-inf,inf],     "",     "Flory-Huggins parameter"],
             ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q,
       phi_p,
       nu,
       b,
       n,
       sldp,
       slds,
       vm,
       vs,
       chi):
    """
    :param q:              Input q-value
    :param phi_p:            Polymer volume fraction
    :param nu:             Excluded volume parameter
    :param b:              Kuhn length
    :param n:              Degree of polymerization
    :param sldp:           Polymer scattering length density
    :param slds:           Solvent scattering length density
    :param vm:             Monomer volume
    :param vs:             Solvent molecule volume
    :param chi:            Flory-Huggins parameter
    :return:               Calculated intensity
    """

    # Excl. Vol. Parameters
    onu  = 1.0/nu
    o2nu = 1.0/2.0/nu
 
    # Polymer chain form factor:
    U = q**2 * b**2 * power(n, 2.0*nu) / 6.0
    Pp = onu * power(U, -o2nu) * sas_gamma(o2nu) * sas_gammainc(o2nu, U) - onu * power(U, -onu) * sas_gamma(onu) * sas_gammainc(onu, U)

    # Form the structure factor according to the RPA:
    Sq = power(n*phi_p*vm*Pp, -1.00) + power(vs * (1.0 - phi_p), -1.00) - 2.00*chi/sqrt(vm*vs)
    inten = (sldp - slds) * (sldp - slds) * 1e-4 * power(Sq, -1.0)

    return inten

Iq.vectorized = True  # Iq accepts an array of q values

def random():
    pars = dict(
        scale=1,
        phi = np.random.uniform(0.001, 0.999),
        nu  = np.random.uniform(0.3,0.6),
        b   = np.random.uniform(7,15),
        n   = np.random.uniform(20,200),
        vm  = np.random.uniform(50,200),
        vs  = np.random.uniform(50,200),
        chi = np.random.uniform(-0.5,1.5) 
    )
    return pars

demo = dict(scale=1, background=0,
            nu = 0.5, b=7, n=40)
