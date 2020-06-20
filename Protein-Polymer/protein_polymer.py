
import numpy as np  # type: ignore
from numpy import pi, inf, power, errstate
from scipy.special import gammainc, gamma, j0, j1

name = "protein_polymer"
title = "Protein-polymer conjugate"
description = """
"""
category = "shape:sphere"

structure_factor = True

#             [ "name",       "units",         default, [lower, upper], "type",   "description"],
parameters = [["sld1",        "1e-6/Ang^2",    1.0,     [-inf, inf],    "sld",    "Protein scattering length density"],
              ["sld2",        "1e-6/Ang^2",    1.0,     [-inf, inf],    "sld",    "Polymer scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",    4.4,     [-inf, inf],    "sld",    "Solvent scattering length density"],
              ["rg1",         "Ang",           40,      [1, inf],       "volume", "Radius of gyration of protein"],
              ["rg2",         "Ang",           40,      [1, inf],       "volume", "Radius of gyration of polymer"],
              ["nu1",         "None",          0.50,    [0.25, 1.0],    "",       "Excluded volume parameter of protein"],
              ["nu2",         "None",          0.50,    [0.25, 1.0],    "",       "Excluded volume parameter of polymer"],
              ["v1",          "Ang^3",         30,      [0, inf]   ,    "volume", "Volume of protein"],
              ["v2",          "Ang^3",         30,      [0, inf]   ,    "volume", "Volume of polymer"],
             ]

def Iq(q,
       sld1,
       sld2,
       sld_solvent,
       rg1,
       rg2,
       nu1,
       nu2,
       v1,
       v2):
    """
    :param q:              Input q-value
    :param sld1:           Chain region 1 scattering length density
    :param sld2:           Chain region 2 scattering length density
    :param sld_solvent:    Solvent scattering length density
    :param rg1:            Radius of gyration of chain in region 1
    :param rg2:            Radius of gyration of chain in region 2
    :param nu1:            Excluded volume parameter of chain in region 1
    :param nu2:            Excluded volume parameter of chain in region 2
    :param v1:             Volume of polymer in region 1
    :param v2:             Volume of polymer in region 2
    :return:               Calculated intensity
    """

    # Volume of core regions:
    Vtotal     = (v1 + v2)

    # One over excl. vol. parm.:
    onu1  = 1.0 / nu1
    o2nu1 = 1.0 / 2.0 / nu1
    onu2  = 1.0 / nu2
    o2nu2 = 1.0 / 2.0 / nu2

    # Propagator function:
    E1 = j0(q*rg1)

    # Polymer size variable
    Usub1 = (q*rg1)**2 * (2.0*nu1 + 1.0) * (2.0*nu1 + 2.0) / 6.0
    Usub2 = (q*rg2)**2 * (2.0*nu2 + 1.0) * (2.0*nu2 + 2.0) / 6.0

    # Form factor amplitude of the polymer:
    with errstate(divide='ignore'):
        Fp1 = o2nu1*power(Usub1, -o2nu1) * gamma(o2nu1) * gammainc(o2nu1, Usub1) 
        Fp2 = o2nu2*power(Usub2, -o2nu2) * gamma(o2nu2) * gammainc(o2nu2, Usub2) 


    # Form factor of the polymer (Pp(q) is not simply Fp(q)^2!!):
    with errstate(divide='ignore'):
        Pp1 = onu1 * power(Usub1, -o2nu1)*gamma(o2nu1)*gammainc(o2nu1, Usub1) - onu1 * power(Usub1, -onu1)*gamma(onu1)*gammainc(onu1,Usub1)
        Pp2 = onu2 * power(Usub2, -o2nu2)*gamma(o2nu2)*gammainc(o2nu2, Usub2) - onu2 * power(Usub2, -onu2)*gamma(onu2)*gammainc(onu2,Usub2)

    # Combine all terms to form intensity:
    #
    # Term 1: Protein
    inten = (sld1 - sld_solvent) * (sld1 - sld_solvent) * v1 * v1 * Pp1

    # Term 2: Polymer
    inten = inten + (sld2 - sld_solvent) * (sld2 - sld_solvent) * v2 * v2 * Pp2

    # Term 3: Cross-term
    inten = inten + 2.0 * v1 * v2 * (sld1 - sld_solvent) * Fp1 * E1 * Fp2

    with errstate(divide='ignore'):
        # Convert SLDs to A^-2, and convert intensity to cm^-1. Normalize by particle volume.
        inten = inten * 1.0e-6 * 1.0e-6 * 1.0e8 / Vtotal

    return inten

Iq.vectorized = True  # Iq accepts an array of q values

# Effective Radius for S(Q). This is only an estimation.
def ER(radius, rc, rg1, rg2, v1, v2):
    """Effective radius of a core-chain-chain sphere."""
    eff_rad = radius + rg1 + rg2
    return eff_rad
