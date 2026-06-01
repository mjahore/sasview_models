r"""
Definition
----------

This model describes a star polymer with f identical triblock arms. Each arm is described by three blocks
with different degrees of polymerization N_i, Flory exponents \nu_{i}, and SLD contrasts.

References
----------

Y. Wei and M. J. A. Hore, "Characterizing Polymer Structure with Small-Angle Neutron Scattering: A Tutorial", J. Appl. Phys. 129, 171101 (2021).
B. Hammouda, "Form Factors for Branched Polymers with Excluded Volume", J. of Research of NIST, 121, 139-164 (2016).

Authorship and Verification
----------------------------

* **Author:** Michael J. A. Hore **Date:** 1 Jun 2026
"""

import numpy as np
from numpy import inf, errstate, power,exp, sqrt
from sasmodels.special import sas_gammainc, sas_gamma

name = "triblock_star"
title = "Triblock Star Polymer"
description = """\
      List of default parameters:
      scale = Scaling factor
      b     = Kuhn length
      sld1  = SLD of inner block 1
      sld2  = SLD of middle block 2
      sld3  = SLD of outer block 3
      slds  = SLD of solvent
      N1    = deg. polym. of block 1
      N2    = deg. polym. of block 2
      N3    = deg. polym. of block 3
      nu1   = Flory exp. of block 1
      nu2   = Flory exp. of block 2
      nu3   = Flory exp. of block 3
      background = Incoherent background"""
category = "shape-independent"

# pylint: disable=bad-whitespace, line-too-long
#             ["name", "units", default,       [lower, upper], "type", "description"],
parameters = [
              ["f",        "",                 4,         [1, inf],       "",     "Number of arms"],
              ["b",        "Ang",              7.0,       [1, inf],       "",     "Kuhn length"],
              ["sld1",     "1e-6/Ang^2",       1.0,       [-inf,inf],     "sld",  "Block 1 SLD"],
              ["sld2",     "1e-6/Ang^2",       1.5,       [-inf,inf],     "sld",  "Block 2 SLD"],
              ["sld3",     "1e-6/Ang^2",       1.0,       [-inf,inf],     "sld",  "Block 3 SLD"],
              ["slds",     "1e-6/Ang^2",       6.3,       [-inf,inf],     "sld",  "Solvent SLD"],
              ["N1",       "",                 50,        [0,inf],        "",     "Deg. Polym. 1"],
              ["N2",       "",                 50,        [0,inf],        "",     "Deg. Polym. 2"],
              ["N3",       "",                 50,        [0,inf],        "",     "Deg. Polym. 3"],
              ["nu1",      "",                 0.5,       [0.25,0.999],   "",     "Flory Exp. 1"],
              ["nu2",      "",                 0.5,       [0.25,0.999],   "",     "Flory Exp. 2"],
              ["nu3",      "",                 0.5,       [0.25,0.999],   "",     "Flory Exp. 3"],
             ]
# pylint: enable=bad-whitespace, line-too-long

def Iq(q,
       f,
       b,
       sld1,
       sld2,
       sld3,
       slds,
       N1,
       N2,
       N3,
       nu1,
       nu2,
       nu3):

    # Compute the contrast terms:
    delta1 = sld1 - slds
    delta2 = sld2 - slds
    delta3 = sld3 - slds

    # Excl. Vol. Parameters for each block:
    onu1  = 1.0/nu1
    o2nu1 = 1.0/2.0/nu1
    onu2  = 1.0/nu2
    o2nu2 = 1.0/2.0/nu2
    onu3  = 1.0/nu3
    o2nu3 = 1.0/2.0/nu3

    # Polymer chain form factor and form factor amplitude for each block:
    ## Block 1:
    U1 = q**2 * b**2 * power(N1, 2.0*nu1) / 6.0
    P1 = onu1 * power(U1, -o2nu1) * sas_gamma(o2nu1) * sas_gammainc(o2nu1, U1) - onu1 * power(U1, -onu1) * sas_gamma(onu1) * sas_gammainc(onu1, U1)
    F1 = 0.50 * onu1 * power(U1, -o2nu1) * sas_gamma(o2nu1) * sas_gammainc(o2nu1, U1)

    ## Block 1, double chain length.
    U12 = q**2 * b**2 * power(2.0*N1, 2.0*nu1) / 6.0
    P12 = onu1 * power(U12, -o2nu1) * sas_gamma(o2nu1) * sas_gammainc(o2nu1, U12) - onu1 * power(U12, -onu1) * sas_gamma(onu1) * sas_gammainc(onu1, U12)
    F12 = 0.50 * onu1 * power(U12, -o2nu1) * sas_gamma(o2nu1) * sas_gammainc(o2nu1, U12)

    ## Block 2:
    U2 = q**2 * b**2 * power(N2, 2.0*nu2) / 6.0
    P2 = onu2 * power(U2, -o2nu2) * sas_gamma(o2nu2) * sas_gammainc(o2nu2, U2) - onu2 * power(U2, -onu2) * sas_gamma(onu2) * sas_gammainc(onu2, U2)
    F2 = 0.50 * onu2 * power(U2, -o2nu2) * sas_gamma(o2nu2) * sas_gammainc(o2nu2, U2)

    ## Block 2, double chain length.
    U22 = q**2 * b**2 * power(2.0*N2, 2.0*nu2) / 6.0
    P22 = onu2 * power(U22, -o2nu2) * sas_gamma(o2nu2) * sas_gammainc(o2nu2, U22) - onu2 * power(U22, -onu2) * sas_gamma(onu2) * sas_gammainc(onu2, U22)
    F22 = 0.50 * onu2 * power(U22, -o2nu2) * sas_gamma(o2nu2) * sas_gammainc(o2nu2, U22)

    ## Block 3:
    U3 = q**2 * b**2 * power(N3, 2.0*nu3) / 6.0
    P3 = onu3 * power(U3, -o2nu3) * sas_gamma(o2nu3) * sas_gammainc(o2nu3, U3) - onu3 * power(U3, -onu3) * sas_gamma(onu3) * sas_gammainc(onu3, U3)
    F3 = 0.50 * onu3 * power(U3, -o2nu3) * sas_gamma(o2nu3) * sas_gammainc(o2nu3, U3)

    ## Block 3, double chain length.
    U32 = q**2 * b**2 * power(2.0*N3, 2.0*nu3) / 6.0
    P32 = onu3 * power(U32, -o2nu3) * sas_gamma(o2nu3) * sas_gammainc(o2nu3, U32) - onu3 * power(U32, -onu3) * sas_gamma(onu3) * sas_gammainc(onu3, U32)
    F3 = 0.50 * onu3 * power(U32, -o2nu3) * sas_gamma(o2nu3) * sas_gammainc(o2nu3, U32)

    # Propagator term (approximate):
    E13  = exp(-U2)
    E132 = exp(-U22)

    # Combine to form interbranch and single branch terms:
    Psb = (N1**2) * (delta1**2) * P1 + (N2**2) * (delta2**2) * P2 + (N3**2) * (delta3**2) * P3
    Psb = Psb + 2.0*N1*N2*delta1*delta2*F1*F2 + 2.0*N2*N3*delta2*delta3*F2*F3
    Psb = Psb + 2.0*N1*N3*delta1*delta3*F1*E13*F3
   
    P2sb = 4.0*(N1**2) * (delta1**2) * P12 + 4.0*(N2**2) * (delta2**2) * P22 + 4.0*(N3**2) * (delta3**2) * P33
    P2sb = P2sb + 8.0*N1*N2*delta1*delta2*F12*F22 + 8.0*N2*N3*delta2*delta3*F22*F32
    P2sb = P2sb + 8.0*N1*N3*delta1*delta3*F12*E132*F32

    # Finally! The interbranch interference term!
    Pib  = 2.0 * P2sb - Psb

    Pq = power(f, -2.0) * (f*Psb - f*(f-1)*Pib)

    # Form the structure factor according to the RPA:
    #Sq = power(n*phi_p*vm*Pp, -1.00) + power(vs * (1.0 - phi_p), -1.00) - 2.00*chi/sqrt(vm*vs)
    #inten = 1e-4 * power(Sq, -1.0)

    inten = 1e-4 * Pq
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
