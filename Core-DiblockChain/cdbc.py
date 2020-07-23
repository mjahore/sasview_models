
import numpy as np  # type: ignore
from numpy import pi, inf, power, errstate

name = "cdbc"
title = "Spherically symmetric core with grafted diblock polymer chains having two different conformations."
description = """       """
category = "shape:sphere"

structure_factor = False
form_factor = True

#             [ "name",       "units",         default, [lower, upper], "type",   "description"],
parameters = [["volf",        "None",          0.02,    [0,1],          "",       "Particle volume fraction"],
              ["sld_c",       "1e-6/Ang^2",    3.47,    [-inf, inf],    "sld",    "Core scattering length density"],
              ["sld_s",       "1e-6/Ang^2",    -0.022,  [-inf, inf],    "sld",    "Initiator scattering length densty"],
              ["sld1",        "1e-6/Ang^2",    0.814,   [-inf, inf],    "sld",    "Chain region 1 scattering length density"],
              ["sld2",        "1e-6/Ang^2",    4.24,    [-inf, inf],    "sld",    "Chain region 2 scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",    6.37,    [-inf, inf],    "sld",    "Solvent scattering length density"],
              ["radius",      "Ang",           75,      [0, inf],       "volume", "Core radius"],
              ["i_shell",     "Ang",           10,      [0, inf],       "volume", "Initiator shell thickness"],
              ["poly_sig",    "chains/nm^2",   0.33,    [0, inf],       "volume", "Grafting Density"],
              ["rc",          "Ang",           100,     [0, inf],       "",       "Transition Point"],
              ["C_infty",     "Ang",           12,      [7, 50],        "",       "Characteristic Ratio"],
	      ["M0",          "None",          113,     [0, inf],       "",       "Monomer molar mass (g/mol)"],
              ["M1",          "None",          8900,    [1, inf],       "",       "Mn, Block 1"],
              ["M2",          "None",          9900,    [1, inf],       "",       "Mn, Block 2"],
              ["nu1",         "None",          0.80,    [0.25, 1.0],    "",       "Flory Exp., Block 1"],
              ["nu2",         "None",          0.80,    [0.25, 1.0],    "",       "Flory Exp., Block 2"],
              ["v",           "Ang^3",         149,     [0, inf]   ,    "volume", "Kuhn Monomer Volume"],
              ["I0",          "None",          0.0,     [0, inf],       "volume", "Intensity of free chains"],
              ["rg3",         "Ang",           25.0,    [0, inf],       "",       "Radius of gyration of free chains"],
              ["nu3",         "None",          0.5,     [0.25, 1.0],    "",       "Flory parm for free chain."],
             ]

radius_effective_modes = ["radius", "outer_radius"]
source = ["lib/sas_3j1x_x.c", "lib/sas_gammainc.c", "lib/sas_gamma.c", "cdbc.c"]

