

import numpy as np  # type: ignore
from numpy import pi, inf, power, errstate

name = "ccc"
title = "Spherically symmetric core with grafted polymer chains having two different conformations. Version 2, May 2020."
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
              ["rc",          "Ang",           150,     [0, inf],       "",       "Cutoff distance between region 1 and region 2"],
              ["poly_sig",    "chains/nm^2",   0.33,    [0, inf],       "volume", "Grafting Density"],
              ["rg1",         "Ang",           163,      [1, inf],       "volume", "Radius of gyration of chain in region 1"],
              ["rg2",         "Ang",           100,      [1, inf],       "volume", "Radius of gyration of chain in region 2"],
              ["nu1",         "None",          0.70,    [0.25, 1.0],    "",       "Excluded volume parameter of chain in region 1"],
              ["nu2",         "None",          0.50,    [0.25, 1.0],    "",       "Excluded volume parameter of chain in region 2"],
              ["v1",          "Ang^3",         12000,   [0, inf]   ,    "volume", "Volume of polymer in region 1"],
              ["v2",          "Ang^3",         12000,   [0, inf]   ,    "volume", "Volume of polymer in region 2"],
              ["I0",          "None",          0.0,     [0, inf],       "volume",        "Intensity of free chains"],
              ["rg3",         "Ang",           25.0,    [0, inf],       "",        "Radius of gyration of free chains"],
             ]

radius_effective_modes = ["radius", "outer_radius"]
source = ["lib/sas_3j1x_x.c", "lib/sas_gammainc.c", "lib/sas_gamma.c", "ccc.c"]

