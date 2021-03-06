

import numpy as np  # type: ignore
from numpy import pi, inf, power, errstate

name = "e_ccc"
title = "Empirical model of polymer-grafted nanosphere."
description = """       """
category = "shape:sphere"

structure_factor = False
form_factor = True

#             [ "name",       "units",         default, [lower, upper], "type",   "description"],
parameters = [["m",           "None",          4.0,     [1,4],          "",       "Core Porod exponent"],
              ["sld_c",       "1e-6/Ang^2",    3.47,    [-inf, inf],    "sld",    "Core scattering length density"],
              ["sld1",        "1e-6/Ang^2",    0.814,   [-inf, inf],    "sld",    "Chain region 1 scattering length density"],
              ["sld2",        "1e-6/Ang^2",    4.24,    [-inf, inf],    "sld",    "Chain region 2 scattering length density"],
              ["sld_solvent", "1e-6/Ang^2",    6.37,    [-inf, inf],    "sld",    "Solvent scattering length density"],
              ["R",           "Ang",           75,      [0, inf],       "volume", "Core radius"],
              ["rc",          "Ang",           150,     [0, inf],       "",       "Cutoff distance between region 1 and region 2"],
              ["poly_sig",    "chains/nm^2",   0.33,    [0, inf],       "volume", "Grafting Density"],
              ["rg1",         "Ang",           163,     [1, inf],       "volume", "Radius of gyration of chain in region 1"],
              ["rg2",         "Ang",           100,     [1, inf],       "volume", "Radius of gyration of chain in region 2"],
              ["nu1",         "None",          0.70,    [0.25, 1.0],    "",       "Excluded volume parameter of chain in region 1"],
              ["nu2",         "None",          0.50,    [0.25, 1.0],    "",       "Excluded volume parameter of chain in region 2"],
              ["v1",          "Ang^3",          12000,   [0, inf],       "volume", "Volume of 1 chain in region 1"],
              ["v2",          "Ang^3",          12000,   [0, inf],       "volume", "Volume of 1 chain in region 2"],
             ]

radius_effective_modes = ["radius", "outer_radius"]
source = ["lib/sas_3j1x_x.c", "lib/sas_gammainc.c", "lib/sas_gamma.c", "e_ccc.c"]

