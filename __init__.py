# -*- coding: utf-8 -*-
#   myplotspec_forcefield.__init__.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
General functions.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
def cmap_ff99SB():
    """
    Generates purple->yellow->black colormap.

    Generates colormap in style of:
      Hornak, Viktor, Abel, Robert, Okur, Asim, Stockbine, Bentley,
      Roitberg, Adrian, Simmerling, Carlos, Comparison of Multiple Amber
      Force Fields and Development of Improved Protein Backbone
      Parameters. Proteins: Structure, Function, and Bioinformatics.
      2006. 65. 712-725.

    Returns:
      cmap (LinearSegmentedColormap): ff99SB-style colormap
    """
    from matplotlib.colors import LinearSegmentedColormap

    cdict = {"red":   [(0, 1, 1)],
             "green": [(0, 1, 1)],
             "blue":  [(0, 0, 0)]}
    for i in range(1, 193, 1):
        red   = 1.0
        green = 1.0 - ((i - 1) / 193)
        blue  = 0.0
        cdict["red"]   += [(i / 384, red,   red)]
        cdict["green"] += [(i / 384, green, green)]
        cdict["blue"]  += [(i / 384, blue,  blue)]
    for i in range(193, 385, 1):
        red   = (1.0 - ((i - 193) / 192)) * ((384 - i) / 192) ** 0.3
        green =  0.0
        blue  = (0.0 + ((i - 193) / 192)) * ((384 - i) / 192) ** 0.2
        cdict["red"]   += [(i / 384, red,   red)]
        cdict["green"] += [(i / 384, green, green)]
        cdict["blue"]  += [(i / 384, blue,  blue)]
    cdict["red"]   = tuple(cdict["red"])
    cdict["green"] = tuple(cdict["green"])
    cdict["blue"]  = tuple(cdict["blue"])

    return LinearSegmentedColormap("test", cdict)
