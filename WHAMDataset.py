# -*- coding: utf-8 -*-
#   myplotspec_forcefield.WHAMDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages WHAM datasets
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class WHAMDataset(object):
    """
    Manages WHAM datasets
    """

    def __init__(self, infile, **kwargs):
        """
        Initializes dataset.

        Arguments:
            infile (str): Path to text infile
        """
        from os.path import expandvars
        import pandas
        import numpy as np

        self.data = pandas.read_csv(expandvars(infile), delim_whitespace=True,
                      header=0, names=["x", "y", "free energy", "probability"],
                      na_values=[9999999.000000])

        self.x_centers = np.unique(self.data["x"])
        self.y_centers = np.unique(self.data["y"])
        self.x_width = np.mean(self.x_centers[1:] - self.x_centers[:-1])
        self.y_width = np.mean(self.y_centers[1:] - self.y_centers[:-1])
        self.x_bins  = np.linspace(self.x_centers[0]  - self.x_width/2,
                                   self.x_centers[-1] + self.x_width/2,
                                   self.x_centers.size + 1)
        self.y_bins  = np.linspace(self.y_centers[0]  - self.y_width/2,
                                   self.y_centers[-1] + self.y_width/2,
                                   self.y_centers.size + 1)

        self.free_energy = np.zeros((self.x_centers.size, self.y_centers.size),
                                    np.float)
        self.probability = np.zeros((self.x_centers.size, self.y_centers.size),
                                    np.float)
        self.free_energy[:] = np.nan
        self.probability[:] = np.nan
        for index, x, y, free_energy, probability in self.data.itertuples():
            if not np.isnan(free_energy):
                y_index = np.where(self.x_centers == x)[0][0]
                x_index = np.where(self.y_centers == y)[0][0]
                self.free_energy[x_index, y_index] = free_energy
                self.probability[x_index, y_index] = probability
