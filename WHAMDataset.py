# -*- coding: utf-8 -*-
#   myplotspec_forcefield.WHAMDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages weighted histogram analysis method datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class WHAMDataset(object):
    """
    Manages weighted histogram analysis method datasets.

    Parses and organizes weighted histogram analysis method Ramachandran
    (backbone Φ/Ψ) distributions from umbrella sampling simulations,
    calculated using:
      Grossfield, Alan, WHAM: the weighted histogram analysis method,
      version 2.0.9, <http://membrane.urmc.rochester.edu/content/wham>_
    """

    def __init__(self, infile, wrap=True, loop_edges=True, max_fe=None,
        verbose=1, debug=1, **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          wrap (bool): Wrap x and y coordinates between 180 and 360 to
            between -180 and 0
          loop_edges (bool):
          max_fe (float): 
          verbose (bool): Enable verbose output
          debug (bool): Enable debug output
        """
        from os.path import expandvars
        import pandas
        import numpy as np

        # Load data
        if verbose > 0:
            print("loading from '{0}'".format(infile))
        self.data = pandas.read_csv(expandvars(infile), delim_whitespace=True,
                      header=0, names=["x", "y", "free energy", "probability"],
                      na_values=[9999999.000000])
        if wrap:
            self.data["x"][self.data["x"] > 180] -= 360
            self.data["y"][self.data["y"] > 180] -= 360

        # Organize data
        self.x_centers = np.unique(self.data["x"])
        self.y_centers = np.unique(self.data["y"])
        self.free_energy = np.zeros((self.x_centers.size, self.y_centers.size),
                                    np.float) * np.nan
        self.probability = np.zeros((self.x_centers.size, self.y_centers.size),
                                    np.float) * np.nan
        for index, x, y, free_energy, probability in self.data.itertuples():
            if not np.isnan(free_energy):
                x_index = np.where(self.x_centers == x)[0][0]
                y_index = np.where(self.y_centers == y)[0][0]
                self.free_energy[x_index, y_index] = free_energy
                self.probability[x_index, y_index] = probability
        self.free_energy -= np.nanmin(self.free_energy)

        self.x_width = np.mean(self.x_centers[1:] - self.x_centers[:-1])
        self.y_width = np.mean(self.y_centers[1:] - self.y_centers[:-1])

        # Loop data to allow contour lines to be drawn to edges
        if loop_edges:
            self.x_centers = np.concatenate(
                               ([self.x_centers[0]-self.x_width],
                                 self.x_centers,
                                [self.x_centers[-1]+self.x_width]))
            self.y_centers = np.concatenate(
                               ([self.y_centers[0]-self.y_width],
                                 self.y_centers,
                                [self.y_centers[-1]+self.y_width]))
            temp = np.zeros((self.x_centers.size, self.y_centers.size))*np.nan
            temp[1:-1,1:-1]  = self.free_energy
            temp[1:-1,-1]    = self.free_energy[:,0]
            temp[-1,1:-1]    = self.free_energy[0,:]
            temp[1:-1,0]     = self.free_energy[:,-1]
            temp[0,1:-1]     = self.free_energy[-1,:]
            temp[0,0]        = self.free_energy[-1,-1]
            temp[-1,-1]      = self.free_energy[0,0]
            temp[0,-1]       = self.free_energy[-1,0]
            temp[-1,0]       = self.free_energy[0,-1]
            self.free_energy = temp

        self.x_bins  = np.linspace(self.x_centers[0]  - self.x_width / 2,
                                   self.x_centers[-1] + self.x_width / 2,
                                   self.x_centers.size + 1)
        self.y_bins  = np.linspace(self.y_centers[0]  - self.y_width / 2,
                                   self.y_centers[-1] + self.y_width / 2,
                                   self.y_centers.size + 1)

        # Prepare mask
        if max_fe is not None:
            self.mask = np.ma.masked_where(np.logical_and(
              self.free_energy <= max_fe,
              np.logical_not(np.isnan(self.free_energy))),
              np.ones_like(self.free_energy))
        else:
            self.mask = np.ma.masked_where(
              np.logical_not(np.isnan(self.free_energy)),
              np.ones_like(self.free_energy))
