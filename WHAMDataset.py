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

    @staticmethod
    def get_cache_key(infile, wrap=True, loop_edges=True, max_fe=None,
        *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.
        """
        from os.path import expandvars

        return (WHAMDataset, expandvars(infile), wrap, loop_edges,
                max_fe)

    @staticmethod
    def get_cache_message(cache_key):
        return "previously loaded from '{0}'".format(cache_key[1])

    def __init__(self, infile, wrap=True, loop_edges=True, max_fe=None,
        verbose=1, debug=0, **kwargs):
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
        dist = pandas.read_csv(expandvars(infile), delim_whitespace=True,
                 header=0, names=["phi", "psi", "free energy",
                 "probability"], na_values=[9999999.000000])
        if wrap:
            dist["phi"][dist["phi"] > 180] -= 360
            dist["psi"][dist["psi"] > 180] -= 360

        # Organize data
        x_centers = np.unique(dist["phi"])
        y_centers = np.unique(dist["psi"])
        free_energy = np.zeros((x_centers.size, y_centers.size),
                                    np.float) * np.nan
        probability = np.zeros((x_centers.size, y_centers.size),
                                    np.float) * np.nan
        for index, phi, psi, fe, p in dist.itertuples():
            if not np.isnan(fe):
                x_index = np.where(x_centers == phi)[0][0]
                y_index = np.where(y_centers == psi)[0][0]
                free_energy[x_index, y_index] = fe
                probability[x_index, y_index] = p
        free_energy -= np.nanmin(free_energy)

        x_width = np.mean(x_centers[1:] - x_centers[:-1])
        y_width = np.mean(y_centers[1:] - y_centers[:-1])

        # Loop data to allow contour lines to be drawn to edges
        if loop_edges:
            x_centers = np.concatenate(([x_centers[0]  - x_width],
                                         x_centers,
                                        [x_centers[-1] + x_width]))
            y_centers = np.concatenate(([y_centers[0]  - y_width],
                                         y_centers,
                                        [y_centers[-1] + y_width]))
            temp = np.zeros((x_centers.size, y_centers.size)) * np.nan
            temp[1:-1,1:-1]  = free_energy
            temp[1:-1,-1]    = free_energy[:,0]
            temp[-1,1:-1]    = free_energy[0,:]
            temp[1:-1,0]     = free_energy[:,-1]
            temp[0,1:-1]     = free_energy[-1,:]
            temp[0,0]        = free_energy[-1,-1]
            temp[-1,-1]      = free_energy[0,0]
            temp[0,-1]       = free_energy[-1,0]
            temp[-1,0]       = free_energy[0,-1]
            free_energy = temp

        self.x_centers = x_centers
        self.y_centers = y_centers
        self.x_width = x_width
        self.y_width = y_width
        self.x_bins  = np.linspace(x_centers[0]  - x_width / 2,
                                   x_centers[-1] + x_width / 2,
                                   x_centers.size + 1)
        self.y_bins  = np.linspace(y_centers[0]  - y_width / 2,
                                   y_centers[-1] + y_width / 2,
                                   y_centers.size + 1)
        self.dist = free_energy

        # Prepare mask
        if max_fe is not None:
            self.mask = np.ma.masked_where(np.logical_and(
              free_energy <= max_fe,
              np.logical_not(np.isnan(free_energy))),
              np.ones_like(free_energy))
        else:
            self.mask = np.ma.masked_where(
              np.logical_not(np.isnan(free_energy)),
              np.ones_like(free_energy))
