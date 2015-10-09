# -*- coding: utf-8 -*-
#   ramaplot.WHAMDataset.py
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
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class WHAMDataset(Dataset):
    """
    Manages weighted histogram analysis method datasets.

    Parses and organizes weighted histogram analysis method Ramachandran
    (backbone Φ/Ψ) distributions from umbrella sampling simulations,
    calculated using:
      Grossfield, Alan, WHAM: the weighted histogram analysis method,
      version 2.0.9, <http://membrane.urmc.rochester.edu/content/wham>_
    """

    @classmethod
    def get_cache_key(cls, infile, zkey="free energy", wrap=True,
        loop_edges=True, mask_cutoff=None, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.
        """
        from os.path import expandvars

        return (cls, expandvars(infile), zkey, wrap, loop_edges, mask_cutoff)

    def __init__(self, zkey="free energy", wrap=True, loop_edges=True,
        mask_cutoff=None, verbose=1, debug=0, **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          zkey (str): Key from which to load distribution; acceptable
            values are  'free energy' and 'probability'
          wrap (bool): Wrap x and y coordinates between 180 and 360 to
            between -180 and 0
          loop_edges (bool): Mirror first and last row and column along
            edges of distribution; enables contours to be plotted
            smoothly to edge
          mask_cutoff (float): Cutoff beyond which distribution is
            masked, if `zkey` is 'free energy', this is a the maximum
            free energy above which the mask will be set, and if `zkey`
            is 'probability', this is the minimum probability below
            which the mask will be set
          verbose (int): Level of verbose output
          debug (int): Level of debug output
          kwargs (dict): Additional keyword arguments
        """
        import pandas
        import numpy as np

        # Manage arguments
        if zkey not in ["free energy", "probability"]:
            raise ValueError("Argument 'zkey' does not support provided " +
              "value '{0}', must be 'free energy or ".format(zkey) +
              "probability.")
        read_csv_kw = dict(delim_whitespace=True, header=0,
          names=("phi", "psi", "free energy", "probability"),
          na_values=(9999999.000000))
        read_csv_kw.update(kwargs.pop("read_csv_kw", {}))

        # Load data
        dist = self.load_dataset(verbose=verbose, read_csv_kw=read_csv_kw,
          **kwargs).data
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
        probability /= np.nansum(probability)

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
            temp[1:-1,1:-1] = free_energy
            temp[1:-1,-1]   = free_energy[:,0]
            temp[-1,1:-1]   = free_energy[0,:]
            temp[1:-1,0]    = free_energy[:,-1]
            temp[0,1:-1]    = free_energy[-1,:]
            temp[0,0]       = free_energy[-1,-1]
            temp[-1,-1]     = free_energy[0,0]
            temp[0,-1]      = free_energy[-1,0]
            temp[-1,0]      = free_energy[0,-1]
            free_energy = temp
            temp = np.zeros((x_centers.size, y_centers.size)) * np.nan
            temp[1:-1,1:-1] = probability
            temp[1:-1,-1]   = probability[:,0]
            temp[-1,1:-1]   = probability[0,:]
            temp[1:-1,0]    = probability[:,-1]
            temp[0,1:-1]    = probability[-1,:]
            temp[0,0]       = probability[-1,-1]
            temp[-1,-1]     = probability[0,0]
            temp[0,-1]      = probability[-1,0]
            temp[-1,0]      = probability[0,-1]
            probability = temp

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

        # Store distribution in instance variable and create mask
        if zkey == "free energy":
            self.dist = free_energy
            if mask_cutoff is not None:
                self.mask = np.ma.masked_where(np.logical_and(
                  free_energy <= mask_cutoff,
                  np.logical_not(np.isnan(free_energy))),
                  np.ones_like(free_energy))
            else:
                self.mask = np.ma.masked_where(
                  np.logical_not(np.isnan(free_energy)),
                  np.ones_like(free_energy))
        elif zkey == "probability":
            self.dist = probability
            if mask_cutoff is not None:
                self.mask = np.ma.masked_where(np.logical_and(
                  probability >= mask_cutoff,
                  np.logical_not(np.isnan(probability))),
                  np.ones_like(probability))
            else:
                self.mask = np.ma.masked_where(
                  np.logical_not(np.isnan(probability)),
                  np.ones_like(probability))
