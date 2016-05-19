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
if __name__ == "__main__":
    __package__ = str("ramaplot")
    import ramaplot
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

    Input data should be provided in the following format::

      #X          Y           Free     Pro
      -178.750000 -178.750000 6.442071 0.000000
      -178.750000 -176.250000 6.559820 0.000000
      ...         ...         ...      ...

    """

    @classmethod
    def get_cache_key(cls, infile, phikey="phi", psikey="psi",
        zkey="free energy", wrap=True, mask_cutoff=None,
        calc_populations=False, plot_populations=False, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documented under :func:`__init__`.
        """
        from os.path import expandvars

        return (cls, expandvars(infile), phikey, psikey, zkey, wrap,
          mask_cutoff, calc_populations, plot_populations)

    def __init__(self, phikey="phi", psikey="psi", zkey="free energy",
        wrap=True, mask_cutoff=None,
        calc_populations=False, plot_populations=False,
        verbose=1, debug=0, **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          phikey (str): Key from which to load Φ
          psikey (str): Key from which to load Ψ
          zkey (str): Key from which to load distribution; acceptable
            values are  'free energy' and 'probability'
          wrap (bool): Wrap x and y coordinates between 180 and 360 to
            between -180 and 0
          mask_cutoff (float): Cutoff beyond which distribution is
            masked, if `zkey` is 'free energy', this is a the maximum
            free energy above which the mask will be set, and if `zkey`
            is 'probability', this is the minimum probability below
            which the mask will be set
          verbose (int): Level of verbose output
          debug (int): Level of debug output
          kwargs (dict): Additional keyword arguments
        """
        import numpy as np
        import pandas as pd

        # Manage arguments
        if zkey not in ["free energy", "probability"]:
            raise ValueError("Argument 'zkey' does not support provided " +
              "value '{0}', must be 'free energy or ".format(zkey) +
              "probability.")

        # Load data
        read_csv_kw = dict(delim_whitespace=True, header=0,
          names=("phi", "psi", "free energy", "probability"),
          na_values=(9999999.000000))
        read_csv_kw.update(kwargs.pop("read_csv_kw", {}))
        dataframe = self.load_dataset(verbose=verbose, read_csv_kw=read_csv_kw,
          **kwargs).dataframe
        if wrap:
            dataframe[phikey][dataframe[phikey] > 180] -= 360
            dataframe[psikey][dataframe[psikey] > 180] -= 360

        # Organize data
        x_centers = np.unique(dataframe[phikey])
        y_centers = np.unique(dataframe[psikey])
        free_energy = np.zeros((x_centers.size, y_centers.size),
                                    np.float) * np.nan
        probability = np.zeros((x_centers.size, y_centers.size),
                                    np.float) * np.nan
        for index, phi, psi, fe, p in dataframe.itertuples():
            if not np.isnan(fe):
                x_index = np.where(x_centers == phi)[0][0]
                y_index = np.where(y_centers == psi)[0][0]
                free_energy[x_index, y_index] = fe
                probability[x_index, y_index] = p
        free_energy -= np.nanmin(free_energy)
        probability /= np.nansum(probability)

        x_width = np.mean(x_centers[1:] - x_centers[:-1])
        y_width = np.mean(y_centers[1:] - y_centers[:-1])

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

        # Calculate state populations
        if calc_populations:
            states = kwargs.get("states",  [
              ("β",       -151,  151),
              ("PPII",     -66,  140),
              ("ξ",       -145,   55),
              ("γ'",       -81,   65),
              ("α",        -70,  -25),
              ("$L_α$",     55,   45),
              ("γ",         73,  -35),
              ("PPII'",     56, -124),
              ("plateau", -100, -130)])
            state_radius = kwargs.get("state_radius", 45)

            distances = np.zeros((len(states), len(x_centers), len(y_centers)))
            xs = []
            ys = []
            for i, (state, x, y) in enumerate(states):
                xs += [x]
                ys += [y]
                # There must be a better way to do this, but this works
                for j, xc in enumerate(x_centers):
                    for k, yc in enumerate(y_centers):
                        dx = (xc - x)
                        if dx <= -180 or dx >= 180:
                            dx = 360 - dx
                        else:
                            dx = dx
                        dy = (yc - y)
                        if dy <= -180 or dy >= 180:
                            dy = 360 - dy
                        else:
                            dy = dy
                        distances[i,j,k] = np.sqrt(dx**2 + dy**2)
            assignments = np.argmin(distances, axis=0)
            assignments[np.min(distances, axis=0) >= state_radius] = \
              len(states) + 1

            index, state_populations = [], []
            for i, (state, x, y) in enumerate(states):
                index += [state]
                state_populations += [(x, y,
                  np.nansum(probability[assignments==i]))]
            state_populations = pd.DataFrame(state_populations, index=index,
              columns=["Φ center", "Ψ center", "population"])
            self.state_populations = state_populations

            if verbose >= 1:
                print(state_populations)

            if plot_populations:
                self.dist = assignments
                self.mask = np.ma.masked_where(
                  np.logical_not(assignments == len(states) + 1),
                  np.ones_like(assignments))
                self.x = np.array(xs)
                self.y = np.array(ys)
                label, label_kw = [], []
                from .myplotspec import multi_get_copy
                default_label_kw = multi_get_copy(["default_label_kw",
                  "label_kw"], kwargs, {})
                for index, row in state_populations.iterrows():
                    label += ["{0}\n{1:2d}%".format(index,
                     int(row["population"]*100))]
                    label_kw += [default_label_kw.copy()]
                    label_kw[-1]["x"] = row["Φ center"]
                    label_kw[-1]["y"] = row["Ψ center"]
                self.label = label
                self.label_kw = label_kw
