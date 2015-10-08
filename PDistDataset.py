# -*- coding: utf-8 -*-
#   ramaplot.PDistDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages probability distribution datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class PDistDataset(object):
    """
    Manages probability distribution datasets.

    Generates probability distribution from series of Φ/Ψ values,
    representing either the probability of Φ/Ψ, or the expectation value
    of a selected measurement (e.g. energy) at that Φ/Ψ.
    """

    @staticmethod
    def get_cache_key(infile, loop_edges=True, mode="hist", bins=72,
        phikey="phi", psikey="psi", max_fe=None, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documentented under :func:`__init__`.
        """
        from os.path import expandvars

        if mode == "hist":
            return (PDistDataset, expandvars(infile), loop_edges, mode, bins,
                    phikey, psikey, max_fe)
        elif mode == "kde":
            kde_kw = kwargs.get("kde_kw", {})
            bandwidth = kde_kw.get("bandwidth", kwargs.get("bandwidth", 5))
            return (PDistDataset, expandvars(infile), loop_edges, mode, bins,
                    phikey, psikey, bandwidth, max_fe)

    @staticmethod
    def get_cache_message(cache_key):
        """
        Generates message to be used when reloading previously-loaded
        dataset.

        Arguments:
            cache_key (tuple): key with which dataset object is stored
              in dataset cache

        Returns:
            cache_message (str): message to be used when reloading
              previously-loaded dataset
        """
        return "previously loaded from '{0}'".format(cache_key[1])

    def __init__(self, infile, loop_edges=True, mode="hist", bins=72,
        phikey="phi", psikey="psi", max_fe=None, verbose=1, debug=0, **kwargs):
        """
        Arguments:
          infile (str): Path to text input file, may contain environment
          loop_edges (bool):
          mode:
          bins:
          phikey:
          psikey:
          max_fe:
          hist_kw:
          bins:
          verbose (int): Level of verbose output
          debug (int): Level of debug output

        .. todo:
            - improve bins support
            - auto-detect phikey and psikey
            - support variable bandwidth
            - PERIODICITY
            - VARIABLE WIDTH
        """
        from copy import copy
        from os.path import expandvars
        import pandas
        import numpy as np

        # Check arguments
        if mode not in ["hist", "kde"]:
            raise ValueError("Argument 'mode' does not support provided " +
              "value '{0}', must be 'hist' or 'kde'".format(mode))

        # Load data
        if verbose > 0:
            print("loading from '{0}'".format(infile))
        dist = pandas.read_csv(expandvars(infile), delim_whitespace=True,
                 index_col=0)


        if mode == "hist":
            hist_kw = copy(kwargs.get("hist_kw", {}))
            if "bins" in hist_kw:
                bins = hist_kw.pop("bins")
            if isinstance(bins, int):
                bins = np.linspace(-180, 180, bins)

            count, x_bins, y_bins = np.histogram2d(
              dist[phikey], dist[psikey], bins=bins, **hist_kw)
            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            y_centers = (y_bins[:-1] + y_bins[1:]) / 2
            x_width = np.mean(x_centers[1:] - x_centers[:-1])
            y_width = np.mean(y_centers[1:] - y_centers[:-1])

            free_energy = -1 * np.log(count / count.sum())
            free_energy[np.isinf(free_energy)] = np.nan
            free_energy -= np.nanmin(free_energy)

        elif mode == "kde":
            from sklearn.neighbors import KernelDensity

            kde_kw = copy(kwargs.get("kde_kw", {}))
            kde_kw["bandwidth"] = kde_kw.get("bandwidth",
              kwargs.get("bandwidth", 5))
            x_bins = np.linspace(-180, 180, bins+1)
            y_bins = np.linspace(-180, 180, bins+1)
            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            y_centers = (y_bins[:-1] + y_bins[1:]) / 2
            x_width = np.mean(x_centers[1:] - x_centers[:-1])
            y_width = np.mean(y_centers[1:] - y_centers[:-1])
            xg, yg = np.meshgrid(x_centers, y_centers)
            xyg = np.vstack([yg.ravel(), xg.ravel()]).T
            samples = np.column_stack((dist[phikey], dist[psikey]))

            kde = KernelDensity(**kde_kw)
            kde.fit(samples)
            pdist_series = np.exp(kde.score_samples(xyg))
            pdist = np.zeros((x_centers.size, y_centers.size),
                      np.float) * np.nan
            for phi, psi, p in np.column_stack((xyg, pdist_series)):
                x_index = np.where(x_centers == phi)[0][0]
                y_index = np.where(y_centers == psi)[0][0]
                pdist[x_index, y_index] = p

            free_energy = -1 * np.log(pdist / pdist.sum())
            free_energy[np.isinf(free_energy)] = np.nan
            free_energy -= np.nanmin(free_energy)
        elif mode == "kde_3d":
            from sklearn.neighbors import KernelDensity

            kde_kw = copy(kwargs.get("kde_kw", {}))
            kde_kw["bandwidth"] = kde_kw.get("bandwidth",
              kwargs.get("bandwidth", 5))
            z_key = kwargs.get("z_key", "energy")
            # Only a single bandwidth is supported; scale z
            #   dimension to span range of 120-240
            scale_range = 340
            z = copy(dist[z_key])
            z -= z.min()        # shift bottom to 0
            z_range = z.max()     # save max

            z *= (scale_range / z_range)
            z += (360 - scale_range) / 2    # Give buffer on top and bottom

            if kwargs.get("left_half", False):
                x_bins = np.linspace(-180, 0, (bins/2)+1)
            else:
                x_bins = np.linspace(-180, 180, bins+1)
            y_bins = np.linspace(-180, 180, bins+1)
            z_bins = np.linspace(   0, 360, bins+1)
            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            y_centers = (y_bins[:-1] + y_bins[1:]) / 2
            z_centers = (z_bins[:-1] + z_bins[1:]) / 2
            x_width = np.mean(x_centers[1:] - x_centers[:-1])
            y_width = np.mean(y_centers[1:] - y_centers[:-1])
            z_width = np.mean(z_centers[1:] - z_centers[:-1])
            xg, yg, zg = np.meshgrid(x_centers, y_centers, z_centers)
            xyzg = np.vstack([xg.ravel(), yg.ravel(), zg.ravel()]).T
            samples = np.column_stack((dist[phikey], dist[psikey], z))

            kde = KernelDensity(**kde_kw)
            kde.fit(samples)
            pdist_series = np.exp(kde.score_samples(xyzg))
            pdist = np.zeros((x_centers.size, y_centers.size, z_centers.size),
                      np.float) * np.nan
            for phi, psi, z, p in np.column_stack((xyzg, pdist_series)):
                x_index = np.where(x_centers == phi)[0][0]
                y_index = np.where(y_centers == psi)[0][0]
                z_index = np.where(z_centers == z)[0][0]
                pdist[x_index, y_index, z_index] = p
            pdist /= np.sum(pdist)                    # Normalize whole thing
            a = np.sum(pdist, axis=2)[:,:,np.newaxis] # Total P in each x/y bin
            b = pdist / a                             # Normalize each x/y bin
            c = np.zeros_like(b) * np.nan
            for i in range(b.shape[0]):
                for j in range(b.shape[1]):
                    c[i,j] = b[i,j] * z_centers # Scale each z bin by weight
            free_energy = np.sum(c, axis=2)     # Weighted average
            free_energy -= (360 - scale_range) / 2  # Shift back down
            free_energy *= (z_range / scale_range) # Back from degrees to E
            free_energy *= 627.503              # Convert to kcal/mol
            free_energy[np.isinf(free_energy)] = np.nan
            free_energy -= np.nanmin(free_energy)
#            print(free_energy, free_energy.shape, free_energy.sum(),
#              free_energy.min(), free_energy.max())
    
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
        self.x = dist[phikey]
        self.y = dist[psikey]

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
