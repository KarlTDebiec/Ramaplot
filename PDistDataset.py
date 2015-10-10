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
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class PDistDataset(Dataset):
    """
    Manages probability distribution datasets.

    Generates probability distribution from series of Φ/Ψ values,
    representing either the probability of Φ/Ψ, or the expectation value
    of a selected measurement (e.g. energy) at that Φ/Ψ.

    Input data should be providied in a whitespace-delimited text file
    including columns for Φ, Ψ, and any additional data, such as this
    output from `cpptraj`'s `multidihedral` command::

      #Frame       phi:2     psi:2    chip:2 ...
             1  -62.1431  144.6768   72.2964 ...
             2  -63.2487  151.6551   71.9101 ...
           ...       ...       ...       ...

    """

    @classmethod
    def get_cache_key(cls, infile, phikey="phi", psikey="psi",
        zkey="free energy", mode="hist", bins=72, bandwidth=5, wrap=True,
        loop_edges=True, mask_cutoff=None, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documentented under :func:`__init__`.
        """
        from os.path import expandvars

        if mode == "hist":
            return (cls, expandvars(infile), phikey, psikey, zkey, mode, bins,
                    wrap, loop_edges, mask_cutoff)
        elif mode == "kde":
            return (cls, expandvars(infile), phikey, psikey, zkey, mode, bins,
                    bandwidth, wrap, loop_edges, mask_cutoff)

    def __init__(self, phikey="phi", psikey="psi", zkey="free energy",
        mode="hist", bins=72, bandwidth=5, wrap=True, loop_edges=True,
        mask_cutoff=None, verbose=1, debug=0, **kwargs):
        """
        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          phikey (str): Key from which to load Φ
          psikey (str): Key from which to load Ψ
          zkey (str): Key from which to load distribution; if 'free
            energy' or 'probability', the 2D probability density of Φ
            and Ψ will be calculated and the selected representation
            returned; for other values a third dimension will be loaded
            from the `zkey` column of `infile`, the 3D probability
            density of Φ, Ψ, and `zkey` will be calculated, and the
            expectation value of `zkey` as a function of Φ and Ψ will be
            returned
          mode (str): Method of calculating probability distribution;
            may be either 'hist', to use a histogram, or 'kde', to use a
            kernel density estimate
          bins (int, list, ndarray): Bins to use for histogram or grid
            to use for kernel density estimate; if int, number of bins
            or gride points between -180° and 180° in Φ and Ψ, if list
            or ndarray, bins or grid directly
          bandwidth (float, optional): Bandwidth to use for kernel
            density estimate
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
          hist_kw: Keyword arguments passed to numpy.histogram2d
          kde_kw: Keyword arguments passed to
            sklearn.neighbors.KernelDensity
          verbose (int): Level of verbose output
          debug (int): Level of debug output
          kwargs (dict): Additional keyword arguments

        .. todo:
            - Auto-detect phikey and psikey
            - Support periodicty
            - Support variable bandwidth
        """
        import numpy as np
        from .myplotspec import multi_get_copy

        # Manage arguments
        if mode not in ["hist", "kde"]:
            raise ValueError("Argument 'mode' does not support provided " +
              "value '{0}', may be 'hist' or 'kde'".format(mode))
        read_csv_kw = dict(delim_whitespace=True, index_col=0)
        read_csv_kw.update(kwargs.pop("read_csv_kw", {}))

        # Load data
        dataframe = self.load_dataset(verbose=verbose, debug=debug,
          read_csv_kw=read_csv_kw, **kwargs).data
        if wrap:
            dataframe[phikey][dataframe[phikey] > 180] -= 360
            dataframe[psikey][dataframe[psikey] > 180] -= 360

        # Analyze
        if mode == "hist" and zkey in ["free energy", "probability"]:
            hist_kw = multi_get_copy("hist_kw", kwargs, {})
            if "bins" in hist_kw:
                bins = hist_kw.pop("bins")
            if isinstance(bins, int):
                bins = np.linspace(-180, 180, bins + 1)

            count, x_bins, y_bins = np.histogram2d(
              dataframe[phikey], dataframe[psikey], bins=bins, **hist_kw)

            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            y_centers = (y_bins[:-1] + y_bins[1:]) / 2
            x_width = np.mean(x_centers[1:] - x_centers[:-1])
            y_width = np.mean(y_centers[1:] - y_centers[:-1])

            probability = count / np.nansum(count)
            probability /= np.nansum(probability)
            free_energy = -1 * np.log(probability)
            free_energy[np.isinf(free_energy)] = np.nan
            free_energy -= np.nanmin(free_energy)

        elif mode == "kde" and zkey in ["free energy", "probability"]:
            from sklearn.neighbors import KernelDensity

            kde_kw = multi_get_copy("kde_kw", kwargs, {})
            kde_kw["bandwidth"] = kde_kw.get("bandwidth", bandwidth)
            if isinstance(bins, int):
                x_bins = y_bins = np.linspace(-180, 180, bins + 1)
            elif isinstance(bins, list):
                x_bins = y_bins = np.array(bins)
            elif isintance(bins, np.ndarray):
                x_bins = y_bins = bins

            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            y_centers = (y_bins[:-1] + y_bins[1:]) / 2
            x_width = np.mean(x_centers[1:] - x_centers[:-1])
            y_width = np.mean(y_centers[1:] - y_centers[:-1])
            xg, yg = np.meshgrid(x_centers, y_centers)
            xyg = np.vstack([yg.ravel(), xg.ravel()]).T
            samples = np.column_stack((dataframe[phikey], dataframe[psikey]))

            kde = KernelDensity(**kde_kw)
            kde.fit(samples)
            probability_series = np.exp(kde.score_samples(xyg))
            probability = np.zeros((x_centers.size, y_centers.size),
                      np.float) * np.nan
            for phi, psi, p in np.column_stack((xyg, probability_series)):
                x_index = np.where(x_centers == phi)[0][0]
                y_index = np.where(y_centers == psi)[0][0]
                probability[x_index, y_index] = p

            probability /= np.nansum(probability)
            free_energy = -1 * np.log(probability)
            free_energy[np.isinf(free_energy)] = np.nan
            free_energy -= np.nanmin(free_energy)

        elif mode == "hist" and zkey not in ["free energy", "probability"]:
            raise TypeError("Support for 3D histogram calculation is not yet "
              "implemented")

        elif mode == "kde" and zkey not in ["free energy", "probability"]:
            from sklearn.neighbors import KernelDensity

            kde_kw = multi_get_copy("kde_kw", kwargs, {})
            kde_kw["bandwidth"] = kde_kw.get("bandwidth", bandwidth)
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
        self.x = dataframe[phikey]
        self.y = dataframe[psikey]

        # Prepare mask
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

        if debug >= 1:
            from .myplotspec.debug import db_s

            db_s("Probability min {0}".format(np.nanmin(probability)))
            db_s("Probability max {0}".format(np.nanmax(probability)))
            db_s("Free energy min {0}".format(np.nanmin(free_energy)))
            db_s("Free energy max {0}".format(np.nanmax(free_energy)))
