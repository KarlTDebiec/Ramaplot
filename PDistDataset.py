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
        mask_cutoff=None, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documented under :func:`__init__`.
        """
        from os.path import expandvars

        if zkey in ["free energy", "probability"]:
            x_bins, y_bins = cls.process_bins_arg(bins, dim=2)
            bins = (tuple(x_bins), tuple(y_bins))
        else:
            x_bins, y_bins, z_bins = cls.process_bins_arg(bins, dim=3)
            bins = (tuple(x_bins), tuple(y_bins), tuple(z_bins))

        if mode == "hist":
            return (cls, expandvars(infile), phikey, psikey, zkey, mode, bins,
                    wrap, mask_cutoff)
        elif mode == "kde":
            return (cls, expandvars(infile), phikey, psikey, zkey, mode, bins,
                    bandwidth, wrap, mask_cutoff)

    @staticmethod
    def process_bins_arg(bins, dim=2):
        """
        Processes bin argument.

        Arguments:
          bins (int, list, ndarray): Bins to use for histogram or grid
            to use for kernel density estimate; if int, number of bins
            or gride points between -180° and 180° in Φ and Ψ, if list
            or ndarray, bins or grid directly

        Returns:
          out_bins (tuple): Processed bins
        """
        import numpy as np
        if dim == 2:
            if isinstance(bins, int):
                x_bins = y_bins = np.linspace(-180, 180, bins + 1)
            elif isinstance(bins, list):
                if len(bins) == 2:
                    if isinstance(bins[0], int):
                        x_bins = np.linspace(-180, 180, bins[0] + 1)
                    elif isinstance(bins[0], list):
                        x_bins = np.array(bins[0])
                    if isinstance(bins[1], int):
                        y_bins = np.linspace(-180, 180, bins[1] + 1)
                    elif isinstance(bins[1], list):
                        y_bins = np.array(bins[1])
                else:
                    x_bins = y_bins = np.array(bins)
            elif isinstance(bins, np.ndarray):
                x_bins = y_bins = bins
            return x_bins, y_bins
        elif dim == 3:
            if isinstance(bins, int):
                x_bins = y_bins = z_bins = np.linspace(-180, 180, bins + 1)
            elif isinstance(bins, list):
                if len(bins) == 2:
                    if isinstance(bins[0], int):
                        x_bins = y_bins = np.linspace(-180, 180, bins[0] + 1)
                    elif (isinstance(bins[0], list)
                    or    isinstance(bins[0], np.ndarray)):
                        x_bins = y_bins = np.array(bins[0])
                    if isinstance(bins[1], int):
                        z_bins = np.linspace(-180, 180, bins[1] + 1)
                    elif (isinstance(bins[1], list)
                    or    isinstance(bins[1], np.ndarray)):
                        z_bins = np.array(bins[1])
                elif len(bins) == 3:
                    if isinstance(bins[0], int):
                        x_bins = np.linspace(-180, 180, bins[0] + 1)
                    elif (isinstance(bins[0], list)
                    or    isinstance(bins[0], np.ndarray)):
                        x_bins = np.array(bins[0])
                    if isinstance(bins[1], int):
                        y_bins = np.linspace(-180, 180, bins[1] + 1)
                    elif (isinstance(bins[1], list)
                    or    isinstance(bins[1], np.ndarray)):
                        y_bins = np.array(bins[1])
                    if isinstance(bins[2], int):
                        z_bins = np.linspace(-180, 180, bins[2] + 1)
                    elif (isinstance(bins[2], list)
                    or    isinstance(bins[2], np.ndarray)):
                        z_bins = np.array(bins[2])
                else:
                    x_bins = y_bins = z_bins = np.array(bins)
            elif isinstance(bins, np.ndarray):
                x_bins = y_bins = z_bins = bins
            return x_bins, y_bins, z_bins
        else:
            raise TypeError()

    def __init__(self, phikey="phi", psikey="psi", zkey="free energy",
        mode="hist", bins=72, bandwidth=5, wrap=True, mask_cutoff=None,
        verbose=1, debug=0, **kwargs):
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
          wrap (bool): Wrap x and y coordinates between 180° and 360° to
            between -180° and 0°
          wrap_z (bool): Wrap z coordinates between -180° and 0 to between 180°
            and 360°; probably only useful for plotting ω
          mask_cutoff (float): Cutoff beyond which distribution is
            masked, if `zkey` is 'free energy', this is a the maximum
            free energy above which the mask will be set, and if `zkey`
            is 'probability', this is the minimum probability below
            which the mask will be set
          hist_kw: Keyword arguments passed to numpy.histogram2d or
            numpy.histogramdd
          kde_kw: Keyword arguments passed to
            sklearn.neighbors.KernelDensity
          verbose (int): Level of verbose output
          debug (int): Level of debug output
          kwargs (dict): Additional keyword arguments

        .. todo:
          - Fix and validate 3D KDE
          - Auto-detect phikey and psikey
          - Support periodicic kernel density estimate
          - Support variable bandwidth kernel density estimate
        """
        from os.path import expandvars
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

        # Option 1: Calculate probability and free energy of Φ, Ψ
        if zkey in ["free energy", "probability"]:
            x_bins, y_bins = self.process_bins_arg(bins, dim=2)
            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            y_centers = (y_bins[:-1] + y_bins[1:]) / 2
            x_width = np.mean(x_centers[1:] - x_centers[:-1])
            y_width = np.mean(y_centers[1:] - y_centers[:-1])

            # Option 1a: Use a histogram (fast but noisy)
            if mode == "hist":
                if verbose >= 1:
                    print("calculating probability distribution of " +
                          "'{0}' and '{1}' using a ".format(phikey, psikey) +
                          "histogram")
                hist_kw = dict(normed=False)
                hist_kw.update(kwargs.get("hist_kw", {}))
                hist_kw["bins"] = hist_kw.get("bins", [x_bins, y_bins])

                probability, _, _ = np.histogram2d(
                  dataframe[phikey], dataframe[psikey], **hist_kw)

            # Option 1b: Use a kernel density estimate (smooth but slow)
            elif mode == "kde":
                if verbose >= 1:
                    print("calculating probability distribution of " +
                          "'{0}' and '{1}' using a ".format(phikey, psikey) +
                          "kernel density estimate")
                from sklearn.neighbors import KernelDensity

                kde_kw = multi_get_copy("kde_kw", kwargs, {})
                kde_kw["bandwidth"] = kde_kw.get("bandwidth", bandwidth)
                xg, yg = np.meshgrid(x_centers, y_centers)
                xyg = np.vstack([yg.ravel(), xg.ravel()]).T
                samples = np.column_stack((dataframe[phikey],
                  dataframe[psikey]))

                kde = KernelDensity(**kde_kw)
                kde.fit(samples)
                probability_series = np.exp(kde.score_samples(xyg))
                probability = np.zeros((x_centers.size, y_centers.size))
                for phi, psi, p in np.column_stack((xyg, probability_series)):
                    x_index = np.where(x_centers == phi)[0][0]
                    y_index = np.where(y_centers == psi)[0][0]
                    probability[x_index, y_index] = p

            # Normalize and calculate free energy
            probability /= np.nansum(probability)
            free_energy = -1 * np.log(probability)
            free_energy[np.isinf(free_energy)] = np.nan
            free_energy -= np.nanmin(free_energy)

        # Option 2: Calculate mean value of a third observable as a
        #   function of Φ, Ψ
        else:
            x_bins, y_bins, z_bins = self.process_bins_arg(bins, dim=3)
            x_centers = (x_bins[:-1] + x_bins[1:]) / 2
            y_centers = (y_bins[:-1] + y_bins[1:]) / 2
            z_centers = (z_bins[:-1] + z_bins[1:]) / 2
            x_width = np.mean(x_centers[1:] - x_centers[:-1])
            y_width = np.mean(y_centers[1:] - y_centers[:-1])

            if kwargs.get("wrap_z"):
                dataframe[zkey][dataframe[zkey] < 0] += 360

            # Option 2a: Use a histogram (fast but noisy)
            if mode == "hist":
                if verbose >= 1:
                    print("calculating mean value  of '{0}'".format(zkey) +
                          "as a function of '{0}' and ".format(phikey) +
                          "'{0}' using a histogram".format(psikey))
                hist_kw = dict(normed=True)
                hist_kw.update(kwargs.get("hist_kw", {}))
                hist_kw["bins"] = hist_kw.get("bins", [x_bins, y_bins, z_bins])

                prob_xyz, _ = np.histogramdd(np.column_stack(
                  (dataframe[phikey], dataframe[psikey], dataframe[zkey])),
                  **hist_kw)
                probability = np.sum(prob_xyz, axis=2)
                prob_z_given_xy = prob_xyz / probability[:,:,np.newaxis]
                weighted_z = prob_z_given_xy*z_centers[np.newaxis,np.newaxis,:]
                mean_z = np.sum(weighted_z, axis=2)

            # Option 2b: Use a kernel density estimate (smooth but slow)
            elif mode == "kde":
                raise()
                from copy import copy
                from sklearn.neighbors import KernelDensity

                kde_kw = multi_get_copy("kde_kw", kwargs, {})
                kde_kw["bandwidth"] = kde_kw.get("bandwidth", bandwidth)
                # Only a single bandwidth is supported; scale z
                #   dimension to span range of 120-240
#                scale_range = 340
                z = copy(dataframe[zkey])
#                z -= z.min()        # shift bottom to 0
#                z_range = z.max()     # save max

#                z *= (scale_range / z_range)
#                z += (360 - scale_range) / 2    # Give buffer on top and bottom

                xg, yg, zg = np.meshgrid(x_centers, y_centers, z_centers)
                xyzg = np.vstack([xg.ravel(), yg.ravel(), zg.ravel()]).T
                samples = np.column_stack((dataframe[phikey],
                  dataframe[psikey], z))

                kde = KernelDensity(**kde_kw)
                kde.fit(samples)
                probability_series = np.exp(kde.score_samples(xyzg))
                prob_xyz = np.zeros((x_centers.size, y_centers.size,
                  z_centers.size), np.float) * np.nan
                for phi,psi,z,p in np.column_stack((xyzg, probability_series)):
                    x_index = np.where(x_centers == phi)[0][0]
                    y_index = np.where(y_centers == psi)[0][0]
                    z_index = np.where(z_centers == z)[0][0]
                    prob_xyz[x_index, y_index, z_index] = p
                prob_xyz /= np.sum(prob_xyz)
                probability = np.sum(prob_xyz, axis=2)
                prob_z_given_xy = prob_xyz / probability[:,:,np.newaxis]
                weighted_z = prob_z_given_xy*z_centers[np.newaxis,np.newaxis,:]
                mean_z = np.sum(weighted_z, axis=2)

#                mean_z -= (360 - scale_range) / 2  # Shift back down
#                mean_z *= (z_range / scale_range) # Back from degrees to E
#                free_energy *= 627.503              # Convert to kcal/mol

            # Normalize and calculate free energy
            probability /= np.nansum(probability)
            free_energy = -1 * np.log(probability)
            free_energy[np.isinf(free_energy)] = np.nan
            free_energy -= np.nanmin(free_energy)

        self.x_centers = x_centers
        self.y_centers = y_centers
        self.x_width = x_width
        self.y_width = y_width
        self.x_bins = x_bins
        self.y_bins = y_bins
        self.x = dataframe[phikey]
        self.y = dataframe[psikey]

        # Prepare mask
        if zkey == "probability":
            self.dist = probability
            if mask_cutoff is not None:
                self.mask = np.ma.masked_where(np.logical_and(
                  probability >= mask_cutoff,
                  np.logical_not(np.isnan(probability))),
                  np.ones_like(probability))
            else:
                self.mask = np.ma.masked_where(
                  np.logical_not(np.isnan(free_energy)),
                  np.ones_like(free_energy))
        elif zkey == "free energy":
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
        else:
            self.dist = mean_z
            if mask_cutoff is not None:
                self.mask = np.ma.masked_where(np.logical_and(
                  free_energy <= mask_cutoff,
                  np.logical_not(np.isnan(free_energy))),
                  np.ones_like(free_energy))
            else:
                self.mask = np.ma.masked_where(
                  np.logical_not(np.isnan(free_energy)),
                  np.ones_like(free_energy))

        # Debug output
        if debug >= 1:
            from .myplotspec.debug import db_s
            db_s("Size {0}".format(probability.size))
            db_s("Probability min {0}".format(np.nanmin(probability)))
            db_s("Probability max {0}".format(np.nanmax(probability)))
            db_s("Free energy min {0}".format(np.nanmin(free_energy)))
            db_s("Free energy max {0}".format(np.nanmax(free_energy)))
            db_s("Free energy nan {0}".format(np.sum(np.isnan(free_energy))))
            db_s("Masked {0}".format(np.sum(self.mask)))
