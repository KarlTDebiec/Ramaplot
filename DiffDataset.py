# -*- coding: utf-8 -*-
#   myplotspec_forcefield.DiffDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages difference datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class DiffDataset(object):
    """
    Manages difference datasets.

    """

    def __init__(self, dataset_1, dataset_2, max_fe=5, verbose=1,
        debug=1, **kwargs):
        """
        Initializes dataset.

        Arguments:
          dataset_1 (Dataset): First dataset
          dataset_2 (Dataset): Second dataset, subtracted from first
          max_fe (float): Maximum free energy; points in the nascent
            free energy difference surface at which the free energy of
            both dataset_1 and dataset_2 are greater than max_fe will be
            set to NaN; used to hide distracting large differences in
            depopulated regions
        """
        from os.path import expandvars
        import pandas
        import numpy as np

        # Validate comparability of datasets
        if not ((dataset_1.x_centers == dataset_2.x_centers).all()
        and     (dataset_1.y_centers == dataset_2.x_centers).all()
        and      dataset_1.x_width   == dataset_2.x_width
        and      dataset_1.y_width   == dataset_2.y_width
        and     (dataset_1.x_bins    == dataset_2.x_bins).all()
        and     (dataset_1.y_bins    == dataset_2.y_bins).all()):
            raise

        # Prepare data
        self.x_centers = dataset_1.x_centers
        self.y_centers = dataset_1.y_centers
        self.x_width = dataset_1.x_width
        self.y_width = dataset_1.y_width
        self.x_bins = dataset_1.x_bins
        self.y_bins = dataset_1.y_bins
        self.free_energy = dataset_1.free_energy - dataset_2.free_energy
        if max_fe is not None:
            self.free_energy[dataset_1.free_energy >= max_fe] = np.nan
        self.probability = dataset_1.probability - dataset_2.probability

        # Format contour settings
        if hasattr(dataset_1, "contour_x_centers"):
            self.contour_x_centers = dataset_1.contour_x_centers
        if hasattr(dataset_1, "contour_y_centers"):
            self.contour_y_centers = dataset_1.contour_y_centers
        if (hasattr(dataset_1, "contour_free_energy")
        and hasattr(dataset_2, "contour_free_energy")):
            self.contour_free_energy = (dataset_1.contour_free_energy
                                     -  dataset_2.contour_free_energy)
            if max_fe is not None:
                self.contour_free_energy[np.logical_and(
                  dataset_1.contour_free_energy >= max_fe,
                  dataset_2.contour_free_energy >= max_fe)] = np.nan

        # Format heatmap settings
        if (hasattr(dataset_1, "imshow_free_energy")
        and hasattr(dataset_2, "imshow_free_energy")):
            self.imshow_free_energy = (dataset_1.imshow_free_energy
                                    -  dataset_2.imshow_free_energy)
            if max_fe is not None:
                self.imshow_free_energy[np.logical_and(
                  dataset_1.imshow_free_energy >= max_fe,
                  dataset_2.imshow_free_energy >= max_fe)] = np.nan
        if hasattr(dataset_1, "imshow_extent"):
            self.imshow_extent = dataset_1.imshow_extent
