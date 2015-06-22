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

    def __init__(self, dataset_1, dataset_2, verbose=1, debug=1, **kwargs):
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
        self.x_centers = dataset_1.x_centers.copy()
        self.y_centers = dataset_1.y_centers.copy()
        self.x_width = dataset_1.x_width.copy()
        self.y_width = dataset_1.y_width.copy()
        self.x_bins = dataset_1.x_bins.copy()
        self.y_bins = dataset_1.y_bins.copy()
        self.dist = dataset_1.dist.copy() - dataset_2.dist

        # Prepare mask
        #   Values that are NaN are all masked, as are values that are
        #   masked in both source datasets; two boolean arrays are
        #   prepared manually; np.logical_and() and np.logical_or()
        #   yield the same result when passed masked arrays, and are
        #   either broken or not intended to be used for this purpose,
        #   and np.ma.mask_or() does not appear to be usable either.
        temp1 = np.zeros_like(dataset_1.dist)
        temp2 = np.zeros_like(dataset_2.dist)
        temp1[dataset_1.mask != 1] = 1
        temp2[dataset_2.mask != 1] = 1
        self.mask = np.ma.masked_where(
          np.logical_and(
            np.logical_not(np.isnan(self.dist)),
            np.logical_or(temp1, temp2)),
          np.ones_like(self.dist))
