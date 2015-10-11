# -*- coding: utf-8 -*-
#   ramaplot.DiffDataset.py
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
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class DiffDataset(Dataset):
    """
    Manages difference datasets.
    """

    @classmethod
    def get_cache_key(cls, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documented under :func:`__init__`.
        """
        dataset_1_kw = kwargs.get("dataset_1_kw")
        dataset_2_kw = kwargs.get("dataset_2_kw")
        dataset_classes = kwargs.get("dataset_classes")
        dataset_1_class = dataset_classes[dataset_1_kw["kind"].lower()]
        dataset_2_class = dataset_classes[dataset_1_kw["kind"].lower()]

        try:
            return (cls, 
                    dataset_1_class.get_cache_key(**dataset_1_kw),
                    dataset_2_class.get_cache_key(**dataset_2_kw))
        except TypeError:
            return None

    def __init__(self, dataset_1_kw=None, dataset_2_kw=None,
        dataset_classes=None, dataset_cache=None, mask_cutoff=None, **kwargs):
        """
        Initializes dataset.

        Arguments:
          dataset_1_kw (dict): Keyword arguments used to generate first
            dataset
          dataset_2_kw (dict): Keyword arguments used to generate second
            dataset
          dataset_cache (dict, optional): Cache of previously-loaded
            datasets
          mask_cutoff (float): Maximum free energy; points in the
            nascent free energy difference surface at which the free
            energy of both dataset_1 and dataset_2 are greater than
            mask_cutoff will be set to NaN; used to hide distracting
            large differences in depopulated regions; Currently only
            supports free energy maximum and not probability minimum
        """
        from copy import copy
        import numpy as np

        # Load datasets
        kind_1 = dataset_1_kw.pop("kind")
        kind_2 = dataset_2_kw.pop("kind")
        dataset_1 = self.load_dataset(dataset_classes[kind_1],
        dataset_cache=dataset_cache, **dataset_1_kw)
        dataset_2 = self.load_dataset(dataset_classes[kind_2],
        dataset_cache=dataset_cache, **dataset_2_kw)

        # Validate comparability of datasets
        if not ((dataset_1.x_centers == dataset_2.x_centers).all()
        and     (dataset_1.y_centers == dataset_2.x_centers).all()
        and      dataset_1.x_width   == dataset_2.x_width
        and      dataset_1.y_width   == dataset_2.y_width
        and     (dataset_1.x_bins    == dataset_2.x_bins).all()
        and     (dataset_1.y_bins    == dataset_2.y_bins).all()):
            raise TypeError()

        # Prepare data
        self.x_centers = copy(dataset_1.x_centers)
        self.y_centers = copy(dataset_1.y_centers)
        self.x_width = copy(dataset_1.x_width)
        self.y_width = copy(dataset_1.y_width)
        self.x_bins = copy(dataset_1.x_bins)
        self.y_bins = copy(dataset_1.y_bins)
        self.dist = copy(dataset_1.dist) - dataset_2.dist

        # Prepare mask
        #   Values that are NaN are all masked, as are values that are
        #   masked in both source datasets; two boolean arrays are
        #   prepared manually; np.logical_and() and np.logical_or()
        #   yield the same result when passed masked arrays, and are
        #   either broken or not intended to be used for this purpose,
        #   and np.ma.mask_or() does not appear to be usable either.
        temp1 = np.zeros_like(dataset_1.dist)
        temp2 = np.zeros_like(dataset_2.dist)
        if mask_cutoff is None:
            dataset_1_mask = dataset_1.mask
            dataset_2_mask = dataset_2.mask
        else:
            dataset_1_mask = np.ma.masked_where(np.logical_and(
              dataset_1.dist <= mask_cutoff,
              np.logical_not(np.isnan(dataset_1.dist))),
              np.ones_like(dataset_1.dist))
            dataset_2_mask = np.ma.masked_where(np.logical_and(
              dataset_2.dist <= mask_cutoff,
              np.logical_not(np.isnan(dataset_2.dist))),
              np.ones_like(dataset_2.dist))
        temp1[dataset_1_mask != 1] = 1
        temp2[dataset_2_mask != 1] = 1
        self.mask = np.ma.masked_where(
          np.logical_and(
            np.logical_not(np.isnan(self.dist)),
            np.logical_or(temp1, temp2)),
          np.ones_like(self.dist))
