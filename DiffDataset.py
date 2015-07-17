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

    @staticmethod
    def get_cache_key(*args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.
        """
        dataset_1_kw = kwargs.get("dataset_1_kw")
        dataset_2_kw = kwargs.get("dataset_2_kw")
        dataset_classes = kwargs.get("dataset_classes")
        dataset_1_class = dataset_classes[dataset_1_kw["kind"].lower()]
        dataset_2_class = dataset_classes[dataset_1_kw["kind"].lower()]

        try:
            return (DiffDataset, 
                    dataset_1_class.get_cache_key(**dataset_1_kw),
                    dataset_2_class.get_cache_key(**dataset_2_kw))
        except TypeError:
            return None

    @staticmethod
    def get_cache_message(cache_key):
        """
        """
        return "previously loaded from '{0}'".format(cache_key[1])

    @staticmethod
    def process_args(**kwargs):
        """
        """

    @staticmethod
    def load_dataset(infile, kind, dataset_cache=None, verbose=1, **kwargs):
        """
        """
        dataset_classes = kwargs.get("dataset_classes")
        dataset_class = dataset_classes[kind.lower()]

        if "dataset_cache" is not None:
            cache_key = dataset_class.get_cache_key(infile=infile, **kwargs)
            if cache_key in dataset_cache:
                if verbose >= 1:
                    print(dataset_class.get_cache_message(cache_key))
                return dataset_cache[cache_key]
            else:
                if verbose >= 1:
                    print("loading from '{0}'".format(infile))
                dataset_cache[cache_key] = dataset_class(infile=infile,
                  verbose=verbose-1, **kwargs)
                return dataset_cache[cache_key]
        else:
            if verbose >= 1:
                print("loading from '{0}'".format(infile))
            return dataset_class(infile=infile, verbose=verbose-1, **kwargs)

    def __init__(self, dataset_1_kw=None, dataset_2_kw=None,
        dataset_classes=None, dataset_cache=None, max_fe=None, **kwargs):
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
        from copy import copy
        import numpy as np

        # Load datasets
        if (dataset_1_kw is None or dataset_2_kw is None or
            "infile" not in dataset_1_kw or "infile" not in dataset_2_kw):
            return None
        dataset_1 = self.load_dataset(dataset_cache=dataset_cache,
                      dataset_classes=dataset_classes, **dataset_1_kw)
        dataset_2 = self.load_dataset(dataset_cache=dataset_cache,
                      dataset_classes=dataset_classes, **dataset_2_kw)

        # Validate comparability of datasets
        if not ((dataset_1.x_centers == dataset_2.x_centers).all()
        and     (dataset_1.y_centers == dataset_2.x_centers).all()
        and      dataset_1.x_width   == dataset_2.x_width
        and      dataset_1.y_width   == dataset_2.y_width
        and     (dataset_1.x_bins    == dataset_2.x_bins).all()
        and     (dataset_1.y_bins    == dataset_2.y_bins).all()):
            raise

        # Prepare data
        self.x_centers = copy(dataset_1.x_centers)
        self.y_centers = copy(dataset_1.y_centers)
        self.x_width   = copy(dataset_1.x_width)
        self.y_width   = copy(dataset_1.y_width)
        self.x_bins    = copy(dataset_1.x_bins)
        self.y_bins    = copy(dataset_1.y_bins)
        self.dist      = copy(dataset_1.dist) - dataset_2.dist

        # Prepare mask
        #   Values that are NaN are all masked, as are values that are
        #   masked in both source datasets; two boolean arrays are
        #   prepared manually; np.logical_and() and np.logical_or()
        #   yield the same result when passed masked arrays, and are
        #   either broken or not intended to be used for this purpose,
        #   and np.ma.mask_or() does not appear to be usable either.
        temp1 = np.zeros_like(dataset_1.dist)
        temp2 = np.zeros_like(dataset_2.dist)
        if max_fe is None:
            dataset_1_mask = dataset_1.mask
            dataset_2_mask = dataset_2.mask
        else:
            dataset_1_mask = np.ma.masked_where(np.logical_and(
              dataset_1.dist <= max_fe,
              np.logical_not(np.isnan(dataset_1.dist))),
              np.ones_like(dataset_1.dist))
            dataset_2_mask = np.ma.masked_where(np.logical_and(
              dataset_2.dist <= max_fe,
              np.logical_not(np.isnan(dataset_2.dist))),
              np.ones_like(dataset_2.dist))
        temp1[dataset_1_mask != 1] = 1
        temp2[dataset_2_mask != 1] = 1
        self.mask = np.ma.masked_where(
          np.logical_and(
            np.logical_not(np.isnan(self.dist)),
            np.logical_or(temp1, temp2)),
          np.ones_like(self.dist))
