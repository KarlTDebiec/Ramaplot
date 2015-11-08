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
    def get_cache_key(cls, dataset_classes=None, mask_cutoff=None, *args,
        **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documented under :func:`__init__`.
        """
        from .myplotspec import multi_get_copy

        minuend_kw = multi_get_copy(["minuend", "minuend_kw"], kwargs, {})
        minuend_class = dataset_classes[minuend_kw["kind"].lower()]
        key = [cls, mask_cutoff, minuend_class.get_cache_key(**minuend_kw)]

        subtrahend_kw = multi_get_copy(["subtrahend", "subtrahend_kw"],
          kwargs, {})
        if isinstance(subtrahend_kw, dict):
            subtrahend_kw = [subtrahend_kw]
        for sh_kw in subtrahend_kw:
            sh_class = dataset_classes[sh_kw.pop("kind").lower()]
            key.append(sh_class.get_cache_key(**sh_kw))

        return tuple(key)

    def __init__(self, dataset_classes=None, mask_cutoff=None, verbose=1,
        debug=0, **kwargs):
        """
        Initializes dataset.

        Arguments:
          minuend_kw (dict): Keyword arguments used to generate first
            dataset
          subtrahend_kw (dict, list): Keyword arguments used to generate
            datasets to subtract; if dict, arguments for single dataset;
            if list, list of dicts of keyword arguments used to generate
            multiple datasets
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
        from .myplotspec import multi_get_copy

        # Load minuend
        minuend_kw = multi_get_copy(["minuend", "minuend_kw"], kwargs, {})
        minuend_class = dataset_classes[minuend_kw.pop("kind").lower()]
        if not "dataset_cache" in minuend_kw:
            minuend_kw["dataset_cache"] = kwargs.get("dataset_cache", {})
        minuend = self.load_dataset(minuend_class, verbose=verbose,
                    debug=debug, **minuend_kw)
        self.x_centers = copy(minuend.x_centers)
        self.y_centers = copy(minuend.y_centers)
        self.x_width = copy(minuend.x_width)
        self.y_width = copy(minuend.y_width)
        self.x_bins = copy(minuend.x_bins)
        self.y_bins = copy(minuend.y_bins)
        self.dist = copy(minuend.dist)

        # Prepare mask
        if mask_cutoff is None:
            self.mask = minuend.mask
        else:
            self.mask = np.ma.masked_where(np.logical_and(
              self.dist <= mask_cutoff,
              np.logical_not(np.isnan(self.dist))),
              np.ones_like(self.dist))

        # Subtract subtrahends
        subtrahend_kw = multi_get_copy(["subtrahend", "subtrahend_kw"],
          kwargs, {})
        if isinstance(subtrahend_kw, dict):
            subtrahend_kw = [subtrahend_kw]
        for sh_kw in subtrahend_kw:
            if not "dataset_cache" in sh_kw:
                sh_kw["dataset_cache"] = kwargs.get("dataset_cache", {})
            sh_class = dataset_classes[sh_kw.pop("kind").lower()]
            sh = self.load_dataset(sh_class, verbose=verbose, debug=debug,
              **sh_kw)

            # Validate comparability of datasets
            if not ((minuend.x_centers == sh.x_centers).all()
            and     (minuend.y_centers == sh.x_centers).all()
            and      minuend.x_width   == sh.x_width
            and      minuend.y_width   == sh.y_width
            and     (minuend.x_bins    == sh.x_bins).all()
            and     (minuend.y_bins    == sh.y_bins).all()):
                raise TypeError()

            # Prepare data
            self.dist -= sh.dist

            # Prepare mask
            #   Values that are NaN are all masked, as are values that
            #   are masked in both source datasets; two boolean arrays
            #   are prepared manually; np.logical_and() and
            #   np.logical_or() yield the same result when passed masked
            #   arrays, and are either broken or not intended to be used
            #   for this purpose, and np.ma.mask_or() does not appear to
            #   be usable either.
            if mask_cutoff is None:
                sh_mask = sh.mask
            else:
                sh_mask = np.ma.masked_where(np.logical_and(
                  sh.dist <= mask_cutoff,
                  np.logical_not(np.isnan(sh.dist))),
                  np.ones_like(sh.dist))
            self.mask = np.ma.masked_where(
              np.logical_or(self.mask != 1, sh_mask  != 1),
              np.ones_like(self.dist))
