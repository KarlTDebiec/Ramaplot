# -*- coding: utf-8 -*-
#   ramaplot.MeanDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages mean datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class MeanDataset(Dataset):
    """
    Manages mean datasets.
    """

    @classmethod
    def get_cache_key(cls, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documented under :func:`__init__`.
        """
        try:
            dataset_classes = kwargs.get("dataset_classes")

            observation_kw = kwargs.get("observation_kw")
            if isinstance(observation_kw, dict):
                observation_kw = [observation_kw]
            for ob_kw in observation_kw:
                ob_class = dataset_classes[ob_kw.pop("kind").lower()]
                key.append(ob_class.get_cache_key(**ob_kw))

            return tuple(key)
        except TypeError:
            return None

    def __init__(self, 
        dataset_classes=None, mask_cutoff=None, **kwargs):
        """
        Initializes dataset.

        Arguments:
          observation_kw (dict, list): Keyword arguments used to
            generate datasets to average ; if dict, arguments for single
            dataset; if list, list of dicts of keyword arguments used to
            generate multiple datasets
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

        # Load dataset arguments
        observation_kw = multi_get_copy(["observation", "observation_kw"],
          kwargs, {})
        if isinstance(observation_kw, dict):
            observation_kw = [observation_kw]

        # Load first dataset
        n_observations = 1
        ob_kw = observation_kw.pop()
        ob_class = dataset_classes[ob_kw.pop("kind").lower()]
        if not "dataset_cache" in ob_kw:
            ob_kw["dataset_cache"] = kwargs.get("dataset_cache", {})
        ob = self.load_dataset(ob_class, **ob_kw)
        self.x_centers = copy(ob.x_centers)
        self.y_centers = copy(ob.y_centers)
        self.x_width   = copy(ob.x_width)
        self.y_width   = copy(ob.y_width)
        self.x_bins    = copy(ob.x_bins)
        self.y_bins    = copy(ob.y_bins)
        self.dist      = copy(ob.dist)

        # Prepare mask
        if mask_cutoff is None:
            self.mask = ob.mask
        else:
            self.mask = np.ma.masked_where(np.logical_and(
              self.dist <= mask_cutoff,
              np.logical_not(np.isnan(self.dist))),
              np.ones_like(self.dist))

        for ob_kw in observation_kw:
            n_observations += 1
            if not "dataset_cache" in ob_kw:
                ob_kw["dataset_cache"] = kwargs.get("dataset_cache", {})
            ob_class = dataset_classes[ob_kw.pop("kind").lower()]
            ob = self.load_dataset(ob_class, **ob_kw)

            # Validate comparability of datasets
            if not ((self.x_centers == ob.x_centers).all()
            and     (self.y_centers == ob.x_centers).all()
            and      self.x_width   == ob.x_width
            and      self.y_width   == ob.y_width
            and     (self.x_bins    == ob.x_bins).all()
            and     (self.y_bins    == ob.y_bins).all()):
                raise TypeError()

            # Prepare data
            self.dist += ob.dist

            # Prepare mask
            #   Values that are NaN are all masked, as are values that
            #   are masked in both source datasets; two boolean arrays
            #   are prepared manually; np.logical_and() and
            #   np.logical_or() yield the same result when passed masked
            #   arrays, and are either broken or not intended to be used
            #   for this purpose, and np.ma.mask_or() does not appear to
            #   be usable either.
            if mask_cutoff is None:
                ob_mask = ob.mask
            else:
                ob_mask = np.ma.masked_where(np.logical_and(
                  ob.dist <= mask_cutoff,
                  np.logical_not(np.isnan(ob.dist))),
                  np.ones_like(ob.dist))
            self.mask = np.ma.masked_where(
              np.logical_or(self.mask != 1, ob_mask  != 1),
              np.ones_like(self.dist))

        # Average
        self.dist /= n_observations
