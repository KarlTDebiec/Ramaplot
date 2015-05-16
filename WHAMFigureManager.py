#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_forcefield.WHAMFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more WHAM figures to specifications provided in a YAML
file.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
import matplotlib
matplotlib.use("agg")
if __name__ == "__main__":
    __package__ = str("myplotspec_forcefield")
    import myplotspec_forcefield
from myplotspec.FigureManager import FigureManager
################################### CLASSES ###################################
class WHAMFigureManager(FigureManager):
    """
    Manages the generation of WHAM figures.
    """
    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs
    from myplotspec.debug import debug_arguments

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          ylabel_kw:
            va:         center
    """

    presets = """
      ramachandran:
        draw_subplot:
          xlabel:           '$\\Phi$'
          xticks:           [-180, -90, 0, 90, 180]
          ylabel:           '$\\Psi$'
          yticks:           [-180, -90, 0, 90, 180]
        draw_dataset:
          contour_kw:
            colors:         '0.25'
            levels:         [0, 1, 2, 3, 4, 5, 6]
          imshow_kw:
            cmap:           'hot_r'
            extent:         [-180, 180, 180, -180]
            interpolation:  'none'
            vmin:           0
            vmax:           6
      poster:
        draw_figure:
          left:         1.2
          sub_width:    3.0
          right:        0.4
          top:          0.3
          sub_height:   3.0
          bottom:       1.0
        draw_subplot:
          title_fp:     36r
          label_fp:     36r
          tick_fp:      24r
          tick_params:
            pad:        10
          lw:           2
          ylabel_kw:
            rotation:   horizontal
      notebook:
        draw_figure:
          left:          0.50
          sub_width:     3.00
          right:         0.20
          top:           0.30
          sub_height:    3.00
          bottom:        0.40
          title_fp:     10b
          label_fp:     10b
          legend_fp:    10b
        draw_subplot:
          title_fp:     10b
          label_fp:     10b
          tick_fp:      8r
          legend:       False
          ylabel_kw:
            rotation:   horizontal
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, infile=None, nan_to_max=True, **kwargs):
        import numpy as np
        from .myplotspec import get_color
        from .WHAMDataset import WHAMDataset

        # Load data; handle missing data gracefully
        if infile is None:
            return
        dataset = WHAMDataset(infile=infile, **kwargs)
        if nan_to_max:
            dataset.free_energy[np.isnan(dataset.free_energy)] = np.nanmax(
              dataset.free_energy)

        # Configure plot settings
        contour_kw = kwargs.get("contour_kw", {})
        if "levels" in kwargs:
            contour_kw["levels"] = kwargs.pop("color")
        elif "levels" not in contour_kw:
            contour_kw["levels"] = range(0,
              int(np.ceil(np.nanmax(dataset.free_energy))))
        imshow_kw = kwargs.get("imshow_kw", {})
        if "extent" in kwargs:
            imshow_kw["extent"] = kwargs.pop("extent")
        elif "extent" not in imshow_kw:
            imshow_kw["extent"] = kwargs.get("extent",
                      [np.min(dataset.x_bins), np.max(dataset.x_bins),
                       np.max(dataset.y_bins), np.min(dataset.y_bins)])

        # Plot
        subplot.contour(dataset.x_centers, dataset.y_centers,
          dataset.free_energy, **contour_kw)
        subplot.imshow(dataset.free_energy, **imshow_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    WHAMFigureManager().main()
