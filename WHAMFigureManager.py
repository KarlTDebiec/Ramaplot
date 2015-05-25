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
        help: Ramachandran (backbone Φ/Ψ) plot, contours at 0-5 kcal/mol
        draw_subplot:
          xlabel:           '$\\Phi$'
          xticks:           [-180, -90, 0, 90, 180]
          ylabel:           '$\\Psi$'
          yticks:           [-180, -90, 0, 90, 180]
        draw_dataset:
          nan_to_max:       True
          contour_kw:
            colors:         '0.25'
            levels:         [0, 1, 2, 3, 4, 5]
          imshow_kw:
            cmap:           'hot_r'
            extent:         [-180, 180, 180, -180]
            interpolation:  'none'
            vmin:           0
            vmax:           5
      poster:
        help: Single plot for poster (width = 4.6", height = 4.3")
        inherits: poster
        draw_figure:
          left:         1.2
          sub_width:    3.0
          right:        0.4
          top:          0.3
          sub_height:   3.0
          bottom:       1.0
        draw_subplot:
          ylabel_kw:
            rotation:   horizontal
      notebook:
        help: Single plot for notebook (width ≤ 6.5", height ≤ 8")
        inherits: notebook
        draw_figure:
          left:          0.50
          sub_width:     3.00
          right:         0.20
          top:           0.30
          sub_height:    3.00
          bottom:        0.40
        draw_subplot:
          legend:       False
          ylabel_kw:
            rotation:   horizontal
      notebook_2:
        help: Two adjacent plots
        extends: notebook
        draw_figure:
          ncols:        2
          wspace:       0.5
        draw_subplot:
          legend:       False
          ylabel_kw:
            rotation:   horizontal
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, infile=None, kind="WHAM",nan_to_max=True,
        **kwargs):
        import numpy as np
        from .myplotspec import get_color
        from .WHAMDataset import WHAMDataset
        from .NDRDDataset import NDRDDataset

        # Load data; handle missing data gracefully
        if infile is None:
            return

        if kind == "WHAM":
            dataset = WHAMDataset(infile=infile, **kwargs)
        elif kind == "NDRD":
            dataset = NDRDDataset(infile=infile, **kwargs)
        if nan_to_max:
            dataset.free_energy[np.isnan(dataset.free_energy)] = np.nanmax(
              dataset.free_energy)

        # Configure contour settings
        contour_kw = kwargs.get("contour_kw", {})
        if hasattr(dataset, "contour_x_centers"):
            contour_x_centers = dataset.contour_x_centers
        else:
            contour_x_centers = dataset.x_centers
        if hasattr(dataset, "contour_y_centers"):
            contour_y_centers = dataset.contour_y_centers
        else:
            contour_y_centers = dataset.y_centers
        if hasattr(dataset, "contour_free_energy"):
            contour_free_energy = dataset.contour_free_energy
        else:
            contour_free_energy = dataset.free_energy
        if "levels" in kwargs:
            contour_kw["levels"] = kwargs.pop("color")
        elif "levels" not in contour_kw:
            contour_kw["levels"] = range(0,
              int(np.ceil(np.nanmax(dataset.free_energy))))

        # Configure heatmap settings
        imshow_kw = kwargs.get("imshow_kw", {})
        if "extent" in kwargs:
            imshow_kw["extent"] = kwargs.pop("extent")
        elif "extent" not in imshow_kw:
            if hasattr(dataset, "imshow_extent"):
                imshow_kw["extent"] = dataset.imshow_extent
            else:
                imshow_kw["extent"] = kwargs.get("extent",
                  [np.min(dataset.x_bins), np.max(dataset.x_bins),
                   np.max(dataset.y_bins), np.min(dataset.y_bins)])
        if hasattr(dataset, "imshow_free_energy"):
            imshow_free_energy = dataset.imshow_free_energy
        else:
            imshow_free_energy = dataset.free_energy

        # Plot
        subplot.contour(contour_x_centers, contour_y_centers,
          contour_free_energy, **contour_kw)
        subplot.imshow(imshow_free_energy, **imshow_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    WHAMFigureManager().main()
