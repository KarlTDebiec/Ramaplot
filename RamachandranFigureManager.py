#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_forcefield.RamachandranFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more Ramachandran figures to specifications provided in
a YAML file.
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
class RamachandranFigureManager(FigureManager):
    """
    Manages the generation of Ramachandran figures.
    """
    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          xlabel: '$\\Phi$'
          xticks: [-180,-90,0,90,180]
          ylabel: '$\\Psi$'
          yticks: [-180,-90,0,90,180]
          ylabel_kw:
            va: center
        draw_dataset:
          contour_kw:
            colors: '0.25'
            levels: [0, 1, 2, 3, 4, 5]
            linestyles: solid
          imshow_kw:
            cmap: 'hot'
            interpolation: 'none'
            vmin: 0
            vmax: 5
    """

    presets = """
      diff:
        help: Difference between two datasets
        draw_dataset:
          contour_kw:
            colors: '0.25'
            levels: [-5,-4,-3,-2,-1,0,1,2,3,4,5]
          imshow_kw:
            cmap: 'seismic'
            vmin: -5
            vmax:  5
          nan_outline: True
      poster:
        help: Single plot for poster (width = 4.6", height = 4.3")
        inherits: poster
        draw_figure:
          left:       1.20
          sub_width:  3.00
          right:      0.40
          top:        0.30
          sub_height: 3.00
          bottom:     1.00
        draw_subplot:
          ylabel_kw:
            rotation: horizontal
      notebook:
        help: Single plot for notebook (width ≤ 6.5", height ≤ 8")
        inherits: notebook
        draw_figure:
          left:       0.50
          sub_width:  2.80
          right:      0.20
          top:        0.30
          sub_height: 2.80
          bottom:     0.40
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
      notebook_2:
        help: Two adjacent plots
        extends: notebook
        draw_figure:
          ncols: 2
          wspace: 0.3
          subplots:
            1:
              ylabel: ""
              yticklabels: []
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
      presentation_wide_6:
        help: Six plots for 16:9 presentation (width = 19.20", height =
              10.80")
        inherits: presentation_wide
        draw_figure:
          ncols: 3
          nrows: 2
          left:       1.50
          sub_width:  3.40
          wspace:     0.30
          sub_height: 3.40
          hspace:     0.30
          bottom:     1.20
          subplots:
            0:
              xlabel: ""
              xticklabels: []
            1:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            2:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            3:
              xticklabels: [-180,-90,0,90]
              yticklabels: [-180,-90,0,90]
            4:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            5:
              ylabel: ""
              yticklabels: []
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
            labelpad: 10
      presentation_wide_6.1:
        help: Right and bottom subplots disabled
        extends: presentation_wide_6
        draw_figure:
          nsubplots: 1
          subplots:
            0:
              xlabel: '$\\Phi$'
              xticklabels: [-180,-90,0,90,180]
      presentation_wide_6.2:
        help: Right and bottom subplots disabled
        extends: presentation_wide_6
        draw_figure:
          nsubplots: 2
          subplots:
            0:
              xlabel: '$\\Phi$'
              xticklabels: [-180,-90,0,90]
            1:
              xlabel: '$\\Phi$'
              xticklabels: [-180,-90,0,90,180]
      presentation_wide_6.3:
        help: Bottom subplots disabled
        extends: presentation_wide_6
        draw_figure:
          nsubplots: 3
          subplots:
            0:
              xlabel: '$\\Phi$'
              xticklabels: [-180,-90,0,90]
            1:
              xlabel: '$\\Phi$'
              xticklabels: [-180,-90,0,90]
            2:
              xlabel: '$\\Phi$'
              xticklabels: [-180,-90,0,90,180]
      presentation_wide_6.4:
        help: Bottom-right subplots disabled
        extends: presentation_wide_6
        draw_figure:
          nsubplots: 4
          subplots:
            1:
              xlabel: '$\\Phi$'
              xticklabels: ["",-90,0,90]
            2:
              xlabel: '$\\Phi$'
              xticklabels: [-180,-90,0,90,180]
      presentation_wide_6.5:
        help: Bottom-right subplot disabled
        extends: presentation_wide_6
        draw_figure:
          nsubplots: 5
          subplots:
            2:
              xlabel: '$\\Phi$'
              xticklabels: ["",-90,0,90,180]
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, infile=None, kind="WHAM",
        nan_to_max=False, nan_outline=False,**kwargs):
        import numpy as np
        import six
        from .myplotspec import get_color
        from .WHAMDataset import WHAMDataset
        from .NDRDDataset import NDRDDataset
        from .DiffDataset import DiffDataset
        parsers = {"WHAM": WHAMDataset, "NDRD": NDRDDataset}

        # Load data
        if infile is None:
            return
        elif isinstance(infile, six.string_types):
            if "load_kw" in kwargs:
                load_kw = kwargs.get("load_kw")
            else:
                load_kw = kwargs
            dataset = parsers[kind](infile=infile, **load_kw)
        elif isinstance(infile, list) and len(infile) == 2:
            load_kw = kwargs.get("load_kw", {})
            if isinstance(load_kw, list) and len(load_kw) == 2:
                load_kw_1 = load_kw[0]
                load_kw_2 = load_kw[1]
            else:
                load_kw_1 = load_kw.copy()
                load_kw_2 = load_kw.copy()
            if isinstance(kind, six.string_types):
                dataset_1 = parsers[kind](infile=infile[0], **load_kw_1)
                dataset_2 = parsers[kind](infile=infile[1], **load_kw_2)
            elif isinstance(kind, list) and len(kind) == 2:
                dataset_1 = parsers[kind[0]](infile=infile[0], **load_kw_1)
                dataset_2 = parsers[kind[1]](infile=infile[1], **load_kw_2)
            dataset = DiffDataset(dataset_1, dataset_2)

#        if nan_to_max:
#            dataset.free_energy[np.isnan(dataset.free_energy)] = np.nanmax(
#              dataset.free_energy)

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

        if nan_outline:
            for x in range(imshow_free_energy.shape[0]):
                for y in range(imshow_free_energy.shape[1]):
                    if np.isnan(imshow_free_energy[x,y]):
                        if     (x != 0
                        and     y != imshow_free_energy.shape[1] - 1
                        and not np.isnan(imshow_free_energy[x-1,y])):
                            subplot.plot(
                              [dataset.y_bins[y], dataset.y_bins[y+1]],
                              [dataset.x_bins[x], dataset.x_bins[x]],
                              color="black")
                        if     (x != imshow_free_energy.shape[0] - 1
                        and     y != imshow_free_energy.shape[1] - 1
                        and not np.isnan(imshow_free_energy[x+1,y])):
                            subplot.plot(
                              [dataset.y_bins[y], dataset.y_bins[y+1]],
                              [dataset.x_bins[x+1], dataset.x_bins[x+1]],
                              color="black")
                        if     (x != imshow_free_energy.shape[0] - 1
                        and     y != 0
                        and not np.isnan(imshow_free_energy[x,y-1])):
                            subplot.plot(
                              [dataset.y_bins[y], dataset.y_bins[y]],
                              [dataset.x_bins[x], dataset.x_bins[x+1]],
                              color="black")
                        if     (x != imshow_free_energy.shape[0] - 1
                        and     y != imshow_free_energy.shape[1] - 1
                        and not np.isnan(imshow_free_energy[x,y+1])):
                            subplot.plot(
                              [dataset.y_bins[y+1], dataset.y_bins[y+1]],
                              [dataset.x_bins[x], dataset.x_bins[x+1]],
                              color="black")

#################################### MAIN #####################################
if __name__ == "__main__":
    RamachandranFigureManager().main()
