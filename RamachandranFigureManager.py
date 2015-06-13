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
          heatmap_kw:
            cmap: hot
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 5
          contour_kw:
            colors: '0.25'
            levels: [0, 1, 2, 3, 4, 5]
            linestyles: solid
          mask_kw:
            cmap: Greys_r
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 1
          outline_kw:
            color: black
    """

    presets = """
      diff:
        help: Difference between two datasets
        draw_dataset:
          heatmap_kw:
            cmap: 'seismic'
            vmin: -5
            vmax:  5
          contour_kw:
            levels: [-5,-4,-3,-2,-1,0,1,2,3,4,5]
          load_kw:
            max_fe: 5
          mask: True
          outline: True
      right_title:
        help: Additional subplot title on right side
        draw_subplot:
          y2ticks: [-180,-90,0,90,180]
          y2ticklabels: []
          y2label_kw:
            labelpad: 20
            rotation: 270
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
    def draw_dataset(self, subplot, infile=None, kind="WHAM", heatmap=True,
        nan_to_max=True, contour=True, mask=False, outline=False,**kwargs):
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

        # Draw heatmap
        if heatmap:
            heatmap_kw = kwargs.get("heatmap_kw", {}).copy()
            heatmap_free_energy = dataset.free_energy.copy()
            if nan_to_max and hasattr(dataset, "mask"):
                heatmap_free_energy[np.isnan(heatmap_free_energy)] = np.nanmax(
                  heatmap_free_energy)
            subplot.pcolormesh(dataset.x_bins, dataset.y_bins,
              heatmap_free_energy.T, zorder=0.1, **heatmap_kw)

        # Draw contour
        if contour:
            contour_kw = kwargs.get("contour_kw", {})
            if "levels" in kwargs:
                contour_kw["levels"] = kwargs.pop("color")
            elif "levels" not in contour_kw:
                contour_kw["levels"] = range(0,
                  int(np.ceil(np.nanmax(dataset.free_energy))))
            subplot.contour(dataset.x_centers, dataset.y_centers,
              dataset.free_energy.T, zorder=0.2, **contour_kw)

        # Draw mask
        if mask and hasattr(dataset, "mask"):
            mask_kw = kwargs.get("mask_kw", {}).copy()
            subplot.pcolormesh(dataset.x_bins, dataset.y_bins,
              dataset.mask.T, zorder=0.3, **mask_kw)

        # Draw outline
        if outline and hasattr(dataset, "mask"):
            outline_kw = kwargs.get("outline_kw", {}).copy()
            for x in range(dataset.free_energy.shape[0]):
                for y in range(dataset.free_energy.shape[1]):
                    if not dataset.mask[x,y]:
                        if (x != 0
                        and y != dataset.mask.shape[1] - 1
                        and dataset.mask[x-1,y]):
                            subplot.plot(
                              [dataset.x_bins[x], dataset.x_bins[x]],
                              [dataset.y_bins[y], dataset.y_bins[y+1]],
                              zorder=0.4, **outline_kw)
                        if (x != dataset.mask.shape[0] - 1
                        and y != dataset.mask.shape[1] - 1
                        and dataset.mask[x+1,y]):
                            subplot.plot(
                              [dataset.x_bins[x+1], dataset.x_bins[x+1]],
                              [dataset.y_bins[y],   dataset.y_bins[y+1]],
                              zorder=0.4, **outline_kw)
                        if (x != dataset.mask.shape[0] - 1
                        and y != 0
                        and dataset.mask[x,y-1]):
                            subplot.plot(
                              [dataset.x_bins[x], dataset.x_bins[x+1]],
                              [dataset.y_bins[y], dataset.y_bins[y]],
                              zorder=0.4, **outline_kw)
                        if (x != dataset.mask.shape[0] - 1
                        and y != dataset.mask.shape[1] - 1
                        and dataset.mask[x,y+1]):
                            subplot.plot(
                              [dataset.x_bins[x], dataset.x_bins[x+1]],
                              [dataset.y_bins[y+1], dataset.y_bins[y+1]],
                              zorder=0.4, **outline_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    RamachandranFigureManager().main()
