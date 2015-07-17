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
          contour: True
          heatmap_kw:
            cmap: hot
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 5
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
          mask_kw:
            cmap: Greys_r
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 1
          outline_kw:
            color: black
          label_kw:
            x: 165
            y: -170
            text_kw:
              ha: right
              va: bottom
              bbox:
                facecolor: white
                alpha: 0.8
                lw: 0
    """

    presets = """
      diff:
        help: Plot difference between two datasets
        draw_dataset:
          kind: diff
          max_fe: 5
          heatmap_kw:
            cmap: seismic
            vmin: -5
            vmax:  5
          contour_kw:
            levels: [-5,-4,-3,-2,-1,0,1,2,3,4,5]
          mask: True
          outline: True
      image:
        help: Plot image of a Ramachandran plot, typically from a publication
        draw_dataset:
          heatmap_kw:
            cmap: hot
            vmin: 0
            vmax: 1
          contour: False
          mask: False
          outline: False
      ff99SB:
        help: Plot heatmap in style of ff99SB paper; omit contours
        draw_dataset:
          heatmap_kw:
            cmap: !!python/object/apply:myplotspec_forcefield.cmap_ff99SB []
          contour: False
      angle:
        help: Plot average value of a backbone angle as a function of φ,ψ
        draw_dataset:
          heatmap_kw:
            cmap: seismic
            vmin: 119
            vmax: 127
          contour: True
          contour_kw:
            levels: [115,116,117,118,119,120,121,122,123,124,125]
          mask: True
          outline: True
      right_title:
        help: Additional subplot title on right side
        draw_subplot:
          y2ticks: [-180,-90,0,90,180]
          y2ticklabels: []
          y2label_kw:
            rotation: 270
      poster:
        help: Single plot for poster (width = 4.6", height = 4.3")
        inherits: poster
        draw_figure:
          left:       1.20
          sub_width:  3.00
          right:      0.40
          bottom:     1.00
          sub_height: 3.00
          top:        0.30
        draw_subplot:
          ylabel_kw:
            rotation: horizontal
      notebook:
        help: Single plot for notebook (width ≤ 6.5", height ≤ 9")
        inherits: notebook
        draw_figure:
          left:       0.50
          sub_width:  2.80
          wspace:     0.30
          right:      0.20
          bottom:     0.40
          sub_height: 2.80
          top:        0.30
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
          y2label_kw:
            labelpad: 6
        draw_dataset:
          label_kw:
            fp: 8r
      notebook_2:
        help: Two adjacent plots
        extends: notebook
        draw_figure:
          ncols: 2
          subplots:
            1:
              ylabel: ""
              yticklabels: []
      notebook_3:
        help: Three adjacent plots for notebook (width ≤ 6.5", height ≤ 9")
        inherits: notebook
        draw_figure:
          ncols: 3
          nrows: 1
          left:       0.50
          sub_width:  1.59
          wspace:     0.10
          right:      0.25
          sub_height: 1.59
          hspace:     0.10
          bottom:     0.40
          top:        0.25
          subplots:
            0:
              xticklabels: [-180,-90,0,90]
            1:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            2:
              ylabel: ""
              yticklabels: []
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
          y2label_kw:
            labelpad: 6
        draw_dataset:
          label_kw:
            fp: 8r
      notebook_6:
        help: Six plots
        extends: notebook_3
        draw_figure:
          nrows: 2
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
      notebook_9:
        help: Nine plots
        extends: notebook_3
        draw_figure:
          nrows: 3
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
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            4:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            5:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            6:
              xticklabels: [-180,-90,0,90]
              yticklabels: [-180,-90,0,90]
            7:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            8:
              ylabel: ""
              yticklabels: []
      notebook_12:
        help: Twelve plots
        extends: notebook_3
        draw_figure:
          nrows: 4
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
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            4:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            5:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            6:
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            7:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            8:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            9:
              xticklabels: [-180,-90,0,90]
              yticklabels: [-180,-90,0,90]
            10:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            11:
              ylabel: ""
              yticklabels: []
      notebook_15:
        help: Fifteen plots (full page)
        extends: notebook_3
        draw_figure:
          nrows: 5
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
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            4:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            5:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            6:
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            7:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            8:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            9:
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            10:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            11:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            12:
              xticklabels: [-180,-90,0,90]
              yticklabels: [-180,-90,0,90]
            13:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            14:
              ylabel: ""
              yticklabels: []
      notebook_18:
        help: Eighteen plots
        extends: notebook_3
        draw_figure:
          nrows: 6
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
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            4:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            5:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            6:
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            7:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            8:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            9:
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            10:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            11:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            12:
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            13:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            14:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            15:
              xticklabels: [-180,-90,0,90]
              yticklabels: [-180,-90,0,90]
            16:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            17:
              ylabel: ""
              yticklabels: []
      notebook_landscape_5:
        help: Five adjacent plots for notebook (width ≤ 6.5", height ≤ 9")
        inherits: notebook
        draw_figure:
          ncols: 5
          nrows: 1
          left:       0.50
          sub_width:  1.59
          wspace:     0.10
          right:      0.25
          sub_height: 1.59
          hspace:     0.10
          bottom:     0.40
          top:        0.25
          subplots:
            0:
              xticklabels: [-180,-90,0,90]
            1:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            2:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            3:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            4:
              ylabel: ""
              yticklabels: []
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
          y2label_kw:
            labelpad: 6
        draw_dataset:
          label_kw:
            fp: 8r
      notebook_landscape_10:
        help: Ten plots
        extends: notebook_landscape_5
        draw_figure:
          nrows: 2
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
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            4:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            5:
              xticklabels: [-180,-90,0,90]
              yticklabels: [-180,-90,0,90]
            6:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            7:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            8:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            9:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
      notebook_landscape_15:
        help: Fifteen plots
        extends: notebook_landscape_5
        draw_figure:
          nrows: 3
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
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            4:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            5:
              xlabel: ""
              xticklabels: []
              yticklabels: [-180,-90,0,90]
            6:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            7:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            8:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            9:
              xlabel: ""
              xticklabels: []
              ylabel: ""
              yticklabels: []
            10:
              xticklabels: [-180,-90,0,90]
              yticklabels: [-180,-90,0,90]
            11:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            12:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            13:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
            14:
              xticklabels: [-180,-90,0,90]
              ylabel: ""
              yticklabels: []
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
          y2label_kw:
            labelpad: 20
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
    def draw_dataset(self, subplot, infile=None, label=None, kind="WHAM",
        nan_to_max=True, heatmap=True, contour=True, mask=False,
        outline=False, **kwargs):
        """
        """
        import numpy as np
        import six
        from .myplotspec import get_color
        from .AnalyticalDataset import AnalyticalDataset
        from .CDLDataset import CDLDataset
        from .DiffDataset import DiffDataset
        from .ImageDataset import ImageDataset
        from .NDRDDataset import NDRDDataset
        from .WHAMDataset import WHAMDataset
        dataset_classes = {"analytical": AnalyticalDataset,
                           "cdl":        CDLDataset,
                           "diff":       DiffDataset,
                           "image":      ImageDataset,
                           "ndrd":       NDRDDataset,
                           "wham":       WHAMDataset}

        # Load data (for DiffDataset, infile is None, and infiles are
        #   instead included in dataset_1_kw and dataset_2_kw
        if  (infile is None
        and ("dataset_1_kw" not in kwargs and "dataset_2_kw" not in kwargs)):
            return
        else:
            dataset = self.load_dataset(dataset_classes[kind.lower()],
                        infile=infile, dataset_classes=dataset_classes,
                        **kwargs)
            if dataset is None or not hasattr(dataset, "dist"):
                return

        # Draw heatmap
        if heatmap:
            heatmap_kw = kwargs.get("heatmap_kw", {}).copy()
            heatmap_dist = dataset.dist.copy()
            if nan_to_max:
                heatmap_dist[np.isnan(heatmap_dist)] = np.nanmax(
                  heatmap_dist)
            subplot.pcolormesh(dataset.x_bins, dataset.y_bins,
              heatmap_dist.T, zorder=0.1, **heatmap_kw)

        # Draw contour
        if contour:
            contour_kw = kwargs.get("contour_kw", {})
            if "levels" in kwargs:
                contour_kw["levels"] = kwargs.pop("color")
            elif "levels" not in contour_kw:
                contour_kw["levels"] = range(0,
                  int(np.ceil(np.nanmax(dataset.dist))))
            subplot.contour(dataset.x_centers, dataset.y_centers,
              dataset.dist.T, zorder=0.2, **contour_kw)

        # Draw mask
        if mask and hasattr(dataset, "mask"):
            mask_kw = kwargs.get("mask_kw", {}).copy()
            subplot.pcolormesh(dataset.x_bins, dataset.y_bins,
              dataset.mask.T, zorder=0.3, **mask_kw)

        # Draw outline
        if outline and hasattr(dataset, "mask"):
            outline_kw = kwargs.get("outline_kw", {}).copy()
            for x in range(dataset.dist.shape[0]):
                for y in range(dataset.dist.shape[1]):
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

        if label is not None:
            from .myplotspec.text import set_text
            label_kw = kwargs.get("label_kw", {}).copy()
            set_text(subplot, s=label, **label_kw)

#################################### MAIN #####################################
if __name__ == "__main__":
    RamachandranFigureManager().main()
