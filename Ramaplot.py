#!/usr/bin/python
# -*- coding: utf-8 -*-
#   ramaplot.Ramaplot.py
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
    __package__ = str("ramaplot")
    import ramaplot
from .myplotspec.FigureManager import FigureManager
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
          y2ticks: [-180,-90,0,90,180]
          y2ticklabels: []
          y2label_kw:
            rotation: 270
        draw_dataset:
          heatmap: True
          heatmap_kw:
            cmap: hot
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 5
          contour: True
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
          mask: True
          mask_kw:
            cmap: Greys_r
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 1
          outline: True
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

    available_presets = """
      ff99SB:
        help: Plot heatmap in style of ff99SB paper; omit contours
        draw_dataset:
          heatmap_kw:
            cmap: !!python/object/apply:ramaplot.cmap_ff99SB []
          contour: False
      free_energy:
        help: Plot free energy as a function of φ,ψ
        draw_dataset:
          heatmap: True
          heatmap_kw:
            cmap: hot
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 5
          contour: True
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
          mask: True
          mask_kw:
            cmap: Greys_r
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 1
          outline: True
          outline_kw:
            color: black
      angle:
        help: Plot average value of a backbone angle as a function of φ,ψ
        draw_dataset:
          heatmap_kw:
            cmap: seismic
            vmin: 110
            vmax: 130
          contour: True
          contour_kw:
            levels: [110,115,120,125,130]
          mask: True
          outline: True
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
            cmap: Greys
            vmin: 0
            vmax: 1
          contour: False
          mask: False
          outline: False
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
          sub_width:  1.59
          wspace:     0.10
          right:      0.25
          bottom:     0.40
          sub_height: 1.59
          top:        0.25
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
            0:
              xticklabels: [-180,-90,0,90]
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
    def draw_dataset(self, subplot, label=None, kind="wham",
        nan_to_max=True, heatmap=True, contour=True, mask=False,
        outline=False, plot=False, verbose=1, debug=0, **kwargs):
        """
        """
        from copy import copy
        from warnings import warn
        import numpy as np
        import six
        from .myplotspec import get_color
        from .AnalyticalDataset import AnalyticalDataset
        from .CDLDataset import CDLDataset
        from .DiffDataset import DiffDataset
        from .ImageDataset import ImageDataset
        from .NDRDDataset import NDRDDataset
        from .PDistDataset import PDistDataset
        from .WHAMDataset import WHAMDataset
        dataset_classes = {"analytical": AnalyticalDataset,
                           "cdl":        CDLDataset,
                           "diff":       DiffDataset,
                           "image":      ImageDataset,
                           "ndrd":       NDRDDataset,
                           "pdist":      PDistDataset,
                           "wham":       WHAMDataset}

        # Load data
        kind = kind.lower()
        dataset_kw = kwargs.get("dataset_kw", kwargs)
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(dataset_classes[kind],
                    dataset_classes=dataset_classes,
                    verbose=verbose, debug=debug, **dataset_kw)
        if dataset is None:
            return

        # Draw heatmap
        if heatmap:
            if not (hasattr(dataset, "dist") and hasattr(dataset, "x_bins")
            and     hasattr(dataset, "y_bins")):
                warn("'heatmap' is enabled but dataset does not have the "
                     "necessary attributes 'mask', 'x_bins', and 'y_bins', "
                     "skipping.")
            else:
                heatmap_kw = copy(kwargs.get("heatmap_kw", {}))
                heatmap_dist = copy(dataset.dist)
                if nan_to_max:
                    heatmap_dist[np.isnan(heatmap_dist)] = np.nanmax(
                      heatmap_dist)
                subplot.pcolormesh(dataset.x_bins, dataset.y_bins,
                  heatmap_dist.T, zorder=0.1, **heatmap_kw)

        # Draw contour
        if contour:
            if not (hasattr(dataset, "dist") and hasattr(dataset, "x_centers")
            and     hasattr(dataset, "y_centers")):
                warn("'contour' is enabled but dataset does not have the "
                     "necessary attributes 'dist', 'x_centers', and "
                     "'y_centers', skipping.")
            else:
                contour_kw = copy(kwargs.get("contour_kw", {}))
                if "levels" in kwargs:
                    contour_kw["levels"] = kwargs.pop("color")
                elif "levels" not in contour_kw:
                    contour_kw["levels"] = range(0,
                      int(np.ceil(np.nanmax(dataset.dist))))
                subplot.contour(dataset.x_centers, dataset.y_centers,
                  dataset.dist.T, zorder=0.2, **contour_kw)

        # Draw mask
        if mask:
            if not (hasattr(dataset, "mask") and hasattr(dataset, "x_bins")
            and     hasattr(dataset, "y_bins")):
                warn("'mask' is enabled but dataset does not have the "
                     "necessary attributes 'mask', 'x_bins', and 'y_bins', "
                     "skipping.")
            else:
                mask_kw = copy(kwargs.get("mask_kw", {}))
                subplot.pcolormesh(dataset.x_bins, dataset.y_bins,
                  dataset.mask.T, zorder=0.3, **mask_kw)

        # Draw outline
        if outline:
            if not (hasattr(dataset, "mask") and hasattr(dataset, "x_bins")
            and     hasattr(dataset, "y_bins")):
                warn("'outline' is enabled but dataset does not have the "
                     "necessary attributes 'mask', 'x_bins', and 'y_bins', "
                     "skipping.")
            else:
                outline_kw = copy(kwargs.get("outline_kw", {}))
                for x in range(dataset.dist.shape[0]):
                    for y in range(dataset.dist.shape[1]):
                        if not dataset.mask[x,y]:
                            if (x != 0
                            and y != dataset.mask.shape[1]
                            and dataset.mask[x-1,y]):
                                subplot.plot(
                                  [dataset.x_bins[x], dataset.x_bins[x]],
                                  [dataset.y_bins[y], dataset.y_bins[y+1]],
                                  zorder=0.4, **outline_kw)
                            if (x != dataset.mask.shape[0] - 1
                            and y != dataset.mask.shape[1]
                            and dataset.mask[x+1,y]):
                                subplot.plot(
                                  [dataset.x_bins[x+1], dataset.x_bins[x+1]],
                                  [dataset.y_bins[y],   dataset.y_bins[y+1]],
                                  zorder=0.4, **outline_kw)
                            if (x != dataset.mask.shape[0]
                            and y != 0
                            and dataset.mask[x,y-1]):
                                subplot.plot(
                                  [dataset.x_bins[x], dataset.x_bins[x+1]],
                                  [dataset.y_bins[y], dataset.y_bins[y]],
                                  zorder=0.4, **outline_kw)
                            if (x != dataset.mask.shape[0]
                            and y != dataset.mask.shape[1] - 1
                            and dataset.mask[x,y+1]):
                                subplot.plot(
                                  [dataset.x_bins[x], dataset.x_bins[x+1]],
                                  [dataset.y_bins[y+1], dataset.y_bins[y+1]],
                                  zorder=0.4, **outline_kw)

        # Draw plot
        if plot:
            if not (hasattr(dataset, "x") and hasattr(dataset, "y")):
                warn("'plot' is enabled but dataset does not have the "
                     "necessary attributes 'x' and 'y', skipping.")
            else:
                plot_kw = copy(kwargs.get("plot_kw", {}))
                subplot.plot(dataset.x, dataset.y, **plot_kw)

        if label is not None:
            from .myplotspec.text import set_text
            label_kw = kwargs.get("label_kw", {}).copy()
            set_text(subplot, s=label, **label_kw)

    def main(self):
        """
        Provides command-line functionality.
        """
        import argparse
        from inspect import getmodule

        parser = argparse.ArgumentParser(
          description     = getmodule(self.__class__).__doc__,
          formatter_class = argparse.RawTextHelpFormatter)

        super(RamachandranFigureManager, self).main(parser=parser)
        
#################################### MAIN #####################################
if __name__ == "__main__":
    RamachandranFigureManager().main()
