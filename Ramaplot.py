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
          multiplot: False
          multi_xticklabels: [-180,-90,0,90,180]
          multi_yticklabels: [-180,-90,0,90,180]
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
          contour: False
          mask: True
          outline: True
      angle_CNA:
        help: Plot average C-N-Cα angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 118
            vmax: 127
          zticks: [119,121,123,125]
          zlabel: C-N-Cα (°)
      angle_NAB:
        help: Plot average N-Cα-Cβ angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 108
            vmax: 115
          zticks: [109,110,111,112,113,114]
          zlabel: N-Cα-Cβ (°)
      angle_NAB_extended:
        help: Plot average N-Cα-Cβ angle as a function of φ,ψ; range
              extended to support proline's smaller angle
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 101
            vmax: 115
          zticks: [102,105,108,111,114]
          zlabel: N-Cα-Cβ (°)
      angle_NAC:
        help: Plot average N-Cα-C angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 106
            vmax: 117
          zticks: [107,109,111,113,115]
          zlabel: N-Cα-C (°)
      angle_BAC:
        help: Plot average B-Cα-C angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 109
            vmax: 118
          zticks: [110,112,114,116]
          zlabel: B-Cα-C (°)
      angle_ACO:
        help: Plot average Cα-C-O angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 117
            vmax: 124
          zticks: [118,119,120,121,122,123]
          zlabel: Cα-C-O (°)
      angle_ACN:
        help: Plot average Cα-C-N angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 113
            vmax: 122
          zticks: [114,116,118,120]
          zlabel: Cα-C-N (°)
      angle_OCN:
        help: Plot average OCN angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 118
            vmax: 126
          zticks: [119,121,123,125]
          zlabel: O-C-N (°)
      angle_OCN:
        help: Plot average OCN angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 120
            vmax: 125
          zticks: [121,122,123,124]
          zlabel: O-C-N (°)
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
          hspace:     0.10
          top:        0.25
          multiplot: True
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
          y2label_kw:
            labelpad: 6
        draw_dataset:
          label_kw:
            fp: 8r
          ztick_fp: 8r
          zlabel_fp: 10b
          ztick_params:
            pad: 2
            bottom: off
            top: off
            left: off
            right: off
      presentation_wide_6:
        help: Six plots for 16:9 presentation (width = 19.20", height =
              10.80")
        inherits: presentation_wide
        draw_figure:
          ncols: 3
          nrows: 2
          multiplot: True
          left:       1.50
          sub_width:  3.40
          wspace:     0.30
          sub_height: 3.40
          hspace:     0.30
          bottom:     1.20
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
            labelpad: 10
          y2label_kw:
            labelpad: 20
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, label=None, kind="wham",
        nan_to_max=True, heatmap=True, colorbar=False, contour=True,
        mask=False, outline=False, plot=False, verbose=1, debug=0, **kwargs):
        """
        """
        from copy import copy
        from warnings import warn
        import numpy as np
        import six
        from .myplotspec import get_color
        from .myplotspec.axes import set_colorbar
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
                pcolormesh=subplot.pcolormesh(dataset.x_bins, dataset.y_bins,
                  heatmap_dist.T, zorder=0.1, **heatmap_kw)
                if colorbar:
                    if not hasattr(subplot, "_mps_partner_subplot"):
                        from .myplotspec.axes import add_partner_subplot
                        add_partner_subplot(subplot, **kwargs)
                    set_colorbar(subplot, pcolormesh, **kwargs)

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
            label_kw = kwargs.get("label_kw", {})
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
