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

    .. todo:
      - Population preset
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
          xlabel: Φ
          xticks: [-180,-90,0,90,180]
          ylabel: Ψ
          yticks: [-180,-90,0,90,180]
          ylabel_kw:
            va: center
          y2ticks: [-180,-90,0,90,180]
          y2ticklabels: []
          y2label_kw:
            rotation: 270
        draw_dataset:
          heatmap_kw:
            cmap: afmhot
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 5
          partner_kw:
            position: right
          colorbar_kw:
            position: right
            orientation: vertical
            zticks: [0,1,2,3,4,5]
            ztick_params:
              bottom: off
              top: off
              left: off
              right: off
            zlabel: 'ΔG (kcal/mol)'
            zlabel_kw:
              rotation: 270
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
          plot_kw:
            marker: 'o'
            ls: None
            mew: 0
            mfc: [0.5,0.5,0.5]
            ms: 0.5
            rasterized: True
          label_kw:
            x: 165
            y: -170
            text_kw:
              ha: right
              va: bottom
            border_lw: 3
    """

    available_presets = """
      ff99SB:
        help: Plot heatmap in style of ff99SB paper; omit contours
        draw_dataset:
          heatmap_kw:
            cmap: !!python/object/apply:ramaplot.cmap_ff99SB []
          contour: False
      potential_energy:
        help: Plot potential energy as a function of φ,ψ
        draw_dataset:
          heatmap: True
          heatmap_kw:
            cmap: bone
            vmin: 0
            vmax: 5
          colorbar_kw:
            zticks: [0,1,2,3,4,5]
            zlabel: 'ΔU (kcal/mol)'
          contour: True
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
          mask: True
          outline: False
      free_energy:
        help: Plot free energy as a function of φ,ψ
        draw_dataset:
          heatmap: True
          heatmap_kw:
            cmap: afmhot
            vmin: 0
            vmax: 5
          colorbar_kw:
            zticks: [0,1,2,3,4,5]
            zlabel: 'ΔG (kcal/mol)'
          contour: True
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
          mask: True
          outline: False
      diff:
        help: Plot difference between two datasets
        draw_dataset:
          kind: diff
          max_fe: 5
          heatmap_kw:
            cmap: RdBu_r
            vmin: -5
            vmax:  5
          colorbar_kw:
            zticks: [-5,-4,-3,-2,-1,0,1,2,3,4,5]
            zlabel: 'ΔΔG (kcal/mol)'
          contour_kw:
            levels: [-5,-4,-3,-2,-1,0,1,2,3,4,5]
          mask: True
          outline: True
      sampling:
        help: Plot sampling as a function of φ,ψ
        draw_dataset:
          heatmap: False
          heatmap_kw:
            cmap: afmhot_r
            vmin: 0
            vmax: 5
          contour: False
          mask: True
          mask_kw:
            cmap: Greys
          outline: False
          plot: True
      bond:
        help: Plot average value of a bond as a function of φ,ψ
        draw_dataset:
          kind: cdl
          heatmap_kw:
            cmap: RdBu
            vmin: 0
            vmax: 3
          contour: False
          mask: True
          outline: True
      bond_CN:
        help: Plot average C-N bond as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 1.31
            vmax: 1.35
          colorbar_kw:
            zticks: [1.32, 1.33, 1.34]
            zlabel: C-N (Å)
      bond_CN_extended:
        help: Plot average C-N bond as a function of φ,ψ; range extended to
              support proline's longer bond
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 1.31
            vmax: 1.40
          colorbar_kw:
            zticks: [1.32, 1.34, 1.36, 1.38]
            zlabel: C-N (Å)
      bond_NA:
        help: Plot average N-Cα bond as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 1.43
            vmax: 1.49
          colorbar_kw:
            zticks: [1.44, 1.45, 1.46, 1.47, 1.48]
            zlabel: N-Cα (Å)
      bond_AB:
        help: Plot average Cα-Cβ bond as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 1.51
            vmax: 1.55
          colorbar_kw:
            zticks: [1.52, 1.53, 1.54]
            zlabel: Cα-Cβ (Å)
      bond_AC:
        help: Plot average Cα-C bond as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 1.50
            vmax: 1.54
          colorbar_kw:
            zticks: [1.51, 1.52, 1.53]
            zlabel: Cα-C (Å)
      bond_CO:
        help: Plot average C-O bond as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 1.22
            vmax: 1.25
          colorbar_kw:
            zticks: [1.23, 1.24]
            zlabel: C-O (Å)
      angle:
        help: Plot average value of an angle as a function of φ,ψ
        draw_dataset:
          kind: cdl
          heatmap_kw:
            cmap: RdBu
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
          colorbar_kw:
            zticks: [119,121,123,125]
            zlabel: C-N-Cα (°)
      angle_NAB:
        help: Plot average N-Cα-Cβ angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 108
            vmax: 115
          colorbar_kw:
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
          colorbar_kw:
            zticks: [102,105,108,111,114]
            zlabel: N-Cα-Cβ (°)
      angle_NAC:
        help: Plot average N-Cα-C angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 106
            vmax: 117
          colorbar_kw:
            zticks: [107,109,111,113,115]
            zlabel: N-Cα-C (°)
      angle_BAC:
        help: Plot average Cβ-Cα-C angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 109
            vmax: 118
          colorbar_kw:
            zticks: [110,112,114,116]
            zlabel: Cβ-Cα-C (°)
      angle_ACO:
        help: Plot average Cα-C-O angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 117
            vmax: 124
          colorbar_kw:
            zticks: [118,119,120,121,122,123]
            zlabel: Cα-C-O (°)
      angle_ACN:
        help: Plot average Cα-C-N angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 113
            vmax: 122
          colorbar_kw:
            zticks: [114,116,118,120]
            zlabel: Cα-C-N (°)
      angle_OCN:
        help: Plot average OCN angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 118
            vmax: 126
          colorbar_kw:
            zticks: [119,121,123,125]
            zlabel: O-C-N (°)
      angle_OCN:
        help: Plot average OCN angle as a function of φ,ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 120
            vmax: 125
          colorbar_kw:
            zticks: [121,122,123,124]
            zlabel: O-C-N (°)
      omega:
        help: Plot average value of omega as a function of φ (i),ψ (i-1)
        draw_subplot:
          xlabel: '$Φ_{i}$'
          ylabel: '$Ψ_{i-1}$'
        draw_dataset:
          kind: cdl
          heatmap_kw:
            cmap: RdBu
            vmin: 170
            vmax: 190
          colorbar_kw:
            zticks: [170,175,180,185,190]
            zlabel: ω
          contour: False
          mask: True
          outline: True
      image:
        help: Plot image of a Ramachandran plot, typically from a publication
        draw_dataset:
          kind: image
          heatmap_kw:
            cmap: Greys_r
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
          colorbar_kw:
            ztick_fp: 8r
            zlabel_fp: 10b
            ztick_params:
              pad: 2
            zlabel_kw:
              labelpad: 14
          label_kw:
            fp: 10b
            border_lw: 4
      presentation_wide:
        inherits: presentation_wide
        draw_figure:
          left:       1.50
          sub_width:  3.40
          wspace:     0.30
          sub_height: 3.40
          hspace:     0.30
          bottom:     1.20
          multiplot: True
        draw_subplot:
          legend: False
          ylabel_kw:
            rotation: horizontal
            labelpad: 10
          y2label_kw:
            labelpad: 20
        draw_dataset:
          contour_kw:
            linewidths: 2
          plot_kw:
            ms: 2
          colorbar_kw:
            ztick_fp:  20r
            zlabel_fp: 24b
            ztick_params:
              pad: 5
            zlw: 3
            zlabel_kw:
              labelpad: 30
          label_kw:
            fp: 24b
            border_lw: 3
      colorbar_right:
        help: Draw colorbar to right of plot
        draw_dataset:
          colorbar: True
          partner_kw:
            position: right
          colorbar_kw:
            position: right
            orientation: vertical
            zlabel_kw:
              rotation: 270
              labelpad: 14
      colorbar_top:
        help: Draw colorbar above plot
        draw_dataset:
          colorbar: True
          partner_kw:
            position: top
          colorbar_kw:
            position: top
            orientation: horizontal
            zlabel_kw:
              rotation: 0
              labelpad: 5
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

        # Draw heatmap and colorbar
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
                pcolormesh = subplot.pcolormesh(dataset.x_bins, dataset.y_bins,
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
                contour = subplot.contour(dataset.x_centers, dataset.y_centers,
                  dataset.dist.T, zorder=0.2, **contour_kw)
                for collection in contour.collections:
                    for path in collection.get_paths():
                        if np.all(path.vertices[0] == path.vertices[-1]):
                            path.vertices = np.append(path.vertices,
                              [path.vertices[1]], axis=0)

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
