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

.. todo:
    - Fix contour closing; worked in some instances when I wrote it but
      does not seem to work in all instances
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
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
          multi_tick_params:
            direction: out
            left: on
            right: off
            bottom: on
            top: off
        draw_subplot:
          title_kw:
            verticalalignment: bottom
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
          tick_params:
            direction: out
            left: on
            right: off
            bottom: on
            top: off
          y2tick_params:
            left: off
            right: off
            bottom: off
            top: off
          grid: True
          grid_kw:
            b: True
            linestyle: '-'
            alpha: 0.1
          label_kw:
            horizontalalignment: right
            verticalalignment: bottom
            zorder: 10
        draw_dataset:
          dataset_kw:
            cls: ramaplot.PDistDataset.PDistDataset
          heatmap_kw:
            cmap: afmhot
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 5
            zorder: 0.1
          colorbar_kw:
            zticks: [0,1,2,3,4,5]
            ztick_params:
              bottom: off
              top: off
              left: off
              right: off
              pad: 1
            zlabel: 'ΔG (kcal/mol)'
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
            zorder: 0.2
          mask_kw:
            cmap: Greys_r
            edgecolors: none
            rasterized: True
            vmin: 0
            vmax: 1
            zorder: 0.3
          outline_kw:
            color: black
            zorder: 0.4
          plot_kw:
            marker: 'o'
            ls: None
            mew: 0
            mfc: [0.5,0.5,0.5]
            ms: 0.5
            rasterized: True
            zorder: 0.5
          label_kw:
            x: 165
            y: -170
            text_kw:
              ha: right
              va: bottom
            border_lw: 1
    """

    available_presets = """
      ff99SB:
        class: appearance
        help: Draw heatmap in style of AMBER ff99SB paper
        draw_subplot:
          grid_kw:
            color: white
            alpha: 0.2
        draw_dataset:
          heatmap_kw:
            cmap: !!python/object/apply:ramaplot.cmap_ff99SB []
          draw_contour: False
      potential_energy:
        class: content
        help: Plot potential energy as a function of Φ,Ψ
        draw_dataset:
          draw_heatmap: True
          heatmap_kw:
            cmap: bone
            vmin: 0
            vmax: 5
          colorbar_kw:
            zticks: [0,1,2,3,4,5]
            zlabel: ΔU (kcal/mol)
          draw_contour: True
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
          draw_mask: True
          draw_outline: False
      free_energy:
        class: content
        help: Plot free energy as a function of Φ,Ψ
        draw_dataset:
          dataset_kw:
            zkey: free energy
            mask_cutoff: 5
          draw_heatmap: True
          heatmap_kw:
            cmap: afmhot
            vmin: 0
            vmax: 5
          colorbar_kw:
            zticks: [0,1,2,3,4,5]
            zlabel: ΔG (kcal/mol)
          draw_contour: True
          contour_kw:
            colors: '0.25'
            levels: [1, 2, 3, 4, 5]
            linestyles: solid
          draw_mask: False
          draw_outline: False
      probability:
        class: content
        help: Plot probability as a function of Φ,Ψ
        draw_dataset:
          dataset_kw:
            zkey: probability
            mask_cutoff: 0.001
          draw_heatmap: True
          heatmap_kw:
            cmap: summer
            vmin: 0
            vmax: 0.02
          colorbar_kw:
            zticks: [0.000,0.005,0.010,0.015,0.02]
            zlabel: Probability
          draw_contour: True
          contour_kw:
            colors: '0.25'
            levels: [0.000,0.005,0.010,0.015,0.020]
            linestyles: solid
          draw_mask: True
          draw_outline: True
      state_populations:
        class: content
        help: Calculate state populations; plot state centers and cutoffs
        draw_dataset:
          dataset_kw:
            calc_populations: True
            state_radius: 45
            plot_populations: True
            default_label_kw:
              border_lw: 1
              text_kw:
                horizontalalignment: center
                verticalalignment: center
                fp: 6b
          heatmap_kw:
            vmax: 9
            cmap: cubehelix
          colorbar_kw:
            zlabel: "Assigned Conformation"
            zticks: [0,1,2,3,4,5,6,7,8]
            zticklabels:
              - "β"
              - "$PP_{II}$"
              - "ξ"
              - "γ'"
              - "α"
              - "$L_α$"
              - "γ"
              - "$PP_{II}$'"
              - "plateau"
          contour_kw:
            levels: [0,1,2,3,4,5,6,7,8]
          draw_contour: False
          draw_outline: True
          draw_mask:    True
          draw_plot:    False
          draw_label:   True
      diff:
        class: content
        help: Plot difference between two datasets
        draw_dataset:
          dataset_kw:
            cls: ramaplot.DiffDataset.DiffDataset
          dataset_1_kw:
            mask_cutoff: 5
          dataset_2_kw:
            mask_cutoff: 5
          heatmap_kw:
            cmap: RdBu_r
            vmin: -5
            vmax:  5
          colorbar_kw:
            zticks: [-5,-4,-3,-2,-1,0,1,2,3,4,5]
            zlabel: ΔΔG (kcal/mol)
          contour_kw:
            levels: [-5,-4,-3,-2,-1,0,1,2,3,4,5]
          draw_mask: True
          draw_outline: True
      sampling:
        class: content
        help: Plot sampling as a function of Φ,Ψ
        draw_dataset:
          dataset_kw:
            cls: ramaplot.PDistDataset.PDistDataset
          draw_heatmap: False
          draw_contour: False
          draw_mask: True
          mask_kw:
            cmap: Greys
          draw_outline: False
          draw_plot: True
      structure:
        class: content
        help: Plot observed Φ,Ψ from a known structure
        draw_dataset:
          dataset_kw:
            cls: ramaplot.PDistDataset.PDistDataset
            mode: none
          draw_heatmap: False
          draw_contour: False
          draw_mask: False
          draw_outline: False
          draw_plot: True
          plot_kw:
            mfc: [0.7,0.7,0.7]
            rasterized: False
            zorder: 9
      bond:
        class: content
        help: Plot average value of a bond as a function of Φ,Ψ
        draw_dataset:
          heatmap_kw:
            cmap: RdBu
            vmin: 0
            vmax: 3
          draw_contour: False
          draw_mask: True
          draw_outline: True
      bond_CN:
        help: Plot average C-N bond as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 1.31
            vmax: 1.35
          colorbar_kw:
            zticks: [1.32, 1.33, 1.34]
            zlabel: C-N (Å)
      bond_CN_extended:
        help: Plot average C-N bond as a function of Φ,Ψ; range extended to
              support proline's longer bond
        extends: bond
        draw_dataset:
          heatmap_kw:
            vmin: 1.31
            vmax: 1.40
          colorbar_kw:
            zticks: [1.32, 1.34, 1.36, 1.38]
            zlabel: C-N (Å)
      bond_NA:
        help: Plot average N-Cα bond as a function of Φ,Ψ
        extends: bond
        draw_dataset:
          heatmap_kw:
            vmin: 1.43
            vmax: 1.49
          colorbar_kw:
            zticks: [1.44, 1.45, 1.46, 1.47, 1.48]
            zlabel: N-Cα (Å)
      bond_AB:
        help: Plot average Cα-Cβ bond as a function of Φ,Ψ
        extends: bond
        draw_dataset:
          heatmap_kw:
            vmin: 1.51
            vmax: 1.55
          colorbar_kw:
            zticks: [1.52, 1.53, 1.54]
            zlabel: Cα-Cβ (Å)
      bond_AC:
        help: Plot average Cα-C bond as a function of Φ,Ψ
        extends: bond
        draw_dataset:
          heatmap_kw:
            vmin: 1.50
            vmax: 1.54
          colorbar_kw:
            zticks: [1.51, 1.52, 1.53]
            zlabel: Cα-C (Å)
      bond_CO:
        help: Plot average C-O bond as a function of Φ,Ψ
        extends: bond
        draw_dataset:
          heatmap_kw:
            vmin: 1.22
            vmax: 1.25
          colorbar_kw:
            zticks: [1.23, 1.24]
            zlabel: C-O (Å)
      angle:
        class: content
        help: Plot average value of an angle as a function of Φ,Ψ
        draw_dataset:
          heatmap_kw:
            cmap: RdBu
            vmin: 110
            vmax: 130
          contour: False
          mask: True
          outline: True
      angle_CNA:
        help: Plot average C-N-Cα angle as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 118
            vmax: 127
          colorbar_kw:
            zticks: [119,121,123,125]
            zlabel: C-N-Cα (°)
      angle_NAB:
        help: Plot average N-Cα-Cβ angle as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 108
            vmax: 115
          colorbar_kw:
            zticks: [109,110,111,112,113,114]
            zlabel: N-Cα-Cβ (°)
      angle_NAB_extended:
        help: Plot average N-Cα-Cβ angle as a function of Φ,Ψ; range
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
        help: Plot average N-Cα-C angle as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 106
            vmax: 117
          colorbar_kw:
            zticks: [107,109,111,113,115]
            zlabel: N-Cα-C (°)
      angle_BAC:
        help: Plot average Cβ-Cα-C angle as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 109
            vmax: 118
          colorbar_kw:
            zticks: [110,112,114,116]
            zlabel: Cβ-Cα-C (°)
      angle_ACO:
        help: Plot average Cα-C-O angle as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 117
            vmax: 124
          colorbar_kw:
            zticks: [118,119,120,121,122,123]
            zlabel: Cα-C-O (°)
      angle_ACN:
        help: Plot average Cα-C-N angle as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 113
            vmax: 122
          colorbar_kw:
            zticks: [114,116,118,120]
            zlabel: Cα-C-N (°)
      angle_OCN:
        help: Plot average OCN angle as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 118
            vmax: 126
          colorbar_kw:
            zticks: [119,121,123,125]
            zlabel: O-C-N (°)
      angle_OCN:
        help: Plot average OCN angle as a function of Φ,Ψ
        extends: angle
        draw_dataset:
          heatmap_kw:
            vmin: 120
            vmax: 125
          colorbar_kw:
            zticks: [121,122,123,124]
            zlabel: O-C-N (°)
      chi:
        class: content
        help: Plot average value of a Χ angle as a function of Φ,Ψ (work in
              progress)
        draw_dataset:
          heatmap_kw:
            cmap: cubehelix
            vmin: -180
            vmax:  180
          contour_kw:
            levels: [-180,-120,-60,0,60,120,180]
      omega:
        class: content
        help: Plot average value of omega as a function of Φ (i),Ψ (i-1)
        draw_subplot:
          xlabel: '$Φ_{i}$'
          ylabel: '$Ψ_{i-1}$'
        draw_dataset:
          heatmap_kw:
            cmap: RdBu
            vmin: 170
            vmax: 190
          colorbar_kw:
            zticks: [170,175,180,185,190]
            zlabel: ω
          draw_contour: False
          draw_mask: True
          draw_outline: True
      image:
        class: content
        help: Plot image of a Ramachandran plot, typically from a publication
        draw_dataset:
          kind: image
          heatmap_kw:
            cmap: Greys_r
            vmin: 0
            vmax: 1
          draw_contour: False
          draw_mask: False
          draw_outline: False
      poster:
        class: target
        help: Poster (width = 4.6", height = 4.3")
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
      analytical:
        class: content
        help: Potential energy calculated directly from force field
        draw_dataset:
          dataset_kw:
            cls: ramaplot.AnalyticalDataset.AnalyticalDataset
            draw_heatmap: True
            draw_contour: True
            contour_kw:
              colors: "0.7"
              levels: [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
            draw_mask: False
            draw_outline: False
            draw_plot: False
      wham:
        class: content
        help: Distribution calculated from umbrella sampling using the Weighted
          Histogram Analysis Method
        draw_dataset:
          dataset_kw:
            cls: ramaplot.WHAMDataset.WHAMDataset
      ndrd:
        class: content
        help: Neighbor-Dependent Ramachandran Dristribution
        draw_dataset:
          dataset_kw:
            cls: ramaplot.NDRDDataset.NDRDDataset
      manuscript:
        class: target
        inherits: manuscript
        draw_figure:
          left:       0.35
          sub_width:  0.89
          wspace:     0.05
          right:      0.12
          bottom:     0.35
          sub_height: 0.89
          hspace:     0.05
          top:        0.25
          multi_xticklabels: [-180,"",0,"",180]
          multi_yticklabels: [-180,"",0,"",180]
          title_kw:
            top: -0.1
            title_fp: 8b
        draw_subplot:
          legend: False
          xlabel_kw:
            labelpad: 3
          ylabel_kw:
            labelpad: 3
            rotation: horizontal
          y2label_kw:
            labelpad: 8
          grid_kw:
            linewidth: 0.5
            alpha: 0.3
            color: [0,0,0]
          draw_label: True
          label_kw:
            fp: 7b
            x: null
            y: null
            xabs: -0.020
            yabs:  0.015
            border_lw: 1
        draw_dataset:
          partner_kw:
            sub_width:  0.05
            wspace:     0.05
            sub_height: 0.05
            hspace:     0.05
          colorbar_kw:
            ztick_fp: 6r
            zlabel_fp: 8b
            zlabel_kw:
              labelpad: 10
          contour_kw:
            linewidths: 0.7
          outline_kw:
            lw: 0.7
          plot_kw:
            ms: 2
          label_kw:
            fp: 6b
      manuscript_tight:
        class: target
        extends: manuscript
        help: Tighter formatting for nine columns rather than seven
        draw_figure:
          sub_width:  0.625
          sub_height: 0.625
        draw_subplot:
          label_kw:
            fp: 6b
            x: null
            y: null
            xabs: -0.020
            yabs:  0.015
        draw_dataset:
          colorbar_kw:
            zlabel_fp: 7b
            zlabel_kw:
              labelpad: 10
          contour_kw:
            linewidths: 0.5
          plot_kw:
            ms: 1.5
      notebook:
        class: target
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
          multi_tick_params:
            left: on
            right: off
            bottom: on
            top: off
          title_kw:
            top: -0.1
        draw_subplot:
          legend: False
          xlabel_kw:
            labelpad: 4
          ylabel_kw:
            labelpad: 4
            rotation: horizontal
          tick_params:
            pad: 4
          y2label_kw:
            labelpad: 8
        draw_dataset:
          colorbar_kw:
            ztick_fp: 8r
            zlabel_fp: 10b
          label_kw:
            fp: 10b
      presentation:
        class: target
        inherits: presentation
        draw_figure:
          fig_width: 
          left:       0.80
          sub_width:  1.50
          wspace:     0.10
          fig_height:
          bottom:     0.60
          sub_height: 1.50
          hspace:     0.10
          top:        0.30
        draw_subplot:
          xlabel_kw:
            horizontalalignment: center
          ylabel_kw:
            rotation: horizontal
            labelpad: 10
          y2label_kw:
            labelpad: 20
          draw_label: True
          label_kw:
            fp: 14r
            x: null
            y: null
            xabs: -0.05
            yabs:  0.05
            border_lw: 2
        draw_dataset:
          contour_kw:
            linewidths: 1
          plot_kw:
            ms: 5
          partner_kw:
            wspace:    0.10
            sub_width: 0.10
          colorbar_kw:
            ztick_fp:  14r
            zlabel_fp: 18r
            ztick_params:
              pad: 5
            zlw: 2
            zlabel_kw:
              labelpad: 25
          label_kw:
            fp: 14r
            border_lw: 3
      presentation_wide:
        class: target
        inherits: presentation_wide
        draw_figure:
          left:       1.50
          sub_width:  3.40
          wspace:     0.30
          sub_height: 3.40
          hspace:     0.30
          bottom:     1.20
        draw_subplot:
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
        class: appearance
        help: Draw colorbar to right of plot
        draw_dataset:
          draw_colorbar: True
          partner_kw:
            position: right
            sub_height: null
          colorbar_kw:
            zlabel_kw:
              rotation: 270
      colorbar_top:
        class: appearance
        help: Draw colorbar above plot
        draw_dataset:
          draw_colorbar: True
          partner_kw:
            position: top
            sub_width: null
          colorbar_kw:
            zlabel_kw:
              rotation: 0
      colorbar_bottom:
        class: appearance
        help: Draw colorbar below plot
        draw_dataset:
          draw_colorbar: True
          partner_kw:
            position: bottom
            sub_width: null
          colorbar_kw:
            zlabel_kw:
              rotation: 0
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, label=None, loop_edges=True,
        nan_to_max=True,
        draw_heatmap=True, draw_colorbar=False, draw_contour=True,
        draw_mask=False, draw_outline=False, draw_plot=False, draw_label=True,
        verbose=1, debug=0, **kwargs):
        """
        Draws a dataset.

        Arguments:
          subplot (Axes): Axes on which to draw
          heatmap (bool): Draw heatmap; requires that dataset has
            attributes 'dist', 'x_bins', and 'y_bins'
          colorbar (bool): Draw colorbar; requires `heatmap` to be True
          contour (bool): Draw contours; requires that dataset has
            attributes 'dist', 'x_centers', and 'y_centers'
          mask (bool): Draw mask; requires that dataset has attributes
            'mask', 'x_bins', and 'y_bins'
          outline (bool): Draw outline; requires that dataset has
            attributes 'mask', 'x_bins', and 'y_bins'
          plot (bool): Draw points; requires that dataset has attributes
            'x' and 'y'
          kwargs (dict): Additional keyword arguments
        """
        from copy import copy
        from warnings import warn
        import numpy as np
        import six
        from .myplotspec import get_colors, multi_get_copy

        # Load data
        dataset_kw = multi_get_copy("dataset_kw", kwargs, {})
        if "infile" in kwargs:
            dataset_kw["infile"] = kwargs["infile"]
        dataset = self.load_dataset(verbose=verbose, debug=debug, **dataset_kw)
        dist      = dataset.dist      if hasattr(dataset,"dist")      else None
        x_bins    = dataset.x_bins    if hasattr(dataset,"x_bins")    else None
        y_bins    = dataset.y_bins    if hasattr(dataset,"y_bins")    else None
        x_width   = dataset.x_width   if hasattr(dataset,"x_width")   else None
        y_width   = dataset.y_width   if hasattr(dataset,"y_width")   else None
        x_centers = dataset.x_centers if hasattr(dataset,"x_centers") else None
        y_centers = dataset.y_centers if hasattr(dataset,"y_centers") else None
        dist_mask = dataset.mask      if hasattr(dataset,"mask")      else None
        x         = dataset.x         if hasattr(dataset,"x")         else None
        y         = dataset.y         if hasattr(dataset,"y")         else None
        if loop_edges:
            if x_centers is not None and y_centers is not None:
                x_centers = np.concatenate(
                              ([x_centers[0]  - x_width], x_centers,
                               [x_centers[-1] + x_width]))
                y_centers = np.concatenate(
                              ([y_centers[0]  - y_width], y_centers,
                               [y_centers[-1] + y_width]))
            if x_bins is not None and y_bins is not None:
                x_bins  = np.linspace(x_centers[0]  - x_width / 2,
                                      x_centers[-1] + x_width / 2,
                                      x_centers.size + 1)
                y_bins  = np.linspace(y_centers[0]  - y_width / 2,
                                      y_centers[-1] + y_width / 2,
                                      y_centers.size + 1)
            if dist is not None:
                temp = np.zeros((x_centers.size, y_centers.size)) * np.nan
                temp[1:-1,1:-1] = dist
                temp[1:-1,-1]   = dist[:,0]
                temp[-1,1:-1]   = dist[0,:]
                temp[1:-1,0]    = dist[:,-1]
                temp[0,1:-1]    = dist[-1,:]
                temp[0,0]       = dist[-1,-1]
                temp[-1,-1]     = dist[0,0]
                temp[0,-1]      = dist[-1,0]
                temp[-1,0]      = dist[0,-1]
                dist = temp

            if dist_mask is not None:
                temp = np.ma.empty(dist.shape)
                temp[1:-1,1:-1] = dist_mask
                temp[1:-1,-1]   = dist_mask[:,0]
                temp[-1,1:-1]   = dist_mask[0,:]
                temp[1:-1,0]    = dist_mask[:,-1]
                temp[0,1:-1]    = dist_mask[-1,:]
                temp[0,0]       = dist_mask[-1,-1]
                temp[-1,-1]     = dist_mask[0,0]
                temp[0,-1]      = dist_mask[-1,0]
                temp[-1,0]      = dist_mask[0,-1]
                dist_mask = temp

        # Draw heatmap and colorbar
        if draw_heatmap:
            if dist is None or x_bins is None or y_bins is None:
                warn("'heatmap' is enabled but dataset does not have the "
                     "necessary attributes 'mask', 'x_bins', and 'y_bins', "
                     "skipping.")
            else:
                heatmap_kw = copy(kwargs.get("heatmap_kw", {}))
                heatmap_dist = copy(dist)
                if nan_to_max:
                    heatmap_dist[np.isnan(heatmap_dist)] = np.nanmax(
                      heatmap_dist)
                pcolormesh = subplot.pcolormesh(x_bins, y_bins, heatmap_dist.T,
                  **heatmap_kw)
                if draw_colorbar:
                    from .myplotspec.axes import set_colorbar
                    if not hasattr(subplot, "_mps_partner_subplot"):
                        from .myplotspec.axes import add_partner_subplot
                        add_partner_subplot(subplot, verbose=verbose,
                          debug=debug, **kwargs)
                    set_colorbar(subplot, pcolormesh, **kwargs)

        # Draw contour
        if draw_contour:
            if dist is None or x_centers is None or y_centers is None:
                warn("'contour' is enabled but dataset does not have the "
                     "necessary attributes 'dist', 'x_centers', and "
                     "'y_centers', skipping.")
            else:
                contour_kw = copy(kwargs.get("contour_kw", {}))
                if "levels" in kwargs:
                    contour_kw["levels"] = kwargs.pop("color")
                elif "levels" not in contour_kw:
                    contour_kw["levels"] = range(0,
                      int(np.ceil(np.nanmax(dist))))
                contour = subplot.contour(x_centers, y_centers, dist.T,
                  **contour_kw)

                # Close bottom of contours
                for collection in contour.collections:
                    for path in collection.get_paths():
                        if np.all(path.vertices[0] == path.vertices[-1]):
                            try:
                                path.vertices = np.append(path.vertices,
                                  [path.vertices[1]], axis=0)
                            except IndexError:
                                continue

        # Draw mask
        if draw_mask:
            if dist_mask is None or x_bins is None or y_bins is None:
                warn("'mask' is enabled but dataset does not have the "
                     "necessary attributes 'mask', 'x_bins', and 'y_bins', "
                     "skipping.")
            else:
                mask_kw = copy(kwargs.get("mask_kw", {}))
                subplot.pcolormesh(x_bins, y_bins, dist_mask.T, **mask_kw)

        # Draw outline
        if draw_outline:
            if dist_mask is None or x_bins is None or y_bins is None:
                warn("'outline' is enabled but dataset does not have the "
                     "necessary attributes 'mask', 'x_bins', and 'y_bins', "
                     "skipping.")
            else:
                outline_kw = copy(kwargs.get("outline_kw", {}))
                for xo in range(dist.shape[0]):
                    for yo in range(dist.shape[1]):
                        if not dist_mask[xo,yo]:
                            if (xo != 0
                            and yo != dist_mask.shape[1]
                            and dist_mask[xo-1,yo]):
                                subplot.plot(
                                  [x_bins[xo], x_bins[xo]],
                                  [y_bins[yo], y_bins[yo+1]],
                                  **outline_kw)
                            if (xo != dist_mask.shape[0] - 1
                            and yo != dist_mask.shape[1]
                            and dist_mask[xo+1,yo]):
                                subplot.plot(
                                  [x_bins[xo+1], x_bins[xo+1]],
                                  [y_bins[yo],   y_bins[yo+1]],
                                  **outline_kw)
                            if (xo != dist_mask.shape[0]
                            and yo != 0
                            and dist_mask[xo,yo-1]):
                                subplot.plot(
                                  [x_bins[xo], x_bins[xo+1]],
                                  [y_bins[yo], y_bins[yo]],
                                  **outline_kw)
                            if (xo != dist_mask.shape[0]
                            and yo != dist_mask.shape[1] - 1
                            and dist_mask[xo,yo+1]):
                                subplot.plot(
                                  [x_bins[xo],   x_bins[xo+1]],
                                  [y_bins[yo+1], y_bins[yo+1]],
                                  **outline_kw)

        # Draw plot
        if draw_plot:
            if x is None or y is None:
                warn("'plot' is enabled but dataset does not have the "
                     "necessary attributes 'x' and 'y', skipping.")
            else:
                plot_kw = copy(kwargs.get("plot_kw", {}))
                get_colors(plot_kw, kwargs)
                subplot.plot(x, y, **plot_kw)

        # Draw label (need to clean this up)
        if draw_label:
            from .myplotspec.text import set_text
            if label is not None:
                label_kw = kwargs.get("label_kw", {})
                if (isinstance(label, six.string_types)
                and isinstance(label_kw, dict)):
                    set_text(subplot, s=label, **label_kw)
                elif (isinstance(label, list) and isinstance(label_kw, list)
                and len(label) == len(label_kw)):
                    for l, kw in zip(label, label_kw):
                        set_text(subplot, s=l, **kw)
                else:
                    raise Exception("bad label arguments")
            elif hasattr(dataset, "label") and dataset.label is not None:
                label = dataset.label
                label_kw = dataset.label_kw
                if (isinstance(label, six.string_types)
                and isinstance(label_kw, dict)):
                    set_text(subplot, s=label, **label_kw)
                elif (isinstance(label, list) and isinstance(label_kw, list)
                and len(label) == len(label_kw)):
                    for l, kw in zip(label, label_kw):
                        set_text(subplot, s=l, **kw)
                else:
                    raise Exception("bad label arguments")


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
