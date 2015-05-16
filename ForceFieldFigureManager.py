#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_forcefield.ForceFieldFigureManager.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Generates one or more force field figures to specifications provided in
a YAML file. Currently hardcoded for dihedrals.
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
class ForceFieldFigureManager(FigureManager):
    """
    Manages the generation of force field figures.
    """
    from .myplotspec.manage_defaults_presets import manage_defaults_presets
    from .myplotspec.manage_kwargs import manage_kwargs
    from myplotspec.debug import debug_arguments

    defaults = """
        draw_figure:
          subplot_kw:
            autoscale_on: False
        draw_subplot:
          xticks:       [-180,-135,-90,-45,0,45,90,135,180]
          xlabel:       Dihedral (°)
          yticks:       [0,1,2,3,4,5]
          ylabel:       'U\n\n(kcal/mol)'
          ylabel_kw:
            va:         center
    """

    presets = """
      poster:
        draw_figure:
          fig_width:    15.00
          left:          3.20
          sub_width:     4.15
          right:         1.00
          bottom:        1.20
          sub_height:    2.40
          top:           0.20
          shared_legend:
            left:       7.50
            sub_width:  3.00
            sub_height: 2.40
            bottom:     1.20
            legend_kw:
              frameon:      False
              legend_fp:    25r
              loc:          6
            legend_lw:  10
        draw_subplot:
          xticks: [-180,-90,0,90,180]
          title_fp:     36r
          label_fp:     36r
          ylabel_kw:
            rotation:   horizontal
            labelpad:   100
          tick_fp:      24r
          tick_params:
            length:     3
            width:      1
            pad:        10
          lw:           2
        draw_dataset:
          plot_kw:
            lw:         2
      presentation:
        draw_figure:
          fig_width:    10.24
          fig_height:    7.68
          left:          1.60
          sub_width:     3.20
          sub_height:    2.00
          bottom:        4.00
          shared_legend:
            left:       4.85
            sub_width:  4.00
            sub_height: 2.00
            bottom:     4.00
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    14r
              loc:          2
            legend_lw:  5
        draw_subplot:
          xticks: [-180,-90,0,90,180]
          title_fp:     18r
          label_fp:     18r
          tick_fp:      14r
          legend:       False
          tick_params:
            length:     3
            width:      1
            pad:        6
          lw:           2
          ylabel_kw:
            rotation:   horizontal
            labelpad:   48
        draw_dataset:
          plot_kw:
            lw:         2
      presentation_wide:
        draw_figure:
          fig_width:    19.2
          fig_height:   10.8
          left:          2.1
          sub_width:     6.0
          right:         2.0
          top:           1.9
          sub_height:    3.5
          bottom:        5.4
          title_fp:     24b
          label_fp:     24b
          legend_fp:    24r
          shared_legend:
            left:       8.2
            sub_width:  4.0
            sub_height: 2.9
            bottom:     5.4
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    24r
              loc:          2
            legend_lw:  5
        draw_subplot:
          title_fp:     24b
          label_fp:     24b
          tick_fp:      16r
          legend:       False
          tick_params:
            length:     5
            width:      2
            pad:        6
          lw:           3
          ylabel_kw:
            rotation: horizontal
            labelpad:      50
        draw_dataset:
          plot_kw:
            lw:         3
      notebook:
        draw_figure:
          left:          1.00
          sub_width:     2.10
          right:         1.00
          top:           0.30
          sub_height:    1.30
          bottom:        0.50
          title_fp:     10b
          label_fp:     10b
          legend_fp:    10b
          shared_legend:
            left:        3.20
            sub_width:   1.20
            sub_height:  1.30
            bottom:      0.50
            legend_lw:  3
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    8r
              loc:          6
        draw_subplot:
          title_fp:     10b
          label_fp:     10b
          tick_fp:      8r
          legend:       False
          ylabel_kw:
            rotation:   horizontal
            labelpad:   30
      notebook_two:
        draw_figure:
          ncols:        2
          left:          1.00
          sub_width:     2.10
          wspace:        0.10
          right:         0.20
          top:           0.30
          sub_height:    1.30
          bottom:        0.90
          title_fp:     10b
          label_fp:     10b
          legend_fp:    10b
          xlabel_kw:
            ha:         center
            bottom:     -0.33
          subplots:
            0:
              xticklabels:  [-180,-135,-90,-45,0,45,90,135]
            1:
              ylabel:       ""
              yticklabels:  []
              xticklabels:  ["",-135,-90,-45,0,45,90,135,180]
          shared_legend:
            left:        1.00
            sub_width:   4.30
            sub_height:  0.50
            bottom:      0.00
            legend_lw:  3
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    8r
              loc:          9
              ncol:         3
        draw_subplot:
          title_fp:     10b
          label_fp:     10b
          tick_fp:      8r
          legend:       False
          ylabel_kw:
            rotation:   horizontal
            labelpad:   30
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, parameters=None, label=None, handles=None,
        **kwargs):
        import numpy as np
        from .myplotspec import get_color

        # Handle missing input gracefully
        if parameters is None:
            return

        # Configure plot settings
        plot_kw = kwargs.get("plot_kw", {})
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw.pop("color"))
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        x = np.linspace(-180, 180, 100)
        y = np.zeros_like(x)
        for height, phase, periodicity in parameters:
            y += height * (1 + np.cos(np.deg2rad(np.abs(periodicity) * x
                            + phase)))
        y -= np.min(y)
        handle = subplot.plot(x, y, **plot_kw)[0]
        if handles is not None and label is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    ForceFieldFigureManager().main()
