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
class ForceFieldFigureManager(FigureManager):
    """
    Manages the generation of ForceFieldFigureManager figures.

    Attributes:
      defaults (str, dict): Default arguments to :func:`draw_report`,
        :func:`draw_figure`, :func:`draw_subplot`, and
        :func:`draw_dataset` functions, in yaml format. Outer level (of
        indentation or keys) provides function names, and inner level
        provides default arguments to each function::

          defaults = \"\"\"
              method_1:
                method_1_arg_1: 1000
                method_1_arg_2: abcd
              method_2
                method_2_arg_1: 2000
                method_2_arg_2: efgh
              ...
          \"\"\"

      presets (str, dict): Available sets of preset arguments to
        :func:`draw_report`, :func:`draw_figure`, :func:`draw_subplot`,
        and :func:`draw_dataset` functions, in yaml format. Outer level
        (of indentation or keys) provides preset names, middle level
        provides function names, and inner level provides arguments to
        pass to each function when preset is active::

          presets = \"\"\"
            preset_1:
              method_1:
                method_1_arg_1: 1001
                method_1_arg_2: abcde
              method_2
                method_2_arg_1: 2001
                method_2_arg_2: efghi
            preset_2:
              method_1:
                method_1_arg_1: 1002
                method_1_arg_2: abcdef
              method_2
                method_2_arg_1: 2002
                method_2_arg_2: efghij
          \"\"\"
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
          fig_width:    4.2
          left:         0.8
          sub_width:    2.2
          top:          0.3
          sub_height:   1.4
          bottom:       0.5
          title_fp:     10b
          label_fp:     10b
          legend_fp:    10b
          shared_legend:
            left:       3.1
            sub_width:  1.2
            sub_height: 1.4
            bottom:     0.5
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
            labelpad:   20
      notebook_two:
        draw_figure:
          ncols:         2
          fig_width:     5.5
          left:          0.8
          sub_width:     2.2
          wspace:        0.1
          top:           0.3
          sub_height:    1.4
          bottom:        1.0
          title_fp:     10b
          label_fp:     10b
          legend_fp:    10b
          xlabel_kw:
            ha:         center
            bottom:     -0.35
          subplots:
            1:
              ylabel: ""
              yticklabels: []
          shared_legend:
            left:       0.9
            sub_width:  4.2
            sub_height: 0.6
            bottom:     0.0
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
            labelpad:   20
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
