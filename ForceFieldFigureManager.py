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

    presets = """
      presentation_wide:
        draw_figure:
          left:         2.0
          sub_width:    6.0
          right:        2.0
          top:          2.0
          sub_height:   3.5
          bottom:       2.0
          title_fp:     24b
          label_fp:     24b
        draw_subplot:
          title_fp:     24b
          label_fp:     24b
          tick_fp:      16r
          legend_fp:    16r
          legend:       False
          tick_params:
            length:     5
            width:      2
            pad:        6
          lw:           2
        draw_dataset:
          plot_kw:
            lw:         2
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, parameters, label=None, handles=None,
        **kwargs):
        import numpy as np
        from .myplotspec import get_color

        # Configure plot settings
        plot_kw = kwargs.get("plot_kw", {})
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw.pop("color"))
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        x = np.linspace(-180, 180, 100)
        y = np.zeros_like(x)
        for height, phase, periodicity in parameters:
            print(height, periodicity, phase)
            y += height * (1 + np.cos(np.deg2rad(np.abs(periodicity) * x
                            + phase)))
        handle = subplot.plot(x, y, **plot_kw)[0]
        if handles is not None and label is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    ForceFieldFigureManager().main()
