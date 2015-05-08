#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_forcefield.MDGXFigureManager.py
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
class MDGXFigureManager(FigureManager):
    """
    Manages the generation of MDGXFigureManager figures.

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
          xticks:       [0, 1]
          xticklabels:  []
          yticks:       [0,1,2,3,4,5]
          ylabel:       '$\\left|{U_{QM}-U_{MM}}\\right|$\n\n(kcal/mol)'
          ylabel_kw:
            va:         center
    """

    presets = """
      presentation_one_of_two:
        draw_figure:
          ncols:         1
          fig_width:    10.24
          fig_height:    7.68
          left:          5.20
          sub_width:     3.30
          sub_height:    2.50
          bottom:        3.00
          shared_legend:
            left:       8.54
            sub_width:  4.00
            sub_height: 2.50
            bottom:     3.00
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    14r
              loc:          6
            legend_lw:  5
        draw_subplot:
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
            labelpad:   50
        draw_dataset:
          plot_kw:
            lw:         2
      presentation_two:
        draw_figure:
          ncols:        2
          fig_width:    10.24
          fig_height:    7.68
          left:          1.70
          sub_width:     3.30
          wspace:        0.20
          sub_height:    2.50
          bottom:        3.00
          subplots:
            1:
              ylabel: ""
              yticklabels: []
          shared_legend:
            left:       8.54
            sub_width:  4.00
            sub_height: 2.50
            bottom:     3.00
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    14r
              loc:          6
            legend_lw:  5
        draw_subplot:
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
            labelpad:   50
        draw_dataset:
          plot_kw:
            lw:         2
      presentation_one_of_four:
        draw_figure:
          ncols:        2
          nrows:        2
          nsubplots:    1
          fig_width:    10.24
          fig_height:    7.68
          left:          1.70
          sub_width:     3.30
          wspace:        0.20
          sub_height:    2.30
          hspace:        0.20
          bottom:        1.00
          shared_legend:
            left:       5.04
            sub_width:  4.00
            sub_height: 2.50
            bottom:     3.50
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    14r
              loc:          6
            legend_lw:  5
        draw_subplot:
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
            labelpad:   50
        draw_dataset:
          plot_kw:
            lw:         2
      presentation_two_of_four:
        draw_figure:
          ncols:        2
          nrows:        2
          nsubplots:    2
          fig_width:    10.24
          fig_height:    7.68
          left:          1.70
          sub_width:     3.30
          wspace:        0.20
          sub_height:    2.30
          hspace:        0.20
          bottom:        1.00
          subplots:
            1:
              ylabel: ""
              yticklabels: []
          shared_legend:
            left:       8.54
            sub_width:  4.00
            sub_height: 2.50
            bottom:     3.50
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    14r
              loc:          6
            legend_lw:  5
        draw_subplot:
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
            labelpad:   50
        draw_dataset:
          plot_kw:
            lw:         2
      presentation_three_of_four:
        draw_figure:
          ncols:        2
          nrows:        2
          nsubplots:    3
          fig_width:    10.24
          fig_height:    7.68
          left:          1.70
          sub_width:     3.30
          wspace:        0.20
          sub_height:    2.30
          hspace:        0.20
          bottom:        1.00
          subplots:
            1:
              ylabel: ""
              yticklabels: []
          shared_legend:
            left:       8.54
            sub_width:  4.00
            sub_height: 2.50
            bottom:     3.50
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    14r
              loc:          6
            legend_lw:  5
        draw_subplot:
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
            labelpad:   50
        draw_dataset:
          plot_kw:
            lw:         2
      presentation_four:
        draw_figure:
          ncols:        2
          nrows:        2
          fig_width:    10.24
          fig_height:    7.68
          left:          1.70
          sub_width:     3.30
          wspace:        0.20
          sub_height:    2.30
          hspace:        0.20
          bottom:        1.00
          subplots:
            1:
              ylabel: ""
              yticklabels: []
            3:
              ylabel: ""
              yticklabels: []
          shared_legend:
            left:       8.54
            sub_width:  4.00
            sub_height: 2.50
            bottom:     3.50
            legend_kw:
              frameon:      False
              labelspacing: 0.5
              legend_fp:    14r
              loc:          6
            legend_lw:  5
        draw_subplot:
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
            labelpad:   50
        draw_dataset:
          plot_kw:
            lw:         2
      notebook:
        draw_figure:
          fig_width:    4.2
          left:         1.0
          sub_width:    2.2
          top:          0.3
          sub_height:   1.4
          bottom:       0.5
          title_fp:     10b
          label_fp:     10b
          legend_fp:    10b
          shared_legend:
            left:       3.3
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
    """

    @manage_defaults_presets()
    @manage_kwargs()
    def draw_dataset(self, subplot, infile=None, bandwidth=0.1, grid=None,
        label=None, handles=None, contains=None, **kwargs):
        import numpy as np
        from sklearn.neighbors import KernelDensity
        from .myplotspec import get_color
        from .MDGXDataset import MDGXDataset

        # Handle missing input gracefully
        if infile is None:
            return

        # Configure analysis settings
        if grid is None:
            grid = np.linspace(0, 10, 1000)

        # Configure plot settings
        plot_kw = kwargs.get("plot_kw", {})
        if "color" in plot_kw:
            plot_kw["color"] = get_color(plot_kw.pop("color"))
        elif "color" in kwargs:
            plot_kw["color"] = get_color(kwargs.pop("color"))

        # Load data
        dataset = MDGXDataset(infile=infile, **kwargs)
        data = dataset.data
        if contains is not None:
            data = dataset.data[dataset.data["restart"].str.contains(
              contains)]
#            print(sorted(list(set(data["restart"]))))
        error = np.abs(data["QM energy"] - data["MM energy"])
        print(np.percentile(error, 25),
        np.percentile(error, 50),
        np.percentile(error, 75),
        np.mean(error))

        # Plot
        x = kwargs.get("x")
        violin = subplot.violinplot(np.array(error), [x],
          points=1000, widths=0.1, showmeans=False, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolors(plot_kw["color"])
            body.set_edgecolors('none')
            body.set_alpha(1)
            body.set_zorder(1)
#        subplot.plot([x-0.05, x+0.05],
#          [np.mean(error), np.mean(error)],
#          color="white", lw=1.5, zorder=1.5)
        subplot.plot([x-0.05, x+0.05],
          [np.percentile(error, 25), np.percentile(error, 25)],
          color="white", lw=1.0, ls="-", zorder=1.5)
        subplot.plot([x-0.05, x+0.05],
          [np.percentile(error, 50), np.percentile(error, 50)],
          color="white", lw=3.0, zorder=1.5)
        subplot.plot([x-0.05, x+0.05],
          [np.percentile(error, 75), np.percentile(error, 75)],
          color="white", lw=1.0, ls="-", zorder=1.5)
        violin = subplot.violinplot(np.array(error), [x],
          points=1000, widths=0.1, showmeans=False, showextrema=False)
        for body in violin["bodies"]:
            body.set_facecolors('none')
            body.set_edgecolors(plot_kw["color"])
            body.set_alpha(1)
            body.set_zorder(2)
        handle = subplot.plot(-1, -1, **plot_kw)[0]
        if handles is not None and label is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    MDGXFigureManager().main()
