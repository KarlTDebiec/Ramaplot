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
Generates one or more MDGX figures to specifications provided in a YAML
file.
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
    Manages the generation of MDGX figures.
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
        draw_dataset:
          violin_kw:
            points:         1000
            widths:         0.1
            showmeans:      False
            showextrema:    False
          median_kw:
            color:          white
            lw:             2.0
          percentile_kw:
            color:          white
            lw:             1.0
    """

    presets = """
      poster_two:
        draw_figure:
          ncols:        2
          fig_width:    15.00
          left:          3.40
          sub_width:     4.15
          wspace:        0.20
          bottom:        0.70
          sub_height:    3.00
          top:           0.20
          subplots:
            1:
              ylabel: ""
              yticklabels: []
          shared_legend:
            left:       12.00
            sub_width:   3.00
            sub_height:  3.00
            bottom:      0.70
            legend_kw:
              frameon:      False
              legend_fp:    24r
              loc:          6
            legend_lw:  10
        draw_subplot:
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
      presentation_one_of_two:
        draw_figure:
          ncols:         1
          left:          5.20
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
      presentation_one_of_four:
        draw_figure:
          nsubplots:    1
          shared_legend:
            left:       5.04
      presentation_two_of_four:
        draw_figure:
          nsubplots:    2
      presentation_three_of_four:
        draw_figure:
          nsubplots:    3
      notebook:
        draw_figure:
          left:          1.00
          sub_width:     2.10
          right:         1.00
          top:           0.30
          sub_height:    1.30
          bottom:        0.20
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
          fig_width:     5.50
          left:          1.00
          sub_width:     2.10
          wspace:        0.10
          top:           0.30
          sub_height:    1.30
          bottom:        0.55
          title_fp:     10b
          label_fp:     10b
          legend_fp:    10b
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
          subplots:
            1:
              ylabel:       ""
              yticklabels:  []
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
    def draw_dataset(self, subplot, infile=None, label=None, handles=None,
        **kwargs):
        import numpy as np
        from .myplotspec import get_color
        from .MDGXDataset import MDGXDataset

        verbose = kwargs.get("verbose", 0)

        # Handle missing input gracefully
        if infile is None:
            return

        # Configure plot settings
        violin_kw = kwargs.get("violin_kw", kwargs.get("plot_kw", {}))
        percentile_kw = kwargs.get("percentile_kw", {})
        median_kw = kwargs.get("median_kw", {})
        if "color" in violin_kw:
            color = get_color(violin_kw.pop("color"))
        elif "color" in kwargs:
            color = get_color(kwargs.pop("color"))
        if "color" in percentile_kw:
            percentile_kw["color"] = get_color(percentile_kw.pop("color"))
        if "color" in median_kw:
            median_kw["color"] = get_color(median_kw.pop("color"))

        # Load data
        dataset = MDGXDataset(infile=infile, **kwargs)
        data = dataset.data
        error = np.abs(data["QM energy"] - data["MM energy"])
        percentiles = {
          "0": np.min(error),
          "25": np.percentile(error, 25),
          "50": np.percentile(error, 50),
          "75": np.percentile(error, 75),
          "100": np.max(error)}
        mean = np.mean(error)
        if verbose >= 1:
            print("min: {0:05.2f} ".format(percentiles["0"]) +
                  "25%: {0:05.2f} ".format(percentiles["25"]) +
                  "50%: {0:05.2f} ".format(percentiles["50"]) +
                  "75%: {0:05.2f} ".format(percentiles["75"]) +
                  "max: {0:05.2f} ".format(percentiles["100"]) +
                  "mean: {0:05.2f}".format(mean))

        # Plot
        x = kwargs.get("x")
        violin = subplot.violinplot(np.array(error), [x], **violin_kw)
        for body in violin["bodies"]:
            body.set_facecolors(color)
            body.set_edgecolors("none")
            body.set_alpha(1)
            body.set_zorder(1)
        subplot.plot([x-0.05, x+0.05], [percentiles["25"], percentiles["25"]],
          zorder=1.5, **percentile_kw)
        subplot.plot([x-0.05, x+0.05], [percentiles["50"], percentiles["50"]],
          zorder=1.5, **median_kw)
        subplot.plot([x-0.05, x+0.05], [percentiles["75"], percentiles["75"]],
          zorder=1.5, **percentile_kw)
        violin = subplot.violinplot(np.array(error), [x], **violin_kw)
        for body in violin["bodies"]:
            body.set_facecolors("none")
            body.set_edgecolors(color)
            body.set_alpha(1)
            body.set_zorder(2)
        handle = subplot.plot(-1, -1, color = color)[0]
        if handles is not None and label is not None:
            handles[label] = handle

#################################### MAIN #####################################
if __name__ == "__main__":
    MDGXFigureManager().main()
