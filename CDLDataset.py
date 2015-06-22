# -*- coding: utf-8 -*-
#   myplotspec_forcefield.CDLDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages Conformation-Dependent Library datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class CDLDataset(object):
    """
    Manages Conformation-dependent library datasets.

    Parses and organizes Conformation-Dependent Library datasets, as
    published in:
      Berkholz, Donald S., Shapovalov, Maxim V., Dunbrack Jr. Roland L.,
      Karplus, P. Andrew. Conformation Dependence of Backbone Geometry
      in Proteins. Structure. 2009. 17. 1316-1325.
    Data for use with this class may be obtained from
    <http://dunbrack.fccc.edu/nmhrcm/>_
    """

    def __init__(self, infile, loop_edges=False,
        **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          selection (str, list): Selected distribution and measurement
          loop_edges (bool): Copy edges to enable plotting to edge of
            plot
          kwargs (dict): additional keyword arguments
        """
        import re
        import pandas
        import numpy as np
        import six

        def format_regex(regex):
            """
            Formats regular expression for parsing selection.

            Arguments:
              regex (str): Regular expression string, for which
                '{distribution}'  and '{whitespace}' will be replaced
                with the appropriate expressions

            Returns:
              regex (Pattern): Compiled regular expression
            """
            return re.compile(regex.format(
                dataset = "(NonPGIV_nonxpro|IleVal_nonxpro"
                          "|Gly_nonxpro|Pro_nonxpro"
                          "|NonPGIV_xpro|IleVal_xpro|IleVal_xpro"
                          "|Gly_xpro|Pro_xpro)",
                measurement = "(CNA|NAB|NAC|BAC|ACO|ACN|OCN"
                              "|CN|NA|AB|AC|CO)",
                whitespace = "\s+"), re.IGNORECASE)

        def load_dataset(infile, selection, verbose=1, debug=0, **kwargs):
            """
            Loads selected dataset from selected infile.

            Arguments:
              infile (str): Path to text input file, may contain
                environment variables
              selection (str): Start of lines containing desired
                dataset
              verbose (int): Enable verbose output
              debug (int): Enable debug output

            Returns:
              dist (DataFrame): Selected dataset
            """
            from os.path import expandvars
            from cStringIO import StringIO

            if verbose > 0:
                print("loading '{0}' from '{1}'".format(selection, infile))
            s = StringIO()

            with open(expandvars(infile)) as f:
                for line in f:
                    if line.startswith(selection):
                        s.write(line)

            s.seek(0)
            dist = pandas.read_csv(s, delim_whitespace=True, header=None,
                     usecols=range(1,29), names=["phi", "psi", "mode", "N",
                     "CNA mean", "CNA sd", "NAB mean", "NAB sd", "NAC mean",
                     "NAC sd", "BAC mean", "BAC sd", "ACO mean", "ACO sd",
                     "ACN mean", "ACN sd", "OCN mean", "OCN sd", "CN mean",
                     "CN sd", "NA mean", "NA sd", "AB mean", "AB sd",
                     "AC mean", "AC sd", "CO mean", "CO sd"])
            return dist

        type_error_text = ("SELECTION NOT UNDERSTOOD")
        sel_dist = format_regex("^(?P<dataset>{dataset})$")
        sel_meas = format_regex("^(?P<measurement>{measurement})$")

        # Unpack list to string
        if isinstance(selection, six.string_types):
            pass
        elif isinstance(selection, list) and len(selection) == 1:
            pass
        elif isinstance(selection, list) and len(selection) == 1:
            sel_1 = selection[0]
            if re.match(sel_dist, selection):
                sel = re.match(sel_dist, selection).groupdict()["dataset"]
                data = load_dataset(infile, sel)
            else:
                raise TypeError(type_error_text)
        else:
            raise TypeError(type_error_text)
        sel = "CNA"
        if sel.endswith("mean"):
            sel = sel.rstrip("mean").strip()

        # Organize data
        x_centers = np.unique(data["phi"])
        y_centers = np.unique(data["psi"])
        x_width = np.mean(x_centers[1:] - x_centers[:-1])
        y_width = np.mean(y_centers[1:] - y_centers[:-1])
        dist = np.zeros((x_centers.size, y_centers.size), np.float) * np.nan
        mask = np.zeros((x_centers.size, y_centers.size), np.bool)
        for index, row in data.iterrows():
            x_index = np.where(x_centers == row["phi"])[0][0]
            y_index = np.where(y_centers == row["psi"])[0][0]
            dist[x_index, y_index] = row["{0} mean".format(sel)]
            if row["mode"] == "B":
                mask[x_index, y_index] = True

        # Loop data to allow contour lines to be drawn to edges
        if loop_edges:
            x_centers = np.concatenate(([x_centers[0] - x_width],
                                         x_centers,
                                        [x_centers[-1] + x_width]))
            y_centers = np.concatenate(([y_centers[0] - y_width],
                                         y_centers,
                                        [y_centers[-1] + y_width]))
            temp = np.zeros((x_centers.size, y_centers.size), np.float)*np.nan
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
            temp = np.zeros((x_centers.size, y_centers.size), np.bool)
            temp[1:-1,1:-1] = mask
            temp[1:-1,-1]   = mask[:,0]
            temp[-1,1:-1]   = mask[0,:]
            temp[1:-1,0]    = mask[:,-1]
            temp[0,1:-1]    = mask[-1,:]
            temp[0,0]       = mask[-1,-1]
            temp[-1,-1]     = mask[0,0]
            temp[0,-1]      = mask[-1,0]
            temp[-1,0]      = mask[0,-1]
            mask = temp

        self.x_centers = x_centers
        self.y_centers = y_centers
        self.x_width = x_width
        self.y_width = y_width
        self.x_bins = np.linspace(x_centers[0]  - x_width / 2,
                                  x_centers[-1] + x_width / 2,
                                  x_centers.size + 1)
        self.y_bins = np.linspace(y_centers[0]  - y_width / 2,
                                  y_centers[-1] + y_width / 2,
                                  y_centers.size + 1)
        self.dist = dist
        self.mask = np.ma.masked_where(mask, np.ones_like(dist))
        print(np.min(dist), np.max(dist))
