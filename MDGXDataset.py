# -*- coding: utf-8 -*-
#   myplotspec_forcefield.MDGXDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages MDGX datasets
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class MDGXDataset(object):
    """
    Manages MDGX datasets
    """

    def __init__(self, infile, contains=None, verbose=0, **kwargs):
        """
        Initializes dataset.

        Arguments:
            infile (str): Path to text infile
        """
        from os.path import expandvars
        import pandas
        if verbose >= 1:
            print("Loading {0}".format(expandvars(infile)))
        data = pandas.read_csv(expandvars(infile), delim_whitespace=True,
          index_col=0, header=None, names=["topology", "restart",
          "QM absolute energy", "QM energy", "MM energy", "weight"])
        if contains is not None:
            import six
            sel_topologies = sel_restarts = None
            if isinstance(contains, six.string_types):
                sel_topologies = [contains]
            elif isinstance(contains, list):
                sel_topologies = contains
            elif isinstance(contains, dict):
                if "topology" in contains:
                    sel_topologies = contains["topology"]
                    if isinstance(sel_topologies, six.string_types):
                        sel_topologies = [sel_topologies]
                if "restart" in contains:
                    sel_restarts = contains["restart"]
                    if isinstance(sel_restarts, six.string_types):
                        sel_restarts = [sel_restarts]
            if sel_topologies is not None:
                if verbose >= 1:
                    print("Selecting topologies containing {0}".format(
                      str(sel_topologies)[1:-1]))
                    data = data[data["topology"].str.contains(
                      "({0})".format("|".join(sel_topologies)))]
            if sel_restarts is not None:
                if verbose >= 1:
                    print("Selecting restarts containing {0}".format(
                      str(sel_restarts)[1:-1]))
                    data = data[data["restart"].str.contains(
                      "({0})".format("|".join(sel_restarts)))]
        if verbose >= 1:
            print("Selected {0} topologies and {1} restarts".format(
              len(set(data["topology"])), len(set(data["restart"]))))
        if verbose >= 2:
            if sel_topologies is not None:
                print("Selected topologies: {0}".format(
                  str(sorted(set(data["topology"])))[1:-1].replace("'","")))
            if sel_restarts is not None:
                print("Selected restarts: {0}".format(
                  str(sorted(set(data["restart"])))[1:-1].replace("'","")))
        self.data = data
