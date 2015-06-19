# -*- coding: utf-8 -*-
#   myplotspec_forcefield.AnalyticalDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages analytical Ramachandran plot datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class AnalyticalDataset(object):
    """
    Manages analytical Ramachandran plot datasets.

      Vn
    ------- * (1 + cos(periodicity * x + phase))
    divider
    """

    def __init__(self, infile, verbose=1, debug=1, **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          verbose (bool): Enable verbose output
          debug (bool): Enable debug output
        """
        from os.path import expandvars
        from warnings import warn
        import six
        import pandas as pd
        import numpy as np
        from .AmberForceField import AmberForceField

        # Load data and initialize
        if verbose > 0:
            print("loading from '{0}'".format(infile))
        ff = AmberForceField(parm=expandvars(infile), verbose=verbose-1,
               debug=debug, **kwargs)
        torsions = ff.parameters["dihedrals"]
        pe = np.zeros((360, 360))
        grid = np.linspace(-180, 180, 360)

        # Apply torsion terms
        for dim in ["x", "y"]:
            
            terms = kwargs.get(dim, [])
            if (isinstance(terms, six.string_types)
            or  isinstance(terms, dict)):
                terms = [terms]

            # Apply each term in selected dimension
            for term in terms:

                # Load offset if provided
                if isinstance(term, dict):
                    offset = float(term.values()[0])
                    term = term.keys()[0]
                elif isinstance(term, six.string_types):
                    offset = 0.0
                type_1, type_2, type_3, type_4 = [t.strip()
                                                   for t in term.split("-")]
                dim_torsions = torsions[(torsions["type_1"] == type_1) &
                                        (torsions["type_2"] == type_2) &
                                        (torsions["type_3"] == type_3) &
                                        (torsions["type_4"] == type_4)]
                if dim_torsions.size == 0:
                    dim_torsions = torsions[(torsions["type_4"] == type_1) &
                                            (torsions["type_3"] == type_2) &
                                            (torsions["type_2"] == type_3) &
                                            (torsions["type_1"] == type_4)]
                    if dim_torsions.size == 0:
                        raise Exception(
                          "Term '{0:2}-{1:2}-{2:2}-{3:2}' ".format(
                            type_1, type_2, type_3, type_4) +
                            "not present in '{0}'".format(expandvars(infile)))
                    elif verbose >= 2:
                        warn( "Term '{0:2}-{1:2}-{2:2}-{3:2}' ".format(
                          type_1, type_2, type_3, type_4) +
                          "not present in '{0}';".format(expandvars(infile)) +
                          "Term '{0:2}-{1:2}-{2:2}-{3:2}' ".format(
                          type_4, type_3, type_2, type_1) +
                          "is present and will be used")
                if verbose >= 1:
                    print(dim_torsions[["type_1", "type_2", "type_3", "type_4",
                      "divider", "barrier", "phase", "periodicity"]])
                for index, torsion in dim_torsions.iterrows():
                    divider     = float(torsion["divider"])
                    barrier     = float(torsion["barrier"])
                    phase       = float(torsion["phase"])
                    periodicity = float(torsion["periodicity"])
                    torsion_pe  = barrier / divider * (1 +
                                    np.cos(
                                      np.deg2rad(
                                        np.abs(periodicity) * grid +
                                          phase + offset)))
                    if dim == "x":  pe += torsion_pe[:,np.newaxis]
                    else:           pe += torsion_pe
        pe -= np.min(pe)

        # Organize data
        self.free_energy = pe
        self.x_centers = grid
        self.y_centers = grid
        self.x_width = np.mean(self.x_centers[1:] - self.x_centers[:-1])
        self.y_width = np.mean(self.y_centers[1:] - self.y_centers[:-1])
        self.x_bins  = np.linspace(self.x_centers[0]  - self.x_width / 2,
                                   self.x_centers[-1] + self.x_width / 2,
                                   self.x_centers.size + 1)
        self.y_bins  = np.linspace(self.y_centers[0]  - self.y_width / 2,
                                   self.y_centers[-1] + self.y_width / 2,
                                   self.y_centers.size + 1)
