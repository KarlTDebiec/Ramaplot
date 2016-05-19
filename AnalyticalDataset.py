# -*- coding: utf-8 -*-
#   ramaplot.AnalyticalDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages analytical Ramachandran plot datasets.

.. todo:
  - Use parmed? Should be able to support many force fields
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

    @staticmethod
    def get_cache_key(infile, phi=None, psi=None, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.
        """
        from os.path import expandvars

        return (AnalyticalDataset, expandvars(infile),
                AnalyticalDataset.process_term_arg(phi),
                AnalyticalDataset.process_term_arg(psi))

    @staticmethod
    def get_cache_message(cache_key):
        """
        Generates message to be used when reloading previously-loaded
        dataset.

        Arguments:
            cache_key (tuple): key with which dataset object is stored
              in dataset cache

        Returns:
            cache_message (str): message to be used when reloading
              previously-loaded dataset
        """
        return "previously loaded from '{0}'".format(cache_key[1])

    @staticmethod
    def load_parm(infile, dataset_cache=None, verbose=1, **kwargs):
        """
        """
        from .AmberForceField import AmberForceField

        if "dataset_cache" is not None:
            cache_key = AmberForceField.get_cache_key(parm=infile, **kwargs)
            if cache_key in dataset_cache:
                if verbose >= 1:
                    print(AmberForceField.get_cache_message(cache_key))
                return dataset_cache[cache_key]
            else:
                if verbose >= 1:
                    print("loading from '{0}'".format(infile))
                dataset_cache[cache_key] = AmberForceField(parm=infile,
                  verbose=verbose-1, **kwargs)
                return dataset_cache[cache_key]
        else:
            if verbose >= 1:
                print("loading from '{0}'".format(infile))
            return AmberForceField(parm=infile, verbose=verbose-1, **kwargs)

    @staticmethod
    def process_term_arg(terms=None):
        """
        Processes torsion term arguments

        Arguments:
          terms (str, list): torsion term(s) to be loaded from parm
            file

        Returns:
          out_terms (tuple): processed terms
        """
        import six

        out_terms = []

        # May be "C -N -CX-C"
        if terms is None:
            pass
        elif isinstance(terms, six.string_types):
            out_terms.append([terms, 0.0])
        elif isinstance(terms, list):
            # May be ["C -N -CX-C "]
            if (len(terms) == 1
            and isinstance(terms[0], six.string_types)):
                out_terms.append([terms[0], 0.0])

            # May be ["C -N -CX-C ", 120]
            elif (len(terms) == 2
            and   isinstance(terms[0], six.string_types)
            and  (isinstance(terms[1], float)
            or    isinstance(terms[1], int))):
                out_terms.append(terms)

            # May be [["C -N -CX-C ", 120], "C -TN-CX -C "]
            else:
                for in_term in terms:
                    if isinstance(in_term, six.string_types):
                        out_terms.append([in_term, 0.0])
                    elif isinstance(in_term, list):
                        if (len(in_term) == 1
                        and isinstance(in_term[0], six.string_types)):
                            out_terms.append([in_term[0], 0.0])
                        elif (len(in_term) == 2
                        and   isinstance(in_term[0], six.string_types)
                        and  (isinstance(in_term[1], float)
                        or    isinstance(in_term[1], int))):
                            out_terms.append(in_term)
                        else:
                            raise()
                    else:
                        raise()
        else:
            raise()
        return tuple(tuple(x) for x in out_terms)

    def __init__(self, infile, verbose=1, **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to Amber parm text file, may contain
            environment variables
          verbose (int): Level of verbose output
        """
        from os.path import expandvars
        from warnings import warn
        import six
        import pandas as pd
        import numpy as np
        from .AmberForceField import AmberForceField

        # Load or reload data
        infile = expandvars(infile)
        ff = self.load_parm(infile, verbose=verbose, **kwargs)

        # Initialize
        torsions = ff.parameters["dihedrals"]
        dist = np.zeros((360, 360))
        grid = np.linspace(-180, 180, 360)

        # Apply torsion terms
        for dim in ["phi", "psi"]:
            
            terms = AnalyticalDataset.process_term_arg(kwargs.get(dim))

            # Apply each term in selected dimension
            for term in terms:

                # Load offset if provided
                term, offset = term
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
                        warn("Term '{0:2}-{1:2}-{2:2}-{3:2}' ".format(
                          type_1, type_2, type_3, type_4) +
                          "not present in '{0}';".format(expandvars(infile)) +
                          "Term '{0:2}-{1:2}-{2:2}-{3:2}' ".format(
                          type_4, type_3, type_2, type_1) +
                          "is present and will be used")
                if verbose >= 2:
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
                    if dim == "phi":
                        dist += torsion_pe[:,np.newaxis]
                    else:
                        dist += torsion_pe
        dist -= np.min(dist)

        # Organize data
        self.x_centers = grid
        self.y_centers = grid
        self.dist = dist
        self.x_width = np.mean(grid[1:] - grid[:-1])
        self.y_width = np.mean(grid[1:] - grid[:-1])
        self.x_bins  = np.linspace(grid[0]  - self.x_width / 2,
                                   grid[-1] + self.x_width / 2,
                                   grid.size + 1)
        self.y_bins  = np.linspace(grid[0]  - self.y_width / 2,
                                   grid[-1] + self.y_width / 2,
                                   grid.size + 1)
