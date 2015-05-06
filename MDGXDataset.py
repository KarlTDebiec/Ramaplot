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

    def __init__(self, infile, **kwargs):
        """
        Initializes dataset.

        Arguments:
            infile (str): Path to text infile
        """
        from os.path import expandvars
        import pandas
        self.data = pandas.read_csv(expandvars(infile), delim_whitespace=True,
          index_col=0, header=None, names=["topology", "restart",
          "QM absolute energy", "QM energy", "MM energy", "weight"])
