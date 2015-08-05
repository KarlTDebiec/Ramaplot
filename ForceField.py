#!/usr/bin/python
# -*- coding: utf-8 -*-
#   ramaplot.ForceField.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Reads and represents force fields
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class ForceField(object):
    """
    """
    def write_hdf5(self, outfilename, verbose=1, debug=0, **kwargs):
        """
        """
        for p in self.parameters.keys():
            print(p)
            print(self.parameters(p))
