# -*- coding: utf-8 -*-
#   myplotspec_forcefield.ImageDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages Ramachandran plot image datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class ImageDataset(object):
    """
    Manages Ramachandran plot image datasets.
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
        from matplotlib.image import imread
        import pandas
        import numpy as np

        # Load data
        if verbose > 0:
            print("loading from '{0}'".format(infile))
        image = np.rot90(imread(expandvars(infile)), 3)
        image = 0.2989*image[:,:,0] + 0.5870*image[:,:,1] + 0.1140*image[:,:,2]
        self.free_energy = image

        # Organize data
        self.x_centers = np.linspace(-180, 180, image.shape[0])
        self.y_centers = np.linspace(-180, 180, image.shape[1])
        self.x_width = np.mean(self.x_centers[1:] - self.x_centers[:-1])
        self.y_width = np.mean(self.y_centers[1:] - self.y_centers[:-1])
        self.x_bins  = np.linspace(self.x_centers[0]  - self.x_width / 2,
                                   self.x_centers[-1] + self.x_width / 2,
                                   self.x_centers.size + 1)
        self.y_bins  = np.linspace(self.y_centers[0]  - self.y_width / 2,
                                   self.y_centers[-1] + self.y_width / 2,
                                   self.y_centers.size + 1)
