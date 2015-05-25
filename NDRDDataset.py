# -*- coding: utf-8 -*-
#   myplotspec_forcefield.NDRDDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages Neighbor-dependent Ramachandran distribution datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class NDRDDataset(object):
    """
    Manages Neighbor-dependent Ramachandran distribution datasets.

    Parses and organizes neighbor-dependent Ramachandran (backbone Φ/Ψ)
    distributions for protein loops, as published in:
      Ting, Daniel, Wang, Guoli, Shapovalov, Maxim, Mitra, Rajib,
      Jordan, Michael I, Dunbrack Jr. Roland L. Neighbor-Dependent
      Ramachandran Probability Distributions of Amino Acids Developed
      from a Hierarchical Dirichlet Process Model. PLoS Computational
      Biology. 2010. 6. e1000763.
    Data for use with this class may be obtained from
    <http://dunbrack.fccc.edu/ndrd>_
    """

    def __init__(self, infile, selection="ALA", **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          selection (str, list): Selected distribution; if string, may
            be an amino acid (e.g. 'CYS', for which the (ALA)-CYS-(ALA)
            distribution will be returned), or a specific distribution
            (e.g. 'CYS left ALA', for which the (ALA)-CYS distribution
            will be returned); if list, may be two amino acids (e.g.
            ['ALA', 'CYS'], for which the (ALA)-CYS distribution will
            be returned), or three amino acids (e.g. ['ALA', 'CYS',
            'ASP'], for which the (ALA)-CYS-(ASP) distribution will be
            returned)
          kwargs (dict): additional keyword arguments
        """
        from cStringIO import StringIO
        from os.path import expandvars
        import re
        import pandas
        import numpy as np
        import six

        def format_regex(regex):
            """
            Formats regular expression for parsing selection.

            Arguments:
              regex (str): Regular expression string, for which
                '{neighbor}', '{central}', '{direction}', and
                '{whitespace}' with the appropriate expressions

            Returns:
              regex (Pattern): Compiled regular expression
            """
            return re.compile(regex.format(
                neighbor   = "(ALA|ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS"
                             "|ILE|LEU|LYS|MET|PHE|PRO|CPR|SER|THR|TRP"
                             "|TYR|VAL|ALL)",
                central    = "(ALA|ALA|ARG|ASN|ASP|CYS|GLN|GLU|GLY|HIS"
                             "|ILE|LEU|LYS|MET|PHE|PRO|CPR|SER|THR|TRP"
                             "|TYR|VAL)",
                direction  = "(left|right)",
                whitespace = "\s+"))

        def load_distribution(infile, selection, verbose=1, **kwargs):
            """
            Loads selected distribution from selected infile.

            Arguments:
              infile (str): Path to text input file, may contain environment
                variables
              selection (str): Start of lines containing desired
                distribution
              verbose (int): Enable verbose output

            Returns:
              dist (DataFrame): Selected distribution
            """
            if verbose > 0:
                print("loading '{0}' from '{1}'".format(selection, infile))
            s = StringIO()

            with open(expandvars(infile)) as f:
                for line in f:
                    if line.startswith(selection):
                        s.write(line)

            s.seek(0)
            dist = pandas.read_csv(s, delim_whitespace=True, header=None,
                     usecols=[3,4,5,6], names=["x", "y", "probability",
                     "free energy"])
            return dist

        def calculate_triplet(dist_1, dist_2, dist_3, **kwargs):
            """
            Calculates distribution for an amino acid triplet.

            Formula used included with NDRD data:
              To calculate probabilies for triplets (center,left,right),
              use:
              log p*(phi,psi | C,L,R) = log p(phi,psi |C,L) +
                                        log p(phi,psi |C,R) -
                                        log p(phi,psi |C,R=ALL)
              Once log p*(phi,psi | C,L,R) is calculated, calculate
              p*(phi,psi |C,L,R) = exp(log(p*(phi,psi | C,L,R)))
              Then sum them up for each Ramachandran map, and normalize
              the probabilities by dividing by the sum:
              p(phi,psi, | C,L,R) = p*(phi,psi | C,L,R) / sum

            Arguments:
              dist_1 (DataFrame): (XX1)>XX2 distribution
              dist_2 (DataFrame): XX2>(XX3) distribution
              dist_3 (DataFrame): XX2>(ALL) distribution

            Returns:
              dist (DataFrame): RES1>RES2>RES3 distribution for RES2
            """
            dist = dist_1.copy()
            dist["free energy"] += dist_2["free energy"]
            dist["free energy"] -= dist_3["free energy"]
            dist["probability"]  = np.exp(dist["free energy"])
            dist["probability"] /= np.sum(dist["probability"])

            return dist

        type_error_text = ("Selection '{0}' ".format(selection) +
                           "not understood. Selection may be a string "
                           "or list; if string, may be an amino acid "
                           "(e.g. 'CYS', for which the (ALA)-CYS-(ALA) "
                           "distribution will be returned), or a "
                           "specific distribution (e.g. 'CYS left "
                           "ALA', for which the (ALA)-CYS distribution "
                           "will be returned); if list, may be two "
                           "amino acids (e.g. ['ALA', 'CYS'], for "
                           "which the (ALA)-CYS distribution will be "
                           "returned), or three amino acids (e.g. "
                           "['ALA', 'CYS', 'ASP'], for which the "
                           "(ALA)-CYS-(ASP) distribution will be "
                           "returned)")
        sel_central  = format_regex("^(?P<central>{central})$")
        sel_neighbor = format_regex("^(?P<neighbor>{neighbor})$")
        sel_dist     = format_regex("^(?P<first>{central}){whitespace}"
                                     "(?P<direction>{direction}){whitespace}"
                                     "(?P<second>{neighbor})$")

        # Unpack list to string
        if isinstance(selection, list) and len(selection) == 1:
            selection = selection[0]
        if isinstance(selection, six.string_types):
            # Selection is a single residue, return ALA>XXX>ALA
            if re.match(sel_central, selection):
                fields = re.match(sel_central, selection).groupdict()
                sel_1   = "{0} left  ALA".format(fields["central"])
                sel_2   = "{0} right ALA".format(fields["central"])
                sel_3   = "{0} right ALL".format(fields["central"])
            # Selection is a single distribution
            elif re.match(sel_dist, selection):
                fields = re.match(sel_dist, selection).groupdict()
                sel_1 = "{0} {1:5} {2}".format(fields["first"],
                          fields["direction"], fields["second"])
                sel_2 = None
                sel_3 = None
            else:
                raise TypeError(type_error_text)
        elif isinstance(selection, list):
            if len(selection) == 2:
                # Selection is a sequence of two resides, return XX1>XX2
                if (re.match(sel_neighbor, selection[0])
                and re.match(sel_central, selection[1])):
                    fields_1 = re.match(sel_neighbor,selection[0]).groupdict()
                    fields_2 = re.match(sel_central,selection[1]).groupdict()
                    res_1 = fields_1["neighbor"]
                    res_2 = fields_2["central"]
                    sel_1 = "{0} left  {1}".format(res_2, res_1)
                    sel_2 = None
                    sel_3 = None
                else:
                    raise TypeError(type_error_text)
            elif len(selection) == 3:
                # Selection is a sequence of three residues
                if (re.match(sel_neighbor, selection[0])
                and re.match(sel_central, selection[1])
                and re.match(sel_neighbor, selection[2])):
                    fields_1 = re.match(sel_neighbor,selection[0]).groupdict()
                    fields_2 = re.match(sel_central,selection[1]).groupdict()
                    fields_3 = re.match(sel_neighbor,selection[2]).groupdict()
                    res_1 = fields_1["neighbor"]
                    res_2 = fields_2["central"]
                    res_3 = fields_3["neighbor"]
                    sel_1 = "{0} left  {1}".format(res_2, res_1)
                    sel_2 = "{0} right {1}".format(res_2, res_3)
                    sel_3 = "{0} right ALL".format(res_2)
                else:
                    raise TypeError(type_error_text)
        if sel_1 is not None and sel_2 is None and sel_3 is None:
            dist = load_distribution(infile, sel_1)
        elif sel_1 is not None and sel_2 is not None and sel_3 is not None:
            dist_1 = load_distribution(infile, sel_1)
            dist_2 = load_distribution(infile, sel_2)
            dist_3 = load_distribution(infile, sel_3)
            dist = calculate_triplet(dist_1, dist_2, dist_3)
        else:
            raise TypeError(type_error_text)

        # Organize data
        self.data = dist
        self.x_centers = np.unique(self.data["x"])
        self.y_centers = np.unique(self.data["y"])
        self.x_width = np.mean(self.x_centers[1:] - self.x_centers[:-1])
        self.y_width = np.mean(self.y_centers[1:] - self.y_centers[:-1])
        self.x_bins  = np.linspace(self.x_centers[0]  - self.x_width/2,
                                   self.x_centers[-1] + self.x_width/2,
                                   self.x_centers.size + 1)
        self.y_bins  = np.linspace(self.y_centers[0]  - self.y_width/2,
                                   self.y_centers[-1] + self.y_width/2,
                                   self.y_centers.size + 1)
        self.free_energy = np.zeros((self.x_centers.size, self.y_centers.size),
                                    np.float)
        self.probability = np.zeros((self.x_centers.size, self.y_centers.size),
                                    np.float)
        self.free_energy[:] = np.nan
        self.probability[:] = np.nan
        for index, x, y, probability, free_energy in self.data.itertuples():
            y_index = np.where(self.x_centers == x)[0][0]
            x_index = np.where(self.y_centers == y)[0][0]
            self.free_energy[x_index, y_index] = free_energy
            self.probability[x_index, y_index] = probability
        self.free_energy -= np.nanmin(self.free_energy)

        # Format contour settings
        self.contour_x_centers = np.zeros(self.x_centers.size + 1, np.float)
        self.contour_x_centers[0:-1] = self.x_centers
        self.contour_x_centers[-1] = self.contour_x_centers[-2] + self.x_width
        self.contour_y_centers = np.zeros(self.y_centers.size + 1, np.float)
        self.contour_y_centers[0:-1] = self.y_centers
        self.contour_y_centers[-1] = self.contour_y_centers[-2] + self.y_width
        self.contour_free_energy = np.zeros((self.x_centers.size + 1,
                                             self.y_centers.size + 1),
                                             np.float)
        self.contour_free_energy[:] = np.nan
        self.contour_free_energy[0:-1,0:-1] = self.free_energy
        self.contour_free_energy[0:-1,-1]   = self.free_energy[:,0]
        self.contour_free_energy[-1,0:-1]   = self.free_energy[0,:]
        self.contour_free_energy[-1,-1]     = self.free_energy[0,0]

        # Format heatmap settings
        self.imshow_free_energy = self.contour_free_energy
        self.imshow_extent = [-182.5, 182.5, -182.5, 182.5]