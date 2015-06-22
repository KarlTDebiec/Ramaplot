# -*- coding: utf-8 -*-
#   myplotspec_forcefield.NDRDDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages Neighbor-Dependent Ramachandran Distribution datasets.
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
################################### CLASSES ###################################
class NDRDDataset(object):
    """
    Manages Neighbor-Dependent Ramachandran Distribution datasets.

    Parses and organizes neighbor-dependent Ramachandran (backbone Φ/Ψ)
    distribution datasets, as published in:
      Ting, Daniel, Wang, Guoli, Shapovalov, Maxim, Mitra, Rajib,
      Jordan, Michael I, Dunbrack Jr. Roland L. Neighbor-Dependent
      Ramachandran Probability Distributions of Amino Acids Developed
      from a Hierarchical Dirichlet Process Model. PLoS Computational
      Biology. 2010. 6. e1000763.
    Data for use with this class may be obtained from
    <http://dunbrack.fccc.edu/ndrd>_
    """

    def __init__(self, infile, selection="ALA", loop_edges=True, max_fe=None,
        **kwargs):
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
                '{whitespace}' will be replaced with the appropriate
                expressions

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
              infile (str): Path to text input file, may contain
                environment variables
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
                     usecols=[3,4,5,6], names=["phi", "psi", "probability",
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
        x_centers = np.unique(dist["phi"])
        y_centers = np.unique(dist["psi"])
        free_energy = np.zeros((x_centers.size, y_centers.size),
                               np.float) * np.nan
        probability = np.zeros((x_centers.size, y_centers.size),
                               np.float) * np.nan
        for index, phi, psi, p, fe in dist.itertuples():
            x_index = np.where(x_centers == phi)[0][0]
            y_index = np.where(y_centers == psi)[0][0]
            free_energy[x_index, y_index] = fe
            probability[x_index, y_index] = p
        free_energy -= np.nanmin(free_energy)

        x_width = np.mean(x_centers[1:] - x_centers[:-1])
        y_width = np.mean(y_centers[1:] - y_centers[:-1])

        # Loop data to allow contour lines to be drawn to edges
        if loop_edges:
            x_centers = np.concatenate(([x_centers[0]  - x_width],
                                         x_centers,
                                        [x_centers[-1] + x_width]))
            y_centers = np.concatenate(([y_centers[0] -  y_width],
                                         y_centers,
                                        [y_centers[-1] + y_width]))
            temp = np.zeros((x_centers.size, y_centers.size)) * np.nan
            temp[1:-1,1:-1]  = free_energy
            temp[1:-1,-1]    = free_energy[:,0]
            temp[-1,1:-1]    = free_energy[0,:]
            temp[1:-1,0]     = free_energy[:,-1]
            temp[0,1:-1]     = free_energy[-1,:]
            temp[0,0]        = free_energy[-1,-1]
            temp[-1,-1]      = free_energy[0,0]
            temp[0,-1]       = free_energy[-1,0]
            temp[-1,0]       = free_energy[0,-1]
            free_energy = temp

        self.x_centers = x_centers
        self.y_centers = y_centers
        self.x_width = x_width
        self.y_width = y_width
        self.x_bins  = np.linspace(x_centers[0]  - x_width/2,
                                   x_centers[-1] + x_width/2,
                                   x_centers.size + 1)
        self.y_bins  = np.linspace(y_centers[0]  - y_width/2,
                                   y_centers[-1] + y_width/2,
                                   y_centers.size + 1)
        self.dist = free_energy

        # Prepare mask
        if max_fe is not None:
            self.mask = np.ma.masked_where(np.logical_and(
              free_energy <= max_fe,
              np.logical_not(np.isnan(free_energy))),
              np.ones_like(free_energy))
        else:
            self.mask = np.ma.masked_where(
              np.logical_not(np.isnan(free_energy)),
              np.ones_like(free_energy))
