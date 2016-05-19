# -*- coding: utf-8 -*-
#   ramaplot.NDRDDataset.py
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

    type_error_text = ("Selection may be a  string or list; if string, "
      "may be an  amino acid (e.g. 'CYS', for which the "
      "(ALA)-CYS-(ALA) distribution will be  returned), or a specific "
      "distribution  (e.g. 'CYS left ALA', for which the  (ALA)-CYS "
      "distribution will be  returned); if list, may be two amino "
      "acids (e.g. ['ALA', 'CYS'], for which  the (ALA)-CYS "
      "distribution will be  returned), or three amino acids (e.g. "
      "['ALA', 'CYS', 'ASP'], for which the  (ALA)-CYS-(ASP) "
      "distribution will be returned.")

    @classmethod
    def get_cache_key(cls, infile, selection="ALA", zkey="free energy",
        mask_cutoff=None,
        calc_populations=False, plot_populations=False,
        *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documented under :func:`__init__`.
        """
        from os.path import expandvars

        return (cls, expandvars(infile), cls.process_selection_arg(selection),
          zkey, mask_cutoff, calc_populations, plot_populations)

    @staticmethod
    def get_cache_message(cache_key):
        """
        Generates message to be used when reloading previously-loaded
        dataset.

        Arguments:
            cache_key (tuple): key with which dataset object is stored
              in dataset cache

        Returns:
            cache_message (str): Message to be used when reloading
              previously-loaded dataset
        """
        return "previously loaded '{0}' from '{1}'".format(
          str(cache_key[2]).replace("u'","").replace("'",""), cache_key[1])

    @staticmethod
    def process_selection_arg(selection):
        """
        Processes selection arguments

        Arguments:
          selection (str, list): Selection(s) to be loaded from file

        Returns:
          out_selection (tuple): Processed selection(s)
        """
        import re
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

        sel_central  = format_regex("^(?P<central>{central})$")
        sel_neighbor = format_regex("^(?P<neighbor>{neighbor})$")
        sel_dataset  = format_regex("^(?P<first>{central}){whitespace}"
                                     "(?P<direction>{direction}){whitespace}"
                                     "(?P<second>{neighbor})$")

        # Unpack list to string
        if isinstance(selection, list) and len(selection) == 1:
            selection = selection[0]
        if isinstance(selection, six.string_types):
            # Selection is a single residue, return ALA>XAA>ALA
            if re.match(sel_central, selection):
                fields = re.match(sel_central, selection).groupdict()
                out_selection = ("{0} left  ALA".format(fields["central"]),
                                 "{0} right ALA".format(fields["central"]),
                                 "{0} right ALL".format(fields["central"]))
            # Selection is a single distribution
            elif re.match(sel_dataset, selection):
                fields = re.match(sel_dataset, selection).groupdict()
                out_selection = ("{0} {1:5} {2}".format(fields["first"],
                  fields["direction"], fields["second"]),)
            else:
                raise TypeError("Selection '{0}' not ".format(selection) +
                                "understood. " + NDRDDataset.type_error_text)
        elif isinstance(selection, list):
            if len(selection) == 2:
                # Selection is a sequence of two resides, return XX1>XX2
                if (re.match(sel_neighbor, selection[0])
                and re.match(sel_central, selection[1])):
                    fields_1 = re.match(sel_neighbor,selection[0]).groupdict()
                    fields_2 = re.match(sel_central,selection[1]).groupdict()
                    out_selection = ("{0} left  {1}".format(
                      fields_1["neighbor"], fields_2["central"]),)
                else:
                    raise TypeError("Selection '{0}' not ".format(selection) +
                                    "understood. " + self.type_error_text)
            elif len(selection) == 3:
                # Selection is a sequence of three residues
                if (re.match(sel_neighbor, selection[0])
                and re.match(sel_central, selection[1])
                and re.match(sel_neighbor, selection[2])):
                    fields_1 = re.match(sel_neighbor,selection[0]).groupdict()
                    fields_2 = re.match(sel_central,selection[1]).groupdict()
                    fields_3 = re.match(sel_neighbor,selection[2]).groupdict()
                    out_selection = (
                      "{0} left  {1}".format( fields_2["central"],
                      fields_1["neighbor"]),
                      "{0} right {1}".format(fields_2["central"],
                      fields_3["neighbor"]),
                      "{0} right ALL".format(fields_2["central"]))
                else:
                    raise TypeError("Selection '{0}' not ".format(selection) +
                                    "understood. " + self.type_error_text)

        return out_selection

    @staticmethod
    def load_distribution(infile, selection, verbose=1, **kwargs):
        """
        Loads selected distribution from selected infile.

        Arguments:
          infile (str): Path to text input file
          selection (str): Start of lines containing desired dataset
          verbose (int): Level of verbose output

        Returns:
          dataset (DataFrame): Selected dataset
        """
        from six import StringIO
        import pandas as pd

        if verbose >= 1:
            print("loading '{0}' from '{1}'".format(selection, infile))

        s = StringIO()
        with open(infile) as open_infile:
            for line in open_infile:
                if line.startswith(selection):
                    s.write(line)
        s.seek(0)

        dataset = pd.read_csv(s, delim_whitespace=True, header=None,
          usecols=[3,4,5,6], names=["phi", "psi", "probability",
          "free energy"])

        return dataset

    @staticmethod
    def calculate_triplet(dataset_1, dataset_2, dataset_3, **kwargs):
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
          dataset_1 (DataFrame): (XX1)>XX2 distribution
          dataset_2 (DataFrame): XX2>(XX3) distribution
          dataset_3 (DataFrame): XX2>(ALL) distribution

        Returns:
          dataset (DataFrame): RES1>RES2>RES3 distribution for RES2
        """
        from copy import copy
        import numpy as np

        dataset = copy(dataset_1)
        dataset["free energy"] += dataset_2["free energy"]
        dataset["free energy"] -= dataset_3["free energy"]
        dataset["probability"]  = np.exp(-1 * dataset["free energy"])
        dataset["probability"] /= np.sum(dataset["probability"])

        return dataset

    def __init__(self, infile, selection="ALA", zkey="free energy",
        loop_edges=True, mask_cutoff=None,
        calc_populations=False, plot_populations=False,
        verbose=1, **kwargs):
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
          kwargs (dict): Additional keyword arguments
        """
        from os.path import expandvars
        import pandas as pd
        import numpy as np

        # Manage arguments
        if zkey not in ["free energy", "probability"]:
            raise ValueError("Argument 'zkey' does not support provided " +
              "value '{0}', must be 'free energy or ".format(zkey) +
              "probability.")

        infile = expandvars(infile)
        selection = self.process_selection_arg(selection)

        if len(selection) == 1:
            dataset = self.load_distribution(infile, selection[0],
              verbose=verbose)
        elif len(selection) == 3:
            dataset = self.calculate_triplet(
                         self.load_distribution(infile, selection[0],
                           verbose=verbose),
                         self.load_distribution(infile, selection[1],
                           verbose=verbose),
                         self.load_distribution(infile, selection[2],
                           verbose=verbose))
        else:
            raise TypeError("Selection '{0}' not ".format(selection) +
                            "understood. " + self.type_error_text)

        # Organize data
        x_centers = np.unique(dataset["phi"])
        y_centers = np.unique(dataset["psi"])
        x_width = np.mean(x_centers[1:] - x_centers[:-1])
        y_width = np.mean(y_centers[1:] - y_centers[:-1])
        free_energy = np.zeros((x_centers.size, y_centers.size),
                               np.float) * np.nan
        probability = np.zeros((x_centers.size, y_centers.size),
                               np.float) * np.nan
        for index, phi, psi, p, fe in dataset.itertuples():
            x_index = np.where(x_centers == phi)[0][0]
            y_index = np.where(y_centers == psi)[0][0]
            free_energy[x_index, y_index] = fe
            probability[x_index, y_index] = p
        free_energy -= np.nanmin(free_energy)
        probability /= np.nansum(probability)

        self.x_centers = x_centers
        self.y_centers = y_centers
        self.x_width = x_width
        self.y_width = y_width
        self.dist = free_energy
        self.x_bins  = np.linspace(x_centers[0]  - x_width / 2,
                                   x_centers[-1] + x_width / 2,
                                   x_centers.size + 1)
        self.y_bins  = np.linspace(y_centers[0]  - y_width / 2,
                                   y_centers[-1] + y_width / 2,
                                   y_centers.size + 1)

        # Store distribution in instance variable and create mask
        if zkey == "free energy":
            self.dist = free_energy
            if mask_cutoff is not None:
                self.mask = np.ma.masked_where(np.logical_and(
                  free_energy <= mask_cutoff,
                  np.logical_not(np.isnan(free_energy))),
                  np.ones_like(free_energy))
            else:
                self.mask = np.ma.masked_where(
                  np.logical_not(np.isnan(free_energy)),
                  np.ones_like(free_energy))
        elif zkey == "probability":
            self.dist = probability
            if mask_cutoff is not None:
                self.mask = np.ma.masked_where(np.logical_and(
                  probability >= mask_cutoff,
                  np.logical_not(np.isnan(probability))),
                  np.ones_like(probability))
            else:
                self.mask = np.ma.masked_where(
                  np.logical_not(np.isnan(probability)),
                  np.ones_like(probability))

#        # Prepare mask
#        if mask_cutoff is not None:
#            self.mask = np.ma.masked_where(np.logical_and(
#              free_energy <= mask_cutoff,
#              np.logical_not(np.isnan(free_energy))),
#              np.ones_like(free_energy))
#        else:
#            self.mask = np.ma.masked_where(
#              np.logical_not(np.isnan(free_energy)),
#              np.ones_like(free_energy))

        # Calculate state populations
        if calc_populations:
            states = kwargs.get("states",  [
              ("β",       -151,  151),
              ("PPII",     -66,  140),
              ("ξ",       -145,   55),
              ("γ'",       -81,   65),
              ("α",        -70,  -25),
              ("$L_α$",     55,   45),
              ("γ",         73,  -35),
              ("PPII'",     56, -124),
              ("plateau", -100, -130)])
            state_radius = kwargs.get("state_radius", 45)

            distances = np.zeros((len(states), len(x_centers), len(y_centers)))
            xs = []
            ys = []
            for i, (state, x, y) in enumerate(states):
                xs += [x]
                ys += [y]
                # There must be a better way to do this, but this works
                for j, xc in enumerate(x_centers):
                    for k, yc in enumerate(y_centers):
                        dx = (xc - x)
                        if dx <= -180 or dx >= 180:
                            dx = 360 - dx
                        else:
                            dx = dx
                        dy = (yc - y)
                        if dy <= -180 or dy >= 180:
                            dy = 360 - dy
                        else:
                            dy = dy
                        distances[i,j,k] = np.sqrt(dx**2 + dy**2)
            assignments = np.argmin(distances, axis=0)
            assignments[np.min(distances, axis=0) >= state_radius] = \
              len(states) + 1

            index, state_populations = [], []
            for i, (state, x, y) in enumerate(states):
                index += [state]
                state_populations += [(x, y,
                  np.nansum(probability[assignments==i]))]
            state_populations = pd.DataFrame(state_populations, index=index,
              columns=["Φ center", "Ψ center", "population"])
            self.state_populations = state_populations

            if verbose >= 1:
                print(state_populations)

            if plot_populations:
                self.dist = assignments
                self.mask = np.ma.masked_where(
                  np.logical_not(assignments == len(states) + 1),
                  np.ones_like(assignments))
                self.x = np.array(xs)
                self.y = np.array(ys)
                label, label_kw = [], []
                from .myplotspec import multi_get_copy
                default_label_kw = multi_get_copy(["default_label_kw",
                  "label_kw"], kwargs, {})
                for index, row in state_populations.iterrows():
                    label += ["{0}\n{1:2d}%".format(index,
                     int(row["population"]*100))]
                    label_kw += [default_label_kw.copy()]
                    label_kw[-1]["x"] = row["Φ center"]
                    label_kw[-1]["y"] = row["Ψ center"]
                self.label = label
                self.label_kw = label_kw
