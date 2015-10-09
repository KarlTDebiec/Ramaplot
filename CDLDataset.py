# -*- coding: utf-8 -*-
#   ramaplot.CDLDataset.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Manages Conformation-Dependent Library datasets.

.. todo
    - Write type_error_text
    - Improve error checking
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
from .myplotspec.Dataset import Dataset
################################### CLASSES ###################################
class CDLDataset(Dataset):
    """
    Manages Conformation-dependent library datasets.

    Parses and organizes Conformation-Dependent Library datasets, as
    published in:
      Berkholz, Donald S., Shapovalov, Maxim V., Dunbrack Jr. Roland L.,
      Karplus, P. Andrew. Conformation Dependence of Backbone Geometry
      in Proteins. Structure. 2009. 17. 1316-1325.
    and:
      Berkholz, Donald S., Driggers, Camden M., Shapovalov, Maxim V., Dunbrack
      Jr., Roland L., Karplus, P. Andrew. Nonplanar Peptide Bonds in Proteins
      are Common and Conserved but not Biased Toward Active Sites. Proceedings
      of the National Acadedemy of the Sciences of the United States of
      America. 2012. 109. 449-453.
    Data for use with this class may be obtained from
    `<http://dunbrack.fccc.edu/nmhrcm>`_ and
    `<http://dunbrack.fccc.edu/omega>`_.
    """

    type_error_text = ("Selection may be a string or list; if string, "
      "may be the name of a class (e.g. 'NonPGIV_nonxpro', for which the "
      "N-CA-C (tau) angle distribution will be returned), or the name of a "
      "field (e.g. 'CN', for which the standard residue C-N bond distribution "
      "will be returned; if list, must be a pair of class and field (e.g. "
      "'['IleVal_xpro', 'BAC']', for which the isoleucine/valine CB-CA-C "
      "angle distribution will be returned.")

    @classmethod
    def get_cache_key(cls, infile, selection="NonPGIV_nonxpro",
        loop_edges=True, *args, **kwargs):
        """
        Generates tuple of arguments to be used as key for dataset
        cache.

        Arguments documentented under :func:`__init__`.
        """
        from os.path import expandvars

        return (cls, expandvars(infile), cls.process_selection_arg(selection),
          loop_edges)

    @staticmethod
    def process_selection_arg(selection):
        """
        Processes selection arguments

        Arguments:
          selection (str, list): selection to be loaded from file

        Returns:
          out_selection (tuple): processed selections

        .. todo:
          - Check that if dataset is All, All_nonxpro, or All_xpro,
            field is omega
        """
        import re
        import six

        def format_regex(regex):
            """
            Formats regular expression for parsing selection.

            Arguments:
              regex (str): Regular expression string, for which
                '{dataset}', '{field}', and '{whitespace}' will be
                replaced with the appropriate expressions

            Returns:
              regex (Pattern): Compiled regular expression
            """
            return re.compile(regex.format(
                dataset = "(All|All_nonxpro|All_xpro"
                          "|NonPGIV_nonxpro|IleVal_nonxpro"
                          "|Gly_nonxpro|Pro_nonxpro"
                          "|NonPGIV_xpro|IleVal_xpro"
                          "|Gly_xpro|Pro_xpro)",
                field = "(CNA|NAB|NAC|BAC|ACO|ACN|OCN"
                              "|CN|NA|AB|AC|CO|W)",
                whitespace = "\s+"), re.IGNORECASE)

        re_dataset = format_regex("^(?P<dataset>{dataset})$")
        re_field = format_regex("^(s|m)?(?P<field>{field})"
                           "({whitespace}(mean|avg|average|sd|std|stdev))?$")
        re_sd   = format_regex(
                    "(^s{field}$|^{field}{whitespace}(sd|std|stdev)$)")

        # Unpack list to string
        if isinstance(selection, list) and len(selection) == 1:
            selection = selection[0]
        if isinstance(selection, six.string_types):
            # Selection is a dataset, default to CNA field
            if re.match(re_dataset, selection):
                dataset, field = selection, "NAC"
            # Selection is a field, default to !PGVI>!P dataset
            elif re.match(re_field, selection):
                dataset, field = "NonPGIV_nonxpro", selection
            else:
                raise TypeError("Selection '{0}' not ".format(selection) +
                                "understood. " + CDLDataset.type_error_text)
        elif isinstance(selection, list) and len(selection) == 2:
            # Selection is a dataset, field pair
            if (re.match(re_dataset, selection[0])
            and re.match(re_field, selection[1])):
                dataset, field = selection
            else:
                raise TypeError("Selection '{0}' not ".format(selection) +
                                "understood. " + CDLDataset.type_error_text)
        else:
            raise TypeError("Selection '{0}' not ".format(selection) +
                            "understood. " + CDLDataset.type_error_text)

        match_sd = re.match(re_sd, field)
        match_field = re.match(re_field, field)
        if match_sd:
            field = "{0} sd".format(match_field.groupdict()["field"])
        else:
            field = "{0} mean".format(match_field.groupdict()["field"])

        return dataset, field

    @staticmethod
    def load_cdl_dataset(infile, selection, verbose=1, **kwargs):
        """
        Loads selected dataset from selected infile.

        Arguments:
          infile (str): Path to text input file
          selection (str): Start of lines containing desired dataset
          verbose (int): Level of verbose output

        Returns:
          dist (DataFrame): Selected dataset
        """
        from cStringIO import StringIO
        import pandas

        if verbose >= 1:
            print("loading '{0}' from '{1}'".format(selection, infile))

        s = StringIO()
        with open(infile) as f:
            omega = True if "omega" in f.readline() else False
            for line in f:
                if line.startswith(selection):
                    s.write(line)
        s.seek(0)

        if omega:
            return pandas.read_csv(s, delim_whitespace=True, header=None,
              usecols=range(1,7), names=["phi", "psi", "mode", "N",
              "W mean", "W sd", ])
        else:
            return pandas.read_csv(s, delim_whitespace=True, header=None,
              usecols=range(1,29), names=["phi", "psi", "mode", "N",
              "CNA mean", "CNA sd", "NAB mean", "NAB sd", "NAC mean", "NAC sd",
              "BAC mean", "BAC sd", "ACO mean", "ACO sd", "ACN mean", "ACN sd",
              "OCN mean", "OCN sd", "CN mean",  "CN sd",  "NA mean",  "NA sd",
              "AB mean",  "AB sd",  "AC mean",  "AC sd",  "CO mean",  "CO sd"])

    def __init__(self, infile, selection="NonPGIV_nonxpro", loop_edges=True,
        **kwargs):
        """
        Initializes dataset.

        Arguments:
          infile (str): Path to text input file, may contain environment
            variables
          selection (str, list): Selected dataset and field
          loop_edges (bool): Copy edges to enable plotting to edge of
            plot
          kwargs (dict): additional keyword arguments
        """
        from os.path import expandvars
        import pandas
        import numpy as np

        infile = expandvars(infile)
        selection = self.process_selection_arg(selection)

        dataset = self.load_cdl_dataset(infile, selection[0])

        # Organize data
        x_centers = np.unique(dataset["phi"])
        y_centers = np.unique(dataset["psi"])
        x_width = np.mean(x_centers[1:] - x_centers[:-1])
        y_width = np.mean(y_centers[1:] - y_centers[:-1])
        dist = np.zeros((x_centers.size, y_centers.size), np.float) * np.nan
        mask = np.zeros((x_centers.size, y_centers.size), np.bool)
        for index, row in dataset.iterrows():
            x_index = np.where(x_centers == row["phi"])[0][0]
            y_index = np.where(y_centers == row["psi"])[0][0]
            dist[x_index, y_index] = row[selection[1]]
            if row["mode"] == "B":
                mask[x_index, y_index] = True

        # Loop distribution to allow contour lines to be drawn to edges
        if loop_edges:
            x_centers = np.concatenate(([x_centers[0] - x_width],
                                         x_centers,
                                        [x_centers[-1] + x_width]))
            y_centers = np.concatenate(([y_centers[0] - y_width],
                                         y_centers,
                                        [y_centers[-1] + y_width]))
            temp = np.zeros((x_centers.size, y_centers.size), np.float)
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
