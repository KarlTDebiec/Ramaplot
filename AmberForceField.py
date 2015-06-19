#!/usr/bin/python
# -*- coding: utf-8 -*-
#   myplotspec_forcefield.AmberForceField.py
#
#   Copyright (C) 2015 Karl T Debiec
#   All rights reserved.
#
#   This software may be modified and distributed under the terms of the
#   BSD license. See the LICENSE file for details.
"""
Reads and represents AMBER-format force fields
"""
################################### MODULES ###################################
from __future__ import absolute_import,division,print_function,unicode_literals
import re
import pandas as pd
if __name__ == "__main__":
    __package__ = str("myplotspec_forcefield")
    import myplotspec_forcefield
from myplotspec_forcefield.ForceField import ForceField
################################### CLASSES ###################################
class AmberForceField(ForceField):
    """
    Represents, reads, and writes AMBER-format force fields

    Warning: the only clear difference between mass and vdw are the
      two spaces in front, should probably implement some
      intelligent warning if the first number is >3
    Warning: need setting to pull out water
    Warning: Not all masses have second number (polarizability?)
    """

    def __init__(self, parm=None, **kwargs):
        """
        Separate init functions to load from hdf5, input files, etc.
        """

        self.parameter_regex   = dict(
          blank           = self.amber_regex("^\s*$"),
          mass            = self.amber_regex("^(?P<type>{t}){w}"
                                              "(?P<mass>{f})"
                                              "(?P<polarizability>{w}{f}|{w})"
                                              "(?P<note>.*$)"),
          atomlist        = self.amber_regex("^({t}{w})*$"),
          bond            = self.amber_regex("^(?P<type_1>{t})-"
                                              "(?P<type_2>{t}){w}"
                                              "(?P<force_constant>{f}){w}"
                                              "(?P<length>{f}){w}"
                                              "(?P<note>.*$)"),
          angle           = self.amber_regex("^(?P<type_1>{t})-"
                                              "(?P<type_2>{t})-"
                                              "(?P<type_3>{t}){w}"
                                              "(?P<force_constant>{f}){w}"
                                              "(?P<angle>{f}){w}"
                                              "(?P<note>.*$)"),
          dihedral        = self.amber_regex("^(?P<type_1>{t})-"
                                              "(?P<type_2>{t})-"
                                              "(?P<type_3>{t})-"
                                              "(?P<type_4>{t}){w}"
                                              "(?P<divider>{i}){w}"
                                              "(?P<barrier>{sf}){w}"
                                              "(?P<phase>{sf}){w}"
                                              "(?P<periodicity>{sf}){w}"
                                              "(?P<note>.*$)"),
          improper        = self.amber_regex("^(?P<type_1>{t})-"
                                              "(?P<type_2>{t})-"
                                              "(?P<type_3>{t})-"
                                              "(?P<type_4>{t}){w}"
                                              "(?P<barrier>{sf}){w}"
                                              "(?P<phase>{sf}){w}"
                                              "(?P<periodicity>{sf}){w}"
                                              "(?P<note>.*$)"),
          hbond        = self.amber_regex("^{w}(?P<type_1>{t}){w}"
                                              "(?P<type_2>{t}){w}"
                                              "(?P<A>{f}){w}"
                                              "(?P<B>{f}){w}"
                                              "(?P<ASOLN>{f}){w}"
                                              "(?P<note>.*$)"),
          vdw_format = self.amber_regex("^.+{w}(?P<vdw_format>SK|RE|AC).*$"),
          vdw          = self.amber_regex("^{w}(?P<type>{t}){w}"
                                              "(?P<radius>{f}){w}"
                                              "(?P<well_depth>{f}){w}"
                                              "(?P<note>.*$)"),
          ljedit_title = self.amber_regex("^LJEDIT$"),
          ljedit       = self.amber_regex("^{w}(?P<type_1>{t}){w}"
                                              "(?P<type_2>{t}){w}"
                                              "(?P<radius_1>{f}){w}"
                                              "(?P<well_depth_1>{f}){w}"
                                              "(?P<radius_2>{f}){w}"
                                              "(?P<well_depth_2>{f}){w}"
                                              "(?P<note>.*$)"),
          end          = self.amber_regex("^END$"))

        self.residue_regex = dict(
          blank      = self.amber_regex("^\s*$"),
          atoms      = self.amber_regex("^\s*\"(?P<name>{a})\"{w}"
                                            "\"(?P<type>{t})\"{w}"
                                              "(?P<type_index>{i}){w}"
                                              "(?P<residue_index>{i}){w}"
                                              "(?P<flags>{i}){w}"
                                              "(?P<atom_index>{i}){w}"
                                              "(?P<element>{i}){w}"
                                              "(?P<charge>{sf}){w}"
                                              "(?P<note>.*$)"),
          atom_edits = self.amber_regex("^\s*\"(?P<name>{a})\"{w}"
                                            "\"(?P<type>{t})\"{w}"
                                              "(?P<type_index>{i}){w}"
                                              "(?P<element>{si}{w})"
                                              "(?P<charge>{sf}{w})"
                                              "(?P<note>.*$)"),
          box          = self.amber_regex("^\s*(?P<box>{sf}){w}"
                                              "(?P<note>.*$)"),
          res_seq      = self.amber_regex("^\s*(?P<childsequence>{i}){w}"
                                              "(?P<note>.*$)"),
          res_connect  = self.amber_regex("^\s*(?P<connect>{i}){w}"
                                              "(?P<note>.*$)"),
          bonds        = self.amber_regex("^\s*(?P<atom_index_1>{i}){w}"
                                              "(?P<atom_index_2>{t}){w}"
                                              "(?P<flag>{i}){w}"
                                              "(?P<note>.*$)"),
          hierarchy  = self.amber_regex("^\s*\"(?P<above_type>U|R|A)\"{w}"
                                              "(?P<above_index>{i}){w}"
                                            "\"(?P<below_type>U|R|A)\"{w}"
                                              "(?P<below_index>{i}){w}"
                                              "(?P<note>.*$)"),
          name       = self.amber_regex("^\s*\"(?P<name>{r})\""
                                              "(?P<note>.*$)"),
          coordinates  = self.amber_regex("^\s*(?P<x>{sfe}){w}"
                                              "(?P<y>{sfe}){w}"
                                              "(?P<z>{sfe}){w}"
                                              "(?P<note>.*$)"),
          res_connect2 = self.amber_regex("^\s*(?P<atom_index_1>{i}){w}"
                                              "(?P<atom_index_2>{i}){w}"
                                              "(?P<atom_index_3>{i}){w}"
                                              "(?P<atom_index_4>{i}){w}"
                                              "(?P<atom_index_5>{i}){w}"
                                              "(?P<atom_index_6>{i}){w}"
                                              "(?P<note>.*$)"),
          residues   = self.amber_regex("^\s*\"(?P<name>{r})\"{w}"
                                              "(?P<residue_index>{i}){w}"
                                              "(?P<child_atom_index>{i}){w}"
                                              "(?P<start_atom_index>{i}){w}"
                                            "\"(?P<residue_type>p|n|w|\?)\"{w}"
                                              "(?P<note>.*$)"),
          pdb_seq      = self.amber_regex("^\s*(?P<residue_index>{i}){w}"
                                              "(?P<note>.*$)"),
          solventcap   = self.amber_regex("^\s*(?P<solventcap>{sf}){w}"
                                              "(?P<note>.*$)"),
          velocities   = self.amber_regex("^\s*(?P<x>{sfe}){w}"
                                              "(?P<y>{sfe}){w}"
                                              "(?P<z>{sfe}){w}"
                                              "(?P<note>.*$)"))

        self.parameters = parameters = {}
        parameters["mass_types"] = pd.DataFrame(columns=
          ["type", "mass", "polarizability", "note"])
        parameters["hydrophobic_types"] = pd.DataFrame(columns=
          ["type"])
        parameters["bonds"] = pd.DataFrame(columns=
          ["type_1", "type_2", "force_constant", "length", "note"])
        parameters["angles"] = pd.DataFrame(columns=
          ["type_1", "type_2", "type_3", "force_constant", "angle", "note"])
        parameters["dihedrals"] = pd.DataFrame(columns=
          ["type_1", "type_2", "type_3", "type_4", "divider", "barrier",
           "phase", "periodicity", "note"])
        parameters["impropers"] = pd.DataFrame(columns=
          ["type_1", "type_2", "type_3", "type_4", "barrier", "phase",
           "periodicity", "note"])
        parameters["hbonds"] = pd.DataFrame(columns=
          ["type_1", "type_2", "A", "B", "ASOLN"])
        parameters["vdw_eq_types"] = pd.DataFrame()
        parameters["vdw_types"] = pd.DataFrame(columns=
          ["type", "radius", "well_depth", "note"])
        parameters["ljedits"] = pd.DataFrame(columns=
          ["type_1", "type_2", "radius_1", "well_depth_1", "radius_2",
           "well_depth_2"])
        self.residues = {}

        if parm is not None:
            self.read_parm(parm, **kwargs)

    @staticmethod
    def amber_regex(regex, title=False):
        """
        Prepares regex for matching AMBER fields

        Arguments:
          regex (string): regular expression

        Returns:
          (string): regular expression
        """
        if title:
            regex = "^!entry\.(?P<residue_name>{r})\.unit\." + regex + "{w}.*$"
        return re.compile(regex.format(
          r   = "[\w][\w][\w][\w]?",            # Residue
          a   = "[\w][\w]?[\w]?[\w]?",          # Atom name
          t   = "[\w][\w \*]?",                 # Atom type
          i   = "\d+",                          # Integer
          si  = "[-]?\d+",                      # Signed Integer
          f   = "\d+\.?\d*?",                   # Float
          sf  = "[-]?\d+\.?\d*?",               # Signed float
          sfe = "[-]?\d+\.?\d*?[E]?[-]?\d*?",   # Signed float in E notation
          w   = "\s+"))                         # Whitespace

    def strip_dict(self, dictionary):
        """
        Strips each string in a dict, and deletes if empty

        Arguements:
          dictionary (dict): dictionary to strip

        Returns:
          (dict): dictionary with each element stripped
        """
        for key, value in dictionary.items():
            value = value.strip()
            if value == "":
                del dictionary[key]
            else:
                dictionary[key] = value
        return dictionary

    def read_parm(self, infilename, verbose=1, debug=0, **kwargs):
        """
        Reads a parm file

        Arguments:
          infilename (str): Path to input parm file
          verbose (int): Enable verbose output
          debug (int): Enable debug output
          kwargs (dict): Additional keyword arguments
        """
        if verbose >= 1:
            print("READING PARM: {0}".format(infilename))

        # References to instance variables and methods to make below
        #   code readable
        strip_dict  = self.strip_dict
        amber_regex = self.amber_regex
        mass_types        = self.parameters["mass_types"]
        hydrophobic_types = self.parameters["hydrophobic_types"]
        bonds             = self.parameters["bonds"]
        angles            = self.parameters["angles"]
        dihedrals         = self.parameters["dihedrals"]
        impropers         = self.parameters["impropers"]
        hbonds            = self.parameters["hbonds"]
        vdw_eq_types      = self.parameters["vdw_eq_types"]
        vdw_types         = self.parameters["vdw_types"]
        ljedits           = self.parameters["ljedits"]
        re_blank        = self.parameter_regex["blank"]
        re_mass         = self.parameter_regex["mass"]
        re_atomlist     = self.parameter_regex["atomlist"]
        re_bond         = self.parameter_regex["bond"]
        re_angle        = self.parameter_regex["angle"]
        re_dihedral     = self.parameter_regex["dihedral"]
        re_improper     = self.parameter_regex["improper"]
        re_hbond        = self.parameter_regex["hbond"]
        re_vdw_format   = self.parameter_regex["vdw_format"]
        re_vdw          = self.parameter_regex["vdw"]
        re_ljedit_title = self.parameter_regex["ljedit_title"]
        re_ljedit       = self.parameter_regex["ljedit"]
        re_end          = self.parameter_regex["end"]

        section = 1

        with open(infilename, "r") as infile:
            line = infile.readline()
            while line:
                # BLANK
                if re.match(re_blank, line):
                    if verbose >= 1:
                        print("BLANK          |{0}".format(line.strip()))
                # 1: TITLE
                elif section <= 1 and not re.match(re_mass, line):
                    if verbose >= 1:
                        print("TITLE          |{0}".format(line.strip()))
                # 2: MASS
                elif section <= 2 and re.match(re_mass, line):
                    section = 2
                    if verbose >= 1:
                        print("MASS           |{0}".format(line.strip()))
                    fields = strip_dict(re.match(re_mass, line).groupdict())
                    mass_types = mass_types.append(fields, ignore_index=True)
                # 3: HYDROPHIC (list of types)
                elif section <= 3 and re.match(re_atomlist, line):
                    section = 3
                    if verbose >= 1:
                        print("HYDROPHOBIC    |{0}".format(line.rstrip()))
                    fields = [{"type": v} for v in
                      amber_regex("{t}").findall(line)]
                    hydrophobic_types = hydrophobic_types.append(fields,
                                        ignore_index=True)
                # 4: BOND
                elif section <= 4 and re.match(re_bond, line):
                    section = 4
                    if verbose >= 1:
                        print("BOND           |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_bond, line).groupdict())
                    bonds = bonds.append(fields, ignore_index=True)
                # 5: ANGLE
                elif section <= 5 and re.match(re_angle, line):
                    section = 5
                    if verbose >= 1:
                        print("ANGLE          |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_angle, line).groupdict())
                    angles = angles.append(fields, ignore_index=True)
                # 6: DIHEDRAL
                elif section <= 6 and re.match(re_dihedral, line):
                    section = 6
                    if verbose >= 1:
                        print("DIHEDRAL       |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_dihedral,
                             line).groupdict())
                    dihedrals = dihedrals.append(fields, ignore_index=True)
                # 7: IMPROPER
                elif section <= 7 and re.match(re_improper, line):
                    section = 7
                    if verbose >= 1:
                        print("IMPROPER       |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_improper,
                             line).groupdict())
                    impropers = impropers.append(fields, ignore_index=True)
                # 8: HBOND
                elif section <= 8 and re.match(re_hbond, line):
                    section = 8
                    if verbose >= 1:
                        print("HBOND          |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_hbond, line).groupdict())
                    hbonds = hbonds.append(fields, ignore_index=True)
                # 9: VDW (equivalent types)
                elif section <= 9 and re.match(re_atomlist, line):
                    section = 9
                    if verbose >= 1:
                        print("VDW EQUIVALENT |{0}".format(line.rstrip()))
                    fields = [{"type_{0}".format(i): v for i, v in
                      enumerate(re.compile(amber_regex("{t}")).findall(line))}]
                    vdw_eq_types = vdw_eq_types.append(fields,
                                           ignore_index=True)
                # 10: VDW (format)
                elif section <= 10.3 and re.match(re_vdw_format, line):
                    if verbose >= 1:
                        print("VDW FORMAT     |{0}".format(line.rstrip()))
                # 10.2: VDW (radius and well depth)
                elif section <= 10.2 and re.match(re_vdw, line):
                    section = 10.2
                    if verbose >= 1:
                        print("VDW            |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_vdw, line).groupdict())
                    vdw_types = vdw_types.append(fields, ignore_index=True)
                # 11: LJEDIT (title)
                elif (section <= 11 and re.match(re_ljedit_title, line)):
                    section = 11
                    if verbose >= 1:
                        print("LJEDIT         |{0}".format(line.rstrip()))
                # 11.1: LJEDIT (atom types, radii, and well depth)
                elif (section <= 11.1 and re.match(re_ljedit, line)):
                    section = 11.1
                    if verbose >= 1:
                        print("LJEDIT         |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_ljedit, line).groupdict())
                    ljedits = ljedits.append(fields, ignore_index=True)
                # END
                elif re.match(re_end, line):
                    if verbose >= 1:
                        print("END            |{0}".format(line.rstrip()))
                    break
                # NO MATCH
                else:
                    if verbose >= 1:
                        print("NOMATCH        |{0}".format(line.rstrip()))
                line = infile.readline()
        if debug >= 1:
            print(mass_types)
            print(hydrophobic_types)
            print(bonds)
            print(angles)
            print(dihedrals)
            print(impropers)
            print(hbonds)
            print(vdw_eq_types)
            print(vdw_types)
            print(ljedits)
        self.parameters["mass_types"]        = mass_types
        self.parameters["hydrophobic_types"] = hydrophobic_types
        self.parameters["bonds"]             = bonds
        self.parameters["angles"]            = angles
        self.parameters["dihedrals"]         = dihedrals
        self.parameters["impropers"]         = impropers
        self.parameters["hbonds"]            = hbonds
        self.parameters["vdw_eq_types"]      = vdw_eq_types
        self.parameters["vdw_types"]         = vdw_types
        self.parameters["ljedits"]           = ljedits

    def read_frcmod(self, infilename, verbose=1, debug=0, **kwargs):
        """
        Arguments:
          infilename (str): Path to input lib file
          verbose (int): Enable verbose output
          debug (int): Enable debug output
          kwargs (dict): Additional keyword arguments
        """
        if verbose >= 1:
            print("READING FRCMOD: {0}".format(infilename))

        are = self.amber_regex
        strip_dict = self.strip_dict
        re_blank = are("^\s*$")

        section = 1
        with open(infilename, "r") as infile:
            line = infile.readline()
            while line:
                # BLANK
                if re.match(re_blank, line):
                    if verbose >= 1:
                        print("BLANK          |{0}".format(line.strip()))
                # 1: TITLE
                elif section <= 1 and not re.match(re_mass, line):
                    if verbose >= 1:
                        print("TITLE          |{0}".format(line.strip()))
                # 2: MASS
                elif section <= 2 and re.match(re_mass, line):
                    section = 2
                    if verbose >= 1:
                        print("MASS           |{0}".format(line.strip()))
                    fields = strip_dict(re.match(re_mass, line).groupdict())
                    mass_types = mass_types.append(fields, ignore_index=True)
                # 3: HYDROPHIC (list of types)
                elif section <= 3 and re.match(re_atomlist, line):
                    section = 3
                    if verbose >= 1:
                        print("HYDROPHOBIC    |{0}".format(line.rstrip()))
                    fields = [{"type": v} for v in are("{t}").findall(line)]
                    hydrophobic_types = hydrophobic_types.append(fields,
                                        ignore_index=True)
                # 4: BOND
                elif section <= 4 and re.match(re_bond, line):
                    section = 4
                    if verbose >= 1:
                        print("BOND           |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_bond, line).groupdict())
                    bonds = bonds.append(fields, ignore_index=True)
                # 5: ANGLE
                elif section <= 5 and re.match(re_angle, line):
                    section = 5
                    if verbose >= 1:
                        print("ANGLE          |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_angle, line).groupdict())
                    angles = angles.append(fields, ignore_index=True)
                # 6: DIHEDRAL
                elif section <= 6 and re.match(re_dihedral, line):
                    section = 6
                    if verbose >= 1:
                        print("DIHEDRAL       |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_dihedral,
                             line).groupdict())
                    dihedrals = dihedrals.append(fields, ignore_index=True)
                # 7: IMPROPER
                elif section <= 7 and re.match(re_improper, line):
                    section = 7
                    if verbose >= 1:
                        print("IMPROPER       |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_improper,
                             line).groupdict())
                    impropers = impropers.append(fields, ignore_index=True)
                # 8: HBOND
                elif section <= 8 and re.match(re_hbond, line):
                    section = 8
                    if verbose >= 1:
                        print("HBOND          |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_hbond, line).groupdict())
                    hbonds = hbonds.append(fields, ignore_index=True)
                # 9: VDW (equivalent types)
                elif section <= 9 and re.match(re_atomlist, line):
                    section = 9
                    if verbose >= 1:
                        print("VDW EQUIVALENT |{0}".format(line.rstrip()))
                    fields = [{"type_{0}".format(i): v for i, v in
                             enumerate(re.compile(are("{t}")).findall(line))}]
                    vdw_eq_types = vdw_eq_types.append(fields,
                      ignore_index=True)
                # 10: VDW (format)
                elif section <= 10.3 and re.match(re_vdw_format, line):
                    if verbose >= 1:
                        print("VDW FORMAT     |{0}".format(line.rstrip()))
                # 10.2: VDW (radius and well depth)
                elif section <= 10.2 and re.match(re_vdw, line):
                    section = 10.2
                    if verbose >= 1:
                        print("VDW            |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_vdw, line).groupdict())
                    vdw_types = vdw_types.append(fields, ignore_index=True)
                # 11: LJEDIT (title)
                elif (section <= 11 and re.match(re_ljedit_title, line)):
                    section = 11
                    if verbose >= 1:
                        print("LJEDIT         |{0}".format(line.rstrip()))
                # 11.1: LJEDIT (atom types, radii, and well depth)
                elif (section <= 11.1 and re.match(re_ljedit, line)):
                    section = 11.1
                    if verbose >= 1:
                        print("LJEDIT         |{0}".format(line.rstrip()))
                    fields = strip_dict(re.match(re_ljedit, line).groupdict())
                    ljedits = ljedits.append(fields, ignore_index=True)
                # END
                elif re.match(re_end, line):
                    if verbose >= 1:
                        print("END            |{0}".format(line.rstrip()))
                    break
                # NO MATCH
                else:
                    if verbose >= 1:
                        print("NOMATCH        |{0}".format(line.rstrip()))
                line = infile.readline()

    def read_lib(self, infilename, verbose=1, debug=0, **kwargs):
        """
        Arguments:
          infilename (str): Path to input lib file
          verbose (int): Enable verbose output
          debug (int): Enable debug output
          kwargs (dict): Additional keyword arguments
        """
        if verbose >= 1:
            print("READING LIB: {0}".format(infilename))

        stripd          = self.strip_dict
        re_blank        = self.residue_regex["blank"]
        re_atoms        = self.residue_regex["atoms"]
        re_atom_edits   = self.residue_regex["atom_edits"]
        re_box          = self.residue_regex["box"]
        re_res_seq      = self.residue_regex["res_seq"]
        re_res_connect  = self.residue_regex["res_connect"]
        re_bonds        = self.residue_regex["bonds"]
        re_hierarchy    = self.residue_regex["hierarchy"]
        re_name         = self.residue_regex["name"]
        re_coordinates  = self.residue_regex["coordinates"]
        re_res_connect2 = self.residue_regex["res_connect2"]
        re_residues     = self.residue_regex["residues"]
        re_pdb_seq      = self.residue_regex["pdb_seq"]
        re_solventcap   = self.residue_regex["solventcap"]
        re_velocities   = self.residue_regex["velocities"]

        # Regular expressions for titles
        re_t_atoms          = self.amber_regex("atoms",          title=True)
        re_t_atom_edits     = self.amber_regex("atomspertinfo",  title=True)
        re_t_box            = self.amber_regex("boundbox",       title=True)
        re_t_res_seq        = self.amber_regex("childsequence",  title=True)
        re_t_res_connect    = self.amber_regex("connect",        title=True)
        re_t_bonds          = self.amber_regex("connectivity",   title=True)
        re_t_hierarchy      = self.amber_regex("hierarchy",      title=True)
        re_t_name           = self.amber_regex("name",           title=True)
        re_t_coordinates    = self.amber_regex("positions",      title=True)
        re_t_res_connect2   = self.amber_regex("residueconnect", title=True)
        re_t_residues       = self.amber_regex("residues",       title=True)
        re_t_pdb_seq        = self.amber_regex("residuesPdbSequenceNumber",
                                                                 title=True)
        re_t_solventcap     = self.amber_regex("solventcap",     title=True)
        re_t_velocities     = self.amber_regex("velocities",     title=True)

        # Regular expressions for contents
        section = 0
        residue = None

        with open(infilename, "r") as infile:
            line = infile.readline()
            while line:
                # BLANK
                if re.match(re_blank, line):
                    if verbose >= 1:
                        print("BLANK          |{0}".format(line.strip()))
                # 1: ATOMS
                elif re.match(re_t_atoms, line):
                    if verbose >= 1:
                        print("ATOMS          |{0}".format(line.strip()))
                    section = 1
                    fields = stripd(re.match(re_t_atoms, line).groupdict())
                    residue = self.residues[fields["residue_name"]] = {}
                    residue["atoms"] = pd.DataFrame(columns=
                      ["name", "type", "type_index", "residue_index", "flags",
                       "atom_index", "element", "charge", "note"])
                elif section == 1 and re.match(re_atoms, line):
                    if verbose >= 1:
                        print("ATOMS          |{0}".format(line.strip()))
                    fields = stripd(re.match(re_atoms, line).groupdict())
                    residue["atoms"] = residue["atoms"].append(
                      fields, ignore_index=True)
                # 2: ATOMSPERTINFO
                elif re.match(re_t_atom_edits, line):
                    if verbose >= 1:
                        print("ATOMSPERTINFO  |{0}".format(line.strip()))
                    section = 2
                    residue["atom_edits"] = pd.DataFrame(columns=
                      ["name", "type", "type_index", "element", "charge",
                       "note"])
                elif section == 2 and re.match(re_atom_edits, line):
                    if verbose >= 1:
                        print("ATOMSPERTINFO  |{0}".format(line.strip()))
                    fields = stripd(re.match(re_atom_edits, line).groupdict())
                    residue["atom_edits"] = residue["atom_edits"].append(
                      fields, ignore_index=True)
                # 3: BOUNDBOX
                elif re.match(re_t_box, line):
                    if verbose >= 1:
                        print("BOUNDBOX       |{0}".format(line.strip()))
                    section = 3
                    box_keys = ["box", "angle", "x_length", "y_length",
                                "z_length"]
                    box_items = []
                elif section == 3 and re.match(re_box, line):
                    if verbose >= 1:
                        print("BOUNDBOX       |{0}".format(line.strip()))
                    fields = stripd(re.match(re_box, line).groupdict())
                    box_items.append(
                      (box_keys.pop(0), [fields["box"]]))
                    if len(box_keys) == 0:
                        residue["box"] = pd.DataFrame.from_items(box_items)
                # 4: CHILDSEQUENCE
                elif re.match(re_t_res_seq, line):
                    if verbose >= 1:
                        print("CHILDSEQUENCE  |{0}".format(line.strip()))
                    section = 4
                    residue["res_seq"] = pd.DataFrame(columns=
                      ["childsequence", "note"])
                elif section == 4 and re.match(re_res_seq, line):
                    if verbose >= 1:
                        print("CHILDSEQUENCE  |{0}".format(line.strip()))
                    fields = stripd(re.match(re_res_seq, line).groupdict())
                    residue["res_seq"] = residue["res_seq"].append(
                      fields, ignore_index=True)
                # 5: CONNECT
                elif re.match(re_t_res_connect, line):
                    if verbose >= 1:
                        print("CONNECT        |{0}".format(line.strip()))
                    section = 5
                    connect_keys = [
                      "connect_atom_index_1", "connect_atom_index_2", "note"]
                    connect_items = []
                elif section == 5 and re.match(re_res_connect, line):
                    if verbose >= 1:
                        print("CONNECT        |{0}".format(line.strip()))
                    fields = stripd(re.match(re_res_connect, line).groupdict())
                    connect_items.append(
                      (connect_keys.pop(0), [fields["connect"]]))
                    if len(connect_keys) == 0:
                        residue["res_connect"] = pd.DataFrame.from_items(
                          connect_items)
                # 6: CONNECTIVITY
                elif re.match(re_t_bonds, line):
                    if verbose >= 1:
                        print("CONNECTIVITY   |{0}".format(line.strip()))
                    section = 6
                    residue["bonds"] = pd.DataFrame(columns=
                      ["atom_index_1", "atom_index_2", "flag", "note"])
                elif section == 6 and re.match(re_bonds, line):
                    if verbose >= 1:
                        print("CONNECTIVITY   |{0}".format(line.strip()))
                    fields = stripd(re.match(re_bonds,
                             line).groupdict())
                    residue["bonds"] = residue["bonds"].append(
                      fields, ignore_index=True)
                # 7: HIERARCHY
                elif re.match(re_t_hierarchy, line):
                    if verbose >= 1:
                        print("HIERARCHY      |{0}".format(line.strip()))
                    section = 7
                    residue["hierarchy"] = pd.DataFrame(columns=
                      ["above_type", "above_index", "below_type",
                       "below_index", "note"])
                elif section == 7 and re.match(re_hierarchy, line):
                    if verbose >= 1:
                        print("HIERARCHY      |{0}".format(line.strip()))
                    fields = stripd(re.match(re_hierarchy,
                             line).groupdict())
                    residue["hierarchy"] = residue["hierarchy"].append(
                      fields, ignore_index=True)
                # 8: NAME
                elif re.match(re_t_name, line):
                    if verbose >= 1:
                        print("NAME           |{0}".format(line.strip()))
                    section = 8
                    residue["name"] = pd.DataFrame(columns=
                      ["childsequence", "note"])
                elif section == 8 and re.match(re_name, line):
                    if verbose >= 1:
                        print("NAME           |{0}".format(line.strip()))
                    fields = stripd(re.match(re_name, line).groupdict())
                    residue["name"] = residue["name"].append(
                      fields, ignore_index=True)
                # 9: POSITIONS
                elif re.match(re_t_coordinates, line):
                    if verbose >= 1:
                        print("POSITIONS      |{0}".format(line.strip()))
                    section = 9
                    residue["coordinates"] = pd.DataFrame(columns=
                      ["x", "y", "z", "note"])
                elif section == 9 and re.match(re_coordinates, line):
                    if verbose >= 1:
                        print("POSITIONS      |{0}".format(line.strip()))
                    fields = stripd(re.match(re_coordinates,
                             line).groupdict())
                    residue["coordinates"] = residue["coordinates"].append(
                      fields, ignore_index=True)
                # 10: RESIDUECONNECT
                elif re.match(re_t_res_connect2, line):
                    if verbose >= 1:
                        print("RESIDUECONNECT |{0}".format(line.strip()))
                    section = 10
                    residue["res_connect2"] = pd.DataFrame(columns=
                      ["atom_index_1", "atom_index_2", "atom_index_3",
                       "atom_index_4", "atom_index_5", "atom_index_6", "note"])
                elif section == 10 and re.match(re_res_connect2, line):
                    if verbose >= 1:
                        print("RESIDUECONNECT |{0}".format(line.strip()))
                    fields = stripd(re.match(re_res_connect2,
                             line).groupdict())
                    residue["res_connect2"] = residue["res_connect2"].append(
                      fields, ignore_index=True)
                # 11: RESIDUES
                elif re.match(re_t_residues, line):
                    if verbose >= 1:
                        print("RESIDUES       |{0}".format(line.strip()))
                    section = 11
                    residue["residues"] = pd.DataFrame(columns=
                      ["name", "residue_index", "child_atom_index",
                       "start_atom_index", "residue_type", "note"])
                elif re.match(re_residues, line):
                    if verbose >= 1:
                        print("RESIDUES       |{0}".format(line.strip()))
                    fields = stripd(re.match(re_residues,
                             line).groupdict())
                    residue["residues"] = residue["residues"].append(
                      fields, ignore_index=True)
                # 12: RESIDUESPDBSEQUENCENUMBER
                elif re.match(re_t_pdb_seq, line):
                    if verbose >= 1:
                        print("PDBSEQUENCENUM |{0}".format(line.strip()))
                    section = 12
                    residue["pdb_seq"] = pd.DataFrame(columns=
                      ["residue_index", "note"])
                elif section == 12 and re.match(re_pdb_seq, line):
                    if verbose >= 1:
                        print("PDBSEQUENCENUM |{0}".format(line.strip()))
                    fields = stripd(re.match(re_pdb_seq, line).groupdict())
                    residue["pdb_seq"] = residue["pdb_seq"].append(
                      fields, ignore_index=True)
                # 13: SOLVENTCAP
                elif re.match(re_t_solventcap, line):
                    if verbose >= 1:
                        print("SOLVENTCAP     |{0}".format(line.strip()))
                    section = 13
                    solventcap_keys = ["solventcap", "angle", "x_length",
                                       "y_length", "z_length"]
                    solventcap_temp = []
                elif section == 13 and re.match(re_solventcap, line):
                    if verbose >= 1:
                        print("SOLVENTCAP     |{0}".format(line.strip()))
                    fields = stripd(re.match(re_solventcap, line).groupdict())
                    solventcap_temp.append(
                      (solventcap_keys.pop(0), [fields["solventcap"]]))
                    if len(solventcap_keys) == 0:
                        residue["solventcap"] = pd.DataFrame.from_items(
                          solventcap_temp)
                # 14: VELOCITIES
                elif re.match(re_t_velocities, line):
                    if verbose >= 1:
                        print("VELOCITIES     |{0}".format(line.strip()))
                    section = 14
                    residue["velocities"] = pd.DataFrame(columns=
                      ["x", "y", "z", "note"])
                elif section == 14 and re.match(re_velocities, line):
                    if verbose >= 1:
                        print("VELOCITIES     |{0}".format(line.strip()))
                    fields = stripd(re.match(re_velocities,
                             line).groupdict())
                    residue["velocities"] = residue["velocities"].append(
                      fields, ignore_index=True)
                # NO MATCH
                else:
                    if verbose >= 1:
                        print("NOMATCH        |{0}".format(line.rstrip()))
                line = infile.readline()
        for name in sorted(self.residues):
            residue = self.residues[name]
            print()
            print(name)
            fields = ["atoms", "atom_edits", "box", "childsequence",
                      "connect", "bonds", "hierarchy", "name",
                      "coordinates", "residueconnect", "residues",
                      "pdbindex", "solventcap", "velocities"]
            for field in fields:
                if field in residue:
                    print(field)
                    print(residue[field])
