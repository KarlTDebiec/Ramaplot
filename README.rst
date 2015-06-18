Introduction
============

MYPlotSpec_ForceField is a Python package used to design matplotlib-based
figures of biomolecular force field components using the simple text format
YAML.

Features
--------

ForceFieldFigureManager
~~~~~~~~~~~~~~~~~~~~~~~

Currently plots only torsion terms.

MDGXFigureManager
~~~~~~~~~~~~~~~~~

Plots debug output from the Amber MDGX program, illustrating the accuracy which
which an MM force field is able to reproduce the relative energies of an
ensemble of comformations obtained from QM calculations.

RamachandranFigureManager
~~~~~~~~~~~~~~~~~~~~~~~~~

Plots Ramachandran distributions. Currently supports input from Umbrella
Sampling simulations analyzed using the `weighted histogram analysis method
<http://membrane.urmc.rochester.edu/content/wham>`_, `neighbor-dependent
Ramachanran distribution datasets <http://dunbrack.fccc.edu/ndrd>`_, and
images.

Dependencies
------------

MYPlotSpec_ForceField supports Python 2.7 and 3.4, and requires the following
packages:

- matplotlib
- numpy
- six
- yaml

MYPlotSpec_ForceField has been tested with Anaconda python 2.2.0 on Arch Linux,
OSX Yosemite, and Windows 8.1.

Installation
------------

Put in your ``$PYTHONPATH``::

    export PYTHONPATH=/path/to/my/python/modules:$PYTHONPATH

where ``/path/to/my/python/modules`` contains ``myplotspec_forcefield``.

Authorship
----------

MYPlotSpec_ForceField is developed by Karl T. Debiec, a graduate student at the
University of Pittsburgh advised by Professors Lillian T. Chong and Angela M.
Gronenborn.

License
-------

Released under a 3-clause BSD license.
