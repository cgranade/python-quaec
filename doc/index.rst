..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.

.. QuaEC documentation master file, created by
   sphinx-quickstart on Fri Mar 23 12:44:13 2012.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

QuaEC: Quantum Error Correction Analysis in Python
==================================================

Contents
--------

.. toctree::
   :maxdepth: 2
   
   groups
   paulicollections
   bsf
   stab
   circuits
   solvers
   predicates
   utils
   exceptions
   bibliography

Introduction
============

QuaEC is a library for working with quantum error correction and
fault-tolerance. In particular, QuaEC provides support for maniuplating Pauli
and Clifford operators, as well as binary symplectic representations of each.

QuaEC is intended to provide easy, automated analysis of error-correcting protocols
based on stabilizer codes. for example, one can define a stabilizer code from a 
pre-existing library, and produce an object representing a circuit to encode data
into that code: 

>>> import qecc as q
>>> perfect_code=q.StabilizerCode.perfect_5q_code()
>>> print perfect_code.encoding_cliffords().next().circuit_decomposition()
    CNOT    1 0
    CNOT    3 0
    CNOT    4 0
    CNOT    2 1
    CNOT    3 2
    CNOT    4 2
    CNOT    4 3
    H       0
    H       1
    H       2
    H       3
    H       4
    CZ      0 1
    CZ      0 2
    CZ      1 2
    CZ      0 3
    CZ      1 4
    CZ      3 4
    H       0
    H       1
    H       2
    H       3
    H       4
    H       4
    SWAP    3 4
    H       4
    CNOT    0 4
    CNOT    0 3
    CNOT    0 2
    CNOT    0 1
    X       1
    X       2
    Z       3
    Y       4

Getting Started with QuaEC
==========================

Obtaining QuaEC
---------------

Currently, QuaEC is hosted `on GitHub <https://github.com/cgranade/python-quaec>`_.
The latest unstable version of QuaEC is available for
`download <https://github.com/cgranade/python-quaec/zipball/master>`_
as a ZIP there. Stable releases can be found on the
`downloads page`_, including installation packages for Windows and common Linux
distributions.

QuaEC is available via `PyPI`_ as well. To obtain it, run ``easy_install quaec``
at the terminal or in the Windows command line.

.. _downloads page: https://github.com/cgranade/python-quaec/downloads
.. _PyPI: http://pypi.python.org/pypi

Installation
------------

Once you have obtained QuaEC, installation is straightforward using the included
``setup.py`` script or the installation packages.

To use ``setup.py`` on Unix-like systems, run the following commands from the
command line::

    $ cd /path/to/quaec/
    $ sudo python setup.py install
    
To use ``setup.py`` on Windows, run ``cmd.exe``, then enter the following
commands::

    C:\> cd C:\path\to\quaec\
    C:\path\to\quaec\> python setup.py install
    
You may be prompted for permission by User Access Control, as the installer
attempts to install QuaEC into the system-wide packages directory.

Once QuaEC has been installed, it is made available as the :mod:`qecc` package:

>>> import qecc as q
>>> print q.Pauli('XYZ', phase=2)
i^2 XYZ

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

