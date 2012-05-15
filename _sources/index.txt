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

Welcome to QuaEC's documentation!
=================================

Contents
--------

.. toctree::
   :maxdepth: 2
   
   groups
   bsf
   circuits
   predicates
   utils

Introduction
============

QuaEC is a library for working with quantum error correction and
fault-tolerance. In particular, QuaEC provides support for maniuplating Pauli
and Clifford operators, as well as binary symplectic representations of each.


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

In the future, QuaEC will be made available via `PyPI`_ as well. At that point,
QuaEC will be obtainable using ``easy_install``.

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
>>> print q.Pauli('XYZ', ph=2)
i^2 XYZ

Next Steps
----------

Until we have a proper introduction written, a good start would be to look at
the :class:`qecc.Pauli` class.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

