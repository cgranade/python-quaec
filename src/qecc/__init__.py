#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# __init__.py: Package definition for qecc.
##
# © 2012 Christopher E. Granade (cgranade@gmail.com) and
#     Ben Criger (bcriger@gmail.com).
# This file is a part of the QuaEC project.
# Licensed under the AGPL version 3.
##
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
##

__version__ = (1, 0, 1)

# All of the modules must be completely imported before we can start importing
# specific names, due to circular dependencies between the various modules.
import singletons as _sing
import PauliClass as _pc
import paulicollections as _pc
import CliffordClass as _cc
import bsf as _bsf
import exceptions as _xpts
import pred as _p
import circuit as _circ
import constraint_solvers as _cs
import stab as _stab

__modules = [_pc, _pc, _cc, _bsf, _xpts, _p, _circ, _cs, _stab]
# Note that we exclude _sing here to prevent changing the id of each singleton.

# Note that the utils module is not exposed, as we wish for that module to
# contain private functions.
#
# Similarly, bsf_decomp is hidden, exposed as a method in
# bsf.BinarySymplecticMatrix.

# So that reload(qecc) does what users expect, we need to reload each module.
map(reload, __modules)

# We now expose the particular names we want to expose, relying on the __all__
# variable in each module to define which names to expose.
from singletons import *
from PauliClass import *
from paulicollections import *
from CliffordClass import *
from bsf import *
from exceptions import *
from pred import *
from circuit import *
from constraint_solvers import *
from stab import *

__all__ = reduce(lambda a, b: a+b, map(lambda mod: mod.__all__, __modules)) + ['__version__']

