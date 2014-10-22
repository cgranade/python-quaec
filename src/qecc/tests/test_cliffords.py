#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# CliffordClass.py: Implementation of qecc.Clifford and related utility
#     functions.
##
# © 2014 Christopher E. Granade (cgranade@gmail.com) and
#        Ben Criger (bcriger@gmail.com).
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

## IMPORTS ###################################################################

import qecc as q

from nose.tools import raises

## TESTS #####################################################################

@raises(IndexError)
def check_clifford_factories_index_errors(factory_fn, n_args):
	args = (4,) * n_args
	factory_fn(1, *args)

def test_clifford_factories_index_errors():
	for factory_fn, n_args in [
			(q.cnot, 2),
			(q.hadamard, 1),
			(q.phase, 1),
			(q.swap, 2),
			(q.cz, 2)
	]:
		yield check_clifford_factories_index_errors, factory_fn, n_args