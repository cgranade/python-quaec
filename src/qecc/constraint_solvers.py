#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# constraint_solvers.py: Solvers for various Pauli and Clifford constraint
#     problems.
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

## ALL ##

__all__ = [
    'solve_commutation_constraints'
]

## IMPORTS ##

from PauliClass import Pauli, elem_gens, com, from_generators
from paulicollections import PauliList
from pred import AllPredicate
from itertools import ifilter

## FUNCTIONS ##

def solve_commutation_constraints(
        commutation_constraints=[],
        anticommutation_constraints=[],
        search_in_gens=None
    ):
    """
    """
    
    # Start by checking that the arguments make sense.
    if len(commutation_constraints) == 0 and len(anticommutation_constraints) == 0:
        raise ValueError("At least one constraint must be specified.")
        
    # Next, normalize our arguments to be PauliLists, so that we can obtain
    # centralizers easily.
    if not isinstance(commutation_constraints, PauliList):
        commutation_constraints = PauliList(commutation_constraints)
    if not isinstance(anticommutation_constraints, PauliList):
        # This is probably not necessary, strictly speaking, but it keeps me
        # slightly more sane to have both constraints represented by the same
        # sequence type.
        anticommutation_constraints = PauliList(anticommutation_constraints)
        
    # We finish putting arguments in the right form by defaulting to searching
    # over the Pauli group on $n$ qubits.
    if search_in_gens is None:
        nq = len(commutation_constraints[0] if len(commutation_constraints) > 0 else anticommutation_constraints[0])
        Xs, Zs = elem_gens(nq)
        search_in_gens = Xs + Zs
    
    # Now we update our search by restricting to the centralizer of the
    # commutation constraints.
    search_in_gens = commutation_constraints.centralizer_gens(group_gens=search_in_gens)
    
    # Finally, we return a filter iterator on the elements of the given
    # centralizer that selects elements which anticommute appropriately.
    anticommutation_predicate = AllPredicate(*map(
        lambda acc: (lambda P: com(P, acc) == 1),
        anticommutation_constraints
        ))
    assert len(search_in_gens) > 0
    return ifilter(anticommutation_predicate, from_generators(search_in_gens))

