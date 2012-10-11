#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# constraint_solvers.py: Solvers for various pc.Pauli and Clifford constraint
#     problems.
##
# © 2012 Christopher E. Granade (cgranade@gmail.pc.com) and
#     Ben Criger (bcriger@gmail.pc.com).
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

import PauliClass as pc
from paulicollections import PauliList
from pred import AllPredicate
from itertools import ifilter

## FUNCTIONS ##

def solve_commutation_constraints(
        commutation_constraints=[],
        anticommutation_constraints=[],
        search_in_gens=None,
        search_in_set=None
    ):
    r"""
    Given commutation constraints on a Pauli operator, yields an iterator onto
    all solutions of those constraints.
    
    :param commutation_constraints: A list of operators :math:`\{A_i\}` such
        that each solution :math:`P` yielded by this function must satisfy
        :math:`[A_i, P] = 0` for all :math:`i`.
    :param commutation_constraints: A list of operators :math:`\{B_i\}` such
        that each solution :math:`P` yielded by this function must satisfy
        :math:`\{B_i, P\} = 0` for all :math:`i`.
    :param search_in_gens: A list of operators :math:`\{N_i\}` that generate
        the group in which to search for solutions. If ``None``, defaults to
        the elementary generators of the pc.Pauli group on :math:`n` qubits, where
        :math:`n` is given by the length of the commutation and anticommutation
        constraints.
    :param search_in_set: An iterable of operators to which the search for 
        satisfying assignments is restricted. This differs from ``search_in_gens``
        in that it specifies the entire set, not a generating set. When this
        parameter is specified, a brute-force search is executed. Use only
        when the search set is small, and cannot be expressed using its generating
        set. 
    :returns: An iterator ``it`` such that ``list(it)`` contains all operators
        within the group :math:`G = \langle N_1, \dots, N_k \rangle\rangle`
        given by ``search_in_gens``, consistent with the commutation and
        anticommutation constraints.
        
    This function is based on finding the generators of the centralizer groups 
    of each commutation constraint, and is thus faster than a predicate-based
    search over the entire group of interest. The resulting iterator can be
    used in conjunction with other filters, however.
    
    >>> import qecc as q
    >>> list(q.solve_commutation_constraints(q.PauliList('XXI', 'IZZ', 'IYI'), q.PauliList('YIY')))
    [i^0 XII, i^0 IIZ, i^0 YYX, i^0 ZYY]
    >>> from itertools import ifilter
    >>> list(ifilter(lambda P: P.wt <= 2, q.solve_commutation_constraints(q.PauliList('XXI', 'IZZ', 'IYI'), q.PauliList('YIY'))))
    [i^0 XII, i^0 IIZ]
    """
        
    # Normalize our arguments to be PauliLists, so that we can obtain
    # centralizers easily.
    if not isinstance(commutation_constraints, PauliList):
        commutation_constraints = PauliList(commutation_constraints)
    if not isinstance(anticommutation_constraints, PauliList):
        # This is probably not necessary, strictly speaking, but it keeps me
        # slightly more sane to have both constraints represented by the same
        # sequence type.
        anticommutation_constraints = PauliList(anticommutation_constraints)

    # Then check that the arguments make sense.
    if len(commutation_constraints) == 0 and len(anticommutation_constraints) == 0:

        raise ValueError("At least one constraint must be specified.")

    #We default to executing a brute-force search if the search set is
    #explicitly specified:
    if search_in_set is not None:
        commutation_predicate = AllPredicate(*map(
            lambda acc: (lambda P: pc.com(P, acc) == 0),
            commutation_constraints
            ))
        commuters = filter(commutation_predicate, search_in_set)
        anticommutation_predicate = AllPredicate(*map(
            lambda acc: (lambda P: pc.com(P, acc) == 1),
            anticommutation_constraints
            ))
        return filter(anticommutation_predicate, commuters)

    # We finish putting arguments in the right form by defaulting to searching
    # over the pc.Pauli group on $n$ qubits.
    if search_in_gens is None:
        nq = len(commutation_constraints[0] if len(commutation_constraints) > 0 else anticommutation_constraints[0])
        Xs, Zs = pc.elem_gens(nq)
        search_in_gens = Xs + Zs
    
    # Now we update our search by restricting to the centralizer of the
    # commutation constraints.
    search_in_gens = commutation_constraints.centralizer_gens(group_gens=search_in_gens)
    
    # Finally, we return a filter iterator on the elements of the given
    # centralizer that selects elements which anticommute appropriately.
    anticommutation_predicate = AllPredicate(*map(
        lambda acc: (lambda P: pc.com(P, acc) == 1),
        anticommutation_constraints
        ))
    assert len(search_in_gens) > 0
    return ifilter(anticommutation_predicate, pc.from_generators(search_in_gens))

