#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# pred.py: Predicate library for use with stabilizer groups.
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

## IMPORTS ##
from sys import version_info
if version_info[0] == 3:
    PY3 = True
    from importlib import reload
elif version_info[0] == 2:
    PY3 = False
else:
    raise EnvironmentError("sys.version_info refers to a version of "
        "Python neither 2 nor 3. This is not permitted. "
        "sys.version_info = {}".format(version_info))

import operator as op
import itertools as it

if PY3:
    from . import PauliClass as pc
else:
    import PauliClass as pc

## ALL ##

# The following names will get imported by
# "from pred import *".
__all__ = [
    'Predicate',
    'AllPredicate', 'AnyPredicate',
    'SetMembershipPredicate', 'PauliMembershipPredicate',
    'commutes_with', 'in_group_generated_by'
]

## CLASSES ##

class Predicate(object):
    """
    Class representing a predicate function on one or more arguments.
    
    >>> from qecc import Predicate
    >>> p = Predicate(lambda x: x > 0)
    >>> p(1)
    True
    >>> p(-1)
    False
    
    Instances can also be constructed by logical operations on existing
    Predicate instances:
    
    >>> q = Predicate(lambda x: x < 3)
    >>> (p & q)(1)
    True
    >>> (p | q)(-1)
    True
    >>> (~p)(2)
    False
    
    """

    def __init__(self, fn):
        self.fn = fn
        
    def __call__(self, *args, **kwargs):
        return self.fn(*args, **kwargs)
        
    def combine(self, other, outer_fn):
        """
        Returns a new :class:`Predicate` that combines this predicate with
        another predicate using a given function to combine the results.
        
        >>> gt_2 = Predicate(lambda x: x > 2)
        >>> even = Predicate(lambda x: x % 2 == 0)
        >>> nand = lambda x, y: not (x and y)
        >>> r = gt_2.combine(even, nand)
        >>> map(r, range(1,5))
        [True, True, True, False]
        
        """
        def new_predicate(*args, **kwargs):
            return outer_fn(self(*args, **kwargs), other(*args, **kwargs))
            
        return Predicate(new_predicate)
        
    def __and__(self, other):
        return AllPredicate(self, other)
        
    def __or__(self, other):
        return AnyPredicate(self, other)
        
    def __invert__(self):
        def not_fn(*args, **kwargs):
            return not self(*args, **kwargs)
            
        return Predicate(not_fn)
 
class AllPredicate(Predicate):
    """
    Predicate class representing the logical AND of several other predicates.
    
    Given two predicates ``p`` and ``q``,
    
    >>> p = lambda: True
    >>> q = lambda: False
    >>> p_and_q = AllPredicate(p, q)
    >>> p_and_q () == p () & q ()
    True
    
    is equivalent to ``p_and_q = p & q``.
    """
    
    def __init__(self, *preds):
        self.preds = preds
        
    def __call__(self, *args, **kwargs):
        return all(p(*args, **kwargs) for p in self.preds)
        
    def __and__(self, other):
        if isinstance(other, AllPredicate):
            return AllPredicate(*(self.preds + other.preds))
        else:
            return AllPredicate(*(self.preds + (other,)))
        
class AnyPredicate(Predicate):
    """
    Predicate class representing the logical OR of several other predicates.
    
    Given two predicates ``p`` and ``q``,
    
    >>> p = lambda: True
    >>> q = lambda: False
    >>> p_or_q = AnyPredicate(p, q)
    >>> p_or_q () == p () | q ()
    True
    
    
    is equivalent to ``p_and_q = p | q``.
    """
    
    def __init__(self, *preds):
        self.preds = preds
        
    def __call__(self, *args, **kwargs):
        return any(p(*args, **kwargs) for p in self.preds)
        
    def __or__(self, other):
        if isinstance(other, AnyPredicate):
            return AnyPredicate(*(self.preds + other.preds))
        else:
            return AnyPredicate(*(self.preds + (other,)))
        
class SetMembershipPredicate(Predicate):
    """
    Given an iterable ``S``, constructs a predicate that returns ``True``
    if and only if its argument is in ``S``.
    
    >>> from qecc import SetMembershipPredicate
    >>> p = SetMembershipPredicate(range(4))
    >>> map(p, range(-1, 5))
    [False, True, True, True, True, False]
    
    """
    
    def __init__(self, S):
        self.S = set(S)
        
    def __call__(self, x):
        return x in self.S
    
    def __repr__(self):
        return "Predicate: f(x) = True if x in {0}".format(repr(self.S))
        

class PauliMembershipPredicate(SetMembershipPredicate):
    """
    Given a set ``S`` of Pauli operators represented as :class:`qecc.Pauli`
    instances, constructs a predicate that returns
    ``True`` for a Pauli ``P`` if and only if ``P`` is in ``S``.
    
    If the keyword argument ``ignore_phase`` is ``True``, then the
    comparison to determine whether ``P`` is in ``S`` only considers the
    operator part of ``P``.
    """
    
    def __init__(self, S, ignore_phase=True):
        super(PauliMembershipPredicate, self).__init__(
            [pc.Pauli(P.op) for P in S] if ignore_phase else S)
        self.ignore_phase = ignore_phase
        
    def __call__(self, P):
        if self.ignore_phase:
            P = pc.Pauli(P.op)
        return P in self.S

## USEFUL PREDICATES ##

def commutes_with(*paulis):
    """
    Returns a predicate that checks whether a Pauli ``P`` commutes with each of
    a given list of Pauli operators.
    """
    paulis = list(map(pc.ensure_pauli, paulis))

    def pred_fn(P):
        # Using imap here instead of map allows all() to short-circuit.
        return all(it.imap(lambda Q: pc.com(P, Q) == 0, paulis))
        
    return Predicate(pred_fn)
        
def in_group_generated_by(*paulis):
    """
    Returns a predicate that selects Pauli operators in the group generated by
    a given list of generators.
    """
    # Warning: This is inefficient for large groups!
    paulis = list(map(pc.ensure_pauli, paulis))
    
    return PauliMembershipPredicate(pc.from_generators(paulis), ignore_phase=True)
        
## TEST ##

if __name__ == "__main__":
    p = Predicate(lambda x: x > 0)
    q = Predicate(lambda x: x < 3)
    
    p_and_q = p & q
    p_or_q  = p | q
    not_p   = ~p
    
    for test in [2, 4, -1]:
        print(test, p(test), q(test), p_and_q(test), p_or_q(test), not_p(test))
        
    print(list(filter(p_and_q, list(range(-4,5)))))
    
    S = set([1, 2, 3])
    in_S = SetMembershipPredicate(S)
    
    print(list(map(in_S, list(range(-1, 5)))))
    
    print(list(filter(
        commutes_with('XX', 'ZZ') & ~in_group_generated_by('XX'),
        pc.pauli_group(2)
        )))
    
