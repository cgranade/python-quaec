#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# paulicollections.py: 
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

import PauliClass as pc
from collections import Sequence
from singletons import Unspecified

import unitary_reps

## ALL ##

__all__ = [
    'PauliList'
    ]
        
## CLASSES ##

class PauliList(list):
    r"""
    Subclass of :obj:`list` offering useful methods for lists of
    :class:`qecc.Pauli` instances.
    
    :param paulis: Instances either of :obj:`str` or :class:`qecc.Pauli`, or
        the special object :obj:`qecc.Unspecified`.
        Strings are passed to the constructor of :class:`qecc.Pauli` for
        convinenence.
    """

    def __init__(self, *paulis):
        if len(paulis) == 1:
            single_arg = paulis[0]
            if isinstance(single_arg, str):
                paulis = [pc.Pauli(single_arg)]
            elif isinstance(single_arg, Sequence) or hasattr(single_arg, '__iter__'):
                paulis = map(pc.ensure_pauli, single_arg)
            else:
                paulis = map(pc.ensure_pauli, paulis)
        else:
            paulis = map(pc.ensure_pauli, paulis)
            
        # FIXME: figure out why super(list, self).__init__ doesn't work.
        list.__init__(self, paulis)
        
    def __getitem__(self, *args):
        item = super(PauliList, self).__getitem__(*args)
        if not isinstance(item, list):
            return item
        else:
            return PauliList(*item)
        
    def __getslice__(self, *args):
        # Note that this must be overrided due to an implementation detail of
        # CPython. See the note at
        #     http://docs.python.org/reference/datamodel.html#additional-methods-for-emulation-of-sequence-types
        return PauliList(*super(PauliList, self).__getslice__(*args))
        
    def __add__(self, other):
        return PauliList(*(super(PauliList, self).__add__(other)))
        
    ## REPRESENTATION AND STRING COVNERSIONS ##
    
    def __repr__(self):
        # For now, he repr and str for a list can be the same. We can improve
        # this in the future.
        return str(self)
    
    def __str__(self):
        return "PauliList({})".format(", ".join(map(repr, self)))
        
    ## OPERATORS ACTING ON PAULI LISTS ##
    
    def __and__(self, other):
        if not isinstance(other, pc.Pauli):
            return NotImplemented
            
        else:
            return PauliList(*[P & other for P in self])
    
    def __rand__(self, other):
        if not isinstance(other, pc.Pauli):
            return NotImplemented
        
        else:
            return PauliList(*[other & P for P in self])
    
    ## OTHER METHODS ##
    def pad(self, extra_bits=0, lower_right=None):
        r"""
        Takes a PauliList, and returns a new PauliList, 
        appending ``extra_bits`` qubits, with stabilizer operators specified by
        ``lower_right``.
        
        :arg pauli_list_in: list of Pauli operators to be padded. 
        :param int extra_bits: Number of extra bits to be appended to the system.
        :param lower_right: list of `qecc.Pauli` operators, acting on `extra_bits` qubits.
        :rtype: list of :class:`qecc.Pauli` objects.
        
        Example:
        
        >>> import qecc as q
        >>> pauli_list = q.PauliList('XXX', 'YIY', 'ZZI')
        >>> pauli_list.pad(extra_bits=2, lower_right=q.PauliList('IX','ZI'))
        PauliList(i^0 XXXII, i^0 YIYII, i^0 ZZIII, i^0 IIIIX, i^0 IIIZI)

        """
        
        len_P = len(self)
        nq_P  = len(self[0]) if len_P > 0 else 0

        if extra_bits == 0 and lower_right is None or len(lower_right) == 0:
            return PauliList(self)
        elif len(lower_right) != 0:
            extra_bits=len(lower_right[0])
                
        setout = PauliList([pc.Pauli(pauli.op + 'I'*extra_bits) for pauli in self])
            
        if lower_right is None:
            setout += [pc.eye_p(nq_P + extra_bits)] * extra_bits
        else:
            setout += [pc.eye_p(nq_P) & P for P in lower_right]
                
        return setout    
    
    def generated_group(self, coset_rep=None):
        """
        Yields an iterator onto the group generated by this list of Pauli
        operators. See also :obj:`qecc.from_generators`.
        """
        return pc.from_generators(self, coset_rep)
    
    def stabilizer_subspace(self):
        r"""
        Returns a :class:`numpy.ndarray` of shape ``(n - k, 2 ** n)`` containing
        an orthonormal basis for the mutual +1 eigenspace of each fully
        specified Pauli in this list. Here, ``n`` is taken to be the number of
        qubits and ``k`` is taken to be the number of independent Pauli
        operators in this list.
        
        Raises a :obj:`RuntimeError` if NumPy cannot be imported.
        
        For example, to find the Bell basis vector :math:`\left|\beta_{00}\right\rangle`
        using the stabilizer formalism:
        
        >>> import qecc as q
        >>> q.PauliList('XX', q.Unspecified, q.Unspecified, 'ZZ').stabilizer_subspace()
        array([[ 0.70710678+0.j,  0.00000000+0.j,  0.00000000+0.j,  0.70710678+0.j]])
        
        Similarly, one can find the codewords of the phase-flip code :math:`S = \langle XXI, IXX \rangle`:
        
        >>> q.PauliList('XXI', 'IXX').stabilizer_subspace()
        array([[ 0.50000000-0.j, -0.00000000-0.j, -0.00000000-0.j,  0.50000000-0.j,
                -0.00000000-0.j,  0.50000000+0.j,  0.50000000-0.j, -0.00000000-0.j],
               [ 0.02229922+0.j,  0.49950250+0.j,  0.49950250+0.j,  0.02229922+0.j,
                 0.49950250+0.j,  0.02229922+0.j,  0.02229922+0.j,  0.49950250+0.j]])
                 
        Note that in this second case, some numerical errors have occured;
        this method does not guarantee that the returned basis
        vectors are exact.
        """
        return unitary_reps.mutual_eigenspace([P.as_unitary() for P in self if P is not Unspecified])
        
    def centralizer_gens(self, group_gens=None):
        r"""
        Returns the generators of the centralizer group
        :math:`\mathrm{C}(P_1, \dots, P_k)`, where :math:`P_i` is the :math:`i^{\text{th}}`
        element of this list. See :meth:`qecc.Pauli.centralizer_gens` for
        more information.
        """
        if group_gens is None:
            # NOTE: Assumes all Paulis contained by self have the same nq.
            Xs, Zs = pc.elem_gens(len(self[0]))
            group_gens = Xs + Zs
            
        if len(self) == 0:
            # C({}) = G
            return PauliList(group_gens)
            
        centralizer_0 = self[0].centralizer_gens(group_gens=group_gens)
            
        if len(self) == 1:
            return centralizer_0
        else:
            return self[1:].centralizer_gens(group_gens=centralizer_0)
        
