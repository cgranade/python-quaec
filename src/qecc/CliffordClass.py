#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# CliffordClass.py: Implementation of qecc.Clifford and related utility
#     functions.
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

## RELOAD FIX ##
# This is a pretty poor way of fixing it, but should work for now.
    
import PauliClass as _pc
import bsf as _bsf

try:
    reload(_pc)
    reload(_bsf)
except:
    # If it fails, we're no worse off, so just ignore the problem.
    pass

## IMPORTS ##

import operator as op

from copy import copy, deepcopy
from itertools import product, chain, combinations
from PauliClass import *
from bsf import *
from numpy import hstack, newaxis

from constraint_solvers import solve_commutation_constraints
from unitary_reps import clifford_as_unitary

from singletons import EmptyClifford, Unspecified

## ALL ##

__all__ = [
    'Clifford',
    'eye_c', 'cnot', 'replace_one_character', 'cz', 'hadamard',
    'phase', 'permutation', 'swap', 'pauli_gate', 'paulify',
    'generic_clifford',
    'gen_cliff', 'transcoding_cliffords', 'min_len_transcoding_clifford',
    'clifford_group'
]

## CONSTANTS ##

VALID_OPS = ['I', 'X', 'Y', 'Z']
VALID_PHS = range(4)

## CLASSES ##

class Clifford(object):
    r"""
    Class representing an element of the Cifford group on :math:`n`
     qubits.
    
    :param xbars: A list of operators :math:`\bar{X}_i` such that the
        represented Clifford operation :math:`C` acts as 
        :math:`C(X_i) = \bar{X}_i`.
    :param zbars: See ``xbars``.
    :type xbars: list of :class:`qecc.Pauli` instances
    :type zbars: list of :class:`qecc.Pauli` instances
    """
    
    def __init__(self, xbars, zbars):
        # TODO: add xbars_in and zbars_in as optional arguments, then pass them
        #       to generic_clifford.
        # TODO: check that at least one output is specified.
        for output_xz in xbars+zbars:
            if (output_xz is not Unspecified) and (not isinstance(output_xz,Pauli)):
                raise TypeError("Output operators must be Paulis.")
        # Prevent fully unspecified operators.
        if all(P is Unspecified for P in xbars + zbars):
            raise ValueError("At least one output must be specified.")
        self.xout=copy(xbars)
        self.zout=copy(zbars)

    def __len__(self):
        """

        Yields the number of qubits on which the Clifford ``self`` acts.

        """
        for P in self.xout + self.zout:
            if P is not Unspecified:
                return len(P)

    def __repr__(self):
        """Prints a Clifford in Pauli notation (yielding a list of 
        input Paulis, and a list of output Paulis.)"""
        left_side_x,left_side_z=elem_gens(len(self))
        right_side=self.xout+self.zout
        return '\n'.join(
                '{gen} |-> {out}'.format(gen=gen, out=out)
                for gen, out in zip(left_side_x + left_side_z, self.xout + self.zout)
            )

    def is_valid(self):
        """
        Checks that the output of the represented Clifford gate obeys
         the proper commutation relations.
        """
        for P in sum(elem_gens(len(self)), []):
            for Q in sum(elem_gens(len(self)), []):
                if com(self.conjugate_pauli(P), self.conjugate_pauli(Q)) != com(P, Q):
                    print P, Q, self.conjugate_pauli(P), self.conjugate_pauli(Q)
                    return False
                    
        return True

    def conjugate_pauli(self,pauli):
        r"""

        Given an instance of :class:`qecc.Pauli` representing the
         operator :math:`P`, calculates the mapping 
         :math:`CPC^{\dagger}`.

        :arg pauli: Representation of the Pauli operator :math:`P`.
        :type pauli: qecc.Pauli
        :returns: Representation of the Pauli operator 
        :math:`CPC^{\dagger}`,
            where :math:`C` is the Clifford operator represented by this
            instance.
        :rtype: qecc.Pauli
        """
        
        if not isinstance(pauli,Pauli):
            # If we don't have a Pauli, maybe we have an iterable.
            try:
                dummy = iter(pauli)
                # Yep. It was an iterable.
                return map(self.conjugate_pauli, pauli)
            except TypeError:
                # Nope. Wasn't iterable. Raise an error.
                raise TypeError("Cliffords conjugate Paulis.")
        #Initialize the output Pauli to the identity:
        rolling_pauli=Pauli('I'*len(pauli))        
        for idx,op in enumerate(pauli.op):
            #For every X/Z the input Pauli contains, multiply by the
            # corresponding output Pauli from self. 
            if op == 'X':
                rolling_pauli=rolling_pauli*self.xout[idx]
            elif op == 'Z':
                rolling_pauli=rolling_pauli*self.zout[idx]
            elif op == 'Y':
                #Y = iXZ:
                rolling_pauli=rolling_pauli*self.xout[idx]*self.zout[idx]
                rolling_pauli.mul_phase(1)
        return rolling_pauli 

    def __eq__(self,other):
        return (self.xout==other.xout)and(self.zout==other.zout)


    def __mul__(self,other):
        """multiplies two Cliffords, self and other, by conjugating the
         output Paulis from other, according to the relations given by
          self, yielding self*other. """
        if not isinstance(other,Clifford):
            return NotImplemented 
        Xs=[]
        Zs=[]
        for ex, zed in zip(other.xout,other.zout):
            Xs.append(self.conjugate_pauli(ex))
            Zs.append(self.conjugate_pauli(zed))
        return Clifford(Xs,Zs)

    def __rand__(self, other):
        if other is EmptyClifford:
            return self
            
        return NotImplemented

    def __and__(self,other):
        """Takes the tensor product of two Cliffords *self* and
         *other*."""
        if other is EmptyClifford:
            return self
            
        if not isinstance(other,Clifford):
            return NotImplemented 
        nq_self=len(self)
        nq_other=len(other)
        id_self_size=eye_p(nq_self)
        id_other_size=eye_p(nq_other)
        """We embed each Clifford into a larger space, and concatenate
         the output lists."""
        exones=[]
        extwos=[]
        zedones=[]
        zedtwos=[]
        for idx in range(nq_self):
            exones.append(self.xout[idx] & id_other_size)
            zedones.append(self.zout[idx] & id_other_size)
        for idx in range(nq_other):
            extwos.append(id_self_size & other.xout[idx])
            zedtwos.append(id_self_size & other.zout[idx])
        return Clifford(exones+extwos,zedones+zedtwos)

    def __call__(self, other):
        if not isinstance(other, Pauli):
            return NotImplemented
        return self.conjugate_pauli(other)
        
    def constraint_completions(self):
        """
        """
        # Note: disregards phases.
        XKIND, ZKIND = range(2)
        
        # Start by finding the first unspecified output.
        nq = len(self)
        X_bars, Z_bars = self.xout, self.zout
        P_bars = [X_bars, Z_bars] # <- Useful for indexing by kinds.
        XZ_pairs = zip(X_bars, Z_bars)
        try:
            unspecified_idx, unspecified_kind = iter(
                (idx, kind)
                for idx, kind
                in product(xrange(nq), range(2))
                if XZ_pairs[idx][kind] is Unspecified
            ).next()
        except StopIteration:
            # If there are no unspecified constraints, then self is the only
            # satisfying completion.
            yield self
            return
        
        # We must always commute with disjoint qubits.
        commutation_constraints = reduce(op.add,
            (XZ_pairs[idx] for idx in xrange(nq) if idx != unspecified_idx),
            tuple()
            )
            
        # On the same qubit, we must anticommute with the opposite operator.
        anticommutation_constraints = [XZ_pairs[unspecified_idx][XKIND if unspecified_kind == ZKIND else ZKIND]]
        
        # Filter out Unspecified constraints.
        specified_pred = lambda P: P is not Unspecified
        commutation_constraints = filter(specified_pred, commutation_constraints)
        anticommutation_constraints = filter(specified_pred, anticommutation_constraints)
        
        # Now we iterate over satisfactions of the constraints, yielding
        # all satisfactions of the remaining constraints recursively.
        Xgs, Zgs = elem_gens(nq)
        for P in solve_commutation_constraints(commutation_constraints, anticommutation_constraints, search_in_gens=Xgs+Zgs):
            P_bars[unspecified_kind][unspecified_idx] = P.mul_phase(-P.ph)
            # I wish I had "yield from" here. Ah, well. We have to recurse
            # manually instead.
            C = Clifford(*P_bars)
            for completion in C.constraint_completions():
                yield completion

    def as_bsm(self):
        """
        Returns a representation of the Clifford operator as a binary
         symplectic matrix.
        
        :rtype: :class:`qecc.BinarySymplecticMatrix`
        """
        def to_col(P):
            v = P.as_bsv()
            out = hstack([v.x, v.z])[..., newaxis]
            return out
        return BinarySymplecticMatrix(hstack(map(to_col, self.xout + self.zout)))
        
    def as_unitary(self):
        return clifford_as_unitary(self)
        
    
## FUNCTIONS ##
def eye_c(nq):
    """
    Yields the identity Clifford, defined to map every generator of the 
    Pauli group to itself.

    :rtype: Clifford
    """
    return Clifford(*elem_gens(nq)) if nq > 0 else EmptyClifford
    
def replace_one_character(string,location,new_character):
    """
    Replaces the character in ``string`` at ``location`` with
     ``new_character``.

    :rtype: str
    """
    return string[:location]+new_character+string[location+1:]
    
def cnot(nq,ctrl,targ):
    """
    Yields the ``nq``-qubit CNOT Clifford controlled on ``ctrl``,
     acting a Pauli :math:`X` on ``targ``.

    :rtype: :class:`qecc.Clifford`
    """
    #Initialize to the identity Clifford:
    cnotto=eye_c(nq)
    #Wherever ctrl has an X, put an X on targ:
    cnotto.xout[ctrl].op=replace_one_character(cnotto.xout[ctrl].op,targ,'X')
    #Wherever targ has a Z, put a Z on ctrl:
    cnotto.zout[targ].op=replace_one_character(cnotto.zout[targ].op,ctrl,'Z')
    return cnotto
    
def cz(nq, q1, q2):
    """
    Yields the ``nq``-qubit C-Z Clifford, acting on qubits ``q1`` and
     ``q2``.

    :rtype: :class:`qecc.Clifford`
    """
    #Initialize to the identity Clifford:
    gate = eye_c(nq)
    #Wherever ctrl or targ get an X, map to XZ:
    gate.xout[q1].op = replace_one_character(gate.xout[q1].op, q2, 'Z')
    gate.xout[q2].op = replace_one_character(gate.xout[q2].op, q1, 'Z')
    return gate
    
def hadamard(nq,q):
    """
    Yields the ``nq``-qubit Clifford, switching :math:`X` and :math:`Z`
     on qubit ``q``, yielding a minus sign on :math:`Y`.

    :rtype: :class:`qecc.Clifford`
    """
    #Switch a Z and an X in the identity Clifford:
    return eye_c(q) & Clifford([Pauli('Z')],[Pauli('X')]) & eye_c(nq-q-1)

def phase(nq,q):
    r"""
    Yields the :math:`\frac{\pi}{4}_z`-rotation Clifford, acting on qubit ``q``.

    :rtype: :class:`qecc.Clifford`
    """
    return eye_c(q) & Clifford([Pauli('Y')],[Pauli('Z')]) & eye_c(nq-q-1)
    
def permutation(lst, p):
    """
    Permutes a list ``lst`` according to a set of indices ``p``.

    :rtype: list
    """
    return [lst[idx] for idx in p]
    
def swap(nq, q1, q2):
    """
    Yields the swap Clifford, on ``nq`` qubits, which swaps the Pauli generators on ``q1`` and ``q2``.

    :rtype: :class:`qecc.Clifford`
    """
    p = range(nq)
    p[q1], p[q2] = p[q2], p[q1]
    
    gate = eye_c(nq)
    gate.xout = permutation(gate.xout, p)
    gate.zout = permutation(gate.zout, p)
    return gate
    
def pauli_gate(pauli):
    """
    Imports an instance of the :class:`qecc.Pauli` class into the :class:`qecc.Clifford` class, representing a Pauli as a series of sign changes.

    :rtype: :class:`qecc.Clifford`
    """
    nq = len(pauli.op)
    return Clifford(*tuple(
        [gen.mul_phase(2*com(pauli,gen))  for gen in gen_set]
        for gen_set in elem_gens(nq)
    ))

def paulify(clinput):
    """
    Tests an input Clifford ``clinput`` to determine if it is, in
    fact, a Pauli. If so, it outputs the Pauli. If not, it
    returns the Clifford. 

    BE WARNED: If you turn a Pauli 
    into a Clifford and back again, the phase will be lost.
    """
    
    nq=len(clinput.xout) #Determine number of qubits.
    test_ex,test_zed=elem_gens(nq) #Get paulis to compare.
    """If the Paulis input to the Clifford are only altered in phase, then the Clifford is also a Pauli."""
    for ex_clif,zed_clif,ex_test,zed_test in zip(clinput.xout, clinput.zout,test_ex,test_zed):
        if ex_clif.op != ex_test.op or zed_clif.op != zed_test.op:
            print "Clifford is not Pauli."
            return clinput
        #If the Clifford is Pauli, determine which by examining operators with altered phases.
        exact=eye_p(nq)
        zedact=eye_p(nq) #Initialize accumulators
        """If a negative sign appears on a given generator, assign a Pauli to that qubit that conjugates the generator to a minus sign, e.g. ZXZ = -X """
        for idx_x in range(nq):
            if clinput.xout[idx_x].ph==2:
                exact.op = replace_one_character(exact.op, idx_x, 'Z')
        for idx_z in range(nq):
            if clinput.zout[idx_z].ph==2:
                zedact.op = replace_one_character(zedact.op, idx_z, 'X')
        return Pauli((exact*zedact).op)

def generic_clifford(paulis_in, paulis_out):
    """
    Given two lists of :class:`qecc.Pauli` instances, ``paulis_in`` and
    ``paulis_out``, produces an instance ``C`` of :class:`qecc.Clifford` such that
    ``C(paulis_in[i]) == paulis_out[i]`` for all ``i`` in ``range(2 * nq)``,
    where ``nq`` is the length of each element of the two lists.
    
    Each of ``paulis_in`` and ``paulis_out`` is assumed to be ordered such that
    the slice ``[0:nq]`` produces a list of logical :math:`X` operators, and
    such that the slice ``[nq:2*nq]`` produces the logical :math:`Z` operators.
    
    :param paulis_in: A list of length ``2 * nq`` logical Pauli operators
        specifying the input constraints for the desired Clifford operation.
    :param paulis_out: A list of length ``2 * nq`` logical Pauli operators
        specifying the output constraints for the desired Clifford operation.
    :return: A Clifford operator mapping the input constraints to the output
        constraints.
    :rtype: qecc.Clifford
    """
    nq=len(paulis_in)/2
    
    xins=paulis_in[0:nq]
    zins=paulis_in[nq:2*nq]
    
    xouts=paulis_out[0:nq]
    zouts=paulis_out[nq:2*nq]
    
    G    = Clifford(xouts,zouts)
    H    = Clifford(xins,zins)    
    Hinv = (H.as_bsm().inv().as_clifford())
    return G*Hinv
    
# For backwards compatibility, we define gen_cliff as an alias.
gen_cliff = generic_clifford

def transcoding_cliffords(paulis_in,paulis_out):
    r"""
    This function produces an iterator onto all Cliffords that take
    the pauli set 'paulis_in' to the set 'paulis_out', looping over all
    possible sets of constraints on the remaining degrees of freedom. 
    """
    nq_in=len(paulis_in[0])
    nq_out=len(paulis_out[0])
    gens_in=len(paulis_in)
    gens_out=len(paulis_out)
    #Set up temporary lists so that the larger Paulis are on the left.
    if nq_in>=nq_out:
        paulis_temp_left=paulis_in
        paulis_temp_right=paulis_out
        nq_tl=nq_in
        nq_tr=nq_out
        gens_tl=gens_in
        gens_tr=gens_out
    elif nq_in<nq_out:
        paulis_temp_left=paulis_out
        paulis_temp_right=paulis_in
        nq_tl=nq_out
        nq_tr=nq_in
        gens_tl=gens_out
        gens_tr=gens_in
    #Construct left_to_2n, a list of Paulis that complete the input Cliffords on the left.
    #First, get a mutually commuting set of nq_in generators, then round up using a
    #single clifford_bottom.
    if gens_tl < nq_tl:
        left_to_n=list(next(mutually_commuting_sets(nq_tl-gens_tl,nq_tl,group_iter=ns_mod_s(*paulis_temp_left))))
        left_to_2n=left_to_n+list(next(clifford_bottoms(paulis_temp_left+left_to_n)))
    else:
        left_to_2n=list(next(clifford_bottoms(paulis_temp_left)))
    #The right side constraints are what we have to iterate over.
    pads_mid_right=mutually_commuting_sets(gens_tl-gens_tr,nq_tl-nq_tr)
    for pad_mr in pads_mid_right:
        paulis_temp_left, paulis_temp_right = pad(paulis_temp_left,paulis_temp_right,lower_right=pad_mr)
        if len(paulis_temp_right) < nq_tr:
            for mcs in mutually_commuting_sets(nq_tr-len(paulis_temp_right),nq_tr,group_iter=ns_mod_s(*paulis_temp_right)):
                paulis_n=paulis_temp_right+mcs
        else:
            paulis_n=paulis_temp_right
        for bottom_half in clifford_bottoms(paulis_n):
            paulis_2n_right=paulis_temp_right+list(bottom_half)
            print list(chain(paulis_temp_left,left_to_2n)),paulis_2n_right
            if nq_in>=nq_out:
                yield gen_cliff(list(chain(paulis_temp_left,left_to_2n)),paulis_2n_right)
            else:
                yield gen_cliff(paulis_2n_right,list(chain(paulis_temp_left,left_to_2n)))

def min_len_transcoding_clifford(paulis_in,paulis_out):
    circuit_iter=map(lambda p: p.as_bsm().circuit_decomposition(), transcoding_cliffords(paulis_in,paulis_out))
    return min(*circuit_iter)
    
def clifford_group(nq, consider_phases=False):
    idx = 0
    for P in pauli_group(nq):
        if P.wt() > 0:
            C = Clifford([P] + [Unspecified]*(nq -1), [Unspecified]*nq)
            for completion in C.constraint_completions():
                if consider_phases:
                    P_bars = [completion.xout, completion.zout]
                    # phase_array is chosen to disregard global phases by
                    # absorbing them into xout[0].
                    for phase_array in product([0, 2], repeat=2*nq):
                        for idx_kind, idx_qubit in product(range(2), range(nq)):
                            P_bars[idx_kind][idx_qubit].ph = phase_array[nq*idx_kind + idx_qubit]
                        yield Clifford(*P_bars)
                else:
                    yield completion
            
