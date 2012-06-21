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
from exceptions import *

from paulicollections import PauliList

from constraint_solvers import solve_commutation_constraints
from unitary_reps import clifford_as_unitary

from singletons import EmptyClifford, Unspecified

import warnings

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
KINDS = ['X', 'Z']
PHASES = [" +", "+i", " -", "-i"]

## CLASSES ##

class Clifford(object):
    r"""
    Class representing an element of the Cifford group on :math:`n`
    qubits.
    
    :param xbars: A list of operators :math:`\bar{X}_i` such that the
        represented Clifford operation :math:`C` acts as 
        :math:`C(X_i) = \bar{X}_i`. Note that in order for the represented
        operator to be an automorphism, each :math:`\bar{X}_i` must have phase
        either 0 or 2. A warning will result if this condition is not met.
    :param zbars: See ``xbars``.
    :type xbars: list of :class:`qecc.Pauli` instances
    :type zbars: list of :class:`qecc.Pauli` instances
    """
    
    def __init__(self, xbars, zbars):
        # Require that all specified operators be Paulis.
        # Moreover, we should warn the caller if the output phase is not either
        # 0 or 2, since such operators are not automorphisms of the Pauli group.
        self.xout = PauliList(*xbars)
        self.zout = PauliList(*zbars)
        
        for output_xz in chain(self.xout, self.zout):
            if output_xz is not Unspecified and output_xz.ph not in [0, 2]:
                warnings.warn(
                    'The output phase of a Clifford operator has been specified as {}, such that the operator is not a valid automorphism.\n'.format(
                        output_xz.ph
                    ) +
                    'To avoid this warning, please choose all output phases to be from the set {0, 2}.'
                )
                
        # Prevent fully unspecified operators.
        if all([P is Unspecified for P in chain(xbars, zbars)]):
            raise ValueError("At least one output must be specified.")
        

    def __len__(self):
        """
        Yields the number of qubits on which the Clifford ``self`` acts.
        """
        for P in self.xout + self.zout:
            if P is not Unspecified:
                return len(P)

    def __repr__(self):
        return "<Clifford operator on {nq} qubit{s} at 0x{id:x}>".format(
            nq=len(self),
            id=id(self),
            s='s' if len(self) > 1 else ''
        )

    def __str__(self):
        """
        Returns a string representing a Clifford in Pauli notation (yielding a
        list of input Paulis, and a list of output Paulis.)
        """
        SPARSE_NQ = 3
        if len(self) > SPARSE_NQ:
            print len(self), SPARSE_NQ
            return self.str_sparse()
        
        left_side_x,left_side_z=elem_gens(len(self))
        right_side=self.xout+self.zout
        return '\n'.join(
                '{gen.op} |-> {outsign}{out}'.format(
                    gen=gen,
                    out=out.op if out is not Unspecified else "Unspecified",
                    outsign=PHASES[out.ph] if out is not Unspecified else ""
                    )
                for gen, out in zip(left_side_x + left_side_z, self.xout + self.zout)
            )
            
    def str_sparse(self):
        out = zip(KINDS, map(enumerate, [self.xout, self.zout]))
        nq = len(self)
        return "\n".join(
            [
                "{}[{}] |-> {}{}".format(kind, idx, PHASES[P.ph], P.str_sparse(incl_ph=False))
                for kind, Ps in out for idx, P in Ps
                if elem_gen(nq, idx, kind) != P
            ])

    def is_valid(self, quiet=True):
        """
        Returns ``True`` if this instance represents a valid automorphism. In
        particular, this method returns ``True`` if all output phase assignments
        are either 0 or 2, and if all of the commutation relations on its
        outputs are obeyed. Unspecified outputs are ignored.
        
        :param bool quiet: If set to ``True``, this method will not print out
            any information, but will return ``True`` or ``False`` as described
            above. Otherwise, if the operator is not a valid Clifford operator,
            diagnostic information will be printed.
        """
        if any(P.ph not in [0, 2] for P in chain(self.xout, self.zout) if P is not Unspecified):
            if not quiet:
                print "At least one output operator has a phase other than 0 or 2."
                
            return False
        
        for P in sum(elem_gens(len(self)), []):
            for Q in sum(elem_gens(len(self)), []):
                UP = self.conjugate_pauli(P)
                UQ = self.conjugate_pauli(Q)
                
                if UP is not Unspecified and UQ is not Unspecified and com(UP, UQ) != com(P, Q):
                    if not quiet:
                        print "c({P}, {Q}) == {cPQ}, but c(U({P}), U({Q})) == c({UP}, {UQ}) == {cUPUQ}.".format(
                                P=P, Q=Q,
                                cPQ=com(P,Q),
                                UP=UP, UQ=UQ,
                                cUPUQ=com(UP, UQ)
                            )
                    return False
                    
        return True
        
    def inv(self):
        print "Phase information will be lost, aBSM notation will fix this."
        return self.as_bsm().inv().as_clifford

    def conjugate_pauli(self,pauli):
        r"""

        Given an instance of :class:`qecc.Pauli` representing the
        operator :math:`P`, calculates the mapping 
        :math:`CPC^{\dagger}`, where :math:`C` is the operator represented by
        this instance.

        :arg pauli: Representation of the Pauli operator :math:`P`.
        :type pauli: qecc.Pauli
        :returns: Representation of the Pauli operator 
            :math:`CPC^{\dagger}`, where :math:`C` is the Clifford operator
            represented by this instance.
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
        rolling_pauli=Pauli('I'*len(pauli), phase=pauli.ph)
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
        if isinstance(other, Clifford):
            return all(P == Q for P, Q in zip(self.xout + self.zout, other.xout + other.zout))
        else:
            return False
            
    def __ne__(self, other):
        return not self.__eq__(other)


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
        Yields an iterator onto possible Clifford operators whose outputs agree
        with this operator for all outputs that are specified. Note that all
        yielded operators assign the phase 0 to all outputs, by convention.
        
        If this operator is fully specified, the iterator will yield exactly one
        element, which will be equal to this operator.
        
        For example:
        
        >>> import qecc as q
        >>> C = q.Clifford([q.Pauli('XI'), q.Pauli('IX')], [q.Unspecified, q.Unspecified])
        >>> it = C.constraint_completions()
        >>> print it.next()
        XI |->  +XI
        IX |->  +IX
        ZI |->  +ZI
        IZ |->  +IZ
        >>> print it.next()
        XI |->  +XI
        IX |->  +IX
        ZI |->  +ZI
        IZ |->  +IY
        >>> print len(list(C.constraint_completions()))
        8
        
        If this operator is not a valid Clifford operator, then this method will
        raise an :class:`qecc.InvalidCliffordError` upon iteraton.
        """
        # Check for validity first.
        if not self.is_valid():
            raise InvalidCliffordError("The specified constraints are invalid or are contradictory.")
        
        # Useful constants.
        XKIND, ZKIND = range(2)
        
        # Start by finding the first unspecified output.
        nq = len(self)
        X_bars, Z_bars = self.xout, self.zout
        P_bars = map(copy, [X_bars, Z_bars]) # <- Useful for indexing by kinds.
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
            P_bars[unspecified_kind][unspecified_idx] = Pauli(P.op, phase=0)
            # I wish I had "yield from" here. Ah, well. We have to recurse
            # manually instead.
            C = Clifford(*P_bars)
            for completion in C.constraint_completions():
                yield completion
                
        return

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
        """
        Returns a :class:`numpy.ndarray` containing a unitary matrix
        representation of this Clifford operator.
        
        Raises a :class:`RuntimeError` if NumPy cannot be imported.
        """
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

def transcoding_cliffords(stab_in,xs_in,zs_in,stab_out,xs_out,zs_out):
    r"""
    
    """
    #Preliminaries:
    nq_in=len(stab_in[0])
    nq_out=len(stab_out[0])
    nq_anc=abs(nq_in-nq_out)

    #Decide left side:
    if nq_in<nq_out:
        stab_left=stab_out
        xs_left=xs_out
        zs_left=zs_out
        stab_right=stab_in
        xs_right=xs_in
        zs_right=zs_in
    else:
        stab_right=stab_out
        xs_right=xs_out
        zs_right=zs_out
        stab_left=stab_in
        xs_left=xs_in
        zs_left=zs_in
        
    cliff_xouts_left=stab_left+xs_left
    cliff_zouts_left=[Unspecified]*len(stab_left)+zs_left
    
    cliff_left=Clifford(cliff_xouts_left,cliff_zouts_left).constraint_completions().next()
    list_left=cliff_left.xout+cliff_left.zout

    for mcset in mutually_commuting_sets(n_gens=len(stab_left)-len(stab_right),n_bits=nq_anc):
        temp_xouts_right=pad(stab_right,lower_right=mcset)+map(lambda p: p&eye_p(nq_anc),xs_right)
        temp_zouts_right=[Unspecified]*len(stab_left)+map(lambda p: p&eye_p(nq_anc),zs_right)
    for completion in Clifford(temp_xouts_right,temp_zouts_right).constraint_completions():
        if nq_in<nq_out:
            yield gen_cliff(completion.xout+completion.zout,list_left)
        else:
            yield gen_cliff(list_left,completion.xout+completion.zout)


def min_len_transcoding_clifford(paulis_in,paulis_out):
    circuit_iter=map(lambda p: p.as_bsm().circuit_decomposition(), transcoding_cliffords(paulis_in,paulis_out))
    return min(*circuit_iter)
    
def clifford_group(nq, consider_phases=False):
    r"""
    Given a number of qubits :math:`n`, returns an iterator that produces all
    elements of :math:`\mathcal{C}_n`, the Clifford group on :math:`n` qubits.
    
    :param int nq: The number of qubits upon which each yielded element will
        act.
    :param bool consider_phases: If ``True``, then Clifford operators whose
        assignments of phases to the generators of the Pauli group differ
        will be treated as distinct. Otherwise, the yielded elements will
        be drawn from the group
        :math:`\hat{\mathcal{C}}_n = \mathrm{Aut}(\hat{\mathcal{P}}_n / \{ i^k I : k \in \mathbb{Z}_4 \})`,
        such that the phases of the outputs are not considered.
    """
    idx = 0
    for P in pauli_group(nq):
        if P.wt > 0:
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

