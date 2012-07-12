#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# stab.py: Classes and methods encapsulating and manipulating stabilizer codes.
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

import operator as op
import itertools as it
import math

import PauliClass as p # Sorry for the confusing notation here.
import CliffordClass as c
import paulicollections as pc

from collections import defaultdict
from singletons import EmptyClifford, Unspecified
from constraint_solvers import *

import warnings

## ALL ##

__all__ = [
    'StabilizerCode'
]

## CLASSES ##

class StabilizerCode(object):
    r"""
    Class representing a stabilizer code specified by the generators of its
    stabilizer group and by representatives for the logical operators acting
    on the code.
    
    :param group_generators: Generators :math:`N_i` such that the stabilizer
        group :math:`S` of the represented code is given by
        :math:`S = \left\langle N_i \right\rangle`.
    :param logical_xs: Representatives for the logical :math:`X` operators
        acting on encoded states.
    :param logical_zs: Representatives for the logical :math:`Z` operators
        acting on encoded states.
    """
    
    ## CONSTRUCTOR ##
    
    def __init__(self, group_generators, logical_xs, logical_zs):
        self.group_generators = pc.PauliList(*group_generators)
        self.logical_xs = pc.PauliList(*logical_xs)
        self.logical_zs = pc.PauliList(*logical_zs)
        
    ## PRETTY PRINTING ##
        
    def __repr__(self):
        return "<[[{n}, {k}, {d}]] StabilizerCode at {id:0x}>".format(
            n=self.nq, k=self.nq_logical,
            d=self.distance if self.nq < 6 else "?",
            id=id(self)
        )
    
    def __str__(self):
        return "S = <{group_generators}>\nXbars = {0.logical_xs}\nZbars = {0.logical_zs}".format(
            self,
            group_generators=", ".join(map(str, self.group_generators)))
       
    ## READ-ONLY PROPERTIES ##
        
    @property
    def nq(self):
        """
        The number of physical qubits into which this code encodes data.
        """
        try:
            return len(iter(
                gen
                for gen in self.group_generators + self.logical_xs + self.logical_zs
                if gen is not Unspecified
            ).next())
        except StopIteration:
            return 0
        
    @property
    def n_constraints(self):
        """
        The number of stabilizer constraints on valid codewords.
        """
        return len(self.group_generators)
        
    @property
    def nq_logical(self):
        """
        The number of logical qubits admitted by this code.
        """
        return self.nq - self.n_constraints

    @property
    def logical_ys(self):
        return pc.PauliList((ex * zed).mul_phase(1) for (ex,zed) in zip(self.logical_xs,self.logical_zs))     

    @property
    def distance(self):
        r"""
        The distance of this code, defined by :math:`\min\text{wt} \{
        P | P \in \text{N}(S) \backslash S \}`, where :math:`S` is the stabilizer group
        for this code.
        
        Warning: this property is currently very slow to compute.
        """
        return min(P.wt for P in self.normalizer_group(mod_s=True))
        
    @property
    def n_correctable(self):
        r"""
        The number of errors :math:`t` correctable by this code, defined by
        :math:`\left\lfloor \frac{d - 1}{2} \right\rfloor`, where :math:`d` is
        the distance of the code, given by the ``distance`` property.
        """
        return math.floor((self.distance - 1) / 2)
 
    ## GROUP ENUMERATION METHODS ##
 
    def stabilizer_group(self, coset_rep=None):
        r"""
        Iterator onto all elements of the stabilizer group :math:`S` describing
        this code, or onto a coset :math:`PS` of the stabilizer group.
        
        :param qecc.Pauli coset_rep: A Pauli operator :math:`P`, so that the
            iterated coset is :math:`PS`. If not specified, defaults to the
            identity.
        :yields: All elements of the coset :math:`PS` of the stabilizer
            group :math:`S`.
        """
        return self.group_generators.generated_group(coset_rep=coset_rep)
        
    def logical_pauli_group(self, incl_identity=True):
        r"""
        Iterator onto the group :math:`\text{N}(S) / S`, where :math:`S` is
        the stabilizer group describing this code. Each member of the group
        is specified by a coset representative drawn from the respective
        elements of :math:`\text{N}(S) / S`. These representatives are
        chosen to be the logical :math:`X` and :math:`Z` operators specified
        as properties of this instance.
        
        :param bool incl_identity: If ``False``, the identity coset :math:`S`
            is excluded from this iterator.
        :yields: A representative for each element of :math:`\text{N}(S) / S`.
        """
        return p.from_generators(self.logical_xs + self.logical_zs, incl_identity=incl_identity)
        
    def normalizer_group(self, mod_s=False):
        r"""
        Returns all elements of the normalizer of the stabilizer group. If
        ``mod_s`` is ``True``, returns the set :math:`N(S)\backslash S`.
        """
        for Pbar in self.logical_pauli_group(incl_identity=not mod_s):
            for normalizer_element in self.stabilizer_group(coset_rep=Pbar):
                yield normalizer_element
        
    ## EN/DE/TRANSCODING METHODS ##
        
    def encoding_cliffords(self):
        r"""
        Returns an iterator onto all Clifford operators that encode into this
        stabilizer code, starting from an input register such that the state to
        be encoded is a state of the first :math:`k` qubits, and such that the
        rest of the qubits in the input register are initialized to
        :math:`\left|0\right\rangle`.
        
        :yields: instances ``C`` of :class:`qecc.Clifford` such that
            ``C(q.StabilizerCode.unencoded_state(k, n - k))`` equals this code.
        """
        C = c.Clifford(
            self.logical_xs + ([Unspecified] * self.n_constraints),
            self.logical_zs + self.group_generators)
        return C.constraint_completions()
        
    def star_decoder(self, for_enc=None):
        r"""
        Returns a tuple of a decoding Clifford and a :class:`qecc.PauliList`
        specifying the recovery operation to perform as a function of the result
        of a :math:`Z^{\otimes{n - k}}` measurement on the ancilla register.
        
        For syndromes corresponding to errors of weight greater than the distance,
        the relevant element of the recovery list will be set to
        :obj:`qecc.Unspecified`.
        
        :param for_enc: If not ``None``, specifies to use a given Clifford
            operator as the encoder, instead of the first element yielded by
            :meth:`encoding_cliffords`.
        """
        def error_to_pauli(error):
            if error == p.I.as_clifford():
                return "I"
            if error == p.X.as_clifford():
                return "X"
            if error == p.Y.as_clifford():
                return "Y"
            if error == p.Z.as_clifford():
                return "Z"
        
        if for_enc is None:
            encoder = self.encoding_cliffords().next()
        else:
            encoder = for_enc
        decoder = encoder.inv()
        
        errors = pc.PauliList(p.eye_p(self.nq)) + pc.PauliList(p.paulis_by_weight(self.nq, self.n_correctable))
        
        syndrome_dict = defaultdict(lambda: Unspecified)
        syndrome_meas = [p.elem_gen(self.nq, idx, 'Z') for idx in range(self.nq_logical, self.nq)]
                
        for error in errors:
            effective_gate = decoder * error.as_clifford() * encoder
            # FIXME: the following line emulates measurement until we have a real
            #        measurement simulation method.
            syndrome = tuple([effective_gate(meas).ph / 2 for meas in syndrome_meas])
            
            recovery = "".join([
                # FIXME: the following is a broken hack to get the phases on the logical qubit register.
                error_to_pauli(c.Clifford([effective_gate.xout[idx][idx]], [effective_gate.zout[idx][idx]]))
                for idx in range(self.nq_logical)
            ])
            
            # For degenerate codes, the syndromes can collide, so long as we
            # correct the same way for each.
            if syndrome in syndrome_dict and syndrome_dict[syndrome] != recovery:
                raise RuntimeError('Syndrome {} has collided.'.format(syndrome))
                
            syndrome_dict[syndrome] = recovery
        
        recovery_list = pc.PauliList(syndrome_dict[syndrome] for syndrome in it.product(range(2), repeat=self.n_constraints))
        
        return decoder, recovery_list

    def minimize_distance_from(self, other, quiet=True):
        """
        Reorders the stabilizer group generators of this code to minimize
        the Hamming distance with the group generators of another code,
        using a greedy heuristic algorithm.
        """
        
        self_gens = self.group_generators
        other_gens = other.group_generators
        
        for idx_generator in range(len(self_gens)):
            min_hdist    = self.nq + 1 # Effectively infinite.
            min_wt       = self.nq + 1
            best_gen = None
            best_gen_decomp = ()
            
            for stab_elems in p.powerset(self_gens[idx_generator:]):
                if len(stab_elems) > 0:
                    stab_elem = reduce(op.mul, stab_elems)
                    hd = stab_elem.hamming_dist(other_gens[idx_generator])
                    
                    if hd <= min_hdist and stab_elem.wt <= min_wt and (hd < min_hdist or stab_elem.wt < min_wt):
                        min_hdist = hd
                        min_wt    = stab_elem.wt
                        best_gen  = stab_elem
                        best_gen_decomp = stab_elems
                    
            assert best_gen is not None, "Powerset iteration failed."
                    
            if best_gen in self_gens:
                # Swap so that it lies at the front.
                idx = self_gens.index(best_gen)
                if not quiet and idx != idx_generator:
                    print 'Swap move: {} <-> {}'.format(idx_generator, idx)
                self_gens[idx_generator], self_gens[idx] = self_gens[idx], self_gens[idx_generator]
                
            else:
                # Set the head element to best_gen, correcting the rest
                # as needed.
                if self_gens[idx_generator] in best_gen_decomp:
                    if not quiet:
                        print 'Set move: {}  =  {}'.format(idx_generator, best_gen)
                    self_gens[idx_generator] = best_gen
                else:
                    if not quiet:
                        print 'Mul move: {} *=  {}'.format(idx_generator, best_gen)
                    self_gens[idx_generator] *= best_gen
                    
        return self

    ## BLOCK CODE METHODS ##

    def block_logical_pauli(self, P):
        r"""
        Given a Pauli operator :math:`P` acting on :math:`k`, finds a Pauli
        operator :math:`\overline{P}` on :math:`nk` qubits that corresponds
        to the logical operator acting across :math:`k` blocks of this code.
        
        Note that this method is only supported for single logical qubit codes.
        """
        
        if self.nq_logical > 1:
            raise NotImplementedError("Mapping of logical Pauli operators is currently only supported for single-qubit codes.")
        
        # TODO: test that phases are handled correctly.
        
        # FIXME: cache this dictionary.
        replace_dict = {
            'I': p.eye_p(self.nq),
            'X': self.logical_xs[0],
            'Y': (self.logical_xs[0] * self.logical_zs[0]).mul_phase(1),
            'Z': self.logical_zs[0]
        }
        
        # FIXME: using eye_p(0) is a hack.
        return reduce(op.and_, 
                (replace_dict[sq_op] for sq_op in P.op),
                p.eye_p(0))

    def mabuchi_operators(stab_code,known_operators=pc.PauliList(),wt2=False):
        """
        Given a stabilizer code, this function yields Pauli operators
        that satisfy all possible commutation relations with the
        stabilizer group generators. This is a strict superset of the 
        detectable errors. To increase the speed of this function, 
        the keyword argument "known_operators" allows the input of the
        detectable errors themselves. 
        """
        n_constraints=stab_code.n_constraints
        enc = stab_code.encoding_cliffords().next()
        dec = enc.inv()
        synd_meas = [p.elem_gen(stab_code.nq, idx, kind) for idx, kind in zip(range(1,stab_code.nq), 'Z'*n_constraints)]
        known_syndromes=[]
        for op in known_operators:
            #Calculate syndromes of known operators, and append to
            #known syndromes
            eff = dec * op.as_clifford() * enc
            known_syndromes.append([eff.conjugate_pauli(meas).ph / 2 for meas in synd_meas])
        code_gens=stab_code.group_generators
        for bits in it.product(range(2),repeat=n_constraints):
            if bits in known_syndromes:
                #If the operator is known, output:
                yield known_operators[known_syndromes.index(bits)]
            else:
                #If the operator is unknown, find a satisfactory operator
                #by solving commutation constraints. 
                anti_coms=list(it.compress(code_gens,bits))
                coms=list(it.compress(code_gens,[1-bit for bit in bits]))
                if wt2:
                    yield min(solve_commutation_constraints(coms,anti_coms),key=lambda j: j.wt)
                else:
                    yield min(solve_commutation_constraints(coms,anti_coms),key=lambda j: j.wt)
                    
    ## OPERATORS ##
    
    def __and__(self, other):
    
        if not isinstance(other, StabilizerCode):
            return NotImplemented
        
        return StabilizerCode(
            (self.group_generators & p.eye_p(other.nq)) +
            (p.eye_p(self.nq) & other.group_generators),
            
            (self.logical_xs & p.eye_p(other.nq)) +
            (p.eye_p(self.nq) & other.logical_xs),
            
            (self.logical_zs & p.eye_p(other.nq)) +
            (p.eye_p(self.nq) & other.logical_zs),
        )
        
    ## CONCATENATION ##
        
    def concatenate(self,other):
        r"""
        Returns the stabilizer for a concatenated code, given the 
        stabilizers for two codes. At this point, it only works for two
        :math:`k=1` codes.
        """
        
        if self.nq_logical > 1 or other.nq_logical > 1:
            raise NotImplementedError("Concatenation is currently only supported for single-qubit codes.")
        
        nq_self = self.nq
        nq_other = other.nq
        nq_new = nq_self * nq_other
        
        # To obtain the new generators, we must apply the stabilizer generators
        # to each block of the inner code (self), as well as the stabilizer
        # generators of the outer code (other), using the inner logical Paulis
        # for the outer stabilizer generators.
        
        # Making the stabilizer generators from the inner (L0) code is straight-
        # forward: we repeat the code other.nq times, once on each block of the
        # outer code. We use that PauliList supports tensor products.
        new_generators = sum(
            (
                p.eye_p(nq_self * k) & self.group_generators & p.eye_p(nq_self * (nq_other - k - 1))
                for k in range(nq_other)
            ),
            pc.PauliList())
                
        # Each of the stabilizer generators due to the outer (L1) code can be
        # found by computing the block-logical operator across multiple L0
        # blocks, as implemented by StabilizerCode.block_logical_pauli.
        new_generators += map(self.block_logical_pauli, other.group_generators)
            
        # In the same way, the logical operators are also found by mapping L1
        # operators onto L0 qubits.
        
        # This completes the definition of the concatenated code, and so we are
        # done.
        
        return StabilizerCode(new_generators,
            logical_xs=map(self.block_logical_pauli, other.logical_xs),
            logical_zs=map(self.block_logical_pauli, other.logical_zs)
        )

    ## COMMON CODES ##

    @staticmethod
    def ancilla_register(nq=1):
        r"""
        Creates an instance of :class:`qecc.StabilizerCode` representing an
        ancilla register of ``nq`` qubits, initialized in the state
        :math:`\left|0\right\rangle^{\otimes \text{nq}}`.
        
        :rtype: qecc.StabilizerCode
        """
        return StabilizerCode(
            p.elem_gens(nq)[1],
            [], []
        )

    @staticmethod
    def unencoded_state(nq_logical=1, nq_ancilla=0):
        """
        Creates an instance of :class:`qecc.StabilizerCode` representing an
        unencoded register of ``nq_logical`` qubits tensored with an ancilla
        register of ``nq_ancilla`` qubits.
        
        :param int nq_logical: Number of qubits to 
        :rtype: qecc.StabilizerCode
        """
        return (
            StabilizerCode([], *p.elem_gens(nq_logical)) &
            StabilizerCode.ancilla_register(nq_ancilla)
        )

    @staticmethod
    def flip_code(n_correctable, stab_kind=p.Z):
        """
        Creates an instance of :class:`qecc.StabilizerCode` representing a
        code that protects against weight-``n_correctable`` flip errors of a
        single kind.
        
        This method generalizes the bit-flip and phase-flip codes, corresponding
        to ``stab_kind=qecc.Z`` and ``stab_kind=qecc.X``, respectively.
        
        :param int n_correctable: Maximum weight of the errors that can be
            corrected by this code.
        :param qecc.Pauli stab_kind: Single-qubit Pauli operator specifying
            which kind of operators to use for the new stabilizer code.
        :rtype: qecc.StabilizerCode
        """
        nq = 2 * n_correctable + 1
        stab_kind = p.ensure_pauli(stab_kind)
        if len(stab_kind) != 1:
            raise ValueError("stab_kind must be single-qubit.")
        
        return StabilizerCode(
            [p.eye_p(j) & stab_kind & stab_kind & p.eye_p(nq-j-2) for j in range(nq-1)],
            ['X'*nq], ['Z'*nq]
        )

    @staticmethod
    def bit_flip_code(n_correctable):
        """
        Creates an instance of :class:`qecc.StabilizerCode` representing a
        code that protects against weight-``n_correctable`` bit-flip errors.
        
        :param int n_correctable: Maximum weight of the bit-flip errors that can
            be corrected by this code.
        :rtype: qecc.StabilizerCode
        """
        return StabilizerCode.flip_code(n_correctable, stab_kind=p.Z)
    @staticmethod
    def phase_flip_code(n_correctable):
        """
        Creates an instance of :class:`qecc.StabilizerCode` representing a
        code that protects against weight-``n_correctable`` phase-flip errors.
        
        :param int n_correctable: Maximum weight of the phase-flip errors that
            can be corrected by this code.
        :rtype: qecc.StabilizerCode
        """
        return StabilizerCode.flip_code(n_correctable, stab_kind=p.X)

    @staticmethod
    def perfect_5q_code():
        """
        Creates an instance of :class:`qecc.StabilizerCode` representing the
        5-qubit perfect code.
        
        :rtype: qecc.StabilizerCode
        """
        return StabilizerCode(
            [
                'XZZXI',
                'IXZZX',
                'XIXZZ',
                'ZXIXZ'
            ],
            ['XXXXX'], ['ZZZZZ']
        )

    @staticmethod
    def steane_code():
        """
        Creates an instance of :class:`qecc.StabilizerCode` representing the
        7-qubit Steane code.
        
        :rtype: qecc.StabilizerCode
        """
        return StabilizerCode(
            [
                'XXXXIII',
                'XXIIXXI',
                'XIXIXIX',
                'ZZZZIII',
                'ZZIIZZI',
                'ZIZIZIZ'                
            ],
            ['XXXXXXX'], ['ZZZZZZZ']
        )

    @staticmethod
    def shor_code():
        """
        Creates an instance of :class:`qecc.StabilizerCode` representing the
        9-qubit Shor code.
        
        :rtype: qecc.StabilizerCode
        """        
        return StabilizerCode.bit_flip_code(1).concatenate(StabilizerCode.phase_flip_code(1))        

    @staticmethod
    def css_code(C1, C2):
        """
        Not yet implemented.
        """
        raise NotImplementedError("Not yet implemented.")

    @staticmethod
    def reed_muller_code(r,t):
        """
        Not yet implemented.
        """
        raise NotImplementedError("Coming Soon: Reed-Muller Codes")
        
    @staticmethod
    def reed_solomon_code(r,t):
        """

        Not yet implemented.

        """
        raise NotImplementedError("Coming Soon: Reed-Solomon Codes")
