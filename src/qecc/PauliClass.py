#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# PauliClass.py: Implementation of qecc.Pauli and related utility functions.
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

from itertools import product, chain, permutations, combinations, ifilter, ifilterfalse, imap, starmap, izip
from copy import copy
import bsf
from operator import mul
import pred

from paulicollections import PauliList
from unitary_reps import pauli_as_unitary

from singletons import Unspecified

import  CliffordClass as cc

## ALL ##

__all__ = [
    'Pauli',
    'ensure_pauli', 'com', 'pauli_group', 'from_generators',
    'is_in_normalizer', 'elem_gen', 'elem_gens', 'eye_p', 'ns_mod_s',
    'pad', 'mutually_commuting_sets',
    'clifford_bottoms',
    'paulis_by_weight'
    ]
        
## CONSTANTS ##

VALID_OPS = ['I', 'X', 'Y', 'Z']
__all__ += VALID_OPS
VALID_PHS = range(4)

MULT_TABLE = {
    ('I', 'I'): (0, 'I'), ('I', 'X'): (0, 'X'), ('I', 'Y'): (0, 'Y'), ('I', 'Z'): (0, 'Z'),
    ('X', 'I'): (0, 'X'), ('X', 'X'): (0, 'I'), ('X', 'Y'): (1, 'Z'), ('X', 'Z'): (3, 'Y'),
    ('Y', 'I'): (0, 'Y'), ('Y', 'X'): (3, 'Z'), ('Y', 'Y'): (0, 'I'), ('Y', 'Z'): (1, 'X'),
    ('Z', 'I'): (0, 'Z'), ('Z', 'X'): (1, 'Y'), ('Z', 'Y'): (3, 'X'), ('Z', 'Z'): (0, 'I')
}

## CLASSES ##

class Pauli(object):
    r"""
    Class representing an element of the Pauli group on :math:`n` qubits.
    
    :param operator: String of I's, X's, Y's and Z's.
    :type operator: str
    :param phase: A phase input as an integer from 0 to 3, interpreted as
        :math:`i^{\mathrm{phase}}`.
    :type phase: int
    """

    def __init__(self, operator, phase=0):

        # NB: operator is a string which contains I, X, Y, Z. We throw an error if it isn't.

        for paulicounter in range(len(operator)):
            if operator[paulicounter] not in VALID_OPS:
                raise ValueError("Input operators I, X, Y or Z.")

        # NB: phase will be defined as an integer from 0 to 3. If it's not integral, we throw an error.
        if not isinstance(phase, int):
            raise ValueError("Input phase must be an integer, preferably 0 to 3.")

        #If it's not in the range, that's cool, we mod 'til it is.
        if not( phase > -1 and phase < 4):
            phase= phase % 4
            
        self.op = operator
        self.ph = phase
        
    def __hash__(self):
        # We need a hash function to store Paulis as dict keys or in sets.
        return hash((self.op, self.ph))
        
    def __len__(self):
        """
        Yields the number of qubits on which the Pauli ``self`` acts.
        """
        return len(self.op)
                
    def __mul__(self, other):
        """
        Multiplies two Paulis, ``self`` and ``other`` symbolically, using tabulated single-Pauli multiplication.
        """
        if not isinstance(other,Pauli):
            return NotImplemented 

        p1 = self
        p2 = other

        if not(len(p1)==len(p2)):
            raise ValueError("These Paulis are not the same length")
            
        newP = Pauli('', p1.ph + p2.ph)
        for paulicounter in range(len(p1)):
            ph, op = MULT_TABLE[(p1.op[paulicounter], p2.op[paulicounter])]
            newP.op=newP.op+op
            newP.ph=newP.ph+ph
            
        newP.ph = newP.ph % 4
            
        return newP

    def __pow__(self,num):
        if num % 2 == 0:
            return eye_p(len(self))
        elif num % 2 == 1:
            return self
        else:
            raise ValueError("Paulis can only be exponentiated with integers")
            
    def __repr__(self):
        """
        Representation for Paulis, printing Paulis in the format ``i^{phase} {op}``.
        """
        return "i^{k} {op}".format(k=self.ph, op=self.op)

    def __str__(self):
        SPARSE_NQ     = 5
        SPARSE_THRESH = 0.3
        
        if len(self) <= SPARSE_NQ or self.wt / len(self) > SPARSE_THRESH:
            return repr(self)
        else:
            return self.str_sparse()
            
    def __getitem__(self, idxs):
        return Pauli(self.op[idxs], phase=self.ph)
            
    def str_sparse(self, incl_ph=True):
        return ("i^{} ".format(self.ph) if incl_ph else "") + (" ".join(
            "{}[{}]".format(P, idx) for idx, P in enumerate(self.op) if P != "I"
        ) if self.wt > 0 else "I")

    def __neg__(self):
        """
        Negates (multiplying by :math:`-1`) a Pauli.
        """
        return copy(self).mul_phase(2)

    def tens(self, other):
        r"""
        Concatenates the op strings of two Paulis, and multiplies their phases, to produce
        the Kronecker product of the two.
        
        :param other: Pauli operator :math:`Q` to be tensored with this instance.
        :type other: qecc.Pauli
        :returns: An instance representing :math:`P\otimes Q`, where :math:`P`
            is the Pauli operator represented by this instance.
        """
        return Pauli(self.op + other.op, self.ph + other.ph)

    def __and__(self,other):
        if not isinstance(other, Pauli):
            return NotImplemented
        return self.tens(other)
        
    def set_phase(self, ph=0):
        self.ph = ph
        return self
        
    def mul_phase(self, ph):
        r"""
        Increments the phase of this Pauli by :math:`i^{\mathrm{ph}}`.
        
        :param ph: Amount the phase is to be incremented by.
        :type ph: int
        :returns: This instance.
        """
        self.ph = (self.ph + ph) % 4
        return self

    def as_gens(self):
        """
        Expresses an input Pauli in terms of the elementary generators :math:`X_j` and :math:`Z_j`, stripping off phases.
        
        :rtype: list of :class:`qecc.Pauli` instances.
        """
        nq=len(self)
        #Initialize generator list to an empty list:
        pauli_gens=[]
        for idx in range(nq):
        #Since phases are being stripped, Y = XZ:
            if self.op[idx]=='Y':
                pauli_gens.append(elem_gen(nq,idx,'X'))
                pauli_gens.append(elem_gen(nq,idx,'Z'))
            elif self.op[idx]=='X':
                pauli_gens.append(elem_gen(nq,idx,'X'))
            elif self.op[idx]=='Z':
                pauli_gens.append(elem_gen(nq,idx,'Z'))
        return pauli_gens

    def as_bsv(self):
        """
        Converts the given Pauli to a binary symplectic vector, discarding phase
        information.
        
        :return: A binary symplectic vector representing this Pauli operator.
        :rtype:  BinarySymplecticVector
        """
        nq = len(self)
        exes=''
        zeds=''
        for idx_q in range(nq):
            if self.op[idx_q]=='X':
                exes=exes+'1'
                zeds=zeds+'0'
            elif self.op[idx_q]=='Y':
                exes=exes+'1'
                zeds=zeds+'1'
            elif self.op[idx_q]=='Z':
                exes=exes+'0'
                zeds=zeds+'1'
            elif self.op[idx_q]=='I':
                exes=exes+'0'
                zeds=zeds+'0' 
        return bsf.BinarySymplecticVector(exes,zeds)

    def __eq__(self,other):
        """
        Tests if two input Paulis, :obj:`self` and :obj:`other`, are equal.

        :rtype: bool
        """
        if self.op==other.op and self.ph==other.ph:
            return True
        else:
            return False

    def __ne__(self,other):
        """
        Tests if two input Paulis, :obj:`self` and :obj:`other`, are not equal.

        :rtype: bool
        """
        if self.op==other.op and self.ph==other.ph:
            return False
        else:
            return True

    def __len__(self):
        """
        Yields the number of qubits on which the Pauli ``self`` acts.
        """
        return len(self.op)
           
    def as_unitary(self):
        """
        Returns a :class:`numpy.ndarray` containing a unitary matrix
        representation of this Pauli operator.
        
        Raises a :obj:`RuntimeError` if NumPy cannot be imported.
        """
        return pauli_as_unitary(self)
              
    def as_clifford(self):
        """
        Returns a :class:`qecc.Clifford` representing conjugation by this Pauli
        operator.
        """
        Xs, Zs = elem_gens(len(self))
        for P in chain(Xs, Zs):
            if com(self, P) == 1:
                P.mul_phase(2)
                
        return cc.Clifford(Xs, Zs)
                
    @property
    def wt(self):
        """
        Measures the weight of a given Pauli.
        
        :rtype: int (between 0 and the number of qubits on which the Pauli is defined)
        :returns: The number of qubits on which the represented Pauli operator
            is is supported.
        """
        return len([op for op in self.op if op != 'I'])
        
    def ct(self):
        """
        The conjugate transpose of this Pauli operator.

        :rtype: an instance of the :class:`qecc.Pauli` class.
        """    
        return Pauli(self.op, 4-self.ph)
        
    def centralizer_gens(self, group_gens=None):
        r"""
        Returns the generators of the centralizer group :math:`\mathrm{C}(P)`,
        where :math:`P` is the Pauli operator represented by this instance.
        If ``group_gens`` is specified, :math:`\mathrm{C}(P)` is taken to be
        a subgroup of the group :math:`G = \langle G_1, \dots, G_k\rangle`,
        where :math:`G_i` is the :math:`i^{\text{th}}` element of
        ``group_gens``.
        
        :param group_gens: Either ``None`` or a list of generators :math:`G_i`.
            If not ``None``, the returned centralizer :math:`\mathrm{C}(P)` is
            a subgroup of the group :math:`\langle G_i \rangle_{i=1}^k`.
        :type group_gens: list of qecc.Pauli instances
        :returns: A list of elements :math:`P_i` of the Pauli group such that
            :math:`\mathrm{C}(P) = \langle P_i \rangle_{i=1}^{n}`, where
            :math:`n` is the number of unique generators of the centralizer.
        """
        if group_gens is None:
            Xs, Zs = elem_gens(len(self))
            group_gens = Xs + Zs
        
        if com(self, group_gens[0]) == 0:
            # That generator commutes, and so we pass it along
            # unmodified.
            return PauliList(group_gens[0], *self.centralizer_gens(group_gens[1:]))
        else:
            # That generator anticommutes, and so we must modify it by
            # multiplication with another anticommuting generator, if one
            # exists.
            found = False
            for idx in range(1, len(group_gens)):
                if com(self, group_gens[idx]) == 1:
                    found = True
                    g_prime = group_gens[idx] * group_gens[0]
                    assert com(self, g_prime) == 0
                    return PauliList(g_prime, *self.centralizer_gens(group_gens[1:]))
            if not found:
                # Generator 0 anticommuted, and so we know two things
                # from our search:
                #     - All other generators commute.
                #     - The anticommuting generator (0) has no match, and must
                #       be excluded.
                return PauliList(*group_gens[1:])

## MORE CONSTANTS ##

I = Pauli('I')
X = Pauli('X')
Y = Pauli('Y')
Z = Pauli('Z')

## FUNCTIONS ##
       
def ensure_pauli(P):
    """
    Given either a string or an instance of :class:`qecc.Pauli`, produces an
    instance of :class:`qecc.Pauli`.
    """
    # TODO: provide a decorator that applies this function to another function's
    #       arguments.
    return P if isinstance(P, Pauli) or P is Unspecified else Pauli(P)
        
# Stolen from itertools cookbook.
def powerset(iterable):
    """
    Sub-function, stolen from itertools cookbook. 

    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return chain.from_iterable(combinations(s, r) for r in range(len(s)+1))
        
def paulis_by_weight(nq, wt):
    def error_by_idxs(idxs, err_string):
        return reduce(mul, starmap(elem_gen, zip([nq]*len(idxs), idxs, err_string)))
        
    if wt == 0:
        return eye_p(nq)
    elif wt == 1:
        return iter(elem_gen(nq, idx, kind) for kind in ['X', 'Y', 'Z'] for idx in range(nq))
    else:
        return chain(
            paulis_by_weight(nq, wt - 1),
            (error_by_idxs(idxs, err_string) for err_string in product('XYZ', repeat=wt) for idxs in combinations(range(nq), wt))
        )
        
        
            
def com(P, Q):
    r"""
    Given two elements *P* and *Q* of a Pauli group, returns 0 if :math:`[P, Q] = 0`
    and returns 1 if :math:`\{P, Q\} = 0`.
    
    :param P: Representation of :math:`P`.
    :type P: qecc.Pauli
    :param Q: Representation of :math:`Q`.
    :type Q: qecc.Pauli
    
    :returns: :math:`c(P, Q)`.
    :rtype: int
    """
    ph1, ph2 = (P*Q).ph, (Q*P).ph
    return 0 if ph1 == ph2 else 1

def pauli_group(nq):
    """
    Generates an iterator onto the Pauli group of :math:`n` qubits,
    where :math:`n` is given as the argument *nq*.
    
    :param nq: The number of qubits acted upon by the returned Pauli group.
    :type nq: int
    
    :returns: An iterator such that ``list(pauli_group(nq))`` produces a list of
        all possible Pauli operators on ``nq`` qubits.
    """
    for op in product(VALID_OPS, repeat=nq):
        yield Pauli("".join(op), 0)
        
def from_generators(gens, coset_rep=None, incl_identity=True):
    """
    Given a list of generators ``gens``, yields an iterator
    onto the group generated by products of elements from
    ``gens``.
    
    If ``coset_rep`` is specified, returns the coset of the group generated by
    ``gens`` represented by ``coset_rep``.
    """
    
    if coset_rep is None:
        coset_rep = eye_p(len(gens[0]))
    
    for prod in powerset(gens):
        if len(prod)>0:
            yield reduce(lambda P, Q: P*Q, prod, coset_rep)
        elif incl_identity:
            yield coset_rep
        # Ignore this element of the powerset if incl_identity is False.
                
        
def is_in_normalizer(pauli, stab):
    """
    Given an element ``pauli`` of a Pauli group and the generators
    ``stab`` of a stabilizer group :math:`S`, returns True if and only if
    ``pauli`` is in the normalizer :math:`N(S)`.
    """
    stab_group = from_generators(stab)
    com_vec = [com(pauli, stab_elem) for stab_elem in stab_group]
    
    return all(map(lambda x: x == 0, com_vec))

def elem_gen(nq, q, op):
    """
    Produces a weight-one Pauli on a specific qubit ``q`` out of ``nq``.
    
    :param int nq: Number of qubits upon which the returned Pauli acts.
    :param int q: Index of the qubit upon which the returned Pauli acts
        non-trivially.
    :param str op: Either ``"I"``, ``"X"``, ``"Y"`` or ``"Z"``, indicating which
        generator to return.
    :returns: An ``nq``-qubit operator that acts as ``op`` on qubit ``q``, and
        that acts as identity on all other qubits.
    :rtype: qecc.Pauli
    """
    if op in VALID_OPS:
        return Pauli('I'*q + op + 'I'*(nq-q-1))
    else:
        raise ValueError('Generators cannot be selected outside I, X, Y, Z.')

def elem_gens(nq):
    """
    Produces all weight-one :math:`X` and :math:`Z` operators on `nq` qubits.
    For example,
    
    >>> import qecc as q
    >>> Xgens, Zgens = q.elem_gens(2)
    >>> print Xgens[1]
    i^0 IX
    
    :param int nq: Number of qubits for each returned operator.
    
    :returns: a tuple of two lists, containing :math:`X` and :math:`Z`
        generators, respectively.
    """
    return tuple(PauliList(*[elem_gen(nq,idx,P) for idx in range(nq)]) for P in ['X','Z'])

def eye_p(nq):
    """
    Given a number of qubits, returns the identity Pauli on that many qubits.

    :param int nq: Number of qubits upon which the returned Pauli acts.
    :rtype: qecc.Pauli
    :returns: A Pauli operator acting as the identity on each of ``nq`` qubits.
    """
    return Pauli('I'*nq)

def ns_mod_s(*stab_gens):
    r"""
    Given the generators of a stabilizer group :math:`S`, returns an iterator
    that yields all elements of the set :math:`\text{N}(S)\\S`.
    
    :param qecc.Pauli stab_gens: Generators of :math:`S`.
    :returns: An iterator such that ``list(ns_mod_s(*stab_gens))`` contains
        each element of :math:`\text{N}(S)\\S`.
    """
    nq = len(stab_gens[0])
    
    return ifilter(
        pred.commutes_with(*stab_gens) & ~pred.in_group_generated_by(*stab_gens),
        pauli_group(nq)
        )

def mutually_commuting_sets(n_elems, n_bits=None, group_gens=None, exclude=None):
    r"""
    Yields an iterator onto tuples representing mutually commuting sets of
    ``n_elems`` independent Pauli operators, excluding the identity Pauli.
    
    :param int n_elems: The number of mutually commuting Pauli operators to
        include in each tuple.
    :param int n_bits: The number of qubits on which each Pauli operator
        considered by this iterator acts. If ``None``, defaults to the number of
        qubits on which the first element of ``group_gens`` acts.
    :param group_gens: The generators of the group in which to search for
        mutually commuting Pauli operators. Defaults to the elementary
        generators of the Pauli group on ``n_bits`` qubits.
    :type group_gens: ``None`` or a sequence of :class:`qecc.Pauli` instances
    :param exclude: If not ``None``, the iterator will omit from its search any
        operators in the group generated by ``exclude``.
    :type exclude: ``None`` or a sequence of :class:`qecc.Pauli` instances
    """

    if n_elems == 0:
        yield ()
        return
    
    if exclude is None:
        if n_bits is not None:
            exclude = [eye_p(n_bits)]
        else:
            if group_gens is not None and len(group_gens) > 0:
                exclude = [eye_p(len(group_gens[0]))]

    assert len(exclude) > 0

    if group_gens is None and n_bits is None:
        raise ValueError('Either a generating set or a number of qubits must be specified.')
    
    if group_gens is None:
        group_gens = add(*elem_gens(n_bits))
    
    assert len(group_gens) > 0
    
    for P in ifilterfalse(lambda q: q.wt==0 or q in from_generators(exclude),from_generators(group_gens)):
        P = Pauli(P.op, phase=0)
        if n_elems==1:
            yield (P,)
        else:
            c_gens=P.centralizer_gens(group_gens)
            for S in mutually_commuting_sets(n_elems-1,group_gens=c_gens, exclude=exclude+[P]):
                yield (P,)+S
        
    
    
def pad(setP, n_eb=0, lower_right=None):
    r"""
    """
    
    # Ensure that we have lists, and make a copy.
    setP= list(setP)
    
    len_P=len(setP)
    n_P=len(setP[0])

    if n_eb == 0 and lower_right is None or len(lower_right)==0:
        return setP
    elif len(lower_right) != 0:
        n_eb=len(lower_right[0])
            
    setout=[]
    for pauli in setP:
        setout.append(Pauli(pauli.op+'I'*n_eb))
    if lower_right is None:
        for jj in n_eb:
            setout.append(eye_p(n_P+n_eb))
    else:
        for opmo in lower_right:
            setout.append(eye_p(n_P)&opmo)
    return setout    

def clifford_bottoms(c_top):
    r"""
    This function yields the next set of nq paulis that mutually 
    commute, and anti-commute with selected elements from c_top. 
    """
    nq=len(c_top)
    possible_zs=[]
    for jj in range(nq):
        applicable_centralizer_gens=PauliList(*(c_top[:jj]+c_top[jj+1:])).centralizer_gens()
        possible_zs.append(filter(lambda a:com(a,c_top[jj])==1,from_generators(applicable_centralizer_gens)))
    for possible_set in product(*possible_zs):
        if all(imap(lambda twolist: pred.commutes_with(twolist[0])(twolist[1]),combinations(possible_set,2))):
            yield possible_set
## MAIN ##

if __name__ == "__main__":
    P = Pauli('XYZ', 0)
    Q = Pauli('YXI', 0)
    R = Pauli('XXI', 0)
    
    print list(pauli_group(2))
    
    print list(from_generators([P, Q, R]))
    print is_in_normalizer(Pauli('ZZX', 0), [P, Q, R])
    
