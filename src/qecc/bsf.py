#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# bsf.py: Implementation of binary symplectic form.
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

import itertools
import string
from exceptions import *

from numpy import s_, array, nonzero, logical_xor, bitwise_xor, bitwise_not, logical_and, bitwise_and, logical_not, all, matrix, hstack, vstack, zeros, eye, dot, empty

from PauliClass import *
from CliffordClass import *

import utils as u

import bsf_decomp
reload(bsf_decomp)

import circuit

## ALL ##

__all__ = [
    'BinarySymplecticVector',
#    'bitstring_to_letterstring',
    'parity', 'bitwise_inner_product', 'all_pauli_bsvs', 'constrained_set',
    'commute', 'xz_switch',
    
    'BinarySymplecticMatrix',
    'is_bsm_valid', 'bsmzeros', 'array_to_pauli', 'directsum'
]

## CLASSES ##

VALID_BITS=range(2)

class BinarySymplecticVector(object):
    """
    Encapsulates a binary symplectic vector representing an element of the Pauli
    group on :math:`n` qubits.
    
    A new :class:`BinarySymplecticVector` can be constructed using either a
    single NumPy array containing both the :math:`X` and :math:`Z` parts of the
    binary symplectic vector. Alternatively, a new vector can be instantiated
    using two NumPy arrays. For example, the following two invocations are
    equivalent:
    
    >>> import qecc
    >>> import numpy as np
    >>> bsv = qecc.BinarySymplecticVector(np.array([1, 0, 0, 0, 0, 0]))
    >>> bsv = qecc.BinarySymplecticVector(np.array([1, 0, 0]), np.array([0, 0, 0]))
    
    The :obj:`len` of a :class:`BinarySymplecticVector` is defined as the number
    of qubits upon which the represented Pauli operator acts, and is thus half
    of the length of a single array containing the same data.
    """

    def __init__(self,*args):
        if len(args) == 1:
            nq = len(args[0])/2
            self._x = array(args[0][0:nq], dtype=int)
            self._z = array(args[0][nq:2*nq], dtype=int)
        elif len(args) == 2:
            xstring,zstring = args
            self._x=array(list(xstring), dtype='int')
            self._z=array(list(zstring), dtype='int')
        else:
            raise ValueError('Wrong number of args.')

    ## MAGIC METHODS ##
    
    def __len__(self):
        return len(self._x)

    def __repr__(self):
        return "( {ex} | {zed} )".format(ex=" ".join(map(str, self._x)), zed=" ".join(map(str, self._z)))
        
    ## PROPERTIES ##

    @property
    def x(self):
        """
        Array containing the :math:`X` part of the binary symplectic vector.
        
        :rtype: :class:`numpy.ndarray`, shape ``(2 * nq, )``.
        """
        return self._x
    
    @property
    def z(self):
        """
        Array containing the :math:`Z` part of the binary symplectic vector.
        
        :rtype: :class:`numpy.ndarray`, shape ``(nq, )``.
        """
        return self._z

    ## OTHER METHODS ##

    def copy(self):
        """
        Returns a copy of the binary symplectic vector such that mutations of
        the copy do not affect this instance. For more details, see the
        :meth:`numpy.ndarray.copy` method.
        """
        return BinarySymplecticVector(self._x.copy(), self._z.copy())

    def as_pauli(self):
        """
        Returns an instance of :class:`qecc.Pauli` representing the same Pauli
        operator as this vector. Note that phase information is not preserved
        by the binary symplectic representation of the Pauli group, and so
        ``P.as_bsv().as_pauli()`` need not equal ``P``.
        """
        exes=bitstring_to_letterstring(self._x,'X')
        zeds=bitstring_to_letterstring(self._z,'Z')
        zedley=Pauli(zeds)
        exeley=Pauli(exes)
        assert type(zedley) is type(exeley)
        pauli=zedley*exeley
        return Pauli(pauli.op,0)
        
    def bsip(self,other):
        r"""
        Returns the binary symplectic inner product :math:`u \odot v` of this
        vector with another vector. Letting :math:`u = (a | b)` and
        :math:`v = (c | d)`, :math:`u\odot v = a \cdot c + b \cdot d`.
        """
        return int(not(commute(self,other)))
        
#FUNCTIONS FOR BSV CLASS     

def bitstring_to_letterstring(bitstring,letter):
    outstring=''
    for idx in range(len(bitstring)):
        if bitstring[idx]==0:
            outstring=outstring+('I')
        elif bitstring[idx]==1:
            outstring=outstring+(letter)
    return outstring

def parity(bitarray):
    """returns True if bitarray is of odd parity, False if it is of even parity."""
    return reduce(logical_xor,bitarray)

def bitwise_inner_product(v1,v2):
    return parity(bitwise_and(v1,v2))

def all_pauli_bsvs(nq):
    r"""
    For a given number of qubits ``nq``, returns an iterator that yields the
    binary symplectic representations of each element of the Pauli group
    :math:`\mathcal{P}_n`.
    """
    for idx_x in itertools.product([0,1],repeat=nq):
        for idx_z in itertools.product([0,1],repeat=nq):
            yield(BinarySymplecticVector(idx_x,idx_z))

def constrained_set(pauli_array_input,logical_array_input):
    r"""
    Given a set of constraints of the form :math:`P_i \odot Q = b_i`, with
    each :math:`P_i` a Pauli operator and each :math:`b_i` a bit, yields an
    iterator onto Pauli operators :math:`Q` such that all constraints are
    satisfied.
    
    :type pauli_array_input: :obj:`list` of :class:`qecc.Pauli` instances.
    :param pauli_array_input: Constraint operators :math:`P_i`.
    :type logical_array_input: :class:`numpy.ndarray` of `dtype=int` and shape
        ``(len(pauli_array_input), )``.
    :param logical_array_input: Constraint values :math:`b_i`.
    """
    
    nq=len(pauli_array_input[0].x)
    #output_array=[]
    logical_bookkeeping_array=[]
    for current_pauli in all_pauli_bsvs(nq):
        logical_bookkeeping_array = map(current_pauli.bsip, pauli_array_input)
        if logical_bookkeeping_array == logical_array_input:
            # output_array.append(current_pauli)
            yield current_pauli
    # return output_array

def commute(bsv1,bsv2):
    """returns True if bsv1 and bsv2 commute by evaluating the symplectic inner product."""
    return logical_not(logical_xor(bitwise_inner_product(bsv1.x,bsv2.z),bitwise_inner_product(bsv1.z,bsv2.x)))

def xz_switch(bsv):
    """
    Given a :class:`qecc.BinarySymplecticVector`, returns a new vector whose
    :math:`X` and :math:`Z` parts have been swapped.
    """
    return BinarySymplecticVector(bsv.z,bsv.x)

## BINARY SYMPLECTIC MATRIX CLASS ##

class BinarySymplecticMatrix(object):
    
    def __init__(self,*args):
        if len(args)==4:
            self._arr=vstack((hstack((args[0],args[1])),hstack((args[2],args[3]))))
        elif len(args)==1:
            self._arr=args[0]
        else:
            raise ValueError("Either one or four arguments must be provided.")
            
    ## PROPERTIES ##
    
    @property
    def nq(self):
        """
        Returns the number of qubits that the binary symplectic matrix acts
        upon.
        """
        return len(self._arr)/2
            
    @property
    def xx(self): 
        nq = self.nq
        return self._arr[0:nq,0:nq]
    @property
    def xz(self):
        nq = self.nq
        return self._arr[0:nq,nq:2*nq]
    @property
    def zx(self):
        nq = self.nq
        return self._arr[nq:2*nq,0:nq]
    @property
    def zz(self):
        nq = self.nq
        return self._arr[nq:2*nq,nq:2*nq]

    @property
    def xc(self):
        nq = self.nq
        return self._arr[:, 0:nq]
    @xc.setter
    def xc(self, newval):
        self._arr[:, 0:self.nq] = newval
    @property
    def zc(self):
        nq = self.nq
        return self._arr[:, nq:2*nq]
    @zc.setter
    def zc(self, newval):
        self._arr[:, self.nq:2*self.nq] = newval

    @property
    def xr(self):
        nq = self.nq
        return self._arr[0:nq, :]
    @xr.setter
    def xr(self, newval):
        self._arr[0:self.nq, :] = newval
    @property
    def zr(self):
        nq = self.nq
        return self._arr[nq:2*nq, :]
    @zr.setter
    def zr(self, newval):
        self._arr[self.nq:2*self.nq, :] = newval


    @xx.setter
    def xx(self, thing):
        nq = self.nq
        self._arr[0:nq,0:nq]=thing
    @xz.setter
    def xz(self, thing):
        nq = self.nq
        self._arr[0:nq,nq:2*nq]=thing
    @zx.setter
    def zx(self, thing):
        nq = self.nq
        self._arr[nq:2*nq,0:nq]=thing
    @zz.setter
    def zz(self, thing):
        nq = self.nq
        self._arr[nq:2*nq,nq:2*nq]=thing

    ## MAGIC METHODS ##
    
    def __getitem__(self, sliceobj):
        return self._arr.__getitem__(sliceobj)

    def __setitem__(self, sliceobj, val):
        self._arr.__setitem__(sliceobj, val)

    def __mul__(self,other):
        return BinarySymplecticMatrix(dot(self._arr, other._arr)%2)

    def __repr__(self):
        # TODO: We could make this a bit
        #       fancier by putting lines between
        #       the blocks, but that'd take a
        #       silly amount of effort, so let's do
        #       that later.
        return repr(self._arr)

    def __eq__(self, other):
        return all(self._arr == other._arr)

    def __and__(self, other):
        """
        Returns the direct sum of this binary symplectic
        matrix with another matrix. Note that resulting
        binary symplectic matrix represents the tensor
        product of the two Clifford gates represented.
        """
        attrs = ['xx', 'xz', 'zx', 'zz']
        blocks = []
        for att in attrs:
            blocks.append(directsum(
                getattr(self, att),
                getattr(other, att)
            ))

        return BinarySymplecticMatrix(*blocks)
        
    ## GATE METHODS ##
    
    def left_H(self,j):
        r"""
        Multiplies on the left by a Hadamard gate on the :math:`j^{\text{th}}`
        qubit. This method acts in-place, as opposed to acting on a copy of the
        binary symplectic matrix. In order to preserve the original matrix,
        use the :meth:`copy` method:
        
        >>> new_bsm = bsm.copy().left_H(idx) # doctest: +SKIP
        """
        u.array_swap(self.zr[j, :], self.xr[j, :])
        return self

    def right_H(self,j):
        r"""
        Multiplies on the right by a Hadamard gate on the :math:`j^{\text{th}}`
        qubit. See :meth:`left_H` for more details.
        """
        u.array_swap(self.zc[:, j], self.xc[:, j])
        return self

    def right_H_all(self):
        r"""
        Multiplies on the right by a Hadamard gate on each
        qubit. See :meth:`left_H` for more details.
        """
        u.array_swap(self.zc, self.xc)
        return self

    def left_SWAP(self,j,k):
        r"""
        Multiplies on the left by a SWAP gate between the :math:`j^{\text{th}}`
        and :math:`k^{\text{th}}` qubits. This method acts in-place, as opposed
        to acting on a copy of the binary symplectic matrix. In order to
        preserve the original matrix, use the :meth:`copy` method:
        
        >>> new_bsm = bsm.copy().left_SWAP(j, k) # doctest: +SKIP
        """
        u.array_swap(self.xr[j, :], self.xr[k, :])
        u.array_swap(self.zr[j, :], self.zr[k, :])
        return self
        
    def right_SWAP(self,j,k):
        r"""
        Multiplies on the right by a SWAP gate between the :math:`j^{\text{th}}`
        and :math:`k^{\text{th}}` qubits. See :meth:`left_SWAP` for more
        details.
        """
        u.array_swap(self.xc[:, j], self.xc[:, k])
        u.array_swap(self.zc[:, j], self.zc[:, k])
        return self
        
    def left_CNOT(self, c, t):
        r"""
        Multiplies on the left by a CNOT gate controlled by the
        :math:`c^{\text{th}}` qubit and targeting the :math:`k^{\text{th}}`
        qubit. This method acts in-place, as opposed to acting on a copy of the
        binary symplectic matrix. In order to preserve the original matrix, use
        the :meth:`copy` method:
        
        >>> new_bsm = bsm.copy().left_CNOT(c, t) # doctest: +SKIP
        """
        if c == t:
            raise ValueError("Control and target qubits cannot be the same.")
        self.xr[t, :] += self.xr[c, :]
        self.xr[t, :] %= 2
        self.zr[c, :] += self.zr[t, :]
        self.zr[c, :] %= 2
        return self
        
    def right_CNOT(self, c, t):
        r"""
        Multiplies on the right by a CNOT gate controlled by the
        :math:`c^{\text{th}}` qubit and targeting the :math:`k^{\text{th}}`
        qubit. For more details, see :meth:`left_CNOT`.
        """
        if c == t:
            raise ValueError("Control and target qubits cannot be the same.")
        self.xc[:, c] += self.xc[:, t]
        self.xc[:, c] %= 2
        self.zc[:, t] += self.zc[:, c]
        self.zc[:, t] %= 2
        return self
        
    def left_R_pi4(self, i):
        r"""
        Multiplies on the left by an :math:`R_{\pi/4}` gate acting on the
        :math:`i^{\text{th}}` qubit. This method acts in-place, as opposed to
        acting on a copy of the binary symplectic matrix. In order to preserve
        the original matrix, use the :meth:`copy` method:
        
        >>> new_bsm = bsm.copy().left_R_pi4(c, t) # doctest: +SKIP
        """
        self.zr[i, :] += self.xr[i, :]
        self.zr[i, :] %= 2
        return self
        
    def right_R_pi4(self, i):
        r"""
        Multiplies on the right by an :math:`R_{\pi/4}` gate acting on the
        :math:`i^{\text{th}}` qubit. For more details, see
        :meth:`left_R_pi4`.
        """
        self.xc[:, i] += self.zc[:, i]
        self.xc[:, i] %= 2
        return self
        
    def left_CZ(self, c1, c2):
        r"""
        Multiplies on the left by an controlled-:math:`Z` gate acting between
        the :math:`c_1^{\text{th}}` and :math:`c_2^{\text{th}}` qubits. This
        method acts in-place, as opposed to acting on a copy of the binary
        symplectic matrix. In order to preserve the original matrix, use the
        :meth:`copy` method:
        
        >>> new_bsm = bsm.copy().left_CZ(c, t) # doctest: +SKIP
        """
        if c1 == c2:
            raise ValueError("Control qubits cannot be the same.")
        self.zr[c2, :] += self.xr[c1, :]
        self.zr[c2, :] %= 2
        self.zr[c1, :] += self.xr[c2, :]
        self.zr[c1, :] %= 2        
        return self
        
    def right_CZ(self, c1, c2):
        r"""
        Multiplies on the right by an controlled-:math:`Z` gate acting between
        the :math:`c_1^{\text{th}}` and :math:`c_2^{\text{th}}` qubits. For more
        details, see :meth:`left_CZ`.
        """
        if c1 == c2:
            raise ValueError("Control qubits cannot be the same.")
        self.xc[:, c2] += self.zc[:, c1]
        self.xc[:, c2] %= 2
        self.xc[:, c1] += self.zc[:, c2]
        self.xc[:, c1] %= 2        
        return self
        
    ## OTHER METHODS ##
        
    def inv(self, check_validity=True):
        """
        Returns the inverse of this binary symplectic matrix, assuming
        that this matrix represents a valid Clifford gate.
        
        Note that if the matrix :math:`H` does not represent a valid Clifford,
        this method will return a matrix :math:`G` such that :math:`H G` is
        not the identity matrix.
        
        :param bool check_validity: If ``True``, then the matrix is first
            checked to ensure that it is a valid Clifford.
        :raises: :class:`qecc.InvalidCliffordError` if ``check_validity`` is ``True``
            and the binary symplectic matrix being inverted does not represent
            a valid Clifford group element.
        """
        if check_validity and not self.is_valid():
            raise InvalidCliffordError('Matrix cannot be inverted.')
        return BinarySymplecticMatrix(self.zz.T,self.xz.T,self.zx.T,self.xx.T)

    def as_clifford(self, check_validity=True):
        """
        Converts this binary symplectic matrix into a Clifford representation.
        
        :param bool check_validity: If ``True``, then the matrix is first
            checked to ensure that it is a valid Clifford.
        
        :rtype: :class:`qecc.Clifford`
        :returns: The same gate as this binary symplectic matrix, represented
            as an instance of :class:`qecc.Clifford`.
        """
        if check_validity and not self.is_valid():
            raise InvalidCliffordError('Matrix cannot be converted.')
        return Clifford(map(array_to_pauli,self.xc.T),map(array_to_pauli,self.zc.T))

    def is_valid(self):
        return is_bsm_valid(self)
        
    def copy(self):
        """
        Returns a copy of this binary symplectic matrix, pointing to a distinct
        location in memory.
        """
        return BinarySymplecticMatrix(self._arr.copy())
        
    def circuit_decomposition(self, validate=True):
        """
        Decomposes the binary symplectic matrix 
        """
        # The circuit decomposition algorithm is long enough that it was moved
        # into another module, bsf_decomp.
        left, right = bsf_decomp.circuit_decomposition_part1(self.copy())
        circ = circuit.Circuit(*(right + list(reversed(left))))
        if validate:
            assert all(circ.as_clifford().as_bsm()._arr == self._arr), "Decomposition failed to produce desired BSM."
        return circ
        
        
## BSM FUNCTIONS ##

def is_bsm_valid(input_bsm):
    xrows = map(BinarySymplecticVector,input_bsm.xc.T)
    zrows = map(BinarySymplecticVector,input_bsm.zc.T)
    for idx_j in range(len(xrows)):
        for idx_k in range(len(zrows)):
            if xrows[idx_j].bsip(zrows[idx_k]) != (idx_j == idx_k):
                return False
            elif xrows[idx_j].bsip(xrows[idx_k]) != 0:
                return False
            elif zrows[idx_j].bsip(zrows[idx_k]) != 0:
                return False
    return True

def bsmzeros(nq):
    """
    Returns a binary symplectic matrix on :math:`n` qubits, initialized to all
    zeros.
    
    :param int nq: Number of qubits that the created matrix will act upon.
    :returns: A binary symplectic matrix containing all zeros.
    :rtype: BinarySymplecticMatrix
    """
    return BinarySymplecticMatrix(
        zeros((2 * (nq),) * 2, dtype=int))

def array_to_pauli(bsv_array):
    return BinarySymplecticVector(bsv_array).as_pauli()    

def directsum(A, B):
    r"""
    Given two matrices :math:`A` and :math:`B` with two indices each,
    returns the direct sum :math:`A \oplus B`.
    
    :type A: ndarray, shape (sA[0], sA[1])
    :type B: ndarray, shape (sB[0], sB[1])
    :rtype: ndarray, shape (sA[0] + sB[0], sA[1] + sB[1])
    :returns: :math:`A \oplus B`
    """
    sA = A.shape
    sB = B.shape
    return hstack([
        vstack([A, zeros((sB[0], sA[1]))]),
        vstack([zeros((sA[0], sB[1])), B])
    ])

## EXAMPLE USAGE ##
        
if __name__ == "__main__":
    bsf_P = Pauli('XYZ').as_bsv()
    print bsf_P
    # Should print:
    #     ( 1 1 0 | 0 1 1 )
    
    print bsf_P.x.astype('int')
    #     array([1, 1, 0])
    
    bsf_cnot = cnot_gt(2, 0, 1).as_bsm()
    print bsf_cnot
    #     ( 1 1 | 0 0
    #       0 1 | 0 0
    #       ----+----
    #       0 0 | 1 1
    #       0 0 | 0 1 )
    
    print int(bsf_cnot.xx[0,1])
    #     1
    
