#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# unitary_reps.py: Internal module for working with unitary representations.
#     Functionality of this module is exposed via methods of other classes.
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

from functools import wraps

import operator as op

import  PauliClass as pc

try:
    import numpy as np
    import numpy.linalg as la
except ImportError:
    import warnings
    warnings.warn("NumPy is missing; some functionality may not be available.")
    np = None
    la = None
    
## DECORATORS ##

def requires_numpy(func):
    @wraps(func)
    def checked_func(*args, **kwargs):
        if np is None:
            raise RuntimeError('{} requires NumPy, but NumPy could not be imported.'.format(func.__name__))
        return func(*args, **kwargs)
        
    return checked_func
    
## CONSTANTS ##

if np is not None:
    PAULIS = {
        'I': np.eye(2),
        'X': np.array([[0, 1], [1, 0]]),
        'Y': np.array([[0, -1j], [1j, 0]]),
        'Z': np.array([[1, 0], [0, -1]])
    }
else:
    PAULIS = {
        'I': None,
        'X': None,
        'Y': None,
        'Z': None
    }
    
## FUNCTIONS ##

@requires_numpy
def mutual_eigenspace(arrays, desired_eigval=1, thresh=1e-7):
    def proj_for_array(array):
        eigvals, eigvecs = la.eig(array)
        good_idxs = np.nonzero(np.abs(eigvals - desired_eigval) < thresh)
        good_eigvecs = eigvecs.T[good_idxs]
        proj = sum(np.dot(vec[..., np.newaxis], vec[..., np.newaxis].T.conj()) for vec in good_eigvecs)
        return proj

    # Assumes square matrices of the same size.
    proj = np.eye(arrays[0].shape[0])
    proj = reduce(np.dot, map(proj_for_array, arrays), proj)
    eigvals, eigvecs = la.eig(proj)
    return eigvecs.T[np.abs(eigvals - 1) < thresh]
        
@requires_numpy
def bitvector(nq, bits):
    v = np.zeros((2**nq, ))
    v[int(''.join(map(str, bits)), 2)] = 1
    return v
    
@requires_numpy
def pauli_as_unitary(P):
    ph = 1j ** P.ph
    return ph * reduce(np.kron, (PAULIS[op] for op in P.op))
    
@requires_numpy
def clifford_as_unitary(C):
    nq = len(C)
    dim = 2**nq
    U = np.zeros((dim,dim), dtype='complex')
    psi_0 = mutual_eigenspace(map(pauli_as_unitary, C.zout)).T
    for b in xrange(dim):
        bits = '{{0:0>{nq}b}}'.format(nq=nq).format(b)
        Xb   = reduce(op.mul, (C.xout[idx] for idx in xrange(nq) if bits[idx] == '1'), pc.eye_p(nq)).as_unitary()
        for a in xrange(dim):
            bra_a = np.zeros((1, dim))
            bra_a[0, a] = 1
            U[a, b] = reduce(np.dot, [bra_a, Xb, psi_0])[0,0]
            
    return U
