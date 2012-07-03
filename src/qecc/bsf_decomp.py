#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# bsf_decomp.py: Decomposes binary symplectic matrices into circuits.
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

import numpy as np

import PauliClass as pc
import CliffordClass as cc
import bsf

import utils as u

## FUNCTIONS ##

def next1(arr):
    idx_1s = np.nonzero(arr)
    if len(idx_1s[0]) == 0:
        return None
    else:
        return np.min(idx_1s)

def circuit_decomposition_part1(bsm):
    
    left_gateseq = []
    right_gateseq = []
    
    for pivot in xrange(bsm.nq):
        
        ## PROCEDURE 6.5 ##
        
        # STEP 1. Use left multiplication by H and SWAP to move a 1 to the
        #         upper-left corner.
        if bsm[pivot, pivot] != 1:
            # Find which row has the 1 in the pivot column.
            idx_pivot1 = next1(bsm[:, pivot])
            
            # If the length of idx_pivot is 0, then we got the empty list,
            # indicating a rank deficiency.
            if idx_pivot1 is None:
                raise RankDeficientError("The input binary symplectic matrix must be full-rank.")
                
            # If we're still here, then we need to move the 1. We start by moving
            # it to the XX block if needed.
            if idx_pivot1 >= bsm.nq:
                idx_pivot1 -= bsm.nq
                left_gateseq.append(("H", idx_pivot1))
                bsm.left_H(idx_pivot1)
                
            # We can now assume that the 1 is in the XX block, so we move it
            # to the right row of that block with a left SWAP.
            if idx_pivot1 != pivot:
                left_gateseq.append(("SWAP", pivot, idx_pivot1))
                bsm.left_SWAP(pivot, idx_pivot1)

        assert bsm[pivot, pivot] == 1, "Pivot element not correctly set."
        
        # STEP 2. Do column reduction on the pivot column of XX, using left CNOT
        #         to elimate any other 1s in that column.
        while True:
            idx_next_1 = next1(bsm.xx[pivot+1:, pivot])
            if idx_next_1 is None:
                break
            else:
                idx_next_1 += pivot + 1
            
            left_gateseq.append(("CNOT", pivot, idx_next_1))
            bsm.left_CNOT(pivot, idx_next_1)
        
        
        # STEP 3. Do row reduction on the pivot row of XX, using right CNOT to
        #         eliminate any other 1s in that row.
        while True:
            idx_next_1 = next1(bsm.xx[pivot, pivot+1:])
            if idx_next_1 is None:
                break
            else:
                idx_next_1 += pivot + 1
                
            right_gateseq.append(("CNOT", idx_next_1, pivot))
            bsm.right_CNOT(idx_next_1, pivot)
        
        
        # STEP 4. Use left multiplication by R_pi/4 and CZ to column eliminate
        #         the pivot column of ZX.
        while True:
            idx_next_1 = next1(bsm.zx[:, pivot])
            
            if idx_next_1 is None:
                break                
            if idx_next_1 == pivot:
                left_gateseq.append(("R_pi4", pivot))
                bsm.left_R_pi4(pivot)
            else:
                left_gateseq.append(("CZ", pivot, idx_next_1))
                bsm.left_CZ(pivot, idx_next_1)
        
        # STEPS 5-6. Repeat Steps 1 to 4 for pivot += 1.
        # This is taken care of by the loop above.
        
    # STEP 8.
    bsm.right_H_all()
    # STEP 9. Right multiply by R_pi/4 to eliminate the diagonal of A, right multiply by  
    # Since B=C=I, and A is symmetric, we can loop over the lower triangular part and the diagonal:

    gateseq_8910=[]
    Hs_on=set([])
    for idx_r in range(bsm.nq):
        if bsm.xx[idx_r,idx_r]==1:
            gateseq_8910.append(("R_pi4", idx_r))
            bsm.right_R_pi4(idx_r)
            Hs_on.add(idx_r)
        for idx_c in range(idx_r):
            if bsm.xx[idx_r,idx_c]==1:
                gateseq_8910.append(("CZ", idx_c, idx_r))
                bsm.right_CZ(idx_c,idx_r)
                Hs_on.update(set([idx_c,idx_r]))

    # STEP 10.
    bsm.right_H_all()

    appropriate_Hs=map(lambda idx: ("H",idx), list(Hs_on))

    gateseq_8910=appropriate_Hs + gateseq_8910 + appropriate_Hs

    right_gateseq += gateseq_8910

    # Do a final check that the identity matrix is obtained. This should always
    # be the case, and if it is not, that indicates a pretty serious bug.
    # Note that this check is pretty quick compared to the above procedure, so
    # we are not too worried about the slowdown.
    # Moreover, this check doesn't occur in any inner loops, hopefully.
    if not np.all(bsm._arr == np.eye(2 * bsm.nq)):
        print bsm._arr
        raise RuntimeError("Internal error in bsf_decomp.py; decomposition did not produce identity.")

    return left_gateseq, right_gateseq
