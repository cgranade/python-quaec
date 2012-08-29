#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# utils.py: Miscellaneous utility functions used internally by QuaEC.
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

from itertools import imap
from functools import wraps
import PauliClass as PC
import warnings

## FUNCTIONS ##

def array_swap(A, B):
    temp = A.copy()
    A[...] = B
    B[...] = temp
    
def inv_dict(d):
    return dict(imap(reversed, d.iteritems()))

## DECORATORS ##

def memoize(func):
    
    @wraps(func)
    def memoized_func(*args):
        if args in memoized_func.__memoize_cache__:
            return memoized_func.__memoize_cache__[args]
        else:
            retval = func(*args)
            memoized_func.__memoize_cache__[args] = retval
            return retval

    memoized_func.__memoize_cache__ = dict()
    return memoized_func

def deprecated(explanation='Deprecated'):
    def decorator(func):
        @wraps(func)
        def decorated(*args, **kwargs):
            warnings.warn(explanation)
            return func(*args, **kwargs)
        return decorated
    return decorator

def ensure_args_pauli(func):
    @wraps(func)
    def pauli_tested_func(*args):
        map(PC.ensure_pauli, args)
        func(*args)
    return pauli_tested_func
            
## TEST ##

if __name__ == "__main__":
    from numpy import array
    
    A = array([1, 2])
    B = array([3, 4])
    array_swap(A, B)
    print A, B
