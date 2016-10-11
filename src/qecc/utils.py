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

## IMPORTS #####################################################################
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

from functools import wraps
import warnings

if PY3:
    from . import PauliClass as PC
else:
    import PauliClass as PC


## FUNCTIONS ###################################################################

def array_swap(A, B):
    temp = A.copy()
    A[...] = B
    B[...] = temp
    del temp
    
def inv_dict(d):
    return dict(map(reversed, iter(d.items())))
    
# FIXME: Once more LaTeX utils have been written, start a new latex.py and
#         put this function in that module.
def latex_array_contents(cells):
    return " \\\\\n            ".join(" & ".join(row) for row in cells)
     
def transpose(lol):
    return list(map(list, list(zip(*lol))))

## DECORATORS ##################################################################

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
        args = list(map(PC.ensure_pauli, args))
        func(*args)
    return pauli_tested_func
            
## TEST ##

if __name__ == "__main__":
    from numpy import array
    
    A = array([1, 2])
    B = array([3, 4])
    array_swap(A, B)
    print(A, B)
