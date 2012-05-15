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

## FUNCTIONS ##

def array_swap(A, B):
    temp = A.copy()
    A[...] = B
    B[...] = temp
        
## TEST ##

if __name__ == "__main__":
    from numpy import array
    
    A = array([1, 2])
    B = array([3, 4])
    array_swap(A, B)
    print A, B
