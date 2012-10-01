#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# exceptions.py: Exception classes used by QuaEC.
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

__all__ = ["InvalidCliffordError", "RankDeficientError"]

class InvalidCliffordError(ValueError):
    """We raise this exception wherever an automated procedure
    has produced a list of output Pauli matrices that does not 
    commute as the result of an automorphism on the Paulis, or
    where some other idiocy has occurred."""
    pass
    
class RankDeficientError(ValueError):
    """
    Indicates that an algorithm failed due to its input being rank-deficient.
    """
    pass
