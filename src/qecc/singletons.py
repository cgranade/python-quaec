#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# singletons.py: Singleton object values used for tracking special conditions
#     in QuaEC.
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

## ALL ##

__all__ = ['EmptyClifford', 'Unspecified']

## METACLASSES ##

# The following metaclass is borrowed from:
# http://stackoverflow.com/questions/6760685/creating-a-singleton-in-python

class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]
        
## CLASSES ##

class EmptyCliffordType(object):
    """
    TODO
    """
    __metaclass__ = Singleton
    
    def __repr__(self):
        return "EmptyClifford"
    
EmptyClifford = EmptyCliffordType()
    
class UnspecifiedType(object):
    """
    Marks that a given constraint is unspecified.
    """
    __metaclass__ = Singleton
    
    def __repr__(self):
        return "Unspecified"
        
        
    # We make it so that Unspecified times any other type is also Unspecified.
    def __mul__(self, other):
        return self
        
    def __rmul__(self, other):
        return self
    
Unspecified = UnspecifiedType()

