#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# circuit.py: Modeling for stabilizer circuits.
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

from abc import ABCMeta, abstractmethod, abstractproperty
from operator import add
from functools import partial

import PauliClass as pc
import CliffordClass as cc

import utils as u

## ALL ##

__all__ = [
    'Location', 'GateLocation', 'CNOTLoc', 'HadaLoc', 'PhaseLoc', 'CHP',
    'WaitLoc', 'PrepLoc', 'MeasLoc', 'Circuit'
]

## CLASSES ##

class Location(object):
    __metaclass__ = ABCMeta
    
    def __init__(self, *qubits):
        self.qubits = qubits
        if len(qubits) != self.n_acts_on:
            raise ValueError(
                "This location acts on {0} qubits, but {1} have been specified."
                .format(self.n_acts_on, len(qubits)))
        
    def __mul__(self, other):
        return Circuit(self, other)
        
    def __and__(self, other):
        return Circuit(self, other.shift_by(1 + max(self.qubits)))
        
    def __repr__(self):
        return "\t{0}\t{1}".format(self.name, " ".join(map(str, self.qubits)))
        
    def named_repr(self, qubit_names):
        return "\t{0}\t{1}".format(self.name, " ".join(map(qubit_names.__getitem__, self.qubits)))
        
    def shift_by(self, nq):
        self.qubits = map(partial(add, nq), self.qubits)
        return self
        
    @abstractproperty
    def n_acts_on(self): pass
    @abstractproperty
    def name(self): pass
        
class GateLocation(Location):
    @abstractmethod
    def as_clifford(self, nq):
        pass
        
    def as_bsm(self, nq):
        return self.as_clifford(nq).as_bsm()
        
class CNOTLoc(GateLocation):
    def as_clifford(self, nq):
        return cc.cnot(nq, *self.qubits)
    
    @property    
    def n_acts_on(self): return 2
    @property
    def name(self): return 'CNOT'
    
class HadaLoc(GateLocation):        
    def as_clifford(self, nq):
        return cc.hadamard(nq, *self.qubits)
        
    @property
    def n_acts_on(self): return 1
    @property
    def name(self): return 'H'
    
class PhaseLoc(GateLocation):        
    def as_clifford(self, nq):
        return cc.phase(nq, *self.qubits)
        
    @property
    def n_acts_on(self): return 1
    @property
    def name(self): return 'P'
    
CHP = (CNOTLoc, HadaLoc, PhaseLoc)
# Useful for the following idiom:
# >>> C, H, P = circuit.CHP
    
class WaitLoc(GateLocation):     
    def as_clifford(self, nq):
        return cc.eye_c(nq)
        
    @property
    def n_acts_on(self): return 1
    @property
    def name(self): return '1'
    
class PrepLoc(Location):
    pass
    
class MeasLoc(Location):
    pass
    
class Circuit(object):
    def __init__(self, *elems):
        self.circuit_elems = list(elems)
        
    def __len__(self):
        return len(self.circuit_elems)
        
    def __and__(self, other):
        pass
        
    def __mul__(self, other):
        if isinstance(other, Location):
            return Circuit(*(self.circuit_elems + [other]))
        else:
            return Circuit(*(self.circuit_elems + other.circuit_elems))
        
    def __repr__(self):
        return "\n".join(map(repr, self.circuit_elems))
        
    @property
    def nq(self):
        return 1 + max(map(lambda elem: max(elem.qubits), reduce(add, self.circuit_elems)))

## EXAMPLE USAGE ##

if __name__ == "__main__":
    pass

    # C = Circuit([Prep('q1', 'X'), Prep('q2', 'X'), Prep('q3', 'Z')], [CNOT('q1', 1), Had('q3')], [CNOT('q1', 'q3'), Wait('q2')])
    # C = C & Circuit([Had('a1')]) # <- Pads with an extra time slice comprised of the empty list.
    # C = C & Circuit([Phase('a2')], [Wait('a3')])
    
    
