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
from copy import copy
from operator import add
from functools import partial

import PauliClass as pc
import CliffordClass as cc

import utils as u

## ALL ##

__all__ = [
    'Location', 'Circuit'
]

## INTERNAL FUNCTIONS ##

def qubits_str(qubits, qubit_names=None):
    if qubit_names is None:
        return ' '.join('q{}'.format(idx + 1) for idx in qubits)
    else:
        return ' '.join(qubit_names[idx] for idx in qubits)

## CLASSES ##

class Location(object):
    """
    Represents a gate, wait, measurement or preparation location in a
    circuit.
    
    Note that currently, only gate locations are implemented.
    
    :param kind: The kind of location to be created. Each kind is an
        abbreviation drawn from ``Location.KIND_NAMES``, or is the index in
        ``Location.KIND_NAMES`` corresponding to the desired location kind.
    :type kind: int or str
    :param qubits: Indicies of the qubits on which this location acts.
    :type qubits: tuple of ints.
    """

    ## CLASS CONSTANTS ##
    KIND_NAMES = [
        'I', 'X', 'Y', 'Z', 'H', 'R_pi4', 'CNOT', 'CZ', 'SWAP'
    ]

    def __init__(self, kind, *qubits):
        if isinstance(kind, int):
            self._kind = kind
        elif isinstance(kind, str):
            self._kind = self.KIND_NAMES.index(kind)
        else:
            raise TypeError("Location kind must be an int or str.")
        
        self._qubits = tuple(qubits)
        
    def __str__(self):
        return "\t{}\t{}".format(self.kind, ' '.join(map(str, self.qubits)))
    def __repr__(self):
        return "<{} Location on qubits {}>".format(self.kind, self.qubits)
    def __hash__(self):
        return hash((self._kind,) + self.qubits)
        
    def as_qcviewer(self, qubit_names=None):
        """
        Returns a representation of this location in a format suitable for
        inclusion in a QCViewer file.
        """
        # FIXME: link to QCViewer in the docstring here.
        # FIXME: map location kinds to QCViewer equivalents.
        return '\t{gatename}\t{gatespec}\n'.format(
            gatename=self.kind,
            gatespec=qubits_str(self.qubits, qubit_names),
            )
        
    @property
    def kind(self):   return self.KIND_NAMES[self._kind]        
    @property
    def qubits(self): return self._qubits
    @property
    def nq(self):     return 1 + max(self.qubits)

def ensure_loc(loc):
    if isinstance(loc, tuple):
        loc = Location(*loc)
    elif not isinstance(loc, Location):
        raise TypeError('Locations must be specified either as Location instances or as tuples.')
    return loc

class Circuit(list):
    def __init__(self, *locs):
        # Circuit(('CNOT', 0, 2), ('H', 1)) works, but
        # Circuit('CNOT', 0, 2) doesn't work.
        list.__init__(self, map(ensure_loc, locs))

    ## SEQUENCE PROTOCOL ##
            
    def append(self, newval):
        super(Circuit, self).append(ensure_loc(newval))
    def insert(self, at, newval):
        super(Circuit, self).insert(at, ensure_loc(newval))
            
    def __repr__(self):
        return "Circuit({})".format(", ".join(map(repr, self)))
    def __str__(self):
        return "\n".join(map(str, self))
        
    def __getitem__(self, *args):
        item = super(Circuit, self).__getitem__(*args)
        if not isinstance(item, list):
            return item
        else:
            return Circuit(*item)
            
    def __getslice__(self, *args):
        return Circuit(*super(Circuit, self).__getslice__(*args))
        
    def __add__(self, other):
        if not isinstance(other, Circuit):
            other = Circuit(*other)
        return Circuit(*super(Circuit, self).__add__(other))
    def __iadd__(self, other):
        if not isinstance(other, Circuit):
            other = Circuit(*other)
        return Circuit(*super(Circuit, self).__iadd__(other))

    ## PROPERTIES ##
                
    @property
    def nq(self):
        return max(loc.nq for loc in self)
        
    @property
    def n_timesteps(self):
        return len(list(self.group_by_time()))

    ## PRETTY PRINTING ##

    def as_qcviewer(self, inputs=(0,), outputs=None, qubit_names=None):
        if outputs is None:
            outputs = (idx for idx in range(self.nq) if idx not in inputs)
            
        acc = '.v ' + qubits_str(range(self.nq), qubit_names) + '\n'
        acc += '.i ' + qubits_str(inputs, qubit_names) + '\n'
        acc += '.o ' + qubits_str(outputs, qubit_names) + '\n'
            
        acc += 'BEGIN\n'
        for loc in self:
            acc += loc.as_qcviewer(qubit_names)

        acc += 'END\n'
        return acc

    ## CIRCUIT SIMPLIFICATION METHODS ##
        
    def cancel_selfinv_gates(self, start_at=0):
        SELFINV_GATES = ['H', 'X', 'Y', 'Z']
        
        if start_at == len(self):
            return self
            
        loc = self[start_at]
        if len(loc.qubits) == 1 and loc.kind in SELFINV_GATES:
            # TODO: add two-qubit gates.
            q = loc.qubits[0]
            
            for idx_future in xrange(start_at + 1, len(self)):
                if q in self[idx_future].qubits:
                    # Check that the kind matches.
                    if self[idx_future].kind == loc.kind:
                        self.pop(idx_future)
                        self.pop(start_at)
                        return self.cancel_selfinv_gates(start_at=start_at)
                    else:
                        # Go on to the next gate, since there's another gate
                        # between here.
                        return self.cancel_selfinv_gates(start_at=start_at+1)
        
        return self.cancel_selfinv_gates(start_at=start_at+1)
        
    def replace_cz_by_cnot(self):
        # FIXME: this is inefficient as hell right now.
        try:
            idx = (idx for idx in range(len(self)) if self[idx].kind == 'CZ').next()
            q = self[idx].qubits
            self[idx] = Location('CNOT', *q)
            self.insert(idx + 1, ('H', q[1]))
            self.insert(idx, ('H', q[1]))
            return self.replace_cz_by_cnot()
        except StopIteration:
            return self
        
    def group_by_time(self, pad_with_waits=False):
        nq = self.nq
        
        found = [False] * nq        
        group_acc = Circuit()
        
        for loc in self:
            if any(found[qubit] for qubit in loc.qubits):
                if pad_with_waits:
                    group_acc += [('I', qubit) for qubit in range(nq) if not found[qubit]]
                yield group_acc
                
                found = [False] * nq
                group_acc = Circuit()
                
            for qubit in loc.qubits:
                found[qubit] = True
                
            group_acc.append(loc)
            
        
        if pad_with_waits:
            group_acc += [('I', qubit) for qubit in range(nq) if not found[qubit]]
        yield group_acc

## EXAMPLE USAGE ##

if __name__ == "__main__":
    pass

    # C = Circuit([Prep('q1', 'X'), Prep('q2', 'X'), Prep('q3', 'Z')], [CNOT('q1', 1), Had('q3')], [CNOT('q1', 'q3'), Wait('q2')])
    # C = C & Circuit([Had('a1')]) # <- Pads with an extra time slice comprised of the empty list.
    # C = C & Circuit([Phase('a2')], [Wait('a3')])
    
    
