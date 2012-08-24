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

from copy import copy
from operator import add, mul

import PauliClass as pc
import CliffordClass as cc

import utils as u

## ALL ##

__all__ = [
    'Location', 'Circuit',
    'ensure_loc'
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
    CLIFFORD_GATE_KINDS = [
        'I', 'X', 'Y', 'Z', 'H', 'R_pi4', 'CNOT', 'CZ', 'SWAP'
    ]
    
    CLIFFORD_GATE_FUNCS = {
        'I': lambda nq, idx: cc.eye_c(nq),
        'X': lambda nq, idx: pc.elem_gen(nq, idx, 'X').as_clifford(),
        'Y': lambda nq, idx: pc.elem_gen(nq, idx, 'Y').as_clifford(),
        'Z': lambda nq, idx: pc.elem_gen(nq, idx, 'Z').as_clifford(),
        'H': cc.hadamard,
        'R_pi4': cc.phase,
        'CNOT': cc.cnot,
        'CZ': cc.cz,
        'SWAP': cc.swap
    }
    
    KIND_NAMES = sum([
        CLIFFORD_GATE_KINDS
    ], [])
    
    QCVIEWER_NAMES = {
        'I': 'I',           # This one is implemented by a gate definition
                            # included by Circuit.as_qcviewer().
        'X': 'X', 'Y': 'Y', 'Z': 'Z',
        'H': 'H',
        'R_pi4': 'P',
        'CNOT': 'tof',
        'CZ': 'Z',
        'SWAP': 'swap'
    }

    def __init__(self, kind, *qubits):
        if isinstance(kind, int):
            self._kind = kind
        elif isinstance(kind, str):
            self._kind = self.KIND_NAMES.index(kind)
        else:
            raise TypeError("Location kind must be an int or str.")
        
        if not all(isinstance(q, int) for q in qubits):
            raise TypeError('Qubit indices must be integers.')
        
        self._qubits = tuple(qubits)
        self._is_clifford = bool(self.kind in self.CLIFFORD_GATE_KINDS)
        
    def __str__(self):
        return "    {:<4}    {}".format(self.kind, ' '.join(map(str, self.qubits)))
    def __repr__(self):
        return "<{} Location on qubits {}>".format(self.kind, self.qubits)
    def __hash__(self):
        return hash((self._kind,) + self.qubits)
        
    ## IMPORT METHODS ##
    
    @staticmethod
    def from_quasm(source):
        """
        Returns a :class:`qecc.Location` initialized from a QuASM-formatted line.
        
        :type str source: A line of QuASM code specifying a location.
        :returns qecc.Location: The location represented by the given QuASM source.
        """
        parts = source.split()
        return Location(parts[0], *map(int, parts[1:]))
        
    ## PROPERTIES ##
        
    @property
    def kind(self):
        """
        Returns a string defining which kind of location this instance
        represents. Guaranteed to be a string that is an element of
        ``Location.KIND_NAMES``.
        """
        return self.KIND_NAMES[self._kind]        
        
    @property
    def qubits(self):
        """
        Returns a tuple of ints describing which qubits this location acts upon.
        """
        return self._qubits
        
    @property
    def nq(self):
        """
        Returns the number of qubits in the smallest circuit that can contain
        this location without relabeling qubits. For a :class:`qecc.Location`
        ``loc``, this property is defined as ``1 + max(loc.nq)``.
        """
        return 1 + max(self.qubits)
        
    @property
    def is_clifford(self):
        """
        Returns ``True`` if and only if this location represents a gate drawn
        from the Clifford group.
        """
        return self._is_clifford
      
    ## SIMULATION METHODS ##
        
    def as_clifford(self, nq=None):
        """
        If this location represents a Clifford gate, returns the action of that
        gate. Otherwise, a :obj:`RuntimeError` is raised.
        
        :param int nq: Specifies how many qubits to represent this location as
            acting upon. If not specified, defaults to the value of the ``nq``
            property.
        :rtype: qecc.Clifford
        """
        if not self.is_clifford:
            raise RuntimeError("Location must be a Clifford gate.")
        else:
            if nq is None:
                nq = self.nq
            elif nq < self.nq:
                raise ValueError('nq must be greater than or equal to the nq property.')
                
            return self.CLIFFORD_GATE_FUNCS[self.kind](nq, *self.qubits)
            
        
    ## EXPORT METHODS ##
        
    def as_qcviewer(self, qubit_names=None):
        """
        Returns a representation of this location in a format suitable for
        inclusion in a QCViewer file. 
            
        :param qubit_names: If specified, the given aliases will be used for the
            qubits involved in this location when exporting to QCViewer.
            Defaults to "q1", "q2", etc.
        :rtype: str
        
        Note that the identity (or "wait") location requires the following to be
        added to QCViewer's ``gateLib``::
        
            NAME wait
            DRAWNAME "1"
            SYMBOL I
            1 , 0
            0 , 1
        """
        # FIXME: link to QCViewer in the docstring here.
        return '    {gatename}    {gatespec}\n'.format(
            gatename=self.QCVIEWER_NAMES[self.kind],
            gatespec=qubits_str(self.qubits, qubit_names),
            )
        
    ## OTHER METHODS ##
    
    def relabel_qubits(self, relabel_dict):
        """
        Returns a new location related to this one by a relabeling of the
        qubits. The relabelings are to be indicated by a dictionary that
        specifies what each qubit index is to be mapped to.
        
        >>> import qecc as q
        >>> loc = q.Location('CNOT', 0, 1)
        >>> print loc
            CNOT    0 1
        >>> print loc.relabel_qubits({1: 2})
            CNOT    0 2
            
        :param dict relabel_dict: If `i` is a key of `relabel_dict`, then qubit
            `i` will be replaced by `relabel_dict[i]` in the returned location.
        :returns qecc.Location: A new location with the qubits relabeled as 
            specified by `relabel_dict`.
        """
        return Location(self.kind, *tuple(relabel_dict[i] if i in relabel_dict else i for i in self.qubits))

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
        
    append.__doc__ = list.append.__doc__
        
    def insert(self, at, newval):
        super(Circuit, self).insert(at, ensure_loc(newval))
        
    insert.__doc__ = list.insert.__doc__
        
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
        """
        Returns the number of qubits on which this circuit acts.
        """
        return max(loc.nq for loc in self)
        
    @property
    def size(self):
        """
        Returns the number of locations in this circuit. Note that this property
        is synonymous with :obj:`len`, in that ``len(circ) == circ.size`` for
        all :class:`qecc.Circuit` instances.
        """
        return len(self)
        
    @property
    def depth(self):
        """
        Returns the minimum number of timesteps required to implement exactly
        this circuit in parallel.
        """
        return len(list(self.group_by_time()))

    ## IMPORT CLASS METHODS ##
    
    @staticmethod
    def from_quasm(source):
        if not isinstance(source, str):
            # Assume source is a file-like, so that iter(source) returns lines
            # in the file.
            it = iter(source)
        else:
            it = iter(source.split('\n'))
            
        return Circuit(*map(Location.from_quasm, it))

    ## PRETTY PRINTING ##
            
    def __repr__(self):
        return "Circuit({})".format(", ".join(map(repr, self)))
    def __str__(self):
        return "\n".join(map(str, self))
        
    def as_quasm(self):
        """
        Returns a representation of the circuit in an assmembler-like format.
        In this format, each location is represented by a single line where
        the first field indicates the kind of location and the remaining fields
        indicate the qubits upon which the location acts.
        
        >>> import qecc as q
        >>> circ = q.Circuit(('CNOT', 0, 2), ('H', 2), ('SWAP', 1, 2), ('I', 0))
        >>> print circ.as_quasm()
            CNOT    0 2
            H       2
            SWAP    1 2
            I       0
        """
        return str(self)

    def as_qcviewer(self, inputs=(0,), outputs=(0,), qubit_names=None):
        """
        Returns a string representing this circuit in the format recognized by
        `QCViewer`_.
        
        :param tuple inputs: Specifies which qubits should be marked as inputs
            in the exported QCViewer circuit.
        :param tuple outputs: Specifies which qubits should be marked as outputs
            in the exported QCViewer circuit.
        :param qubit_names: Names to be used for each qubit when exporting to
            QCViewer.
        
        .. _QCViewer: http://qcirc.iqc.uwaterloo.ca/index.php?n=Projects.QCViewer
        """
            
        header = '.v ' + qubits_str(range(self.nq), qubit_names) + '\n'
        header += '.i ' + qubits_str(inputs, qubit_names) + '\n'
        header += '.o ' + qubits_str(outputs, qubit_names) + '\n'
            
        circ_text = 'BEGIN\n'
        for loc in self:
            circ_text += loc.as_qcviewer(qubit_names)
        circ_text += 'END\n'
        
        return header + circ_text

    ## CIRCUIT SIMULATION METHODS ##
    
    def as_clifford(self):
        """
        If this circuit is composed entirely of Clifford operators, converts it
        to a :class:`qecc.Clifford` instance representing the action of the
        entire circuit. If the circuit is not entirely Clifford gates, this method
        raises a :obj:`RuntimeError`.
        """
        if not all(loc.is_clifford for loc in self):
            raise RuntimeError('All locations must be Clifford gates in order to represent a circuit as a Clifford operator.')
            
        nq = self.nq
        return reduce(mul, (loc.as_clifford(nq) for loc in reversed(self)), cc.eye_c(nq))
        

    ## CIRCUIT SIMPLIFICATION METHODS ##
        
    def cancel_selfinv_gates(self, start_at=0):
        """
        Transforms the circuit, removing any self-inverse gates from the circuit
        if possible. Note that not all self-inverse gates are currently
        supported by this method.
        
        :param int start_at: Specifies which location to consider first. Any
            locations before ``start_at`` are not considered for cancelation by
            this method.
        """
        
        SELFINV_GATES = ['H', 'X', 'Y', 'Z', 'CNOT']
        
        if start_at == len(self):
            return self
            
        loc = self[start_at]
        if loc.kind in SELFINV_GATES:
            if len(loc.qubits) == 1:
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
        """
        Changes all controlled-:math:`Z` gates in this circuit to
        controlled-NOT gates, adding Hadamard locations as required.
        """
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
        """
        Returns an iterator onto subcircuits of this circuit, each of depth 1.
        
        :param bool pad_with_waits: If ``True``, each subcircuit will have
            wait locations added such that every qubit is acted upon in every
            subcircuit.
            
        :yields: each depth-1 subcircuit, corresponding to time steps of the
            circuit
        """
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

    ## OTHER METHODS ##
    
    def relabel_qubits(self, relabel_dict):
        """
        Returns a new circuit related to this one by a relabeling of the
        qubits. The relabelings are to be indicated by a dictionary that
        specifies what each qubit index is to be mapped to.
        
        >>> import qecc as q
        >>> loc = q.Location('CNOT', 0, 1)
        >>> print loc
            CNOT    0 1
        >>> print loc.relabel_qubits({1: 2})
            CNOT    0 2
            
        :param dict relabel_dict: If `i` is a key of `relabel_dict`, then qubit
            `i` will be replaced by `relabel_dict[i]` in the returned circuit.
        :returns qecc.Location: A new circuit with the qubits relabeled as 
            specified by `relabel_dict`.
        """
        return Circuit(*[
            loc.relabel_qubits(relabel_dict) for loc in self
        ])
