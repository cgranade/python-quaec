..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.

===================================
Circuit Manipulation and Simulation
===================================

Introduction
============

Quantum circuits are modeled in QuaEC by a sequence type, :class:`qecc.Circuit`,
that stores zero or more circuit elements, known as `locations`. Each location
has a `kind` that indicates if it is a gate, measurement or preparation
location, as well as which gate, which measurement or which preparation is
indicated.

Creating a :class:`qecc.Location` instance consists of specifying the kind of
location along with a sequence of indices indicating which qubits that location
acts upon.

>>> import qecc as q
>>> loc = q.Location('CNOT', 0, 2)

The :meth:`qecc.Location.as_clifford` method allows converting gate locations
back into a :class:`qecc.Clifford` representation if applicable.

>>> print loc.as_clifford()
XII |->  +XIX
IXI |->  +IXI
IIX |->  +IIX
ZII |->  +ZII
IZI |->  +IZI
IIZ |->  +ZIZ

When creating a :class:`qecc.Circuit`, you may specify each location either as
an instance of :class:`qecc.Location` or as a tuple of arguments to :class:`qecc.Location`'s constructor.

>>> circ = q.Circuit(('CNOT', 0, 2), ('H', 1), ('X', 0))

Printing a circuit or location results in that instance being represented in
the QuASM format, a plaintext representation of quantum circuits.

>>> print loc
    CNOT    0 2
>>> print circ
    CNOT    0 2
    H       1
    X       0
        
The number of qubits, depth and size of each location and circuit can be found
by querying the appropriate properties of a :class:`qecc.Location` or :class:`qecc.Circuit`:

>>> print loc.nq
3
>>> print circ.nq, circ.depth, circ.size, len(circ)
3 2 3 3

Once constructed, a :class:`qecc.Circuit` can be transformed in several ways,
including simplifications and representations in terms of depth-1 subcircuits.

>>> circ = q.Circuit(('CNOT', 0, 2), ('H', 1), ('X', 0), ('H', 1))
>>> print circ
    CNOT    0 2
    H       1
    X       0
    H       1
>>> print circ.cancel_selfinv_gates()
    CNOT    0 2
    X       0
>>> circ = q.Circuit(('CZ', 0, 2), ('H', 1), ('X', 0))
>>> print circ.replace_cz_by_cnot()
    H       2
    CNOT    0 2
    H       2
    H       1
    X       0
>>> print "\n    --\n".join(map(str, circ.group_by_time()))
    H       2
    --
    CNOT    0 2
    --
    H       2
    H       1
    X       0
        
Note that, except for :meth:`qecc.Circuit.group_by_time`, each of these
transformations mutates the circuit, so that the original circuit is lost.

>>> print circ
    H       2
    CNOT    0 2
    H       2
    H       1
    X       0

If a circuit consists entirely of Clifford gate locations, then its entire
action may be represented as a :class:`qecc.Clifford` instance:

>>> circ = q.Circuit(('CZ', 0, 2), ('H', 1), ('X', 0))
>>> print circ.as_clifford()
XII |->  +XIZ
IXI |->  +IZI
IIX |->  -ZIX
ZII |->  -ZII
IZI |->  +IXI
IIZ |->  +IIZ

Finally, circuits can be exported to `QCViewer`_ files (``*.qcv``) for easy
integration with QCViewer's functionality.

>>> print circ.as_qcviewer()
.v q1 q2 q3
.i q1
.o q1
BEGIN
    Z    q1 q3
    H    q2
    X    q1
END
<BLANKLINE>

Note that, by default, qubits in the QCViewer export are named "q1", "q2" and so
on. This may be overriden by passing a sequence of strings as the
``qubit_names`` argument. Which qubits get assigned to the ``.i`` and ``.o``
headers in the QCViewer file are controlled by the ``inputs`` and ``outputs``
arguments, respectively.

>>> print circ.as_qcviewer(inputs=(0,), outputs=(0,), qubit_names=["in1", "anc1", "anc2"])
.v in1 anc1 anc2
.i in1
.o in1
BEGIN
    Z    in1 anc2
    H    anc1
    X    in1
END
<BLANKLINE>

.. _QCViewer: http://qcirc.iqc.uwaterloo.ca/index.php?n=Projects.QCViewer

:class:`qecc.Location` - Class representing locations in a circuit
==================================================================

Class Reference
---------------

.. autoclass:: qecc.Location
    :members:
    :undoc-members:
    
:class:`qecc.Circuit` - Class modeling arrangements of locations
================================================================

Class Reference
---------------

.. autoclass:: qecc.Circuit
    :members:
    :undoc-members:
    
