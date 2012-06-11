..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.

Pauli and Clifford Groups
=========================

:class:`qecc.Pauli` - Class representing Pauli group elements
-------------------------------------------------------------

.. autoclass:: qecc.Pauli
    :members:
    :undoc-members:

The :class:`qecc.Pauli` class supports multiplication, tensor products and
negation by the ``*``, ``&`` and ``-`` operators, respectively.

>>> P = qecc.Pauli('X')
>>> Q = qecc.Pauli('Y')
>>> P * Q
i^1 Z
>>> P & Q
i^0 XY
>>> -P * Q
i^3 Z

Additionally, instances of :class:`qecc.Pauli` can be tested for equality.

>>> -P * Q == P * -Q
True
>>> P * Q != Q * P
True

The length of a :class:`qecc.Pauli` is defined as the number of qubits it acts
upon.

>>> len(qecc.Pauli('XYZI'))
4

Iterating Over Groups and Subgroups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: qecc.pauli_group

.. autofunction:: qecc.from_generators

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: qecc.com

.. autofunction:: qecc.elem_gens

.. autofunction:: qecc.eye_p

.. autofunction:: qecc.pad

Searching Over Pauli Group Elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

QuaEC provides useful tools for searching over elements of the Pauli group.
A few particlar searches are provided built-in, while other searches can be
efficiently built using the predicates described in :doc:`predicates`.

.. autofunction:: qecc.is_in_normalizer

.. autofunction:: qecc.mutually_commuting_sets



:class:`qecc.Clifford` - Class representing Clifford group elements
-------------------------------------------------------------------

.. autoclass:: qecc.Clifford
    :members:
    :undoc-members:
    
.. todo::
    Document ``len``, ``&``, etc.

Alternate Constructors
~~~~~~~~~~~~~~~~~~~~~~

In addition to specifying the outputs of a Clifford operator acting on the
elementary generators of the Pauli group, one can also create a :class:`Clifford`
instance by specifying the ouput of an operator on an arbitrary generating set.
In particlar, the function :func:`qecc.generic_clifford` takes the inputs and outputs
of a given Clifford operator in order to create a :class:`qecc.Clifford`
instance.

.. autofunction:: qecc.generic_clifford

Iterators onto the Clifford Group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: qecc.clifford_group

Common Clifford Gates
~~~~~~~~~~~~~~~~~~~~~

The :mod:`qecc` package provides support for several common Clifford operators.
These functions can be used to quickly analyze small circuits. For more
extensive circuit support, please see :doc:`circuits`.

.. autofunction:: qecc.eye_c

.. autofunction:: qecc.cnot

.. autofunction:: qecc.hadamard

.. autofunction:: qecc.phase

.. autofunction:: qecc.swap

.. autofunction:: qecc.cz

.. autofunction:: qecc.pauli_gate

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: qecc.paulify

