..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.

.. todo::
    Add documentation on Pauli indexing.
    Mention how Clifford can use Unspecified.

Pauli and Clifford Groups
=========================

:class:`qecc.Pauli` - Class representing Pauli group elements
-------------------------------------------------------------

The class :class:`qecc.Pauli` is used to represent elements of the Pauli group
:math:`\mathcal{P}_n` on :math:`n` qubits. Instances can be constructed by
specifying strings of ``I``, ``X``, ``Y`` and ``Z``, corresponding to the
specification of an operator in the Pauli group.

>>> import qecc as q
>>> P = q.Pauli('X')
>>> print P
i^0 X
>>> Q = q.Pauli('XZZXI')
>>> print Q
i^0 XZZXI
>>> R = q.Pauli('XYZ')
>>> print R
i^0 XYZ

Additionaly, a phase can be provided. Since only integer powers of :math:`i` are
allowed as phases, the phase of a :class:`qecc.Pauli` instance is represented
by an integer in ``range(4)``. Any other integer is converted to an integer
in that range that is equivalent mod 4.

>>> print q.Pauli('X', 2)
i^2 X

The :class:`qecc.Pauli` class supports multiplication, tensor products and
negation by the ``*``, ``&`` and ``-`` operators, respectively.

>>> import qecc
>>> P = qecc.Pauli('X')
>>> Q = qecc.Pauli('Y')
>>> P * Q
i^1 Z
>>> P & Q
i^0 XY
>>> -P * Q
i^3 Z

Using these operators, it is straightforward to construct instances of
:class:`qecc.Pauli` from existing instances. To make this easier, QuaEC provides
single-qubit operators ``I``, ``X``, ``Y`` and ``Z``.

>>> from qecc import I, X, Y, Z
>>> print q.Pauli('XZZXI') & I
i^0 XZZXII

Additionally, instances of :class:`qecc.Pauli` can be tested for equality.

>>> -P * Q == P * -Q
True
>>> P * Q != Q * P
True

The length of a :class:`qecc.Pauli` is defined as the number of qubits it acts
upon.

>>> print len(qecc.Pauli('XYZI'))
4

This information is also exposed as the property ``nq``.

>>> print qecc.Pauli('XYZI').nq
4

Class Reference
~~~~~~~~~~~~~~~
.. autoclass:: qecc.Pauli
    :members:
    :undoc-members:


Iterating Over Groups and Subgroups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autofunction:: qecc.pauli_group

.. autofunction:: qecc.from_generators

Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: qecc.com

.. autofunction:: qecc.elem_gens

.. autofunction:: qecc.eye_p

Searching Over Pauli Group Elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

QuaEC provides useful tools for searching over elements of the Pauli group.
A few particlar searches are provided built-in, while other searches can be
efficiently built using the predicates described in :doc:`predicates`.

.. autofunction:: qecc.is_in_normalizer

.. autofunction:: qecc.mutually_commuting_sets

.. autofunction:: qecc.solve_commutation_constraints


:class:`qecc.Clifford` - Class representing Clifford group elements
-------------------------------------------------------------------

Elements of the automorphism group of the Pauli group (known as the Clifford
group) are represented by the class :class:`qecc.Clifford`. Instances of
``Clifford`` are constructed by specifying the mappings of the generators of
the Pauli group, such that the action of a ``Clifford`` instance is defined
for all input Pauli group elements.

>>> import qecc as q
>>> C = q.Clifford(['XX', 'IX'], ['ZI', 'ZZ'])
>>> print C
XI |->  +XX
IX |->  +IX
ZI |->  +ZI
IZ |->  +ZZ

Once an instance of :class:`qecc.Clifford` has been constructed in this way,
its action on elements of the Pauli group can be calculated by calling the
``Clifford`` instance as a function.

>>> from qecc import I, X, Y, Z
>>> C(X & Y)
i^0 YZ
>>> map(C, ['XI', 'IX', 'YI', 'IY', 'ZI', 'IZ'])
[i^0 XX, i^0 IX, i^0 YX, i^0 ZY, i^0 ZI, i^0 ZZ]

Note that in this example, ``C`` has converted strings to :class:`qecc.Pauli`
instances. This is done automatically by :class:`qecc.Clifford`.

Instances of ``Clifford`` can be combined by multiplication (``*``) and by
tensor products (``&``). Multiplication of two ``Clifford`` instances returns
a new instance representing their composition, while the tensor product returns
a new instance that acts on each register independently.

>>> import qecc as q
>>> C = q.Clifford(['XX', 'IX'], ['ZI', 'ZZ'])
>>> D = q.Clifford(['XI', 'IZ'], ['ZI', 'IX'])
>>> print C * D
XI |->  +XX
IX |->  +ZZ
ZI |->  +ZI
IZ |->  +IX
>>> print C & D
X[0] |->  +X[0] X[1]
X[3] |->  +Z[3]
Z[1] |->  +Z[0] Z[1]
Z[3] |->  +X[3]

Note that in the second example, the printing of the Clifford operator has
switched to a sparse format that suppresses printing lines for qubits that are
not acted upon by the operator (in this case, qubits ``1`` and ``2`` are
trivially acted upon by ``C & D``).

As with :class:`qecc.Pauli`, the length of a :class:`qecc.Clifford` instance
is defined as the number of qubits on which that instance acts. This information
is also exposed as the property ``nq``.

>>> import qecc as q
>>> C = q.Clifford(['XX', 'IX'], ['ZI', 'ZZ'])
>>> print len(C)
2
>>> print C.nq
2

Class Reference
~~~~~~~~~~~~~~~

.. autoclass:: qecc.Clifford
    :members:
    :undoc-members:

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

