..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.

Collections of Pauli Operators
==============================

:class:`qecc.PauliList`: Sequence type for Pauli operators
----------------------------------------------------------

For convinenence, the :mod:`qecc` package provides a subclass of :obj:`list`
intended for use with Pauli operators. ``PauliList`` instances can be created
either by converting an existing instance of a sequence type, or by providing
the elements of the new ``PauliList``.

>>> import qecc as q
>>> L = ['I', 'X', 'Y', 'Z']
>>> print q.PauliList(L)
PauliList(i^0 I, i^0 X, i^0 Y, i^0 Z)
>>> print q.PauliList('XYZ', 'YZX', 'ZXY')
PauliList(i^0 XYZ, i^0 YZX, i^0 ZXY)

Tensor products of a :class:`qecc.Pauli`` with a ``PauliList`` result in
tensoring the given Pauli group element onto each element of the list.

>>> from qecc import X
>>> print q.PauliList(L) & X
PauliList(i^0 IX, i^0 XX, i^0 YX, i^0 ZX)

In general, a :class:`qecc.PauliList` can be used anywhere that a list of
:class:`qecc.Pauli` instances is appropriate. For example, the constructor of
:class:`qecc.Clifford` accepts :class:`qecc.PauliList` instances:

>>> import qecc as q
>>> C = q.Clifford(q.PauliList('XX', q.Unspecified), q.PauliList(q.Unspecified, q.Pauli('ZZ', phase=2)))
>>> print C
XI |->  +XX
IX |-> Unspecified
ZI |-> Unspecified
IZ |->  -ZZ

Class Reference
~~~~~~~~~~~~~~~
.. autoclass:: qecc.PauliList
    :members:
    :undoc-members:
    
