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
intended for use with Pauli operators.

.. autoclass:: qecc.PauliList
    :members:
    :undoc-members:

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
