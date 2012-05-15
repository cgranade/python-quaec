..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.

Predicates and Filters
======================

:class:`qecc.Predicate`
-----------------------

The :mod:`qecc` package provides a class :class:`Predicate` to represent a
predicate function; that is, a function which returns a :obj:`bool`.

.. autoclass:: qecc.Predicate
    :members:
    :undoc-members:
    
Specific Predicate Functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
Several useful predefined predicates are provided by :mod:`qecc`.

.. autoclass:: qecc.SetMembershipPredicate
    :members:
    :undoc-members:

.. autoclass:: qecc.PauliMembershipPredicate
    :members:
    :undoc-members:
    
In addition, utility functions are provided for constructing predicates based
on commutation properties of the Pauli group.

.. autofunction:: qecc.commutes_with
.. autofunction:: qecc.in_group_generated_by

Usage Examples
~~~~~~~~~~~~~~

Predicate functions can be used to quickly generate collections of
:class:`qecc.Pauli` operators having a given set of properties.

>>> print filter(
...     commutes_with('XX', 'ZZ') & ~in_group_generated_by('XX'),
...     pauli_group(2)
...     )
[i^0 YY, i^0 ZZ]


