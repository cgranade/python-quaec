..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.

Binary Symplectic Form
======================

The :mod:`qecc` package provides support for elements of the Pauli and Clifford groups in 
binary symplectic form, including support for algorithms acting on these representations.
Note that all classes and functions documented here depend on the :mod:`numpy` package. For
more information on the binary symplectic representation, read [CRSS96]_, Section 2.

Class Reference
===============

:class:`qecc.BinarySymplecticVector` - Binary symplectic representation of Pauli group elements
-----------------------------------------------------------------------------------------------

.. autoclass:: qecc.BinarySymplecticVector
    :members:
    :undoc-members:

.. todo::
    Need to fill out more documentation and examples here.

Vectors in binary symplectic form can be initialized from a single list of integers of
length :math:`2n`, or from a pair of integer lists of length :math:`n`:

>>> import qecc as q
>>> pauli_from_2n = q.BinarySymplecticVector([1, 0, 0, 0, 1 ,0])
>>> pauli_from_2n
( 1 0 0 | 0 1 0 )

>>> import qecc as q
>>> pauli_from_2_lists = q.BinarySymplecticVector([1, 0, 0],[0, 1 ,0])
>>> pauli_from_2_lists
( 1 0 0 | 0 1 0 )


Utility Functions
~~~~~~~~~~~~~~~~~

.. autofunction:: qecc.all_pauli_bsvs

.. autofunction:: qecc.constrained_set

.. autofunction:: qecc.commute

.. autofunction:: qecc.xz_switch

Class Reference
===============

:class:`qecc.BinarySymplecticMatrix` - Binary symplectic representation of Clifford group elements
--------------------------------------------------------------------------------------------------

.. autoclass:: qecc.BinarySymplecticMatrix
    :members:
    :undoc-members:

.. autofunction:: qecc.is_bsm_valid

.. autofunction:: qecc.bsmzeros

.. autofunction:: qecc.array_to_pauli



