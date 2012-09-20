..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.


Binary Symplectic Form
======================

Introduction
~~~~~~~~~~~~

The :mod:`qecc` package provides support for elements of the Pauli and Clifford groups in 
binary symplectic form, including support for algorithms acting on these representations.
Note that all classes and functions documented here depend on the :mod:`numpy` package. For
more information on the binary symplectic representation, read [CRSS96]_, Section 2.

:class:`qecc.BinarySymplecticVector` - Binary symplectic representation of Pauli group elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The class :class:`qecc.BinarySymplecticVector` provides a means of representing elements of the
Pauli group (neglecting global phases) using binary vectors :math:`a` and :math:`b` such that an
element :math:`P` of the Pauli group acting on :math:`n` qubits is :math:`X^{a}Z^{b} = X^{a_1}Z^{b_1}
\otimes \ldots \otimes X^{a_n}Z^{b_n}`. Binary symplectic vectors can be obtained from a single binary
list, two binary lists, or converted from another Pauli instance (removing the phase):

>>> import qecc as q
>>> a=[1, 0, 1]; b=[0, 1, 1]
>>> q.BinarySymplecticVector(a,b)==q.BinarySymplecticVector(a+b)
True

>>> import qecc as q
>>> a=[1, 0, 1]; b=[0, 1, 1]
>>> q.BinarySymplecticVector(a,b)
( 1 0 1 | 0 1 1 )

>>> import qecc as q
>>> q.Pauli('XYIYIIZ',2).as_bsv()
( 1 1 0 1 0 0 0 | 0 1 0 1 0 0 1 )


Class Reference
---------------

.. autoclass:: qecc.BinarySymplecticVector
    :members:
    :undoc-members:

.. todo::
    Need to fill out more documentation and examples here.

Utility Functions
-----------------

.. autofunction:: qecc.all_pauli_bsvs

.. autofunction:: qecc.constrained_set

.. autofunction:: qecc.commute

.. autofunction:: qecc.xz_switch

:class:`qecc.BinarySymplecticMatrix` - Binary symplectic representation of Clifford group elements
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Class Reference
---------------

.. autoclass:: qecc.BinarySymplecticMatrix
    :members:
    :undoc-members:

.. autofunction:: qecc.is_bsm_valid

.. autofunction:: qecc.bsmzeros

.. autofunction:: qecc.array_to_pauli



