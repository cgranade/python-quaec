..
    This work is licensed under the Creative Commons Attribution-
    NonCommercial-ShareAlike 3.0 Unported License. To view a copy of this
    license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send a
    letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
    California, 94041, USA.

================
Stabilizer Codes
================

Introduction
============

QuaEC includes a class, :class:`qecc.StabilizerCode`, that represents
error-correcting codes specified using the stabilizer formalism. To construct
a stabilizer code in QuaEC, the generators of a stabilizer group must be
specified along with a particular assignment of logical operators acting on
states encoded in the stabilizer code.

>>> import qecc as q
>>> stab = q.StabilizerCode(['ZZI', 'IZZ'], ['XXX'], ['ZZZ'])
>>> print stab
S = <i^0 ZZI, i^0 IZZ>
Xbars = PauliList(i^0 XXX)
Zbars = PauliList(i^0 ZZZ)

For convienence, several static methods are provided to create instances for
well-known stabilizer codes.

>>> stab = q.StabilizerCode.perfect_5q_code()
>>> print stab
5-qubit perfect code
S = <i^0 XZZXI, i^0 IXZZX, i^0 XIXZZ, i^0 ZXIXZ>
Xbars = PauliList(i^0 XXXXX)
Zbars = PauliList(i^0 ZZZZZ)

Once constructed, an instance of :class:`qecc.StabilizerCode` exposes properties
that describe the number of physical and logical qubits, as well as the distance
of the code. (Please note that calculating the distance can be extremely slow
for large codes.)

>>> print (stab.nq, stab.nq_logical, stab.distance)
(5, 1, 3)

Encoders and decoders for stabilizer codes can be found in a straightforward
manner using :class:`qecc.StabilizerCode`.

>>> enc = stab.encoding_cliffords().next()
>>> print enc
X[0] |->  +X[0] X[1] X[2] X[3] X[4]
X[1] |->  +X[0] X[2] X[3] X[4]
X[2] |->  +X[1] X[2]
X[3] |->  +Y[0] X[1] X[3] Y[4]
X[4] |->  +X[0] X[1] Y[3] Y[4]
Z[0] |->  +Z[0] Z[1] Z[2] Z[3] Z[4]
Z[1] |->  +X[0] Z[1] Z[2] X[3]
Z[2] |->  +X[1] Z[2] Z[3] X[4]
Z[3] |->  +X[0] X[2] Z[3] Z[4]
Z[4] |->  +Z[0] X[1] X[3] Z[4]
>>> print enc.inv()
X[0] |->  -X[0] Z[3] X[4]
X[1] |->  +X[0] X[1]
X[2] |->  +X[0] X[1] X[2]
X[3] |->  -X[0] X[2] X[3] Z[4]
X[4] |->  +X[0] Y[3] Y[4]
Z[0] |->  -Z[0] Y[1] Y[3] Z[4]
Z[1] |->  -Z[0] Y[2] Z[3] Y[4]
Z[2] |->  +Z[0] Z[1] Z[2] X[3]
Z[3] |->  -Z[0] Y[1] Z[3] Y[4]
Z[4] |->  -Z[0] Z[1] X[2] Z[3] Z[4]

Stabilizer codes may be combined by the tensor product (reprsented in QuaEC by
``&``), or by concatenation:

>>> print stab & stab
S = <i^0 XZZXIIIIII, i^0 IXZZXIIIII, i^0 XIXZZIIIII, i^0 ZXIXZIIIII, i^0 IIIIIXZZXI, i^0 IIIIIIXZZX, i^0 IIIIIXIXZZ, i^0 IIIIIZXIXZ>
Xbars = PauliList(i^0 XXXXXIIIII, i^0 IIIIIXXXXX)
Zbars = PauliList(i^0 ZZZZZIIIII, i^0 IIIIIZZZZZ)
>>> print q.StabilizerCode.bit_flip_code(1).concatenate(q.StabilizerCode.phase_flip_code(1))
S = <i^0 Z[0] Z[1], i^0 Z[1] Z[2], i^0 Z[3] Z[4], i^0 Z[4] Z[5], i^0 Z[6] Z[7], i^0 Z[7] Z[8], i^0 XXXXXXIII, i^0 IIIXXXXXX>
Xbars = PauliList(i^0 XXXXXXXXX)
Zbars = PauliList(i^0 ZZZZZZZZZ)

.. todo::
    This introduction needs to be finished.

Class Reference
===============

:class:`qecc.StabilizerCode`
----------------------------

.. autoclass:: qecc.StabilizerCode
    :members:
    :undoc-members:
    
