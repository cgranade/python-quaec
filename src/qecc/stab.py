#!/usr/bin/python
# -*- coding: utf-8 -*-
##
# stab.py: Classes and methods encapsulating and manipulating stabilizer codes.
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

import operator as op

import PauliClass as p # Sorry for the confusing notation here.
import CliffordClass as c
import paulicollections as pc


from singletons import EmptyClifford, Unspecified

import warnings

## ALL ##

__all__ = [
    'StabilizerCode'
]

## CLASSES ##

class StabilizerCode(object):
    r"""
    TODO
    """
    
    def __init__(self, group_generators, logical_xs, logical_zs):
        self.group_generators = pc.PauliList(*group_generators)
        self.logical_xs = pc.PauliList(*logical_xs)
        self.logical_zs = pc.PauliList(*logical_zs)
        
    @property
    def nq(self):
        """
        TODO
        """
        return len(iter(gen for gen in self.group_generators if gen is not Unspecified).next())
        
    @property
    def n_constraints(self):
        return len(self.group_generators)
        
    @property
    def n_logical(self):
        return self.nq - self.n_constraints
        
    @property
    def distance(self):
        raise NotImplementedError("Not yet implemented.")
        
    def encoding_cliffords(self):
        C = c.Clifford(self.logical_xs + self.group_generators, self.logical_zs + ([Unspecified] * self.n_constraints))
        return C.constraint_completions()

    def concatenate(self,other):
        """
        Returns the stabilizer for a concatenated code, given the 
        stabilizers for two codes. At this point, it only works for two
        k=1 codes.
        """
        nq_self=self.nq
        nq_other=other.nq
        nq_new=nq_self*nq_other
        stab_new=StabilizerCode([],[],[])
        for k in range(nq_other):
            for low_level_gen in self.group_generators:
                stab_new.group_generators.append(p.Pauli('I'*k*nq_self+low_level_gen.op+'I'*(nq_new-nq_self*(k+1))))
        for phys_gen in other.group_generators:
            log_gen_op=''
            log_gen_ph=0
            for letter in phys_gen.op:
                if letter == 'X':
                    log_gen_op+=self.logical_xs[0].op
                elif letter == 'Z':
                    log_gen_op+=self.logical_zs[0].op
                elif letter == 'Y':
                    log_gen_op+=(self.logical_xs[0]*self.logical_zs[0]).op
                    log_gen_ph+=1
                elif letter == 'I':
                    log_gen_op+='I'*nq_self
            stab_new.group_generators.append(p.Pauli(log_gen_op,log_gen_ph))
        for phys_log in other.logical_xs:
            log_log_op=''
            log_log_ph=0
            for letter in phys_log.op:
                if letter == 'X':
                    log_log_op+=self.logical_xs[0].op
                elif letter == 'Z':
                    log_log_op+=self.logical_zs[0].op
                elif letter == 'Y':
                    log_log_op+=(self.logical_xs[0]*self.logical_zs[0]).op
                    log_log_ph+=1
                elif letter == 'I':
                    log_log_op+='I'*nq_self
            stab_new.logical_xs.append(p.Pauli(log_log_op,log_log_ph))
        for phys_log in other.logical_zs:
            log_log_op=''
            log_log_ph=0
            for letter in phys_log.op:
                if letter == 'X':
                    log_log_op+=self.logical_xs[0].op
                elif letter == 'Z':
                    log_log_op+=self.logical_zs[0].op
                elif letter == 'Y':
                    log_log_op+=(self.logical_xs[0]*self.logical_zs[0]).op
                    log_log_ph+=1
                elif letter == 'I':
                    log_log_op+='I'*nq_self
            stab_new.logical_zs.append(p.Pauli(log_log_op,log_log_ph))
        return stab_new

    @staticmethod
    def bit_flip_code(x_dist):
        nq=2*x_dist+1
        return StabilizerCode(
            ['I'*j+'ZZ'+'I'*(nq-j-2) for j in range(nq-1)],
            ['X'*nq], ['Z'*nq]

        )
    @staticmethod
    def phase_flip_code(z_dist):
        nq=2*z_dist+1
        return StabilizerCode(
            ['I'*j+'XX'+'I'*(nq-j-2) for j in range(nq-1)],
            ['X'*nq], ['Z'*nq]
        )

    @staticmethod
    def perfect_5q_code():
        return StabilizerCode(
            [
                'XZZXI',
                'IXZZX',
                'XIXZZ',
                'ZXIXZ'
            ],
            ['XXXXX'], ['ZZZZZ']
        )

    @staticmethod
    def steane_code():
        return StabilizerCode(
            [
                'XXXXIII',
                'XXIIXXI',
                'XIXIXIX',
                'ZZZZIII',
                'ZZIIZZI',
                'ZIZIZIZ'                
            ],
            ['XXXXXXX'], ['ZZZZZZZ']
        )

    @staticmethod
    def shor_code():
        return StabilizerCode.bit_flip_code(1).concatenate(StabilizerCode.phase_flip_code(1))        

    @staticmethod
    def reed_muller_code(r,t):
        raise NotImplementedError("Coming Soon: Reed-Muller Codes")
    @staticmethod
    def reed_solomon_code(r,t):
        raise NotImplementedError("Coming Soon: Reed-Solomon Codes")
