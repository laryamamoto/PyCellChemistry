#---------------------------------------------------------------------------
#
# Reaction.py: explicit representation of chemical reactions using multisets
#
# class Reaction: a reaction made of an educt multiset, a product
#                 multiset and a rate coefficient k
# class ReactionQueue: a FIFO queue of reactions
#
# by Lidia Yamamoto, Univ. Basel, Switzerland, January 2010
# June 2013: adapted to the PyCellChemistry package
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#   Copyright (C) 2015 Lidia A. R. Yamamoto
#   Contact: http://www.artificial-chemistries.org/
#
#   This file is part of PyCellChemistry.
#
#   PyCellChemistry is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   version 3, as published by the Free Software Foundation.
#
#   PyCellChemistry is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with PyCellChemistry, see file LICENSE. If not, see
#   http://www.gnu.org/licenses/
#

import sys
from collections import deque

#---------------------------------------------------------------------------
#
# Reaction: a chemical reaction is composed of an educt multiset, a
# product multiset, and a kinetic coefficient k that contributes to
# the speed of the reaction
#
#                       k
# for example: 2 A + B ----> 3 C
#
# in this reaction, the educt multiset is { A:2, B:1 }, and the
# product multiset is { C:3 }

class Reaction():
    def __init__( self, edmset, prmset, k=1.0 ):
        """ create a new reaction with given educt multiset edmset,
            product multiset prmset, and macroscopic kinetic coefficient k
        """
	self.educts = edmset
	self.products = prmset
        self.ratek = k

    def add_educt( self, mol, mult=1 ):
        """ add 'mult' copies of molecule 'mol' to the educt multiset """
	self.educts.inject(mol, mult)

    def add_product( self, mol, mult=1 ):
        """ add 'mult' copies of molecule 'mol' to the product multiset """
	self.products.inject(mol, mult)

    def del_educt( self, mol='', mult=1 ):
        """ delete 'mult' copies of molecule 'mol' from the educt multiset """
	self.educts.expel(mol, mult)

    def del_product( self, mol='', mult=1 ):
        """ delete 'mult' copies of molecule 'mol' from the product multiset """
	self.products.expel(mol, mult)

    def get_educts( self ):
        """ returns (a reference to) the educt multiset """
	return self.educts

    def get_products( self ):
        """ returns (a reference to) the product multiset """
	return self.products

    def set_coefficient( self, k ):
        """ set the kinetic coefficient of the reaction to k """
	self.ratek = k

    def get_coefficient( self ):
        """ get the kinetic coefficient of the reaction """
	return self.ratek

    def trace( self ):
        """ print reaction, for debug purposes """
	self.educts.trace()
	print >> sys.stderr, " --> ",
	self.products.trace()
	print >> sys.stderr, " , k=", self.ratek

#---------------------------------------------------------------------------
#
# ReactionQueue: an ordered list of reactions
#

class ReactionQueue( deque ):

    def add( self, reaction ):
        """ add reaction to list; returns position of reaction in queue """
	self.append(reaction)
        return len(self) - 1 

    def remove( self ):
        """ remove first reaction from list """
	return self.popleft()

    def peek( self, i=0 ):
        """ returns (a reference to) the reaction at position i """
	if (len(self) <= i):
            return None
	return self[i]

    def update( self, i, reaction ):
        """ update the reaction at position i with the new reaction
            object provided
        """
	if (len(self) <= i):
            return
        self[i] = reaction

    def empty( self ):
        """ True if the reaction list is empty """
	return len(self) == 0

    def nreactions( self ):
        """ total number of reactions in the list """
	return len(self)

    def trace( self ):
        """ print reaction list, for debugging """
	print >> sys.stderr, 'start rqueue'
	for r in self: r.trace()
	print >> sys.stderr, 'end rqueue'

