#---------------------------------------------------------------------------
#
# Multiset.py: a multiset class
#
# a multiset is a bag of objects (molecules in our context), like a
# set where every object may occur more than once
#
# Example: { A:2, B:3 } represents a bag containing 2 molecules of
# species (type) A and 3 molecules of type B
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
import numpy as np

class Multiset():
    def __init__( self ):
        """ initialize empty multiset """
	self.mset = {} # multiset implemented as a python dictionary
	self.total = 0 # total amount of molecules in the multiset

    def inject( self, mol, mult=1 ):
        """ inject a given amount of a molecule in the multiset """
	if (mol == '' or mult < 1): return
	if (mol in self.mset): self.mset[mol] += mult
	else: self.mset[mol] = mult
	self.total += mult

    def expel( self, mol, mult=1 ):
        """ expel a given amount of a molecule from the multiset;
            returns the number of molecules actually expelled
        """
	if (mol == '' or mult < 1 or mol not in self.mset): return 0
	if (mult > self.mset[mol]): mult = self.mset[mol]
	self.mset[mol] -= mult
	self.total -= mult
	if (self.mset[mol] <= 0): del self.mset[mol]
        return mult

    def clear (self ):
        """ delete all molecules from the multiset """
        self.mset.clear()
        self.total = 0

    def rndmol( self ):
        """ peek at a random molecule from the multiset, without
            removing it
        """
        if (self.total <= 0):
            return ''
        molid = np.random.randint(self.total)
        for mol in self.keys():
            m = self.mult(mol)
            if (molid < m):
                return mol
            molid -= m
        return ''

    def expelrnd( self ):
        """ expel a random molecule from the multiset """
        mol = self.rndmol()
        self.expel(mol)
        return mol

    def keys( self ):
        """ list containing the set of different objects in the multiset """
	return self.mset.keys()

    def empty( self ):
        """ true if this multiset is empty """
	return len(self.mset) == 0

    def mult( self, mol='' ):
        """ multiplicity: number of molecules of type 'mol' present
            in the multiset
        """
	if (mol == ''): return self.total
	elif (mol in self.mset): return self.mset[mol]
	else: return 0

    def nmolecules( self ):
        """ total number of molecules in the multiset """
	return self.total

    def nspecies( self ):
        """ number of distinct molecular species in the multiset """
	return len(self.mset)

    def absorb( self, victim ):
        """ absorb victim multiset into my own (victim will return empty) """
        k = victim.keys()
        for mol in k:
            n = victim.mult(mol)
            victim.expel(mol, n)
            self.inject(mol, n)

    def divide( self ):
        """ divide multiset in two parts; returns daughter multiset created """
        newmset = Multiset()
        n = self.total / 2  # to do: choose n randomly?
        for i in range(n):
            mol = self.expelrnd()
            newmset.inject(mol)
        return newmset

    def topkeys( self ):
        """ list of molecules in the multiset, ordered by decreasing
            multiplicity
        """
        return sorted(self.mset, key=self.mset.get, reverse=True)

    def trace( self ):
        """ print the multiset object (typically for debug purposes) """
	#print "%s m=%d n=%d" % ( self.mset, self.mult(), self.nspecies() ),
	print >> sys.stderr, "%s" % ( self.mset )

    def trace_frequencies( self ):
        """ print a histogram of molecule multiplicities """
        for k in self.topkeys():
            print >> sys.stderr, k, self.mult(k)

def rndpolymer(alphabet, length):
    """ generate a random string of given length with characters from the
        given alphabet
    """
    s = ''
    for i in range(length):
        n = np.random.randint(len(alphabet))
        s += alphabet[n]
    return s

