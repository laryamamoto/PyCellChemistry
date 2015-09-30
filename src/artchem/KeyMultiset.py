#---------------------------------------------------------------------------
#
# KeyMultiset: a multiset indexed by a given key.
#
# the key is usually a substring of the molecule, akin to a binding site,
# but it can also be a shorter tag that is computed from the information
# within the molecule
#
# by Lidia Yamamoto, Belgium, October 2013
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

from Multiset import *

class KeyMultiset():
    def __init__( self ):
        """ initialize empty multiset """
	self.keymset = {} # key multiset is dictionary of multisets
	self.total = 0 # total amount of molecules in all multisets

    def inject( self, key, mol, mult=1 ):
        """ inject a given amount of a molecule in the multiset,
            indexed by the provided key
        """
	if (key == '' or mol == '' or mult < 1): return
	if (key in self.keymset):
            self.keymset[key].inject(mol, mult)
	else:
            self.keymset[key] = Multiset()
            self.keymset[key].inject(mol, mult)
	self.total += mult

    def expel( self, key, mol, mult=1 ):
        """ expel a given amount of a molecule from the multiset """
	if (key == '' or mol == '' or mult < 1): return
        if (key not in self.keymset): return
        mult = self.keymset[key].expel(mol, mult)
	self.total -= mult
	if (self.keymset[key].mult() <= 0): del self.keymset[key]

    def rndmol( self, key ):
        """ peek at a random molecule with given key, without removing it """
        if (key not in self.keymset): return ''
        return self.keymset[key].rndmol()

    def expelrnd( self, key ):
        """ expel a random molecule with given key """
        mol = self.rndmol(key)
        if (mol != ''): self.expel(key, mol)
        return mol

    def keys( self ):
        """ list of keys used to index molecules in the key multiset """
	return self.keymset.keys()

    def empty( self ):
        """ true if this multiset is empty """
	return len(self.keymset) == 0

    def mult( self, mol='' ):
        """ multiplicity: number of molecules of species 'mol' present
            in the multiset (regardless of key)
        """
	if (mol == ''): return self.total
        # CAUTION: no key provided, so must loop through keys, inefficient
        for k in self.keys():
            m = self.keymset[k].mult()
            if (m > 0): return m
        return 0

    def multk( self, key ):
        """ number of molecules (multiplicity) with the given key """
	if (key == '' or key not in self.keymset): return 0
        return self.keymset[key].mult()

    def nmolecules( self ):
        """ total number of molecules in the multiset """
	return self.total

    def nspecies( self ):
        """ number of distinct molecular species in the multiset """
        n = 0
        for k in self.keys():
            n += self.keymset[k].nspecies()
        return n

    def trace( self ):
        """ print the multiset object (typically for debug purposes) """
        for k in self.keys():
            print >> sys.stderr, ("key[%s]=" % k),
            self.keymset[k].trace()

