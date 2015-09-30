#---------------------------------------------------------------------------
#
# Cell.py: hierarchical compartments containing molecules that react,
# and possibly other compartments too, like in P Systems
#
# by Lidia Yamamoto, July 2013
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

class Cell():

    def __init__( self ):
        """ initialize empty cell """
        self.mycells = [] # list of internal compartments within this cell
        self.cyto = {}    # molecules in the cytoplasm of this cell
        self.parent = {}  # parent cell
        self.prop = [] # list of propensities, from cytoplast to subcells
        self.wt = 0.0  # sum of propensities
        self.vt = 0.0  # virtual time

    def set_cytoplasm( self, c):
        """ assign a cytoplasm object to this cell; a cytoplasm represents
            the set of molecules floating inside the cell, which are not
            inside another subcell; it is typically implemented as a Multiset
            or a more elaborated class derived from it.
        """
        self.cyto = c

    def get_cytoplasm( self ):
        """ get cytoplasm object for this cell """
        return self.cyto

    def add( self, c ):
        """ add a new internal compartment c to this cell """
        self.mycells.append(c)

    def get_cells():
        """ get list of internal compartments for this cell """
        return self.mycells

    def clear_cells():
        """ delete all the internal compartments from this cell """
        self.mycells = []

    def dissolve( self, i ):
        """ dissolve internal compartment [i], and release its content to the
            cytoplasm
        """
        icell = self.mycells[i]
        icyto = icell.get_cytoplasm()
        if icyto != {}:
            self.cyto.absorb(icyto)
        icomp = icell.get_cells()
        if icomp != []:
            for c in icomp:
                self.add(c)
            icell.clear_cells()

    def divide( self ):
        """ cell division: approximately half of the objects (chosen randomly)
            goes to the daughter cell; returns the newly created daughter cell
        """
        newcell = Cell() # create daughter cell
        if self.cyto != {}:
            # divide cytoplasm molecules evenly among the two cells
            newcyto = self.cyto.divide()
            newcell.set_cytoplasm(newcyto)
        if self.mycells != []:
            # send half of the internal compartments to daughter cell
            n = len(self.mycells) / 2
            for i in range(n):
                p = np.random.randint(len(self.mycells))
                c = self.mycells.pop(p)
                newcell.add(c)

    def inert( self ):
        """ true if none of the vessels in this cell (including the
            cytoplasm and all the subcells, recursively) has reactions to
            perform
        """
        if self.cyto != {}:
            if not self.cyto.inert(): return False
        for cell in self.mycells:
            if not cell.inert(): return False
        return True

    def propensity( self ):
        """ calculate propensities for the elements within this cell """
        self.wt = 0.0
        self.prop = []
        if self.cyto != {}:
            w = self.cyto.propensity()
            self.prop.append(w)
            self.wt += w
        for cell in self.mycells:
            w = cell.propensity()
            self.prop.append(w)
            self.wt += w
        return self.wt

    def gillespie( self ):
        """ hierarchical gillespie for cell and its subcells """
        if self.wt <= 0: return
        w = np.random.random() * self.wt
        i = 0
        if self.cyto != {}:
            if self.prop[0] > 0 and w < self.prop[0]:
                self.cyto.react(w)
                return
            w -= self.prop[0]
            i += 1
        for cell in self.mycells:
            if (self.prop[i] > 0 and w < self.prop[i]):
                cell.react(w)
                return
            w -= self.prop[i]
            i += 1
        # pending: self.vt += ... IN TOP CELL ONLY!!

    def run( self, niter ):
        """ run 'niter' iterations of gillespie for this cell, or until cell
            becomes inert
        """
        for i in range(niter):
            print >> sys.stderr, "ITER=", i
            self.propensity()
            if self.inert(): return
            self.gillespie()

