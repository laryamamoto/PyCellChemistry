#---------------------------------------------------------------------------
#
# Chameleons.py: colored chameleons, acbook ch. 2
#
# by Lidia Yamamoto, Belgium, July 2013
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

from artchem.Multiset import *

class Chameleons:

    def __init__( self, N=2700 ):
        """ initial population of N individuals of 2 colors """
        self.mset = Multiset()
        self.mset.inject(1, N)
        self.mset.inject(2, N)

    def trace_mult( self, time ):
        """ print molecule counts per species at the given time """
        n1 = self.mset.mult(1)
        n2 = self.mset.mult(2)
        n3 = self.mset.mult(3)
        print "%g\t%d\t%d\t%d" % (time, n1, n2, n3)

    def run( self, niter=12000 ):
        """ run the Chameleons AC for niter iterations """
        print "time\tR\tG\tB"
        self.trace_mult(0)
        for i in range(niter):
            m1 = self.mset.expelrnd()
            m2 = self.mset.expelrnd()
            if (m1 != m2): # effective collision, e.g. r + g -> 2b
                m3 = 6 - (m1 + m2)
                self.mset.inject(m3, 2)
            else: # elastic collision
                self.mset.inject(m1)
                self.mset.inject(m2)
            self.trace_mult(i+1)

if __name__ == '__main__':
    chams = Chameleons()
    chams.run()

