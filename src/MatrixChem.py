#---------------------------------------------------------------------------
#
# MatrixChem.py: matrix chemistry, acbook ch. 3
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

import artchem.BinaryStrings as bs
from artchem.Multiset import *

class MatrixChem:

    def __init__( self, nmols=1000 ):
        """ initialize matrix chemistry with 'nmols' random molecules """
        self.mset = Multiset() # "soup" where the molecules will float
        self.dimxy = 2 # dimxy=sqrt(N) <= 5 otherwise integer overflow
        self.nbits = self.dimxy ** 2 # nbits=N
        self.maxnspec = 2**self.nbits # number of possible molecular species
        self.debug = False # set this to True to print reaction traces

        for i in range(nmols): # initialize population with random molecules
            mol = self.randmol()
            self.mset.inject(mol)

    def randmol( self ):
        """ generate a random molecule (a molecule is an N-bit integer) """
        return bs.randbin(self.nbits)

    def squash( self, n ):
        """ squash function sigma:
            returns one bit out of the given integer n
        """
        # set threshold to half of the possible scalar product sum,
        # such that half of the non-binary values are assigned bit
        # value 0, the other half are assigned bit value 1
        if (n < self.dimxy / 2):
            bit = 0
        else:
            bit = 1
        return bit

    def fold( self, mol, transpose=False ):
        """ fold binary string 'mol' into matrix operator
            transpose=False: folds horizontally,
                             e.g. [ [s1 s2] [s3 s4] ] for N=4
            transpose=True:  folds vertically,
                             e.g. [ [s1 s3] [s2 s4] ] for N=4
        """
        foldm = np.matrix(np.zeros((self.dimxy, self.dimxy)), int)
        for x in range(self.dimxy):
            for y in range(self.dimxy):
                bit = (mol % 2)  # read least significant bit
                mol = (mol >> 1) # consume it
                if (transpose):
                    foldm[y,x] = bit
                else:
                    foldm[x,y] = bit
        return foldm

    def unfold( self, fmat, transpose=False ):
        """ unfold integer matrix into binary string, applying squash function
            transpose=False: unfolds horizontally
            transpose=True:  unfolds vertically
        """
        mol = 0
        # we read the most significant bits first (from the end of the
        # matrix), and then shift them to the left to make room for
        # the lower bits
        for x in range(self.dimxy - 1, -1, -1): # for (dimxy - 1) downto 0
            for y in range(self.dimxy - 1, -1, -1):
                if (transpose):
                    num = fmat[y,x]
                else:
                    num = fmat[x,y]
                bit = self.squash(num)
                mol = (mol << 1) + bit
        return mol

    def react( self, m1, m2 ):
        """ chemical reaction: apply operator m1 to operand m2:
                               P(m1) x m2 => m3
            method:
            - fold the operator horizontally obtaining matrix fm
            - fold the operand vertically obtaining matrix op
            - multiply both matrices: the resulting matrix contains
              the result in a vertically-folded form
            - unfold the result matrix to obtain the product of the reaction

           the method adopted here looks different from the original,
           but actually it is fully equivalent and it leads to exactly
           the same result
        """
        fm = self.fold(m1)       # fold m1 as [ [s1 s2] [s3 s4] ]
        op = self.fold(m2, True) # fold m2 as [ [s1 s3] [s2 s4] ]
        res = fm * op            # apply standard matrix multiplication
        m3 = self.unfold(res, True) # result is in transposed form, unfold it
        #print "m1=", m1, "m2=", m2, "m3=", m3
        #print "fm=", fm, "op=", op, "res=", res
        return m3

    def reacttable( self ):
        """ print all-to-all reaction table """
        print "   |",
        for j in range (self.maxnspec):
            print "%2d" % j,
        print
        print "------------------------------------------------------------"
        for i in range (self.maxnspec):
            print "%2d |" % i,
            for j in range (self.maxnspec):
                m = self.react(i, j)
                print "%2d" % m,
            print

    def trace( self ):
        """ print debug traces with multiset information """
        if self.debug:
            print "nspec=", self.mset.nspecies(),
            print "nmol=", self.mset.nmolecules()
            print "mset=",
            self.mset.trace()

    def trace_title( self, prefix='' ):
        """ output tab-separated title line for plotting molecule counts;
            if prefix is provided, preceeds all molecule names with it
        """
        tline = []
        for mol in range(self.maxnspec):
            tline.append("%s" % mol)
        tline.insert(0, "time")
        sep = ("\t%s" % prefix)
        print sep.join(tline)

    def trace_mult( self, time ):
        """ output tab-separated line with all molecule counts """
        tline = [ "%g" % time ]
        for mol in range(self.maxnspec):
            tline.append("%g" % self.mset.mult(mol))
        print "\t".join(tline)

    def run( self, niter=10000 ):
        """ run random molecular collisions for 'niter' iterations """
        neffect = 0   # count number of effective collisions
        npoints = 100 # max number of data points to be plotted (per species)
        plotinterv = niter / npoints # plot interval

        if self.debug:
            self.reacttable()
            self.trace()
        self.trace_title('M')
        for i in range(niter):

            # expel two random molecules from the soup (although these
            # are actually not consumed and will therefore be
            # reinjected later, expelling them prevents us from
            # accidentally picking the same molecule twice)
            m2 = self.mset.expelrnd()
            m1 = self.mset.expelrnd()
            if self.debug:
                print "iter=%d: %d + %d =>" % (i, m1, m2),

            if (m1 == 0 or m2 == 0):
                # zero molecules can quickly become invaders if left
                # to react, so reactions with them are made elastic
                if self.debug:
                    print "elastic"
            else:
                m3 = self.react(m1, m2)   # perform the reaction
                m4 = self.mset.expelrnd() # discard a random molecule
                self.mset.inject(m3)      # inject result in soup
                neffect += 1
                if self.debug:
                    print "%d , del %d " % (m3, m4)

            self.mset.inject(m1)
            self.mset.inject(m2)
            if (i % plotinterv == 0):
                self.trace_mult(i)

        self.trace()
        if self.debug:
            print "neffective=", neffect

if __name__ == '__main__':
    matchem = MatrixChem()
    matchem.run()

