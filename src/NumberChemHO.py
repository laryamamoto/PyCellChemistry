#---------------------------------------------------------------------------
#
# NumberChemHO.py:
# a reimplementation of NumChem.py using HighOrderChem.py
#
# by Lidia Yamamoto, Belgium, May 2014
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

from HighOrderChem import *

def isprime( n ):
    """ check if n is a prime, adapted from:
       http://www.daniweb.com/software-development/python/code/216880/
              check-if-a-number-is-a-prime-number-python
    """
    if (n <= 0 or n % 2 == 0):
        return False
    if (n < 3):
        return True
    for i in range(3, int(n**0.5)+1, 2):
        if (n % i == 0):
            return False
    return True

class NumberChemHO( HighOrderChem ):
    def __init__( self, popsize=100, minn=2, maxn=1000 ):
        """ create a random soup of 'popsize' numbers ranging from 'minn'
            to 'maxn'
        """
	HighOrderChem.__init__( self )
        for i in range(popsize):
            dice = np.random.randint(minn, maxn+1)
            self.mset.inject(dice)
        #rule = 'r[%d,%d] self.divrule'
        rule = 'self.divrule(%d,%d)'
        self.rset.inject(rule, 4)
        #self.rset.trace()

    def divrule(self, m1, m2):
        """ division rule from NumberChem.py, but now as a separate function """
        p = [m1, m2]
        if m1 > m2 and m1 % m2 == 0:
            p[0] = m1 / m2
        elif m2 > m1 and m2 % m1 == 0:
            p[1] = m2 / m1
        return p

    def nprimes( self ):
        """ count the number of prime numbers in the soup """
        np = 0
        for mol in self.mset.keys():
            m = self.mset.mult(mol)
            if (isprime(mol)):
                np += m
        return np

    def trace_title( self ):
        """ print a tab-separated title line for plotting """
        print "iter\tspecies\teffect\tprimes\tcollision"

    def trace( self, i, reaction ):
        """ print a tab-separated data line for plotting """
        print "%d\t%d\t%d\t%d\t%s" % (i, self.mset.nspecies(), \
              self.neffect, self.nprimes(), reaction)

    def run( self, niter=10000 ):
        """ run 'niter' iterations of the number chemistry,
            by binding as many molecules as needed at once
        """
        self.neffect = 0  # number of effective collisions
        self.trace_title()
        self.trace(0, '')
        for i in range(niter):
            (rule, educts, products) = self.iterate()
            reaction = str(educts).strip('[ ]')
            if self.is_effective(educts, products):
                self.neffect += 1
                reaction += "-->" + str(products).strip('[ ]')
            self.trace(i+1, reaction)

if __name__ == '__main__':
    numchem2 = NumberChemHO()
    numchem2.run()

