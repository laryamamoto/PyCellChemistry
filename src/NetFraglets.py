#---------------------------------------------------------------------------
#
# NetFraglets: networked examples of Fraglets, a chemically-inspired
# programming language for computer networks
#
# so far only one example available:
# CDP: confirmed delivery protocol
#
# Reference:
#
# C. Tschudin. Fraglets: A metabolistic execution model for communication
# protocols. Proc. 2nd Annual Symposium on Autonomous Intelligent Networks
# and Systems (AINS), July 2003.
#
# python implementation by Lidia Yamamoto, Belgium, October 2013
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

from Fraglets import *
from artchem.Cell import *

def cdp():
    """ CDP: confirmed delivery protocol """
    niter = 10
    n1 = Fraglets('a')
    n2 = Fraglets('b')
    n1.add_cnx('b', n2)
    n2.add_cnx('a', n1)
    n1.inject(n1.parse('matchp c send b split send a k *'))
    #n1.inject('cd')
    net = Cell()
    net.add(n1)
    net.add(n2)
    #net.run(10)

    print  >> sys.stderr, "INIT:"
    n1.trace()
    n2.trace()
    cnt = 0
    for i in range(niter):
        print >> sys.stderr, "ITER=", i
        if i % 4 == 0:
            mol = 'cd' + str(cnt)
            n1.inject(mol)
            cnt = (cnt + 1) % 10
        net.propensity()
        net.gillespie()

    print  >> sys.stderr, "END:"
    n1.trace()
    n2.trace()

if __name__ == '__main__':
    cdp()

