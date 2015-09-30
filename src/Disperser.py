#---------------------------------------------------------------------------
#
# Disperser.py: implementation of Thomas Meyer's disperser protocol on
# a fixed (example) network topology.
#
# This protocol implements a discretized mass-conserving stochastic
# diffusion process over an arbitrary network of interconnected computers.
# Here an example network with 4 nodes is shown.
#
# References:
#
#   T. Meyer and C. Tschudin. Chemical networking protocols. Proc. 8th
#   ACM Workshop on Hot Topics in Networks (HotNets-VIII), Oct. 2009.
#
#   T. Meyer, L. Yamamoto, and C. Tschudin. An artificial chemistry
#   for networking. Bio-Inspired Computing and Communication, volume
#   5151 of LNCS, pages 45-57. Springer, 2008.
#
# Python implementation by Lidia Yamamoto, Belgium, June 2013
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

from artchem.ReactionVessel import *

class Disperser:

    def __init__( self ):
        """ create networked topology inside single reaction vessel:
            multiple vessels are mimicked via different groups of molecular
            species
        """
        topo = [ [1, 2], [1, 3], [1, 4], [2, 3] ]
        (catalysts, reactionstrs) = self.mkreactions(topo)
 
        # create reactor, choosing N_A*V such that c=k=1 overall
	self.reactor = GillespieVessel(nav=1)

        self.reactor.parse(reactionstrs)
        self.reactor.inject('X4', 1000)
        for c in catalysts:
            self.reactor.inject(c, 1)

    def run( self, finalvt ):
        """ run the disperser protocol until simulation time finalvt """
        plotall = False
        self.reactor.trace_title('X')
        self.reactor.trace_mult('X')
        ptime = 0.0
        injected = False
        expelled = False
        stchange = False
        while (self.reactor.vtime() <= finalvt):
            t = self.reactor.vtime()
            self.reactor.iterate()
            dt = t - ptime
            if (t >= 10 and not injected):
                self.reactor.inject('X2', 600)
                injected = True
                stchange = True
            if (t >= 20 and not expelled):
                self.reactor.expel('X2', 300)
                expelled = True
                stchange = True
            if (dt >= 0.1 or plotall or stchange):
                self.reactor.trace_mult('X')
                ptime = t
                stchange = False

    def mkreactions( self, topo ):
        """ constructs the list of reactions (with their catalysts)
            automatically for a given topology 'topo';
            'topo' is a list of links between 2 nodes (reactors);
            each link is described by a pair [id1, id2] of node ids
            identifying the 2 nodes to be interconnected.
            CAUTION: only works with node ids made of single digits!!
        """
        clist = []
        rlist = []
        for lnk in topo:
            n0 = lnk[0]
            n1 = lnk[1]
            c0 = "C%d%d" % (n0, n1)
            c1 = "C%d%d" % (n1, n0)
            r0 = "%s + X%d --> %s + X%d" % (c0, n0, c0, n1)
            r1 = "%s + X%d --> %s + X%d" % (c1, n1, c1, n0)
            clist.append(c0)
            clist.append(c1)
            rlist.append(r0)
            rlist.append(r1)
        return (clist, rlist)

if __name__ == '__main__':
    disperser = Disperser()
    disperser.run(30.0)

