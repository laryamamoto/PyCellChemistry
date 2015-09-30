#---------------------------------------------------------------------------
#
# Repressilator.py: stochastic simulation of the repressilator
#
# the repressilator is a genetic regulatory network made of 3 genes that
# repress each other in a cycle, leading to oscillations in gene expression
# dynamics
#
# by Lidia Yamamoto, Belgium, June 2013
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

class Repressilator:
    def __init__( self ):
        """ initialize repressilator simulation (stochastic simulation by
            default): initially one gene is activated and the other two are
            repressed
        """
        usestoch = True
        reactionstrs = [
            # Protein represses other gene's expression
            "2 pa + gb --> cb , k=1",
            "2 pb + gc --> cc , k=1",
            "2 pc + ga --> ca , k=1",
            "cb --> 2 pa + gb , k=1",
            "cc --> 2 pb + gc , k=1",
            "ca --> 2 pc + ga , k=1",
            
            # Transcription (gene -> mRNA):
            "ga --> ga + ma  , k=5",
            "gb --> gb + mb  , k=5",
            "gc --> gc + mc  , k=5",

            # Leaky transcription
            "ca --> ca + ma , k=0",
            "cb --> cb + mb , k=0",
            "cc --> cc + mc , k=0",

            # Translation (mRNA -> Protein):
            "ma --> ma + pa , k=1",
            "mb --> mb + pb , k=1",
            "mc --> mc + pc , k=1",

            # Decay of mRNA and Protein:
            "ma -->  ,  k=0.5",
            "mb -->  ,  k=0.5",
            "mc -->  ,  k=0.5",
            "pa -->  ,  k=0.1",
            "pb -->  ,  k=0.1",
            "pc -->  ,  k=0.1",
        ]

	self.reactor = GillespieVessel()
	self.reactor.parse(reactionstrs)

        nav = 100
        self.reactor.set_nav(nav)

        # initial conditions: one gene is activated, the other two are repressed
        self.reactor.inject('ga',  1 * nav)
        self.reactor.inject('cb', 1 * nav)
        self.reactor.inject('cc', 1 * nav)
        #self.reactor.trace()

    def run( self, finalvt=200 ):
        """ run the repressilator until the simulation time reaches finalvt """
        plotall = False
        self.reactor.trace_title()
        self.reactor.trace_mult()
        ptime = 0.0
        while (self.reactor.vtime() <= finalvt):
            t = self.reactor.vtime()
            self.reactor.iterate()
            dt = t - ptime
            if (dt >= 0.2 or plotall):
                self.reactor.trace_mult()
                ptime = t

if __name__ == '__main__':
    repress = Repressilator()
    repress.run()
