#---------------------------------------------------------------------------
#
# DimerStoch.py: an extended version of Dimer.py that allows the user to
# choose between stochastic or deterministic simulation algorithms
# (the default is a stochastic simulation using Gillespie SSA)
#
# A + B <--> A.B
#
# by Lidia Yamamoto, Belgium, August 2013
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

class Dimerization:

    def __init__( self, stoch ):
        """ reversible dimerization, stochastic (stoch=True) or
            deterministic (stoch=False)
        """
        self.stoch = stoch
        reactionstrs = [
            "A + B --> C , k=2",
            "C --> A + B , k=0.2"
        ]
        if self.stoch:
            self.reactor = GillespieVessel(nav=10)
        else:
            self.reactor = WellStirredVessel() 

        self.reactor.parse(reactionstrs)
        self.reactor.deposit('A', 2.0)
        self.reactor.deposit('B', 1.4)

    def run( self, finalvt=10.0, dt=0.1 ):
        """ execute dimerization reactions, until simulation time finalvt is
            reached; a plot line is generated every interval dt
        """
        self.reactor.trace_title()
        while (self.reactor.vtime() <= finalvt):
            self.reactor.trace_conc()
            self.reactor.integrate(dt)

if __name__ == '__main__':
    dimer = Dimerization(True)
    dimer.run()

