#---------------------------------------------------------------------------
#
# Logistic.py: chemical implementation of the logistic equation,
# comparison between ODE integration and Gillespie SSA
#
# logistic equation: dx/dt = r * x * (1 - x/K)
# where K is the maximum carrying capacity of the system
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

class LogisticGrowth:
    def __init__( self ):
        """ initialize logistic growth example """
        x0 = 0.01 # initial concentration for ODE
        popsize = 100 # max population size (in number of molecules) for SSA
        n0 = int(x0 * popsize) # initial number of molecules

        self.rsys = WellStirredVessel() # ODE vessel

        # SSA with single initial molecule n0=1, and Km=100, K=Km/NAV
        mynav = 1.0 * n0 / x0 # set NAV=NA*V such that x0 = n0 / NAV
        self.rssa = GillespieVessel(nav=mynav)

        reactionstrs = [
            # logistic equation for r=1, K = r/d = 1 => d=r
            "X --> 2 X , k=1", # kf=r
            "2 X --> X , k=1", # kr=d
        ]

        self.rsys.parse(reactionstrs)
        self.rssa.parse(reactionstrs)
        self.rsys.deposit('X', x0)
        self.rssa.inject('X', n0)
        print >> sys.stderr, "n0=", n0, "popsize=", popsize

    def run( self, finalvt ):
        """ run logistic growth example until simulation time finalvt """
        plotall = False
        ptime = 0.0
        niter = 0
        print "time\tODE\tSSA"
        self.trace_conc(ptime)
        while (self.rssa.vtime() <= finalvt):
            self.rssa.iterate()
            niter += 1
            dt = self.rssa.vtime() - ptime
            if (dt >= 0.1 or plotall):
                self.rsys.integrate(dt)
                ptime = self.rssa.vtime()
                self.trace_conc(ptime)
        km = self.rssa.get_nav() # Km=number of molecules for K=1.0
        print >> sys.stderr, "niter=", niter, "Km=", km

    def trace_conc( self, time ):
        """ produce tab-separated time vs. concentrations for plotting """
        xo = self.rsys.get_conc('X')
        xs = self.rssa.get_conc('X')
        print "%g\t%g\t%g" % (time, xo, xs)

if __name__ == '__main__':
    logist = LogisticGrowth()
    logist.run(50.0)
