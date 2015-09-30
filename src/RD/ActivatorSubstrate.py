#---------------------------------------------------------------------------
#
# ActivatorSubstrate.py: activator-depleted substrate model
#
# from p. 39 of:
# Hans Meinhardt, "Models of Biological Pattern Formation", 
# Academic Press, London, UK, 1982
#
# originally implemented on top of the breve simulator, www.spiderland.org
# by Lidia Yamamoto, Univ. Basel, Switzerland, January 2010
# 20150914: removed breve dependencies to run within PyCellChemistry
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

from ReactionDiffusion import *

class ActivatorSubstrateDemo( ):
    def __init__( self ):
        self.gsize = 32

        reactionstrs = [ "S + 2 A --> 3 A" ,
                         "S + 2 A --> 2 A" ,
                        #"        --> A  " , # k=0
                         "        --> S  " ,
                         "      A -->    " ,
                         "      S -->    "   ]

        self.rsys = ReactionDiffusionSystem(self.gsize, self.gsize)
        self.rsys.parse(reactionstrs)

        # Koch's parameters

        # diffusion coefficients from ActivatorInhibitor.py:
        DS = 0.2
        #DA = 0.005
        DA = 0.008

        self.rsys.set_coefficient(0, 0.01) # 0.01 rho_a
        self.rsys.set_coefficient(1, 0.01) # k = rho_s - rho_a, rho_s=0.02
        self.rsys.set_coefficient(2, 0.02) # sigma_s
        self.rsys.set_coefficient(3, 0.01) # mu_a
        self.rsys.set_coefficient(4, 0)    # mu_s

        # actually it seems that S=0 is sufficient to start
        #self.rsys.deposit("S", 1.0, self.gsize/2, self.gsize/2)
        #self.rsys.resetAll("S", 1.0)
        self.rsys.resetAll("A", 1.0)

        self.rsys.set_diffcoef("S", DS )
        self.rsys.set_diffcoef("A", DA )
        #self.rsys.perturbAll("S", 0.333, 0.2) # no need to perturb S
        self.rsys.perturbAll("A", 0.333, 0.2)

	self.rsys.set_color( "A", (1, 0, 0) )
	self.rsys.set_color( "S", (0, 0, 1) )

    def run( self, finalvt=3000.0, dt=0.1 ):
        i = 0
        self.rsys.animate(sleeptime=0.01, transparent=False, blur=True)
        while (self.rsys.vtime() <= finalvt):
            self.rsys.integrate(dt)
            i += 1
            if (i % 50 == 0):
                self.rsys.animate(sleeptime=0.01, transparent=False, blur=True)

if __name__ == '__main__':
    actsubs = ActivatorSubstrateDemo()
    actsubs.run()

