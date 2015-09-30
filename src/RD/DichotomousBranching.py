#---------------------------------------------------------------------------
#
# DichotomousBranching.py: dichotomous branching pattern based on the
# combination of an activator-depleted substrate model and a vein
# formation substance Y.
#
# Y formation is triggered by the activator A. The vascular system to
# be formed drains the substrate S, such that the catalysis of
# activator is inhibited due to depletion of S, thus inhibiting
# further vein formation in this area.
#
# from p. 170-172 of:
# Hans Meinhardt, "Models of Biological Pattern Formation", 
# Academic Press, London, UK, 1982
#
# originally implemented on top of the breve simulator, www.spiderland.org
# by Lidia Yamamoto, Univ. Basel, Switzerland, February 2010
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

class DichotomousBranchingDemo( ):
    def __init__( self ):
        self.gsize = 64

        reactionstrs = [ "S + 2 A --> 3 A" ,
                         "S + 2 A --> 2 A" ,
                        #"        --> A  " , # k=0
                         "        --> S  " ,
                         "      A -->    " ,
                         "      S -->    " ,

                         "  S + Y --> Y  " ,     # vein drains substrate
                         "E + 2 Y --> Z  " ,     # vein signal is produced
                         "      Z --> E + 3 Y" , # by enzyme E
                         "  A + Y --> A + 2 Y" , # act. produces vein sig.
                         #"  A --> A + Y" , # act. produces vein sig.
                         "      Y -->        "   # vein sig. decays
                       ]


        self.rsys = ReactionDiffusionSystem(self.gsize, self.gsize)
        self.rsys.parse(reactionstrs)

        # parameters from ActivatorSubstrate.py: 
        DS = 0.2
        DA = 0.005
        self.rsys.set_coefficient(0, 0.01)  # rho_a
        self.rsys.set_coefficient(1, 0.01)  # k = rho_s - rho_a, rho_s=0.02
        self.rsys.set_coefficient(2, 0.02)  # sigma_s
        self.rsys.set_coefficient(3, 0.01)  # mu_a
        self.rsys.set_coefficient(4, 0.001) # mu_s

        self.rsys.set_coefficient(5, 0.03)  # epsilon
        self.rsys.set_coefficient(6, 0.02)  # k1
        self.rsys.set_coefficient(7, 0.01)  # k2 << k1
        self.rsys.set_coefficient(8, 0.02)  # d (differentiation rate)
        self.rsys.set_coefficient(9, 0.01)  # e (decay of y)

        self.rsys.resetAll("A", 1.0)
        #self.rsys.deposit("E", 0.01, self.gsize-1, self.gsize-1)
        #self.rsys.deposit("Y", 0.01, self.gsize-1, self.gsize-1)
        self.rsys.resetAll("E", 0.01)
        self.rsys.resetAll("Y", 0.01)

        self.rsys.set_diffcoef("S", DS)
        self.rsys.set_diffcoef("A", DA)
        self.rsys.perturbAll("A", 0.333, 0.2)
        #self.rsys.set_patch_at("A", 1.0, self.gsize/2, self.gsize/2)
        #self.rsys.set_patch_at("A", 1.0, self.gsize/2-3, self.gsize/2-3)
        #self.rsys.set_patch_at("A", 1.0, self.gsize/2-3, self.gsize/2+3)
        #self.rsys.set_patch_at("A", 1.0, self.gsize/2+3, self.gsize/2-3)
        #self.rsys.set_patch_at("A", 1.0, self.gsize/2+3, self.gsize/2+3)

	self.rsys.set_color( "A", (1, 0, 0) )
	self.rsys.set_color( "S", (0, 0, 1) )
	self.rsys.set_color( "Y", (0, 1, 0) )

    def run( self, finalvt=3000.0, dt=0.1 ):
        i = 0
        self.rsys.animate(sleeptime=0.01, transparent=False, blur=True)
        while (self.rsys.vtime() <= finalvt):
            self.rsys.integrate(dt)
            i += 1
            if (i % 50 == 0):
                self.rsys.animate(sleeptime=0.01, transparent=False, blur=True)

if __name__ == '__main__':
    dicbranch = DichotomousBranchingDemo()
    dicbranch.run()

