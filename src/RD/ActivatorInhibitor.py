#---------------------------------------------------------------------------
#
# ActivatorInhibitor.py: Gierer & Meinhardt's activator-inhibitor system
#
# Reference:
#
# [Koch & Meinhardt 1994] A. J. Koch and H. Meinhardt,
# "Biological pattern formation: from basic mechanisms to complex structures",
# Reviews of Modern Physics, volume 66, number 4, pages 1481-1508, October 1994.
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

class ActivatorInhibitorDemo( ):
    def __init__( self ):
        self.gsize = 32

        # following [Koch & Meinhardt 1994]:
        reactionstrs = [ "C + 2 A --> 3 A + C , k=0.01" ,
                         "  C + I --> I       , k=0.1", # ???
                         "    2 A --> 2 A + I , k=0.02",
                         #"        --> A       , k=0",
                         #"        --> I       , k=0",
                         "        --> C       , k=0.1", # ???
                         "      A -->         , k=0.01",
                         "      I -->         , k=0.02"
                       ]

        # inhibitor produces catalyst too (i.e. catalyst cannot grow
        # where there's no inhibitor):
        xxreactionstrs = [ "C + 2 A --> 3 A + C , k=0.01" ,
                         "  C + 2 I --> 2 I    , k=0.1", # ???
                         "    2 A --> 2 A + I , k=0.02",
                         #"        --> A       , k=0",
                         #"        --> I       , k=0",
                         "      I --> I + C    , k=0.1", # ???
                         "      A -->         , k=0.01",
                         "      I -->         , k=0.02"
                       ]

        self.rsys = ReactionDiffusionSystem(self.gsize, self.gsize)
        self.rsys.parse(reactionstrs)

        self.rsys.set_diffcoef("A", 0.005)
        self.rsys.set_diffcoef("C", 0.0)
        self.rsys.set_diffcoef("I", 0.2)

        # initial concentrations: A and I perturbed around homogeneous
        # steady state (A0=I0=C0=1), C constant at C=C0.

        self.rsys.resetAll("A", 1.0)
        self.rsys.resetAll("C", 1.0) # set it to 0 if C produced by I
        self.rsys.resetAll("I", 1.0)
        self.rsys.perturbAll("A", 0.33, 0.2)
        self.rsys.perturbAll("I", 0.33, 0.2)

        #self.rsys.set_colorscale(4.5)
	self.rsys.set_color( "A", (1, 0, 0) )
	self.rsys.set_color( "C", (0, 0, 1) )
	#self.rsys.set_color( "D", (0, 0, 1) )
	self.rsys.set_color( "I", (0, 1, 0) )

    def run( self, finalvt=3000.0, dt=0.1 ):
        i = 0
        self.rsys.animate(sleeptime=0.01, transparent=False, blur=True)
        while (self.rsys.vtime() <= finalvt):
            self.rsys.integrate(dt)
            i += 1
            if (i % 50 == 0):
                self.rsys.animate(sleeptime=0.01, transparent=False, blur=True)

if __name__ == '__main__':
    actinh = ActivatorInhibitorDemo()
    actinh.run()

