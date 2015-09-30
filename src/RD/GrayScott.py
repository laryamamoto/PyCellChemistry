#---------------------------------------------------------------------------
#
# GrayScott.py: Gray Scott Model of Reaction Diffusion
#
# based on breve/demos/Chemistry/GrayScott.tz by Jon Klein, www.spiderland.org
# and on http://groups.csail.mit.edu/mac/projects/amorphous/GrayScott/
#
# by Lidia Yamamoto, Univ. Basel, Switzerland, January 2010
# 20150910: ported from breve to numpy+vpython for use within PyCellChemistry
#
# References:
#
# [Pearson1993] John E. Pearson, "Complex Patterns in a Simple System",
# Science, volume 261, number 5118, July 1993, pages 189-192.
# doi=10.1126/science.261.5118.189
#
# [Mazin1996] W. Mazin, K. E. Rasmussen, E. Mosekilde, P. Borckmans,
# and G. Dewel, "Pattern formation in the bistable Gray-Scott model",
# Mathematics and Computers in Simulation, Volume 40, 1996, pages 371-396.
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

class GrayScottDemo():
    def __init__( self ):
        self.gsize = 128
        self.dx = 2.5 / 256 # parameter from [Pearson1993]

        reactionstrs = [ "U + 2 V --> 3 V",
                         "      V -->    ",
                         "        --> U  ",
                         "      U -->    "  ]

        # conditions for pattern formation:
        # sigma should be between 2 and 6 [Mazin1996]
        # F should be < 0.25, and K should be < 1/16 = 0.0625, but
        # only a very narrow range of F & K combinations lead to
        # patterns, and this range becomes narrower with decreasing sigma,
        # (see [Mazin1996] p. 379)

        sigma = 4.0
        DU = 2e-5
        DV = DU / sigma
        #F = 0.01
        #K = 0.04
        F = 0.04
        K = 0.06

        # spiral waves (!!) for gsize=256 (can't see anything for gsize=32)
        #sigma = 2.0
        #DU = 2e-5 #0.1
        #DV = DU / sigma
        #F=0.01
        #K=0.04

        self.rsys = ReactionDiffusionSystem(self.gsize, self.gsize, self.dx)
        self.rsys.parse(reactionstrs)
 
        self.rsys.set_coefficient(1, F+K)  # V decays with rate F+K
        self.rsys.set_coefficient(2, F)    # U is injected with rate F
        self.rsys.set_coefficient(3, F)    # U also decays with rate F

        self.rsys.set_diffcoef('U', DU)
        self.rsys.set_diffcoef('V', DV)

        # color conventions as in [Pearson1993]:
	self.rsys.set_color( 'U', (0.5, 0, 0) ) # substrate U is dark red
	self.rsys.set_color( 'V', (0, 0, 1) )   # autocatalyst V is blue

        # initial conditions around U=0.5, V=0.5*alpha for delta ~= 0
        #alpha = F / (F+K)
        #self.rsys.resetAll('U', 0.5)
        #self.rsys.resetAll('V', 0.5 * alpha)
        #self.rsys.perturbAll('V', 0.333, 0.2)

        # spot that divides in 4 then makes a circle at the center
        # for F=0.04 and K=0.06, also for F=0.01 and K=0.04 (but no circle):
        self.rsys.resetAll('U', 0.5)
        initpos = self.gsize/2
        self.rsys.deposit('V', 0.5 * 2, initpos, initpos)

        #self.startcond0(F, K)

        print "dt should be smaller than", self.dx ** 2 / (2 * max(DU, DV))

    def startcond0( self, F, K ):
        # starting conditions from breve's Gray-Scott demo:
        self.rsys.resetAll("U", 1.0)
        for x in range(self.gsize):
            for y in range(self.gsize):
                rnd = np.random.random() - 0.5
                kf = 1.0 + K/F
                c = 0.5 + np.sqrt(np.abs(0.25 - F*kf*kf)) + 0.02*rnd
                self.rsys.set_conc('U', c, x, y)
                rnd = np.random.random() - 0.5
                c = (1.0 - c)/kf + 0.02*rnd
                self.rsys.set_conc('V', c, x, y)

    def run( self, finalvt=3000.0, dt=0.1 ):
        #self.rsys.trace()
        i = 0
        n = 10 # want n snapshots at most
        m = int(finalvt / (dt * n))
        pngsize = 512 # generate 512 x 512 PNG files
        self.rsys.writepng('gs0.png', pngsize, pngsize, transparent=False)
        while (self.rsys.vtime() <= finalvt):
            self.rsys.integrate(dt)
            i += 1
            if (i % m == 0):
                fname = ('gs%d.png' % (i / m) )
                self.rsys.writepng(fname, pngsize, pngsize, transparent=False)
            if (i % 50 == 0):
                self.rsys.animate(sleeptime=0.01, transparent=False)

if __name__ == '__main__':
    demo = GrayScottDemo()
    demo.run()

