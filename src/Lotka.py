#---------------------------------------------------------------------------
#
# Lotka.py: Lotka-Volterra demo, ODE vs. stochastic
#
# usage: python Lotka.py [ stoch ]
#
#   if the string 'stoch' is passed as argument, a stochastic simulation
#   is performed, else a deterministic ODE integration is performed
#
# by Lidia Yamamoto, Univ. Basel, Switzerland, January 2010
#
# 20130621: ported from an initial version running on top of breve
# 20150621: easier switch between stochastic and deterministic modes
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

class LotkaVolterra:

    def __init__( self, usestoch ):
        """ Lotka-Volterra predator-prey system, with option for stochastic
            (usestoch=True) or deterministic simulation (usestoch=False)
        """
        self.stoch = usestoch
        reactionstrs = [
            "rabbit + grass --> 2 rabbit + grass , k=1",
            "fox + rabbit   --> 2 fox            , k=1",
            "fox            -->                  , k=1" ]

        if self.stoch:
            self.reactor = GillespieVessel(nav=40)
        else:
            self.reactor = WellStirredVessel()

        self.reactor.parse(reactionstrs)
        self.reactor.deposit('grass', 1.0)
        self.reactor.deposit('rabbit', 5.0)
        self.reactor.deposit('fox', 2.0)

    def extinct( self ):
        """ true if either predator or prey speces have gone extinct """
        if not self.stoch: return False
        nr = self.reactor.mult('rabbit') 
        nf = self.reactor.mult('fox') 
        return ( nr == 0 or nf == 0)

    def exploded( self ):
        """ true if either predator or prey exploded to infinity """
        cr = self.reactor.get_conc('rabbit') 
        cf = self.reactor.get_conc('fox') 
        return ( cr > 1e6 or cf > 1e6 )

    def run( self, finalvt=40.0, dt=0.001 ):
        """ run Lotka-Volterra system until simulation time 'finalvt' is
            reached; a plot line is generated every integration interval dt
        """
        self.reactor.trace_title()
        while (not self.extinct() and not self.exploded() and \
               self.reactor.vtime() <= finalvt):
            self.reactor.trace_conc()
            self.reactor.integrate(dt)

if __name__ == '__main__':
    args = sys.argv[1:]
    usestoch = False
    if (len(args) > 0 and args[0] == 'stoch'):
        usestoch = True

    lotka = LotkaVolterra(usestoch)
    lotka.run()
