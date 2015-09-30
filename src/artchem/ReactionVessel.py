#---------------------------------------------------------------------------
#
# ReactionVessel.py: a reaction vessel for well-stirred explicit
# Artificial Chemistries (S, R, A) in which molecules and reactions
# are represented explicitly, and for which the set of possible
# molecular species S and the set of reactions R does not change.
#
# Two algorithms A are implemented as specializations of the basic
# ReactionVessel class:
# - WellStirredVessel: a simple ODE integration algorithm based on
#   Euler's method
# - GillespieVessel: a simple implementation of Gillespie's stochastic
#   simulation algorithm (SSA) using the direct method
#
# by Lidia Yamamoto, Univ. Basel, Switzerland, January 2010
# June 2013: adapted to the PyCellChemistry package
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

import sys
import numpy as np

from Multiset import *
from ReactionParser import *

class ReactionVessel( ):
    def __init__( self ):
        """ initialize empty reaction vessel """
	self.species = []
	self.reactions = ReactionQueue()
        self.time = 0.0
        self.ns = 0
        self.nr = 0

    def add( self, reaction ):
        """ add a new chemical reaction to the vessel's set of reactions;
            'reaction' must be an object of the Reaction class
        """
	idx = self.reactions.add(reaction)
        self.nr += 1
	for mol in reaction.get_educts().keys():
	    if (self.species.count(mol) <= 0):
                self.species.append(mol)
                self.ns += 1
	for mol in reaction.get_products().keys():
	    if (self.species.count(mol) <= 0):
                self.species.append(mol)
                self.ns += 1
        return idx

    def parse(self, reactlist ):
        """ invoke the ReactionParser to parse a list of reactions """
	parser = ReactionParser()
	for rs in reactlist:
	    reaction = parser.parse_line(rs)
            self.add(reaction)
        self.close() # invoke close method of the subclass

    def parse_input( self, infile ):
        """ invoke the ReactionParser to parse an already open input file """
        parser = ReactionParser()
        reactions = parser.parse_input(infile)
	for reaction in reactions:
            self.add(reaction)
        self.close() # invoke close method of the subclass
        
    def parse_stdin( self ):
        """ invoke the ReactionParser to parse the standard input (stdin) """
        return self.parse_input(sys.stdin)

    def parse_file( self, fname ):
        """ invoke the ReactionParser to parse the file named 'fname' """
        parser = ReactionParser()
        reactions = parser.parse_file(fname)
	for reaction in reactions:
            self.add(reaction)
        self.close() # invoke close method of the subclass

    def get_names( self ):
        """ list of all possible molecular species (set S of S,R,A) """
        return self.species

    def set_coefficient( self, i, k ):
        """ update kinetic coefficient for reaction """
        reaction = self.reactions.peek(i)
        if reaction != None:
            reaction.set_coefficient(k)

    def get_coefficient( self, i ):
        """ get kinetic coefficient for reaction """
        reaction = self.reactions.peek(i)
        if reaction == None: return 0.0
        return reaction.get_coefficient()

    def nspecies( self ):
        """ total number of molecular species: |S| in S,R,A """
        return self.ns

    def nreactions( self ):
        """ total number of reactions: |R| in S,R,A """
        return self.nr

    def vtime( self ):
        """ current simulation time """
        return self.time

    def trace( self ):
        """ print internal variables, for debug purposes """
	self.reactions.trace()
        print >> sys.stderr, "ns =", self.ns
        print >> sys.stderr, "nr =", self.nr
        print >> sys.stderr, "species = ", self.species

    def trace_title( self, prefix='' ):
        """ output tab-separated title line for plotting;
            an optional prefix can be given to print only molecules
            starting with that character
        """
        if (prefix == ''):
            tline = list(self.species)
        else:
            tline = []
            for mol in self.species:
                if (mol[0] == prefix):
                    tline.append(mol)
        tline.insert(0, "time")
        print "\t".join(tline)

class ODESystem( ReactionVessel ):
    def __init__( self ):
        """ ODE system variables: stoichiometric matrices and rate vector """
	ReactionVessel.__init__( self )
        self.me = None # stoichiometric matrix for educts
        self.mp = None # stoichiometric matrix for products
        self.ms = None # net stoichiometric matrix for each reaction
        self.kv = None # kinetic coefficients

    def close( self ):
        """ produce stoichiometric matrices and rate vector """
        if self.me != None: return # already closed
        ns = self.nspecies()
        nr = self.nreactions()
        if (ns <= 0 or nr <= 0):
            return

        self.me = np.matrix(np.zeros((ns, nr)), int)
        self.mp = np.matrix(np.zeros((ns, nr)), int)
        self.kv = np.array(np.zeros(nr))             

        i = 0
        for mol in self.species:
            j = 0
	    for reaction in self.reactions:
                self.me[i, j] = reaction.get_educts().mult(mol)
                self.mp[i, j] = reaction.get_products().mult(mol)
                self.kv[j] = reaction.get_coefficient()
                j += 1
            i += 1
        self.ms = self.mp - self.me  # net stoichiometric matrix

    def set_coefficient( self, i, k ):
        """ update kinetic coefficient for reaction """
        ReactionVessel.set_coefficient(self, i, k)
        self.kv[i] = k

    def get_coefficient( self, i ):
        """ get kinetic coefficient for reaction (overrides default
            ReactionVessel method in order to use vector kv directly;
            the result of both methods should be the same though)
        """
        return self.kv[i]

class WellStirredVessel( ODESystem ):
    def __init__( self ):
        """ initialize empty reaction vessel for ODE integration """
	ODESystem.__init__( self )
        self.conc = None
        self.rate = None
        self.dcdt = None
        self.dilution = 0.0

    def close( self ):
        """ after all reactions have been parsed, this method should
            be called in order to create all the data structures
            needed for numeric integration
        """
        if self.me != None: return # already closed
        ODESystem.close(self)

        ns = self.nspecies()
        nr = self.nreactions()
        self.conc = np.matrix(np.zeros((ns, 1)))  # concentration vector
        self.rate = np.matrix(np.zeros((nr, 1)))  # rate vector
        self.dcdt = np.matrix(np.zeros((ns, 1)))  # concentration increment

    def get_conc( self, mol ):
        """ get the concentration of molecule 'mol' """
        if (mol == '' or self.species.count(mol) < 1):
            return 0.0
        i = self.species.index(mol)
        return self.conc[i,0]

    def set_conc( self, mol, conc ):
        """ set the concentration of molecule 'mol' to 'conc' """
        if (mol == '' or self.species.count(mol) < 1):
            return
        i = self.species.index(mol)
        if (conc < 0):
            conc = 0.0
        self.conc[i,0] = conc

    def reset( self, mol='', conc=0.0 ):
        """ reset the concentration of molecule 'mol' to 'conc';
            if 'mol' is not specified (mol='') then set the
            concentration of all molecules to the value 'conc';
            hence a reset() without parameters will reset all the
            concentration values to zero
        """
        if (mol == '' ):
            self.conc.fill(conc)
            return
        self.set_conc(mol, conc)

    def deposit( self, mol, conc ):
        """ increase the concentration of molecule 'mol' by an amount
            'conc'
        """
        if (mol == '' or self.species.count(mol) < 1):
            return
        i = self.species.index(mol)
        self.conc[i,0] += conc
        if (self.conc[i,0] < 0):
            self.conc[i,0] = 0.0

    def get_dilution( self ):
        """ get vessel capacity (measured in concentration) """
        return self.dilution

    def set_dilution( self, cap ):
        """ set dilution to vessel capacity, where 'cap' is the
            maximum concentration that the vessel should contain
        """
        self.dilution = cap

    def apply_dilution( self ):
        """ apply dilution to vessel capacity: normalize the
            concentrations to the vessel capacity, such that the sum
            of concentrations of all species in the system does not
            exceed the vessel capacity
        """
        if (self.dilution > 0.0):
            totc = self.conc.sum()
            overflow = totc - self.dilution
            if (overflow > 0):
                self.conc *= self.dilution / totc

    def integrate( self, dt=1.0 ):
        """ ODE integration for one timestep of size dt """
        for j in range(self.nr):
            concprod = 1.0
            for i in range(self.ns):
                concprod *= self.conc[i,0] ** self.me[i,j]
            self.rate[j,0] = self.kv[j] * concprod

        self.conc += dt * (self.ms * self.rate)
        self.apply_dilution()
        for i in range(self.ns):
            if (self.conc[i,0] < 0):
                self.conc[i,0] = 0.0
        self.time += dt

    def trace( self ):
        """ print internal variables, for debug purposes """
	ReactionVessel.trace(self)
        print >> sys.stderr, "me=\n", self.me
        print >> sys.stderr, "mp=\n", self.mp
        print >> sys.stderr, "ms=\n", self.ms
        print >> sys.stderr, "kv=\n", self.kv

    def trace_conc( self ):
        """ trace concentrations to tab-separated line for plotting """
        tline = [ "%g" % self.time ]
        for i in range(self.ns):
            tline.append("%g" % self.conc[i,0])
        print "\t".join(tline)

class GillespieVessel( ReactionVessel, Multiset ):
    def __init__( self, nav=1000 ):
        """ initialize empty reaction vessel for stochastic simulation;
            'nav' is the initial N_A * V quantity (Avogadro constant
            times volume of the vessel); typically a smaller nav
            (hence smaller volume) leads to more prominent stochastic
            fluctuations
        """
	ReactionVessel.__init__( self )
	Multiset.__init__( self )
        self.na = 6.02214e23 # Avogadro constant (N_A)
        self.nav = nav       # NAV = N_A * V
        self.volume = 1.0 * self.nav / self.na # volume such that N_A * V = NAV
        self.cv = None # microscopic kinetic coefficients
        self.wv = None
        self.wt = 0.0

    def get_volume( self ):
        """ volume of the vessel """
        return self.volume

    def get_nav( self ):
        """ current N_A * V value """
        return self.nav

    def set_volume( self, newvol ):
        """ sets the volume of the vessel, recalculating NAV = N_A*V
            and the mesoscopic kinetic constants accordingly
        """
        self.volume = newvol
        self.nav = 1.0 * self.volume * self.na
        if self.cv != None:
            for j in range(self.nr):
                self.wolkenhauer(j)

    def set_nav( self, newnav ):
        """ sets the value of N_A*V = NAV:

            this is useful when we have the parameters from the ODE
            (initial concentrations, kinetic constants k) and want to
            perform a stochastic simulation based on the same
            parameters
    
            in order to adjust NAV we use the relation:

            [S] = ns / NAV => ns = [S]*NAV

            thus by setting NAV small, we make the reactor small for
            the same concentration, that is, the number of molecules
            in the reactor will be smaller, hence the stochastic
            effects can be more easily seen; increase NAV to simulate
            larger particle numbers in the reactor, for a given
            concentration.
        """
        self.nav = newnav
        self.volume = 1.0 * self.nav / self.na
        if self.cv != None:
            for j in range(self.nr):
                self.wolkenhauer(j)

    def set_coefficient( self, i, k ):
        """ update kinetic coefficient for reaction, together with its
            corresponding mesoscopic value
        """
        ReactionVessel.set_coefficient(self, i, k)
        if self.cv != None:
            self.wolkenhauer(i)

    def get_conc( self, mol ):
        """ calculate concentration of molecule 'mol': [S] = ns / NAV """
        s = 1.0 * self.mult(mol) / self.nav
        return s

    def deposit( self, mol, conc ):
        """ deposit a concentration 'conc' of molecule 'mol';
            first must convert 'conc' [S] to number of molecules (ns)
            to be injected: [S] = ns / NAV => ns = [S] * NAV
        """
        if (mol == '' or self.species.count(mol) < 1):
            return
        self.inject(mol, int(conc * self.nav))

    def nspecies( self ):
        """ total number of possible molecular species (static) """
        return ReactionVessel.nspecies( self )

    def nspecies_mset( self ):
        """ number of species currently present in the multiset (dynamic) """
        return Multiset.nspecies( self )

    def close( self ):
        """ create data structures necessary for stochastic simulation;
            convert the macroscopic kinetic constants (k) to
            microscopic equivalents (c) using Wolkenhauer's relation
        """
        self.cv = np.array(np.zeros(self.nr))
        self.wv = np.array(np.zeros(self.nr))
        for j in range(self.nr):
            self.wolkenhauer(j)

    def factorial( self, n, m=1 ):
        """ recursively computes the factorial of n, stopping at m:
            factorial(n, m) = n*(n-1)*...*m for 1 < m < n
        """
        if n < 2 or n < m: return 1
        else: return (n * self.factorial(n-1, m))

    def nchoosek( self, n, k ):
        """ computes the binomial coefficient or "n choose k" function:
            C(n,k) for two integers n & k, where n > k
            C(n,k) = n!/((k!)(n-k)!) = n(n-1)(n-2)...(n-k+1)/k!
            efficient when either k or (n-k) is small (n can be big)
        """
        if (k > n): return 0 # ignore invalid input
        k = min(k, n - k)
        return self.factorial(n, n-k+1) / self.factorial(k)

    def wolkenhauer( self, j ):
        """ calculate the microscopic stochastic rate constants (c)
            from the macroscopic ones (k), according to [Wolkenhauer2004]:

            c = k / (NAV)^(m - 1) * prod_i=1;L(li!)

            where:
            . m is the total number of colliding molecules (size of
              educt multiset)
            . L is the number of species in the educt multiset
            . li is the educt stoichiometric coefficient of substance i
            . NAV = N_A * V (CAUTION!!! must be refreshed if NAV changes!!)
        """

        if self.cv == None: return
        reaction = self.reactions.peek(j)
        if reaction == None: return
        k = reaction.get_coefficient()
        educts = reaction.get_educts()
        m = educts.nmolecules()
        prod = 1
        for mol in educts.keys():
            n = educts.mult(mol)
            prod *= self.factorial(n)
        self.cv[j] = (k / ((self.na * self.volume) ** (m - 1))) * prod

    def perform_reaction( self, rno ):
        """ perform a given reaction, extracting the needed educts
            from the multiset and injecting the resulting products
        """
        reaction = self.reactions.peek(rno)
        if reaction == None: return # this shouldn't happen though
        educts = reaction.get_educts()
        prods = reaction.get_products()
        for mol in educts.keys():
            n = educts.mult(mol)
            self.expel(mol, n)
        for mol in prods.keys():
            n = prods.mult(mol)
            self.inject(mol, n)

    def propensity( self ):
        """ calculate reaction propensities, store them in wv[j], total=wt """
        self.wt = 0.0
        for j in range(self.nr):
            reaction = self.reactions.peek(j)
            educts = reaction.get_educts()
            prod = 1
            for mol in educts.keys():
                n = educts.mult(mol) # number of molecules needed
                m = self.mult(mol)   # number of molecules available
                if (m < n):
                    prod = 0         # no reaction if not enough molecules
                    break
                prod *= self.nchoosek(m, n)
            self.wv[j] = self.cv[j] * prod
            #print "j=", j, "wvj=", wv[j]
            self.wt += self.wv[j]
        return self.wt

    def react( self, dice ):
        """ perform the reaction pointed by the given dice
            relies on pre-calculated propensity values stored in wv[j]
        """
        for j in range(self.nr):
            if (dice < self.wv[j]):
                #print "performing reaction j =", j
                self.perform_reaction(j)
                return
            dice -= self.wv[j]

    def iterate( self ):
        """ one iteration of the Gillespie stochastic simulation algorithm """
        wt = self.propensity()
        if (wt <= 0):
            return

        # choose a random reaction with prob. proportional to propensity
        dice = wt * np.random.random()
        self.react(dice)

        # update simulation time
        dt = - np.log(np.random.random()) / wt
        self.time += dt

    def integrate( self, dt=1.0 ):
        """ iterate Gillespie SSA for an interval of size at least dt """
        t0 = self.time
        while (self.time - t0 < dt):
            self.iterate()

    def trace_mult( self, prefix='' ):
        """ print molecule counts to tab-separated line for plotting.
            a one-letter prefix can be given to print only molecules
            starting with that letter
        """
        tline = [ "%g" % self.time ]
        for mol in self.species:
            if (prefix == '' or mol[0] == prefix):
                tline.append("%g" % self.mult(mol))
        print "\t".join(tline)

    def trace_conc( self, prefix='' ):
        """ print concentrations, like trace_mult but concentrations
            are calculated as [S] = ns / NAV
        """
        tline = [ "%g" % self.time ]
        for mol in self.species:
            if (prefix == '' or mol[0] == prefix):
                tline.append("%g" % self.get_conc(mol))
        print "\t".join(tline)

    def trace( self ):
        """ print internal variables, for debug purposes """
	ReactionVessel.trace( self )
        print >> sys.stderr, "MULTISET=",
        Multiset.trace( self )
        print >> sys.stderr, 'NAV=', self.nav, "VOL=", self.volume
        for j in range(self.nr):
            reaction = self.reactions.peek(j)
            print >> sys.stderr, "REACTION j=", j, "cv=", self.cv[j], ":"
            reaction.trace()

