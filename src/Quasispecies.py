#---------------------------------------------------------------------------
#
# Quasispecies.py: demonstration of quasispecies evolutionary dynamics
# using binary strings, book chapter 7
#
# usage: python quasispecies.py [ <mp> <histfile> ]
# where:
#   mp is the proportion of mutations with respect to the error threshold
#   histfile is the name of an output file that will contain the histogram
#     of species at the end of the simulation, grouped by their distance to
#     the optimum (number of bits that differ from the optimum sequence
#     (here '1111111111')
#
# by Lidia Yamamoto, Kraainem, Belgium, July 2013
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

from Evolution import *

class Quasispecies(Evolution):

    def __init__( self, mp ):
        """ create a random initial population of molecules with 
            intentionally bad fitness; set the mutation probability
            per bit as a proportion 'mp' of the error threshold
        """
	Evolution.__init__( self )
        self.errth = 1.0 / self.nbits # error threshold
        self.qsmutp = mp * self.errth # mutation prob. per bit
        self.killparent = False # kill parent
        self.hillclimb  = False # use hill climb

        print >> sys.stderr, "mp=", mp, "qsmutp=", self.qsmutp

    def propensity( self, mol ):
        """ propensity of a given reaction """
        m = self.soup.mult(mol)
        f = self.fitness(mol)
        p = 1.0 * m * f
        return p

    def set_propensities( self ):
        """ initial calculation of all propensities """
        self.prop = {} # propensities
        self.sump = 0.0 # sum of propensities
        for mol in self.soup.keys():
            p = self.propensity(mol)
            if p > 0.0:
                self.prop[mol] = p
                self.sump += p

    def update_propensity( self, mol, dm ):
        """ update propensity of molecule, given that dm units of it were
            added (dm > 0) or removed (dm < 0)
        """
        if (mol == '' or dm == 0):
            return
        dp = 1.0 * self.fitness(mol) * dm
        if dp == 0.0:
            return
        if (mol in self.prop):
            self.prop[mol] += dp
            if self.prop[mol] <= 0.0:
                del self.prop[mol]
        else:
            self.prop[mol] = dp
        self.sump += dp

    def react( self ):
        """ one iteration of Gillespie's SSA on the quasi-species equation """
        if self.sump <= 0.0:
            return 0.0
        dice = bs.randprob() * self.sump

        for mol in self.prop.keys(): # pick reaction among existing candidates
            if dice <= self.prop[mol]:
                # mutate each bit with probability qsmutp
                mut = bs.mutatep(mol, self.nbits, self.qsmutp)
                if self.hillclimb and self.fitness(mut) < self.fitness(mol):
                    # hill climb: use only fitter mutants
                    mut = mol
                if (self.killparent): # parent dies at childbirth
                    die = mol
                    self.soup.expel(mol)
                else: # a random molecule dies to make room for new one
                    die = self.soup.expelrnd()
                self.soup.inject(mut)
                dt = - np.log(bs.randprob()) / self.sump # time increment

                if (mut != die):
                    self.update_propensity(mut, +1)
                    self.update_propensity(die, -1)
                return dt
            dice -= self.prop[mol]
        print >> sys.stderr, "ERROR! dice =", dice, "prop=", self.prop
        exit(-1)
        #return 0.0

    def dist_histogram( self ):
        """ histogram of number of species at hamming distance d=[0;nbits]
            from master sequence
        """
        opt = self.optimum()
        neigh = {}
        for d in range(self.nbits + 1):
            neigh[d] = 0
        for mol in self.soup.keys():
            m = self.soup.mult(mol)
            d = bs.hamming(mol, opt)
            neigh[d] += m
        return neigh

    def plot_dist_histogram( self, fname ):
        """ plot distance histogram to file 'fname' """
        hdist = self.dist_histogram()
        try:
            outfd = open(fname, 'w')
        except IOError:
            print >> sys.stderr, "Error creating output file", fname
            exit(-1)
        print >> outfd, "distance\tcount"
        for d in range(self.nbits + 1):
            n = hdist[d]
            print >> outfd, ("%d\t%d" % (d, n))
        outfd.close()

    def plot( self, vtime ):
        """ plot fitness information """
        ns = 1.0 * self.soup.nspecies() / self.popsize
        (best, worst) = self.bestworstfit(self.soup)
        fb = self.fitness(best)
        fw = self.fitness(worst)
        avg = self.avgfitness()
        #tot = self.soup.mult() / self.popsize
        #print "%g\t%g\t%g\t%g\t%g\t%g" % (vtime, ns, fb, fw, avg, tot)
        print "%g\t%g\t%g\t%g\t%g" % (vtime, ns, fb, fw, avg)

    def run( self, histfile ):
        """ stochastic implementation of quasi-species equation using
            Gillespie's SSA
        """
        self.set_propensities()
        finalvt = 100.0 # finish time
        nlines = 500    # max number of lines to be plotted
        plotinterv = 1.0 * finalvt / nlines # plot interval
        vtime = 0.0
        lastvt = 0.0
        #print "time\tnspecies\tbest\tworst\tavg\ttotal"
        print "time\tspecies\tbest\tworst\taverage"
        self.plot(vtime)
        while (vtime < finalvt):
            dt = self.react()
            vtime += dt
            if (vtime <= finalvt and vtime - lastvt >= plotinterv):
                self.plot(vtime)
                lastvt = vtime

        (best, worst) = self.bestworstfit(self.soup)
        print >> sys.stderr, "best=", best, "=", bin(best), hex(best), \
            "worst=", worst, "=", bin(worst)
        print >> sys.stderr, "freq="
        self.soup.trace_frequencies()
        if histfile != '':
            self.plot_dist_histogram(histfile)

if __name__ == '__main__':
    args = sys.argv[1:]
    mp = 1.0
    histfile = ''
    if (len(args) > 0):
        mp = float(args[0])
    if (len(args) > 1):
        histfile = args[1]

    qs = Quasispecies(mp)
    qs.run(histfile)

