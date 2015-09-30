#---------------------------------------------------------------------------
#
# Tournament.py: GA with tournament selection in a chemistry, similar
# to the algorithm presented in the catalytic search papers below,
# adapted for chapter 7 of the book.
#
# References:
#
#   Lidia Yamamoto and Wolfgang Banzhaf, Catalytic search in dynamic 
#   environments. Artificial Life XII. MIT Press 2010, pages 277-285.
#
#   Lidia Yamamoto, Evaluation of a catalytic search algorithm,
#   Studies in Computational Intelligence volume 284, Springer 2010,
#   pages 75-87.
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

class Tournament(Evolution):

    def __init__( self ):
        """ start the tournament experiment with an intentionally bad
            initial population
        """
	Evolution.__init__( self )
        self.mutprob = 0.1      # mutation probability
        self.crossprob  = 0.9   # crossover probability

    def tournament1(self, toursize):
        """ tournament competition among individuals, variant 1:
            looser dies, winner reproduces asexually, the child gets mutated
            with a probability
        """
        competitors = Multiset()
        if (self.soup.mult() < toursize):
            return
        for i in range(toursize):
            competitors.inject(self.soup.expelrnd())
        (best, worst) = self.bestworstfit(competitors)
        competitors.expel(worst) # worst individual dies
        #if bs.dice(self.mutprob):
        #    child = bs.mutate1(best, self.nbits) # mutated child
        #else:
        #    child = best # cloned child
        child = bs.mutatep(best, self.nbits, self.mutprob)
        self.soup.inject(child)
        self.soup.absorb(competitors)
        s1 = bin(best)
        s2 = bin(worst)
        s3 = bin(child)
        #print "%s + %s -> %s + %s" % (s1, s2, s1, s3)

    def tournament2(self, toursize):
        """ tournament variant 2: two winners recombine with a probability,
            their children replace two loosers
        """
        competitors = Multiset()
        if (self.soup.mult() < toursize):
            return
        for i in range(toursize):
            competitors.inject(self.soup.expelrnd())
        (b1, w1) = self.bestworstfit(competitors)
        competitors.expel(b1)
        competitors.expel(w1)
        (b2, w2) = self.bestworstfit(competitors)
        competitors.expel(w2)
        if bs.dice(self.crossprob):
            (c1, c2) = bs.crossover(b1, b2, self.nbits)
        else:
            c1 = b1
            c2 = b2
        self.soup.inject(b1)
        self.soup.inject(c1)
        self.soup.inject(c2)
        self.soup.absorb(competitors)

    def tournament3(self, toursize):
        """ tournament variant 3: with crossover and mutation combined """
        competitors = Multiset()
        if (self.soup.mult() < toursize):
            return
        for i in range(toursize):
            competitors.inject(self.soup.expelrnd())
        (b1, w1) = self.bestworstfit(competitors)
        competitors.expel(b1)
        competitors.expel(w1)
        (b2, w2) = self.bestworstfit(competitors)
        competitors.expel(w2)
        c1 = b1
        c2 = b2
        if bs.dice(self.crossprob):
            (c1, c2) = bs.crossover(b1, b2, self.nbits)
        #if bs.dice(self.mutprob):
        #    c1 = bs.mutate1(c1, self.nbits)
        #if bs.dice(self.mutprob):
        #    c2 = bs.mutate1(c2, self.nbits)
        c1 = bs.mutatep(c1, self.nbits, self.mutprob)
        c2 = bs.mutatep(c2, self.nbits, self.mutprob)
        self.soup.inject(b1)
        self.soup.inject(c1)
        self.soup.inject(c2)
        self.soup.absorb(competitors)

    def run( self, niter=1000 ):
        """ run tournaments for 'niter' iterations """
        toursize = 4
        nlines = 500 # max number of lines to be plotted
        plotinterv = niter / nlines # plot interval
        print "time\tnspecies\tbest\tworst\ttotal"
        for i in range(niter):
            # several tournament variants for testing, pick one:
            #self.tournament1(toursize)
            #self.tournament2(toursize)
            self.tournament3(toursize)
            if (i % plotinterv == 0):
                ns = 1.0 * self.soup.nspecies() / self.popsize
                (best, worst) = self.bestworstfit(self.soup)
                fb = self.fitness(best)
                fw = self.fitness(worst)
                tot = self.soup.mult() / self.popsize
                print "%d\t%g\t%g\t%g\t%g" % (i, ns, fb, fw, tot)
        print >> sys.stderr, "best=", bin(best), "worst=", bin(worst)

if __name__ == '__main__':
    tour = Tournament()
    tour.run()
