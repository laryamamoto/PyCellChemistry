#---------------------------------------------------------------------------
#
#  MolecularTSP.py: a reimplementation of the Molecular TSP algorithm by
#  [Banzhaf1990] W. Banzhaf, The "molecular" traveling salesman,
#  Bio. Cybern. 64, 7-14 (1990)
#
#  Usage: python MolecularTSP.py [ ncities ]  # default 10 cities
#
#  by Lidia Yamamoto, Belgium, September 2014
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

import numpy as np
from HighOrderChem import *

class Point:
    def __init__( self, x, y):
        """ a 2D point with coordinates (x,y) """
        self.x = x
        self.y = y

class TSPgraph:
    def __init__( self, nc, ring=True ):
        """ create a random map containing 'nc' cities;
            if the 'ring' flag is set, puts all the cities around a ring;
            else scatters them randomly on a 2D surface
        """
        self.maxncities = nc  # maximum (target) number of cities
        gs = 2 * nc # cities will be placed on a grid of size gs X gs
        self.gridsize = gs
        self.mindist = gs * 1.0 / nc # minimum distance between two cities
        self.ncities = 0   # current number of cities created so far
        self.coord = []    # list of city coordinates
        self.roads = []    # list of all roads leaving a given city

        if ring:
            self.createRingTopo()
            self.buildFullMesh()
        else:
            self.createRndTopo()
            #self.buildRoads()
            self.buildFullMesh()

    def distance( self, p1, p2 ):
        """ Euclidean distance between two points """
        return np.sqrt((p1.x - p2.x)**2 + (p1.y - p2.y)**2)

    def cityDistance( self, city1, city2 ):
        """ Euclidean distance between two cities """
        return self.distance(self.coord[city1], self.coord[city2])

    def occupied( self, x, y ):
        """ true if point (x,y) is too close to another city """
        c0 = Point(x, y)
        for i in range(self.ncities):
            c = self.coord[i]
            if self.distance(c, c0) < self.mindist:
                return True
        return False

    def createRingTopo( self ):
        """ ring topology: cities are placed on a ring with radius
            proportional to the total number of cities
        """
        radius = self.gridsize / 2
        angle = 2 * np.pi / self.maxncities
        theta = 0.0
        for i in range(self.maxncities):
            x = np.cos(theta) * radius + radius
            y = np.sin(theta) * radius + radius
            self.coord.append(Point(x, y))
            theta += angle
            self.ncities += 1

    def createRndTopo( self ):
        """ random topology: cities placed at random coordinates """
        while (self.ncities < self.maxncities):
            x = np.random.randint(0, self.gridsize)
            y = np.random.randint(0, self.gridsize)
            if not self.occupied(x, y):
                self.coord.append(Point(x, y))
                self.ncities += 1

    def buildFullMesh( self ):
        """ fully-meshed graph with roads from all cities to all others """
        for i in range(self.ncities):
            k = 0
            self.roads.append([])
            for j in range(self.ncities):
                if j == i:
                    continue
                self.roads[i].append(j)

    #def buildRoads( self ):
    #    """ random graph of roads (pending, must ensure graph is connected) """

    def randomTour( self, ns ):
        """ walk ns steps randomly on the graph, without going back to the
            same city; returns the tour generated
        """
        c0 = np.random.randint(0, self.ncities)
        tour = [ c0 ]
        c1 = c0
        while (len(tour) < ns):
            roads = list(self.roads[c1])
            while c1 in tour and len(roads) > 0:
                i = np.random.randint(0, len(roads))
                c1 = roads[i]
                roads.remove(c1)
            if c1 in tour:
                return tour
            tour.append(c1)
        return tour

    def trace( self ):
        """ print city number and coordinates for plotting """
        for i in range(self.ncities):
            c = self.coord[i]
            print ("%d\t%g\t%g\t#" % (i, c.x, c.y) ), self.roads[i]

    def traceTour( self, tour, fname ):
        """ prints tour for plotting """
        # PENDING: check if there is a road between two cities
        if type(tour) is not list or tour == []:
            return
        try:
            outfd = open(fname, 'w')
        except IOError:
            print >> sys.stderr, "Error creating output file", fname
            exit(-1)
        print >> outfd, "x\ty"
        for city in tour:
            c = self.coord[city]
            print >> outfd, ( "%g\t%g" % (c.x, c.y) )
        c = self.coord[tour[0]] # closes the tour
        print >> outfd, ( "%g\t%g" % (c.x, c.y) )
        outfd.close()

class MolecularTSP( HighOrderChem ):
    def __init__( self, nc, ring=True ):
        """ initialize molecular TSP with a graph of cities and a reactor
            with random candidate tours and rules that operate on tours;
            the graph contains 'nc' cities forming a ring (when ring=True),
            or scattered randomly (when ring=False) 
        """
	HighOrderChem.__init__( self )
        self.ncities = nc
        self.popsize = 9
        self.mscale = 100 # machine scale = 1/tR on p. 11 of [Banzhaf1990]
        self.maxgen = 1000 # max number of generations

        self.tsp = TSPgraph(self.ncities, ring)
        gs = self.tsp.gridsize
        self.toofar = self.tsp.distance(Point(0,0), Point(gs, gs))
        self.worstfitness = self.ncities * self.toofar
        self.bestfitness = np.pi * gs

        self.tsp.trace()
        print "best=", self.bestfitness, "worst=", self.worstfitness

        # start with a population of random molecules
        for i in range(self.popsize):
            mol = self.randomMolecule()
            self.mset.inject(mol)
        #self.mset.trace()

        # inject a number of rules proportional to their execution probability
        rule = 'self.exchangeMachine(%s)'
        self.rset.inject(rule, self.mscale)
        rule = 'self.cutMachine(%s)'
        self.rset.inject(rule, self.mscale)
        rule = 'self.invertMachine(%s)'
        self.rset.inject(rule, self.mscale)
        rule = 'self.recombinationMachine(%s,%s)'
        self.rset.inject(rule, 1)
        #self.rset.trace()

    def fitness( self, tour ):
        """ calculate the fitness of a tour """
        visited = []
        fit = 0
        for i in range(len(tour)):
            j = (i + 1) % len(tour)
            c1 = tour[i]
            c2 = tour[j]
            if c2 in self.tsp.roads[c1]:
                fit += self.tsp.cityDistance(c1, c2)
            else: # penalize for inexistent road
                fit += 2 * self.toofar
            if c1 in visited: # penalize for city visited more than once
                fit += self.toofar
            visited.append(c1)
        return fit

    def fitter( self, fit1, fit2 ):
        """ true if the fitness value fit1 is better than value fit2 """
        return fit1 < fit2

    def newMolecule( self, fit, tour ):
        """ constructs a molecule with syntax:
            fitness [ city1, city2, ... ]
        """
        mol = ( "%g %s" % (fit, tour))
        return mol

    def parseMolecule( self, mol ):
        """ parses the molecule, returning its fitness and tour components """
        (sfit, stour) = mol.split('[')
        fit = float(sfit)
        ltour = stour.strip('] ').split(',')
        tour = []
        for s in ltour:
            tour.append(int(s))
        return (fit, tour)

    def bestMolecule( self ):
        """ returns the best molecule in the data multiset """
        bfit = self.worstfitness
        bmol = ''
        for mol in self.mset.keys():
            (fit, tour) = self.parseMolecule(mol)
            if self.fitter(fit, bfit):
                bfit = fit
                bmol = mol
        return (bfit, bmol)

    def randomMolecule( self ):
        """ create a random candidate molecule (fit, c1, c2, ...) """
        tour = self.tsp.randomTour(self.ncities)
        fit = self.fitness(tour)
        mol = self.newMolecule(fit, tour)
        return mol

    def randomPair( self, tour ):
        """ returns a random pair of positions within a tour """
        p1 = np.random.randint(len(tour))
        p2 = p1
        while (p1 == p2):
            p2 = np.random.randint(len(tour))
        return (p1, p2)

    def splitTour( self, tour, p1, p2 ):
        """splits the tour in two segments:
            1. the segment starting at position p1, until position p2-1
            2. the segment starting at position p2, until position p1-1
            in a circular way, for example:

            1 2 3 4 5 6 7 ==> [2 3], [4 5 6 7 1]
              ^   ^
              p1  p2
            1 2 3 4 5 6 7 ==> [4 5 6 7 1], [2, 3]
              ^   ^
              p2  p1

            if the positions p1 and p2 coincide, the first segment is
            equal to the tour, and the second segment is empty
        """

        if p1 % len(tour) == p2 % len(tour):
            return (tour, [])
        i = p1
        seg1 = []
        while (i != p2):
            seg1.append(tour[i])
            i = (i+1) % len(tour)
        i = p2
        seg2 = []
        while (i != p1):
            seg2.append(tour[i])
            i = (i+1) % len(tour)
        return (seg1, seg2)

    def exchangeOperator( self, tour ):
        """ E-Machine operator:
            chooses two random cities in a tour and inverts their order
            e.g. 1 2 3 4 5 6 7 ==> 1 2 3 6 5 4 7
                       ^   ^             ^   ^
        """
        (p1, p2) = self.randomPair(tour)
        newtour = list(tour)
        newtour[p1] = tour[p2]
        newtour[p2] = tour[p1]
        # pending: check validity of tour; repeat op if invalid
        return newtour

    def exchangeMachine( self, mol ):
        """E-Machine implementation:
            invokes the E-Machine operator, evaluates the fitness of
            the new molecule, compares it with the original one, and
            returns the fittest molecule
        """
        print "E-",
        (fit, tour) = self.parseMolecule(mol)
        newtour = self.exchangeOperator(tour)
        newfit = self.fitness(newtour)
        if self.fitter(newfit, fit):
            return [ self.newMolecule(newfit, newtour) ]
        else:
            return [ mol ]

    def cutOrInvertOperator( self, tour, invert=False ):
        """ implements the C-Machine operator when invert=False
            implements the I-Machine operator when invert=True
        """
        newtour = list(tour)
        if len(tour) < 3: # too small to do anything
            return newtour
        while (newtour == tour):
            (p1, p2) = self.randomPair(tour)
            (seg1, seg2) = self.splitTour(tour, p1, p2)
            if len(seg2) > 1:
                p3 = np.random.randint(len(seg2) - 1) + 1
                (seg3, seg4) = self.splitTour(seg2, 0, p3)
            else:
                seg3 = seg2
                seg4 = []
            if invert:
                seg1.reverse()
            newtour = seg3 + seg1 + seg4
        return newtour

    def cutMachine( self, mol ):
        """ C-Machine: transposes the segment lying between two cities
            behind a third city
            e.g. 1 2 3 4 5 6 7 ==> 1 4 5 2 3 6 7
                   ^ ^    |         x   |^ ^
        """
        print "C-",
        (fit, tour) = self.parseMolecule(mol)
        newtour = self.cutOrInvertOperator(tour, False)
        newfit = self.fitness(newtour)
        if self.fitter(newfit, fit):
            return [ self.newMolecule(newfit, newtour) ]
        else:
            return [ mol ]

    def invertMachine( self, mol ):
        """ I-Machine: same as the C-Machine but the segment cut is
            reversed before insertion
            e.g. 1 2 3 4 5 6 7 ==> 1 4 5 3 2 6 7
                   ^ ^    |         x   |^ ^
        """
        print "I-",
        (fit, tour) = self.parseMolecule(mol)
        newtour = self.cutOrInvertOperator(tour, True)
        newfit = self.fitness(newtour)
        if self.fitter(newfit, fit):
            return [ self.newMolecule(newfit, newtour) ]
        else:
            return [ mol ]

    def recombinationOperator( self, tour1, tour2 ):
        """ R-Machine operator: recombines two tours
            example:
            [1 2 3 4 5 6 7], [3 6 1 4 5 2 7] ==>
             ^                    |-----|
            [1 4 5 2 2 3 4 5 6 7] ==> [1 4 5 2 3 6 7]
             |-----|
        """

        # extract a random segment from tour1
        (p1, p2) = self.randomPair(tour1)
        (seg1, notused) = self.splitTour(tour1, p1, p2)
        city = seg1.pop(0) # remove overlapping city

        newtour = list(tour2)
        for c in seg1:  # remove duplicated cities
            newtour.remove(c)

        p3 = newtour.index(city) # graft point is overlapping city
        (seg2, seg3) = self.splitTour(newtour, 0, p3+1)
        newtour = seg2 + seg1 + seg3 # insert segment at graft point
        return newtour

    def recombinationMachine( self, mol1, mol2 ):
        """ R-Machine: recombines 2 tours and selects the 2 fitter tours
            out of the multiset of 4 tours formed by the 2 parents plus the
            2 generated offsprings
        """

        print "R-",
        (fit1, tour1) = self.parseMolecule(mol1)
        (fit2, tour2) = self.parseMolecule(mol2)
        tour3 = self.recombinationOperator(tour1, tour2)
        fit3 = self.fitness(tour3)
        mol3 = self.newMolecule(fit3, tour3)
        if self.fitter(fit1, fit3) and self.fitter(fit2, fit3):
            return [ mol1, mol2 ]
        elif self.fitter(fit1, fit2) and self.fitter(fit3, fit2):
            return [ mol1, mol3 ]
        else:
            return [ mol2, mol3 ]

    def traceAll( self, gen ):
        """ print all molecules in the data multiset """
        i = 0
        for mol in self.mset.keys():
            m = self.mset.mult(mol)
            (fit, tour) = self.parseMolecule(mol)
            fname = ( "gen%d-mol%d-fit%g-x%d.gr" % (gen, i, fit, m))
            self.tsp.traceTour(tour, fname)
            i += 1

    def run( self ):
        """ run the Molecular TSP algorithm for up to a maximum number of
            generations, or until a sufficiently optimal solution is
            found;
            produces a set of output files, one for each individual
            in the initial population, and one for each individual in
            the final population
        """

        # number of machine operations in one generation:
        # genops = c = ceil(M / p_eff), on p. 9 of [Banzhaf1990]
        #   where M = self.popsize = number of data strings
        #         p_eff = sum(t_j * m_j) = |rset| / mscale
        n = self.popsize * self.mscale * 1.0 / self.rset.mult()
        genops = int(np.ceil(n))

        print "popsize=", self.popsize, "scale=", self.mscale, \
              "mult=", self.rset.mult(), "genops=", genops
        bmol = ''

        self.traceAll(0)
        gen = 0
        (bfit, bmol) = self.bestMolecule()
        while gen <= self.maxgen and self.fitter(self.bestfitness, bfit):
            print "gen=", gen, "pop=", self.mset.mult(), "bmol=", bmol,
            if gen < self.maxgen:
                for j in range(genops):
                    self.iterate()
            print
            (bfit, bmol) = self.bestMolecule()
            gen += 1
        self.traceAll(gen - 1)

if __name__ == '__main__':
    args = sys.argv[1:]
    n = 10 # default number of cities
    if (len(args) > 0):
        try:
            n = int(args[0])
        except:
            print >> sys.stderr, "usage: python MolecularTSP.py [ ncities ]"
            exit(-1)

    moltsp = MolecularTSP(n, ring=True)
    moltsp.run()

