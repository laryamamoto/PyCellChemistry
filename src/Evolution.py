#---------------------------------------------------------------------------
#
# Evolution.py: basics of evolutionary dynamics, Evolution chapter
#
# see Quasispecies.py and Tournament.py for application examples
#
# by Lidia Yamamoto, Belgium, July 2013
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
from artchem.Multiset import *
import artchem.BinaryStrings as bs

class Evolution:

    def __init__( self ):
        """ create a random initial population of molecules with 
            intentionally bad fitness
        """
        self.soup = Multiset()  # population container
        self.nbits    = 10      # molecule (binary string) length in bits
        self.popsize  = 100     # population size
        self.fitfunct = '1'     # fitness function (see below)
        self.target   = 1023    # numeric target for fitfunct = 'T' or 'N'

        if (self.fitfunct == 'N'): # [Nowak&Schuster1989]
            # target initially present, just see how it grows/survives
            self.soup.inject(self.target)

        while self.soup.mult() < self.popsize:
            mol = self.randmol()
            f = self.fitness(mol)
            if f < 0.4:
                # force a 'bad' random initial population
                self.soup.inject(mol)

    def randmol( self ):
        """ generate a random molecule in the form of an N-bit integer """
        return bs.randbin(self.nbits)

    def fitness( self, binstr ):
        """ calculate the fitness of an individual (normalized to one) """
        if self.fitfunct == 'E': # minimize the entropy of the string
            return 1.0 - bs.entropy(binstr, self.nbits)
        if self.fitfunct == '1': # maximize the number of bits set to one
            return 1.0 * bs.count_ones(binstr) / self.nbits
        if self.fitfunct == '0': # maximize the number of bits set to zero
            return 1.0 * (self.nbits - bs.count_ones(binstr)) / self.nbits
        if self.fitfunct == 'M': # maximize the numeric value of the string
            return 1.0 * binstr / (2**self.nbits)
        if self.fitfunct == 'T': # minimize the distance to a given target
            return 1.0 - 1.0 * abs(self.target - binstr) / (2**self.nbits)
        if self.fitfunct == 'N': # [Nowak&Schuster1989] simplest possible
            if (binstr == self.target):
                return 1.0 # target has maximum fitness
            else:
                return 0.2 # other sequence have equal lower fitness
        return 0.0

    def optimum( self ):
        """ produce an optimum individual for the desired fitness function """
        if self.fitfunct == 'E' or self.fitfunct == '0':
            return 0 # another solution for 'E': 2**self.nbits - 1
        if self.fitfunct == '1' or self.fitfunct == 'M':
            return 2**self.nbits - 1
        if self.fitfunct == 'T' or self.fitfunct == 'N':
            return self.target
        return None

    def avgfitness( self ):
        """ compute the average fitness of the population """
        avg = 0.0
        for mol in self.soup.keys():
            f = self.fitness(mol)
            m = self.soup.mult(mol)
            avg += f * m
        avg = avg / self.soup.mult()
        return avg

    def bestworstfit( self, mset ):
        """ find the best and worst individuals in a given multiset """
        fmax = 0.0
        fmin = 1.0
        best = ''
        worst = ''
        for mol in mset.keys():
            f = self.fitness(mol)
            if f > fmax or best == '':
                best = mol
                fmax = f
            if f < fmin or worst == '':
                worst = mol
                fmin = f
        return (best, worst)
