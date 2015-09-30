#---------------------------------------------------------------------------
#
# NKlandscape.py:
#
# a simple implementation of Kauffman's NK model of rugged fitness landscapes
#
# by Lidia Yamamoto, Kraainem, Belgium, November 2013
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
import sys
import artchem.BinaryStrings as bs

class NKlandscape():
    def __init__( self, n, k ):
        """ create an NK landscape with given n and k, where:
            n is the total number of genes, and 
            k is the number of epistatic interactions per gene
        """
        self.N = n
	self.K = k
        self.origk = True # use original Kauffman's model
        if (k >= n):
            print >> sys.stderr, \
                "K too large, must be integer in the range 0 <= K < N"
            exit(-1)
        if self.origk:
            # use original model as proposed by Kauffman
            self.epi = self.mk_epistatic() # epistatic interactions
        else:
            # attempt to reduce memory and computation load,
            # but doesn't work well
            self.epi = self.mk_epistatic_ALT() # epistatic interactions
            self.fit = self.mk_fitness_ALT()   # fitness without epistasis

    def mk_epistatic( self ):
        """ build matrix of epistatic interactions:
            influence of K other genes on the fitness of each of the N genes
            epi[i,j] = fitness of gene i when the state of K+1 genes
            (gene i itself plus the K other genes that interact with it) is j
            CAUTION!! needs big matrix of size N * 2^(K+1), not scalable!!
        """
        epi = np.matrix(np.zeros((self.N, 2**(self.K + 1))))
        for i in range(self.N):
            for j in range(2**(self.K + 1)):
                #influence of K other genes on gene i
                epi[i,j] = bs.randprob()
                epi[i,j] = 1.0 * int(10.0 * epi[i,j]) / 10.0
        return epi

    def mk_epistatic_ALT( self ):
        """ alternative way to build a matrix of epistatic interactions:
            epi[i,j,v] = amount of influence of gene i on gene j when the
            value of gene i is v (this leads to a smaller matrix, but
            doesn't work very well, disabled by default)
        """
        epi = np.array(np.zeros((self.N, self.N, 2)))
        for i in range(self.N):
            for j in range(self.N):
                for v in [ 0, 1 ]:
                    #amount of influence of gene i on gene j
                    epi[i,j,v] = bs.randprob()
        return epi

    def mk_fitness_ALT( self ):
        """ make a fitness matrix in order to reduce computation load
            (goes with the mk_epistatic_ALT attempt, doesn't work very
            well, disabled by default)
        """
        fm = np.matrix(np.zeros((self.N, 2)))
        for i in range(self.N):
            for j in [ 0, 1 ]:
                fm[i,j] = bs.randprob()
        return fm
                
    def gene_fitness( self, pos, genome ):
        """ returns the fitness of the gene at position 'pos' in the given
            genome; CAUTION!! uses random epistatic matrix, not scalable!!!
        """
        genev = bs.getbitvalue(genome, pos)
        #print genev,
        infl = genev # genes that influence fitness of the gene at pos
        for ki in range(self.K):
            # assume genes pos+1 to pos+K+1 have an influence on gene at pos
            idx = (pos + ki + 1) % self.N
            #print "(", idx, ")",
            gk = bs.getbitvalue(genome, idx)
            infl = infl | (gk << (ki+1))
        wi = self.epi[pos, infl]
        #print "[", bin(infl), "w=", wi, "]",
        return wi

    def gene_fitness_ALT( self, pos, genome ):
        """ returns the fitness of the gene at position 'pos' in the given
            genome, using alternative pre-built fitness matrix and smaller
            matrix of epistatic interactions, different from the original
            method by Kauffman; (doesn't work very well, disabled by default)
        """
        genev = bs.getbitvalue(genome, pos)
        wgene = self.fit[pos, genev]
        for ki in range(self.K):
            # assume genes pos+1 to pos+K+1 have an influence on gene at pos
            idx = (pos + ki + 1) % self.N
            #print "(", idx, ")",
            gk = bs.getbitvalue(genome, idx)
            if gk: # positive interaction
                wgene += self.epi[idx, pos, gk] / self.N
            else:  # negative interaction
                wgene -= self.epi[idx, pos, gk] / self.N
        #wgene = (int(wgene * 10000.0) % 10000) / 10000.0
        return wgene

    def fitness( self, genome ):
        """ compute the fitness value for the given genome: the average of
            the individual fitness values of each gene in the genome
        """
        wavg = 0.0
        for i in range(self.N):
            # fitness of gene i given current genome state
            if self.origk:
                wi = self.gene_fitness(i, genome)
            else:
                wi = self.gene_fitness_ALT(i, genome)                
            wavg += wi
        wavg = wavg / self.N # fitness of genome = avg fit of all genes
        return wavg

    def dist2neighbors(self, genome):
        """ fitness distance between a given genome and all its 1-mutant
            neighbors
        """
        dist = np.array(np.zeros(self.N))
        w0 = self.fitness(genome)
        for pos in range(self.N):
            neigh = bs.flipbit(genome, pos, self.N)
            w = self.fitness(neigh)
            dist[pos] = abs(w - w0)
        return dist

def sweepNKdistance(n):
    """ run the NK model for given N and varying K """
    print "K\tAVGDIST"
    for k in range(n):
        #print '==================== N=', n, 'K=', k, ' ===================='
        nk = NKlandscape(n, k)
        avg = 0.0
        #print "epi=", nk.epi
        for i in range(10):
            genome = bs.randbin(n)
            #print "genome=", genome, 
            d = nk.dist2neighbors(genome)
            fit = nk.fitness(genome)
            #print "dist=", d, "avg=", np.average(d), "fit=", fit
            avg += np.sum(d)
        avg = avg / (n * 10.0)
        #print "K=", k, "AVGDIST=", avg
        print ("%d\t%g" % (k, avg))

if __name__ == '__main__':
    args = sys.argv[1:]
    n = 4
    if (len(args) > 0):
        try:
            n = int(args[0])
        except:
            print >> sys.stderr, "usage: python NKlandscape.py [ N ]"
            exit(-1)

    sweepNKdistance(n)
