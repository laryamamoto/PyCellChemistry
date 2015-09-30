#---------------------------------------------------------------------------
#
# BinaryStrings.py: a module to manipulate binary strings as integers
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

import numpy as np
import sys
import time

def randprob():
    """ random floating point number between 0 and 1 """
    return np.random.random()

def dice(prob):
    """ throw a dice: true if the dice falls within the given probability """
    if prob == 1.0:
        return True
    if prob <= 0.0:
        return False
    d = randprob()
    return (d <= prob)

def randbin( nbits, p1=0.5 ):
    """ random integer with a given number of bits 'nbits'
        bits are uniformly distributed, with a possible bias:
        p1 = probability of a one
        p0 = 1 - p1 =  probability of a 0
    """
    if p1 == 0.5:
        # unbiased coin
        if (nbits <= 16):
            return np.random.randint(2**nbits)
        n = 0
        while (nbits > 0):
            b = min(nbits, 16)
            x = np.random.randint(2**b)
            n = (n << b) + x
            nbits -= b
    else:
        # biased coin: CAUTION! can be over 10 times slower than unbiased!
        n = 0
        for i in range(nbits):
            r = np.random.random()
            if (r < p1):
                b = 1
            else:
                b = 0
            n = (n << 1) + b
    return n

def getbitvalue(binstr, pos):
    """ get value of the bit at position 'pos' in binary string 'binstr' """
    return (binstr >> pos) & 1

def flipbit( binstr, pos, nbits ):
    """ flip a given bit in the string """
    mask = (1 << pos)
    return (binstr ^ mask)

def getsegment(binstr, p0, p1):
    """ get value of segment from p0 to p1 in binary string 'binstr' """
    seg = 0
    for p in range(p1 - 1, p0 - 1, -1):
        b = getbitvalue(binstr, p)
        seg = (seg << 1) | b
    return seg

def mutate1( binstr, nbits ):
    """ mutate a random bit in the string """
    pos = np.random.randint(nbits)
    return flipbit(binstr, pos, nbits)

def mutatep( binstr, nbits, prob ):
    """ mutate each bit in the string with a given probability
        caution: this method is computationally much slower than mutate1()
    """
    if prob <= 0 or nbits <= 0:
        return binstr
    mask = 1
    newb = binstr
    for i in range(nbits):
        if dice(prob):
            newb = (newb ^ mask)
        mask = (mask << 1)
    return newb

def crossover( s1, s2, nbits ):
    """ crossover between two binary strings s1 and s2 at a random position;
        returns two children c1 and c2
    """
    pos = np.random.randint(1, nbits) # cut in 2 segments at position pos
    m2 = 2**nbits - 1
    m1 = 2**pos - 1    # mask for lower segment [0;pos-1]
    m2 = (m2 & ~ m1)   # mask for upper segment [pos;len-1]
    #print "pos=", pos, "m1=", bin(m1), "m2=", bin(m2)
    c1 = ((s1 & m2) | (s2 & m1)) # child1 = upper(s1) + lower(s2) 
    c2 = ((s2 & m2) | (s1 & m1)) # child2 = upper(s2) + lower(s1)
    return (c1, c2)

def count_ones( binstr ):
    """ number of ones in a binary string """
    n = 0
    while (binstr != 0):
        binstr = (binstr & (binstr - 1))
        n += 1
    return n

def hamming( x, y ):
  """ Hamming distance between two integers x and y """
  n = x ^ y
  return count_ones(n)

def entropy( binstr, nbits ):
    """ entropy of a binary string
        CAUTION: for a string of nbits, there are only nbits/2 + 1
        different entropy levels 0 <= entropy <= 1 including 0.0
        (for p0 = 0 or p0 = 1) and 1.0 (for p0=p1=0.5)
    """
    n1 = count_ones(binstr)
    if n1 > nbits:
        raise ValueError('bad bit count for entropy calculation')
    p1 = 1.0 * n1 / nbits
    p0 = 1 - p1
    if p0 == 0 or p1 == 0:
        return 0.0
    else:
        return - p0 * np.log2(p0) - p1 * np.log2(p1)

def int2str( binstr, nbits ):
    """ converts an integer to a string containing the binary number
        in ascii form, without prefix '0b' and filled with leading
        zeros up to nbits
    """
    return str(bin(binstr))[2:].zfill(nbits)

