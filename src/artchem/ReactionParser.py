#---------------------------------------------------------------------------
#
# ReactionParser.py: parser for chemical reactions in text format, such as:
#
#    A + 2 B --> 3 C , k=2.49
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
from Multiset import *
from Reaction import *

class ReactionParser():
    def parse_molecules( self, mollist ):
        """ parse educt or product multiset given as a list of strings
            containing the molecule name and an optional stoichiometric
            coefficient, e.g. [ 'a', '2 fox', '4 b1' ]
        """
	newmols = Multiset()
	for mol in mollist:
            mol = mol.strip()
	    if mol.find(' ') < 0:
		name = mol
		n = 1
	    else:
                (num, name) = mol.split(' ', 2)
                name = name.strip()
                n = int(num)
	    if n > 0 and name != '':
		newmols.inject(name, n)
	return newmols

    def parse_line( self, line ):
        """ parse string containing chemical reaction """
        line2 = line.split('#') # skip comments
        line = line2[0]
        try:
	    if line.find(',') < 0:
                reactstr = line
                k = 1.0
            else:
                (reactstr, kstr) = line.split(',', 2)
                reactstr = reactstr.strip()
                (kstr, kvar) = kstr.split('=', 2)
                if (kstr.strip() != 'k'):
                    raise ValueError
                k = float(kvar)

	    (educts, products) = reactstr.split('-->', 2)
            edlist = educts.split('+')
            prlist = products.split('+')
            edmset = ReactionParser.parse_molecules( self, edlist )
            prmset = ReactionParser.parse_molecules( self, prlist )
        except ValueError:
            print >> sys.stderr, "syntax error on line=", line
            exit(-1)
        if (edmset.empty() and prmset.empty()): return
	newreaction = Reaction(edmset, prmset, k)
        return newreaction

    def parse_input( self, infile ):
        """ parse input file containing multiple chemical reactions,
            one per line.
            'infile' is the file descriptor for the open input file
        """
	reactions = ReactionQueue()
        for line in infile.readlines():
            line = line.strip()
            if (line != '' and line[0] != '#'):
                reaction = self.parse_line(line)
                reactions.add(reaction)
        return reactions

    def parse_stdin( self ):
        """ parse standard input (stdin) """
        return self.parse_input(sys.stdin)

    def parse_file( self, fname ):
        """ open and parse input file fname """
        reactions = None
        try:
            infd = open(fname, 'r')
            reactions = self.parse_input(infd)
            infd.close()
        except IOError:
            print >> sys.stderr, "Error opening file", fname
            exit(-1)
        return reactions

