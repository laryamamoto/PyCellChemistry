#---------------------------------------------------------------------------
#
# HighOrderChem.py:
#
# a simplified high-order chemistry in which the rules are written as
# strings containing a python method call. these strings are inserted
# into a multiset of rules that can in principle be operated upon like
# other multisets of molecules, although at this stage this is not
# done yet.
#
# - rules must have the form function(args), where function is a python
#   function call and args are the parameters for the function
# - rules must return a list of products of the chemical reaction
#
# by Lidia Yamamoto, Belgium, May 2014
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

import re
from artchem.Multiset import *

class HighOrderChem:
    def __init__( self ):
        """ initialize empty data and rule multisets """
        self.mset = Multiset() # multiset of data molecules
        self.rset = Multiset() # multiset of rules

    def isrule( self, mol ):
        """ true if molecule 'mol' contains a reaction rule of the form
            function(args), where function is a python function call
            and args are the parameters for the function
        """
        return type(mol) is str and re.match('[\.\w]+\(.*\)', mol) != None

    def is_effective(self, educts, products):
        """ true if the reaction defined by educts --> products is
            effective, that is the product multiset is different from the
            educt multiset
        """
        if educts == [] and products == []: return False
        return sorted(educts) != sorted(products)

    def iterate( self ):
        """ run one iteration of high-order algorithm: pick a random
            reaction rule, fill its binding sites with as many molecules
            as needed, and fire it.
            returns a triple (r, e, p) where r is the rule that has
            been fired, 'e' is the educt list and 'p' is the product
            list
        """
        rule = self.rset.expelrnd()
        if rule == '':
            return ('', [], []) # no more rules to be processed
        match = re.match('([\.\w]+)\((.*)\)', rule)
        if (match == None): # malformed rule
            return (rule, [], [])
        funct = match.group(1)
        bsite = match.group(2)
        nsites = len(bsite.split(','))
        educts = []
        products = []
        if self.mset.mult() >= nsites:
            for j in range(nsites):
                mol = self.mset.expelrnd()
                educts.append(mol)
            edstr = str(educts).strip('[ ]')
            cmd = ("products = %s(%s)" % (funct, edstr))
            exec(cmd)
            for mol in products:
                self.mset.inject(mol)
        self.rset.inject(rule) # reinject rule for reuse
        return (rule, educts, products)

