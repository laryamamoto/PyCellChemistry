#---------------------------------------------------------------------------
#
# Fraglets.py: a python implementation of Fraglets, a chemically-inspired
# programming language for computer networks
#
# Reference:
#
# C. Tschudin. Fraglets: A metabolistic execution model for communication
# protocols. Proc. 2nd Annual Symposium on Autonomous Intelligent Networks
# and Systems (AINS), July 2003.
#
# python implementation by Lidia Yamamoto, Belgium, October 2013
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

from artchem.KeyMultiset import *

class Fraglets():
    def __init__( self, nid='' ):
        """ create a Fraglets interpreter and reaction vessel with node id
            'nid'
        """
        self.unimol = Multiset()
	self.active = KeyMultiset()
	self.passive = KeyMultiset()
        self.instr = { # instruction set
            'match'  : 'M',
            'matchp' : 'Z',
            'dup'    : 'D',
            'exch'   : 'E',
            'pop'    : 'P',
            'nop'    : 'N',
            'nul'    : 'U',
            'split'  : 'S',
            'send'   : 'X',
            'fork'   : 'F'
        }
        self.op = { # implementation of instruction set
            'M' : [ self.r_match,  'match' ],
            'Z' : [ self.r_matchp, 'matchp' ],
            'D' : [ self.r_dup,    'dup' ],
            'E' : [ self.r_exch,   'exch' ],
            'P' : [ self.r_pop,    'pop' ],
            'N' : [ self.r_nop,    'nop' ],
            'U' : [ self.r_nul,    'nul' ],
            'S' : [ self.r_split,  'split' ],
            'X' : [ self.r_send,   'send' ],
            'F' : [ self.r_fork,   'fork' ]
        }
        self.prop = {}
        self.wt = 0.0
        self.cnx = {} # list of remote connections, identified by a tag
        self.nodeid = nid # tag that identifies this node
        self.idle = True

    def ismatchp(self, mol):
        """ true if fraglet 'mol' starts with a matchp instruction """
        if mol == '': return False
        return mol[0] == self.instr['matchp']

    def isbimol(self, mol):
        """ true if fraglet 'mol' starts with a bimolecular reaction rule """
        if mol == '': return False
        return mol[0] == self.instr['match'] or mol[0] == self.instr['matchp']

    def isunimol(self, mol):
        """ true if fraglet 'mol' starts with a unimolecular reaction rule """
        if mol == '': return False
        return mol[0] in self.op and not self.isbimol(mol)

    def ispassive(self, mol):
        """ true if fraglet 'mol' does not start with a reaction rule """
        if mol == '': return False
        return mol[0] not in self.op

    def getmethod(self, mol):
        """ get method that implements the instruction at the head
            (first symbol) of fraglet 'mol'
        """
        k = mol[0]
        info = self.op[k]
        f = info[0]
        return f

    def getname(self, op):
        """ get human-friendly name for an instruction in character-encoded
            format
        """
        if op not in self.op: return op
        info = self.op[op]
        name = info[1]
        return name

    def set_nodeid(self, tag):
        """ assign nodeid 'tag' to this vessel """
        self.nodeid = tag

    def react1(self, mol):
        """ fire unimolecular reaction involving molecule 'mol' """
        f = self.getmethod(mol)
        result = f(mol)
        self.trace_reaction(mol, result)
        return result

    def react2(self, mol1, mol2):
        """ fire bimolecular reaction between mol1 and mol2 """
        f = self.getmethod(mol1)
        result = f(mol1, mol2)
        self.trace_reaction([mol1, mol2], result)
        return result

    def react(self, w):
        """ perform the selected reaction pointed to by the dice position w
            (typically involked from the hierarchical Gillespie SSA
            implementation in Cell.py)
        """
        if self.wt <= 0: return
        for k in self.active.keys():
            if (k in self.prop):
                if self.prop[k] > 0 and w < self.prop[k]:
                    mol1 = self.active.expelrnd(k)
                    mol2 = self.passive.expelrnd(k)
                    res = self.react2(mol1, mol2)
                    self.inject_list(res)
                    return
                w -= self.prop[k]

    def r_match(self, mol1, mol2):
        """ fire a 'match' fraglet: merge 2 fraglets """
        return mol1[2:] + mol2[1:]

    def r_matchp(self, mol1, mol2):
        """ fire a 'matchp' fraglet: merge 2 fraglets, while keeping a
            copy of the original matchp fraglet
        """
        return [ mol1, self.r_match(mol1, mol2) ]

    def r_dup(self, mol):
        """ fire a 'dup' fraglet: duplicate 3rd symbol """
        if len(mol) < 2: return ''
        if len(mol) < 3: return mol[1]
        return mol[1] + mol[2] + mol[2:]

    def r_exch(self, mol):
        """ fire an 'exch' fraglet: exchange symbols n. 3 and 4 """
        if len(mol) < 2: return ''
        if len(mol) < 4: return mol[1:]
        if len(mol) < 5: return mol[1] + mol[3] + mol[2]
        return mol[1] + mol[3] + mol[2] + mol[4:]

    def r_pop(self, mol):
        """ fire a 'pop' fraglet: consume 'nop' symbol plus 3rd symbol """
        if len(mol) < 2: return ''
        if len(mol) < 4: return mol[1]
        return mol[1] + mol[3:]

    def r_nop(self, mol):
        """ fire a 'nop' fraglet: consume 'nop' symbol """
        if len(mol) < 2:
            return ''
        return mol[1:]

    def r_nul(self, mol):
        """ fire a 'nul' fraglet: delete the fraglet """
        return ''

    def r_split(self, mol):
        """ fire a 'split' fraglet: split at the first occurrence of a '*'
            symbol
        """
        if len(mol) < 2: return ''
        return mol[1:].split('*', 1)

    def r_send(self, mol):
        """ fire a 'send' fraglet by consuming the header symbols and
            injecting the remaining tail in the destination vessel
        """
        # PENDING: send could be implemented as a reaction: if link not
        # connected, wait for connection (so send is a reaction between
        # the fraglet and the link to where the message is sent
        if len(mol) < 3: return ''
        dst = mol[1]
        if dst in self.cnx:
            addr = self.cnx[dst]
            addr.inject(mol[2:])
        return ''

    def r_fork(self, mol):
        """ fire a 'fork' fraglet that duplicates its tail """
        if len(mol) < 2: return ''
        if len(mol) < 3: return mol[1:]
        if len(mol) < 4:
            m1 = mol[1]
            m2 = mol[2]
        else:
            m1 = mol[1] + mol[3:]
            m2 = mol[2:]
        return [m1, m2]

    def add_cnx(self, dst, addr):
        """ add a link between this vessel and another vessel: 'dst' is the
            tag that identifies the link; 'addr' is the reference to the
            Fraglets object corresponding to the destination vessel
        """
        self.cnx[dst] = addr

    def del_cnx(self, dst):
        """ delete a link given its tag id 'dst' """
        del self.cnx[dst]

    def inert(self):
        """ true if reactor is inert: there are no more reactions to fire """
        return self.idle

    def inject( self, mol, mult=1 ):
        """ inject 'mult' copies of fraglet 'mol' in the reactor """
	if (mol == '' or mult < 1): return
        if (self.isbimol(mol)):
            if len(mol) > 1:
                key = mol[1]
                #if key in self.op: # invalid fraglet
                #    while len(mol) > 1 and self.isbimol(mol):
                #        # eliminate [match match match ...] sequences
                #        mol = mol[1:]
                #    #pending: clean fraglet
                self.active.inject(key, mol, mult)
                self.idle = False
            # else discard invalid fraglet
        elif (self.isunimol(mol)):
            if len(mol) > 1:
                self.unimol.inject(mol, mult)
                self.idle = False
            # else discard invalid fraglet
        else:
            key = mol[0]
            self.passive.inject(key, mol, mult)
            self.idle = False

    def inject_list( self, mlist ):
        """ inject the list of fraglets 'mlist' in the reactor """
        if type(mlist) is list:
            for m in mlist:
                self.inject(m)
        else:
            self.inject(mlist)

    def run_unimol(self):
        """ run all unimolecular transformations at once """
        n = 0
        while self.unimol.mult() > 0:
            mol = self.unimol.expelrnd()
            res = self.react1(mol)
            self.inject_list(res)
            n += 1
        return n
        
    def propensity(self):
        """ calculate all propensities of bimolecular reactions """
        self.run_unimol()
        self.prop = {}
        self.wt = 0.0
        for k in self.active.keys():
            m = self.active.multk(k)
            p = self.passive.multk(k)
            w = m * p
            if w > 0:
                self.prop[k] = w
            self.wt += w
        if self.wt <= 0: self.idle = True
        return self.wt

    def run_bimol(self):
        """ pick one bimolecular reaction using Gillespie's SSA;
            assumes propensities are up to date;
            for single reactor only:
            if more than one reactor, use react() within Cell instead
        """
        if self.wt <= 0: return
        w = np.random.random() * self.wt
        self.react(w)

    def iterate(self):
        """ one iteration of fraglets for single-vessel configuration """
        self.propensity()
        if not self.inert(): self.run_bimol()

    def run(self, niter):
        """ run for 'niter' iterations, or until the reactor is inert """
        for i in range(niter):
            print >> sys.stderr, "ITER=", i
            self.iterate()
            if self.inert(): return

    def mol2fraglet(self, mol, mult=1):
        """ convert a 'character-encoded' molecule into a human-readable
            fraglet
        """
        frag = self.nodeid + '['
        for i in range(len(mol)):
            op = mol[i]
            name = self.getname(op)
            frag += ' ' + name
        frag += ' ]'
        if mult > 1:
            frag += str(mult)
        return frag

    def parse(self, frag):
        """ parse fraglet, converting it to a condensed 'character-encoded'
            string with one character per symbol
        """
        # PENDING: deal with '[ ]mult' syntax for multiplicity
        frag = frag.strip('[')
        frag = frag.strip(']')
        slist = frag.split(' ')
        mol = ''
        for s in slist:
            if len(s) <= 0: continue
            if s in self.instr:
                mol += self.instr[s]
            else:
                # tag # TMP: take 1st symbol; TO DO: build symbol table
                mol += s[0]
        return mol

    def trace_mol(self, mol, mult=1):
        """ print a given fraglet, in human-readable format """
        frag = self.mol2fraglet(mol, mult)
        print >> sys.stderr, frag,

    def trace_mlist(self, mlist):
        """ print a list of fraglets """
        if type(mlist) is list:
            prev = False
            for mol in mlist:
                if prev: print >> sys.stderr, ', ',
                self.trace_mol(mol)
                prev = True
        else:
            self.trace_mol(mlist)

    def trace_reaction(self, mlist1, mlist2):
        """ print reaction in the form [frag1], [frag2] --> [frag3]
            (doesn't work for the send reaction)
        """
        self.trace_mlist(mlist1)
        print >> sys.stderr, ' --> ',
        self.trace_mlist(mlist2)
        print >> sys.stderr

    def trace_all_msets(self):
        """ print all fraglets in the vessel, in 'compiled' character-code
            format (for debugging purposes)
        """
        print >> sys.stderr, "UNIMOL=",
        self.unimol.trace()
        print  >> sys.stderr, "ACTIVE=",
        self.active.trace()
        print  >> sys.stderr, "PASSIVE=",
        self.passive.trace()

    def trace(self):
        """ print all fraglets in the vessel, in human-readable format """
        for mol in self.unimol.keys():
            mult = self.unimol.mult(mol)
            self.trace_mol(mol, mult)
            print >> sys.stderr
        for k in self.active.keys():
            mset = self.active.keymset[k]
            for mol in mset.keys():
                mult = mset.mult(mol)
                self.trace_mol(mol, mult)
                print >> sys.stderr
        for k in self.passive.keys():
            mset = self.passive.keymset[k]
            for mol in mset.keys():
                mult = mset.mult(mol)
                self.trace_mol(mol, mult)
                print >> sys.stderr

    def interpret(self, fname=''):
        """ interpret a file containing fraglet code for a single vessel
        """
        # PENDING: extend to multiple networked vessels in NetFraglets.py
        if fname == '':
            infile = sys.stdin
        else:
            try:
                infile = open(fname)
            except:
                print "file not found:", fname
                return
        for line in infile:
            line = line.strip(' \t\n')
            if line == '' or line[0] == '#':
                continue # skip comments and blank lines
            frag = self.parse(line)
            self.inject(frag)
        if fname != '':
            infile.close()

#---------------------------------------------------------------------------

# A few non-networked (single reactor) examples
# for networked (or multiple reactor) examples, see NetFraglets.py

def rndsoup():
    """ generate a random soup of fraglets and run it for a few iterations """
    tags = '12345'
    instrs = 'MDEPNU*' # random instructions except matchp & send
    alphabet = instrs + tags
    probm = 0.5 # probability of a matchp at the beginning
    maxlen = 10
    vessel = Fraglets()
    for i in range(10):
        # create passive fraglet [ tag ... ]
        tag = tags[np.random.randint(len(tags))]
        rndlen = np.random.randint(1, maxlen)
        mol = tag + rndpolymer(alphabet, rndlen)
        p = np.random.random()
        if p < probm: # create active fraglet [ matchp tag ...]
            mol = 'Z' + mol
        vessel.inject(mol)
    print  >> sys.stderr, "INIT:",
    vessel.trace()
    vessel.run(10)
    print  >> sys.stderr, "END:",
    vessel.trace()

def codegrowth():
    """ elongating fraglets """
    n = Fraglets()
    n.inject('SZaDa*aa')
    n.run(10)
    n.trace()

def quine():
    """ quine: a self-replicating fraglet """
    n = Fraglets()
    n.inject('FNbMbFNb')
    n.run(10)
    n.trace()

def test_interpreter(filename):
    """ read fraglets program from a file and run it for a few iterations """
    n = Fraglets()
    n.interpret(filename)
    n.run(100)
    n.trace()

if __name__ == '__main__':
    #rndsoup()
    #codegrowth()
    quine()
    #test_interpreter('test.fra')
