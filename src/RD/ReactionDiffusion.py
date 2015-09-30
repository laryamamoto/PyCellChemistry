#---------------------------------------------------------------------------
#
# ReactionDiffusion.py: implementation of reaction-diffusion chemical systems
#
# originally based on the breve Hypercycle.[tz/py] demo by jon klein
# <jk@spiderland.org>, www.spiderland.org
#
# by Lidia Yamamoto, Univ. Basel, Switzerland, January 2010
# 20150910: removed breve dependencies to run within PyCellChemistry
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
sys.path.insert(0, '..')
from artchem.ReactionVessel import *
import WritePNG as png

class ReactionDiffusionSystem( ODESystem ):
    """ a PDE integrator for a reaction-diffusion system """
    def __init__( self, x, y, dx=1.0, periodic=True, nneighs=4 ):
	ODESystem.__init__( self )
        self.sizex = x # grid size (x, y)
        self.sizey = y
        self.dx = dx   # square grid cells of size (dx, dx)
        self.periodic = periodic # boundary conditions
        self.nneighs = self.set_neighbors(nneighs) # neighborhood config
	self.conc = None   # concentrations of all species on grid
	self.dcdt = None   # spatial dcdt for PDE integration
	self.prod = None
	self.rate = None
	self.diff = None
        self.diffcoef = None  # diffusion coefficients
        self.color = {}       # color of a molecule

    def close( self ):
        """ closes the PDE system such that it can be integrated """
        if self.me != None: return # already closed
        ODESystem.close(self)
        ns = self.nspecies()
        x = self.sizex
        y = self.sizey
	self.conc = np.zeros((ns, x, y))
	self.dcdt = np.zeros((x, y))
	self.prod = np.zeros((x, y)) 
	self.rate = np.zeros((x, y))
	self.diff = np.zeros((x, y))
        self.diffcoef = np.zeros(ns)

    def set_diffcoef( self, mol, value ):
        """ set the diffusion coefficient of molecule mol """
        if (mol == '' or self.species.count(mol) < 1):
            return
        i = self.species.index(mol)
        self.diffcoef[i] = value

    def get_diffcoef( self, mol ):
        """ returns the diffusion coefficient of molecule mol """
        if (mol == '' or self.species.count(mol) < 1):
            return 0.0
        i = self.species.index(mol)
        return self.diffcoef[i]

    def set_color( self, mol, color):
        """ set color of molecule to color=(red, green, blue) """
        self.color[mol] = color

    def get_color( self, mol ):
        """ returns the color that has been assigned to molecule mol """
        if mol == '' or mol not in self.species or mol not in self.color:
            return (0, 0, 0)
        return self.color[mol]

    def set_periodic( self, flag ):
        """ periodic (toroidal) vs. non-periodic boundary conditions """
        self.periodic = flag

    def get_periodic( self ):
        """ True if the current setup is periodic boundary """
        return self.periodic

    def set_neighbors( self, nn ):
        """ set neighborhood configuration: nn = 
            2: simple (north and right neighbors)
            4: Von Neumann
            6: hexagonal lattice
            8: Moore
        """
        if nn not in [ 2, 4, 6, 8 ]:
            print >> sys.stderr, "ns =", self.ns
            exit -1
        self.nneighs = nn

    def get_neighbors( self ):
        """ returns the current neighborhood configuration """
        return self.nneighs

    def get_pos( self, x, y ):
        """ make sure coordinates always fall within boundaries """
        if (self.periodic):
            x = (x + self.sizex) % self.sizex
            y = (y + self.sizey) % self.sizey
        else:
            x = min(max(x, 0), self.sizex)
            y = min(max(y, 0), self.sizey)
        return (x, y)

    def get_conc_by_index( self, idx, x, y ):
        """ get concentration of molecule by index, respecting boundaries """
        if (self.periodic):
            (x, y) = self.get_pos(x, y)
            return self.conc[idx, x, y]
        elif (x >= 0 and x < self.sizex and y >= 0 and y < self.sizey):
            return self.conc[idx, x, y]
        else:
            return 0.0

    def get_conc( self, mol, x, y ):
        """ get concentration of molecule by name, respecting boundaries """
        if (mol == '' or self.species.count(mol) < 1):
            return 0.0
        i = self.species.index(mol)
        return self.get_conc_by_index(i, x, y)

    def set_conc( self, mol, conc, x, y ):
        """ set the concentration of a molecule at position x,y to a
            given value, respecting boundary conditions
        """
        if (mol == '' or self.species.count(mol) < 1):
            return
        i = self.species.index(mol)
        if (conc < 0):
            conc = 0.0
        if (self.periodic):
            (x, y) = self.get_pos(x, y)
            self.conc[i, x, y] = conc
        elif (x >= 0 and x < self.sizex and y >= 0 and y < self.sizey):
            self.conc[i, x, y] = conc

    def deposit( self, mol, conc, x, y ):
        """ deposit a given amount of a molecule at a precise location """
        c = self.get_conc(mol, x, y)
        self.set_conc(mol, c + conc, x, y)

    def rnd_deposit( self, npoints, mol, conc, ampl=0.0):
        """ deposit a random amount of a molecule at random locations """
        if (mol == '' or self.species.count(mol) < 1):
            return
        c = conc
        for i in range(npoints):
            x = np.random.randint(self.sizex)
            y = np.random.randint(self.sizey)
            if (ampl > 0.0):
                c = conc + ampl * np.random.random() - ampl / 2
            self.deposit(mol, c, x, y)

    def reset( self, mol, x, y ):
        """ reset the concentration of a given molecule to zero """
        self.set_conc(mol, 0.0, x, y)

    def resetAll( self, mol='', conc=0.0 ):
        """ reset the concentration of a given molecule to a given
            value, overall on the grid; if no molecule is specified,
            reset the concentrations of all molecules
        """
        if (conc < 0):
            conc = 0.0
        if (mol == ''):
            self.conc.fill(conc)
            return
        if (self.species.count(mol) < 1):
            return
        i = self.species.index(mol)
        self.conc[i].fill(conc)

    def set_patch_at( self, mol, conc, x, y ):
        """ create a patch of chemical at a given location """
        self.set_conc( mol, conc, x, y )
        self.set_conc( mol, conc, x, ( y + 1 ) )
        self.set_conc( mol, conc, x, ( y - 1 ) )
        self.set_conc( mol, conc, ( x + 1 ), y )
        self.set_conc( mol, conc, ( x - 1 ), y )
        self.set_conc( mol, conc, (x - 1), ( y - 1 ) )
        self.set_conc( mol, conc, (x - 1), ( y + 1 ) )
        self.set_conc( mol, conc, ( x + 1 ), (y - 1) )
        self.set_conc( mol, conc, ( x + 1 ), (y + 1) )

    def set_patch( self, mol, conc ):
        """ create a patch of chemical at a random location """
        if (self.sizex < 3 or self.sizey < 3):
            return
        x = 1 + np.random.randint( self.sizex - 2 )
        y = 1 + np.random.randint( self.sizey - 2 )
        self.set_patch_at( mol, conc, x, y )

    def set_patches( self, npatches, mol, initconc ):
        """ create some initial random patches of chemicals """
        m = mol
        for i in range(npatches):
            if (mol == ''):
                c = np.random.randint( self.ns )
                m = self.species[c]
            self.set_patch(m, initconc)


    def add_patch_at( self, mol, conc, x, y ):
        """ add some concentration to a patch of chemical at a given
            location
        """
        self.deposit( mol, conc, x, y )
        self.deposit( mol, conc, x, ( y + 1 ) )
        self.deposit( mol, conc, x, ( y - 1 ) )
        self.deposit( mol, conc, ( x + 1 ), y )
        self.deposit( mol, conc, ( x - 1 ), y )
        self.deposit( mol, conc, (x - 1), ( y - 1 ) )
        self.deposit( mol, conc, (x - 1), ( y + 1 ) )
        self.deposit( mol, conc, ( x + 1 ), (y - 1) )
        self.deposit( mol, conc, ( x + 1 ), (y + 1) )

    def add_patch( self, mol, conc ):
        """ add a patch of chemical at a random location """
        if (self.sizex < 3 or self.sizey < 3):
            return
        x = 1 + np.random.randint( self.sizex - 2 )
        y = 1 + np.random.randint( self.sizey - 2 )
        self.add_patch_at( mol, conc, x, y )

    def add_patches( self, npatches, mol, initconc ):
        """ add some random patches of chemicals """
        m = mol
        for i in range(npatches):
            if (mol == ''):
                c = np.random.randint( self.ns )
                m = self.species[c]
            self.add_patch(m, initconc)

    def reset_region( self, mol, val, x0, y0, x1, y1 ):
        """ set the concentration of substances in this region to the
            given value
        """
        y = y0
        while (y < y1):
            x = x0
            while (x < x1):
                self.set_conc(mol, val, x, y)
                x += 1
            y += 1
        
    def disrupt( self, mol, val, x0, y0, x1, y1 ):
        """ disrupt the concentration of a chemical within a given
            rectangle: clamp the concentr to at most val
        """
        y = y0
        while (y < y1):
            x = x0
            while (x < x1):
                c = self.get_conc(mol, x, y)
                c = min(c, val)
                self.set_conc(mol, c, x, y)
                x += 1
            y += 1

    def perturb( self, mol, prob, ampl, x, y ):
        """ perturb concentration at the given point """
        c0 = self.get_conc(mol, x, y)
        if (c0 > 0):
            dice = np.random.random()
            if dice < prob:
                c = c0 * ( ampl * np.random.random() - ampl / 2)
                self.deposit(mol, c, x, y)

    def perturbAll( self, mol, prob, ampl ):
        """ perturb the given chemical with a given probability and
            amplitude: a random concentration value in the range
            [-amp/2 ; ampl/2] will be deposited at each node with
            probability prob (if the chemical is not in the node,
            nothing happens)
        """
        for y in range(self.sizey):
            for x in range(self.sizex):
                self.perturb(mol, prob, ampl, x, y)

    def perturb_region( self, mol, prob, ampl, x0, y0, x1, y1 ):
        """ perturb concentrations in the given region """
        y = y0
        while (y < y1):
            x = x0
            while (x < x1):
                self.perturb(mol, prob, ampl, x, y)
                x += 1
            y += 1
        
    def diffusion_term_NOTUSED( self, n ):
        """ calculate diffusion term for molecule with index n, leave
            result on self.diff (calculate one by one, far too slow,
            but I leave the code here to illustrate conceptually how
            it works, it's more readable than the matrix form below)
        """
        self.diff.fill(0.0)
        for y in range(self.sizey):
            for x in range(self.sizex):
                c0 = self.get_conc_by_index(n, x, y);
                c1 = self.get_conc_by_index(n, x - 1, y);
                c2 = self.get_conc_by_index(n, x + 1, y);
                c3 = self.get_conc_by_index(n, x, y - 1);
                c4 = self.get_conc_by_index(n, x, y + 1);
                dc = c1 + c2 + c3 + c4 - 4 * c0;
                if self.nneighs == 6:
                    sign = 2 * (x % 2) - 1;
                    c1 = self.get_conc_by_index(n, x + sign, y - 1);
                    c2 = self.get_conc_by_index(n, x + sign, y + 1);
                    dc += c1 + c2 - 2 * c0;
                if self.nneighs == 8:
                    c1 = self.get_conc_by_index(n, x - 1, y - 1);
                    c2 = self.get_conc_by_index(n, x - 1, y + 1);
                    c3 = self.get_conc_by_index(n, x + 1, y - 1);
                    c4 = self.get_conc_by_index(n, x + 1, y + 1);
                    dc += c1 + c2 + c3 + c4 - 4 * c0;
                self.diff[x, y] = dc / (self.dx ** 2)

    def diffusion_term( self, n ):
        """ calculate diffusion term for molecule with index n, for
            whole grid at once using matrix operations; leave result
            on self.diff
        """
        c0 = self.conc[n]
        if self.periodic:
            # rotate whole grid left, right, up and down
            c1 = np.roll(c0, 1, axis=0)
            c2 = np.roll(c0, -1, axis=0)
            c3 = np.roll(c0, 1, axis=1)
            c4 = np.roll(c0, -1, axis=1)
            self.diff = ( c1 + c2 + c3 + c4 - 4 * c0 ) / (self.dx ** 2)
        else:
            # shift whole grid, duplicating edges as if chemicals bounce back
            y = self.sizey
            c1 = np.vstack([c0[1:,:], c0[-1]])
            c2 = np.vstack([c0[0], c0[:-1,:]])
            c3 = np.hstack((c0[:,0:1], c0[:,:-1]))
            c4 = np.hstack((c0[:,1:], c0[:,y-1:y]))
            self.diff = ( c1 + c2 + c3 + c4 - 4 * c0 ) / (self.dx ** 2)
        # PENDING
        #if self.nneighs == 6:
        #if self.nneighs == 8:
        #self.diff = dc

    def reaction_term( self, n):
        """ calculate reaction term for molecule with index n, leaving
            the result on self.dcdt
        """
        self.dcdt.fill(0.0)
        for j in range(self.nreactions()):
            k = self.ms[n,j] * self.kv[j]
            if (k == 0.0): continue
            self.rate.fill(k)
            for i in range(self.nspecies()):
                # calculate [conc]^me
                self.prod = self.conc[i] ** self.me[i,j]
                # calculate r = k * [c1]^e1 * [c2]^p2 ...
                self.rate = np.multiply(self.rate, self.prod)
            self.dcdt += self.rate

    def integrate( self, dt=1.0 ):
        """ PDE integration of reaction-diffusion system with timestep dt """

        for n in range(self.nspecies()):
            self.reaction_term(n)
            if (self.diffcoef[n] != 0.0):
                self.diffusion_term(n)
                # dc/dt = sum(r[j]) + diffcoef*diffterm
                self.dcdt += self.diffcoef[n] * self.diff
            # conc += dc/dt * dt
            self.conc[n] += dt * self.dcdt
        #self.apply_dilution()
        #minc = self.conc.min()
        #if minc < 0.0:
        #    self.conc += minc0
        self.time += dt

    def trace( self ):
        """ print internal variables, for debug purposes """
	ReactionVessel.trace( self )
        print "diffcoef=", self.diffcoef, "\n"

    def trace_title( self ):
        """ output tab-separated title line for plotting """
        print "time",
        for i in range(self.ns):
            print"%s" % self.species[i],
        print''

    def trace_conc( self ):
        """ print the full grid of concentrations for each chemical,
            producing a tab-separated matrix for each chemical
        """
        for mol in self.species:
            i = self.species.index(mol)
            tot = self.conc[i].sum()
            avg = tot / (self.sizex * self.sizey)
            minv = self.conc[i].min()
            print "t=", self.time, " mol=", mol,
            print "sum=", tot, "avg=", avg, "min=", minv
            for y in range(self.sizey):
                for x in range(self.sizex):
                    print "\t", self.get_conc(mol, x, y),
                print ''

    def trace_conc_xy( self, x, y ):
        """ print the concentations of all chemicals at position (x,y) """
        print "%g" % self.time,
        for i in range(self.ns):
            print"%g" % self.conc[i].getValue(x, y),
        print ''

    def computeTextureMatrix( self, transparent ):
        """ produce a texture matrix to be used by animate() """
        if not hasattr(self, 'texmat'):
            # 4D texture matrix with one (r,g,b,a) color per grid position
            self.texmat = np.empty((self.sizex, self.sizey, 4), np.ubyte)
        maxc = self.conc.max()
        if maxc > 1.0:
            norm = self.conc / maxc
        else:
            norm = self.conc
        self.texmat.fill(0)
        if not transparent:
            self.texmat[:,:,3] = 0xFF
        for mol in self.species:
            n = self.species.index(mol)
            (r0, g0, b0) = self.get_color(mol)
            self.texmat[:,:,0] += r0 * 0xFF * norm[n]
            self.texmat[:,:,1] += g0 * 0xFF * norm[n]
            self.texmat[:,:,2] += b0 * 0xFF * norm[n]
            if transparent: self.texmat[:,:,3] += 0xFF * norm[n]

    def animate( self, sleeptime=1.0, transparent=True, blur=False ):
        """ display the reaction-diffusion grid using vpytyon
            transparent: use transparent background instead of black
            blur: produce fuzzy borders between cells on the grid
        """
        import visual as vpy
        self.computeTextureMatrix(transparent)
        self.texture = vpy.materials.texture(data=self.texmat, \
                       mapping='sign', interpolate=blur)
        self.plate = vpy.box(axis=(0,0,1), width=2, height=2, \
                     length=0.1, material=self.texture)
        # TODO: how to reuse existing box and texture???
        print "animate t=", self.time
        vpy.sleep(sleeptime)

    def conc2img(self, cellsizex, cellsizey, transparent):
        """ convert the concentration matrix to the ARGB format
            accepted by png.saveAsPNG()

            if the transparent flag is set, the resulting png image
            will have a transparent background in the spots where
            there are no chemicals;
            else the image will have a black background in those spots

            TODO: implement hexagonal grid
        """
        maxc = self.conc.max()
        img = []
        for y in range(self.sizey):
            row = []
            for x in range(self.sizex):
                if transparent: A = 0
                else: A = 255
                (R, G, B) = (0, 0, 0)
                for mol in self.species:
                    n = self.species.index(mol)
                    (r0, g0, b0) = self.get_color(mol)
                    if maxc > 1.0:
                        norm = self.conc[n,x,y] / maxc
                    else:
                        norm = self.conc[n,x,y]
                    if transparent: A += int(round(norm * 0xFF))
                    R += int(round(norm * r0 * 0xFF))
                    G += int(round(norm * g0 * 0xFF))
                    B += int(round(norm * b0 * 0xFF))
                if transparent: A = min(max(0, A), 255)
                R = min(max(0, R), 255)
                G = min(max(0, G), 255)
                B = min(max(0, B), 255)
                val = (A << 24) | (R << 16) | (G << 8) | B
                for j in range(cellsizex): row.append(val)
            for i in range(cellsizey): img.append(row)
        return img

    def conc2img_DISABLED(self, cellsizex, cellsizey, transparent):
        """ convert the concentration matrix to the ARGB format
            accepted by png.saveAsPNG()

            tried to use matrix ops for speedup, but it turned out
            slower than the element-by-element implementation above
        """
        if not hasattr(self, 'A'):
            self.A = np.empty((self.sizex, self.sizey), np.uint32)
            self.R = np.empty((self.sizex, self.sizey), np.uint32)
            self.G = np.empty((self.sizex, self.sizey), np.uint32)
            self.B = np.empty((self.sizex, self.sizey), np.uint32)
        if transparent:
            self.A.fill(0x00)
        else:
            self.A.fill(0xFF)
        self.R.fill(0x00)
        self.G.fill(0x00)
        self.B.fill(0x00)
        maxc = self.conc.max()
        if maxc > 1.0:
            norm = self.conc / maxc
        else:
            norm = self.conc
        for mol in self.species:
            n = self.species.index(mol)
            (r0, g0, b0) = self.get_color(mol)
            if transparent: self.A += 0xFF * norm[n]
            self.R += r0 * 0xFF * norm[n]
            self.G += g0 * 0xFF * norm[n]
            self.B += b0 * 0xFF * norm[n]
        np.clip(self.A, 0, 0xFF, out=self.A)
        np.clip(self.R, 0, 0xFF, out=self.R)
        np.clip(self.G, 0, 0xFF, out=self.G)
        np.clip(self.B, 0, 0xFF, out=self.B)
        tot = (self.A << 24) | (self.R << 16) | (self.G << 8) | self.B
        img = np.repeat(tot, cellsizex, axis=0)
        img = np.repeat(img, cellsizey, axis=1)
        return img

    def writepng( self, fname, pngsizex=256, pngsizey=256, transparent=True ):
        """ write molecule concentrations to a png image; each grid
            position will have a size of cellsize per cellsize pixels
        """
        cellsizex = max(1, pngsizex / self.sizex)
        cellsizey = max(1, pngsizey / self.sizey)
        img = self.conc2img(cellsizex, cellsizey, transparent)
        data = png.saveAsPNG(img, fname)

#---------------------------------------------------------------------------

if __name__ == '__main__':
    # simple diffusion test
    rsys = ReactionDiffusionSystem(32, 32)
    rsys.parse([ "A --> A", "B --> B" ])
    rsys.set_diffcoef('A', 0.8)
    rsys.set_diffcoef('B', 0.5)
    rsys.set_color('A', (1, 0, 0))
    rsys.set_color('B', (0, 0, 1))
    rsys.deposit('A', 2.0, rsys.sizex / 2, rsys.sizey / 2)
    #rsys.deposit('A', 2.0, 2, 2)
    #rsys.deposit('B', 2.0, 8, 8)
    #rsys.trace_conc()
    rsys.writepng('testAB0.png')
    for i in range(100):
        rsys.integrate(0.1)
        #rsys.trace_conc()
        if (i > 0 and i % 10 == 0):
            #fname = ('testAB%d.png' % i )
            #rsys.writepng(fname)
            rsys.animate()
