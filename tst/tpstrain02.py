########################################################################
# MechSys - Open Library for Mechanical Systems                        #
# Copyright (C) 2005 Dorival M. Pedroso, Raul D. D. Farfan             #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# any later version.                                                   #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY without even the implied warranty of        #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http:#www.gnu.org/licenses/>   #
########################################################################

import mechsys as m

# Constants
L  = 2.0   # Length
H  = 2.0   # Height
E  = 207.0 # Young
nu = 0.3   # Poisson
q  = 1.0   # Load
nx = 20    # Divisions along x (must be even)
ny = 20    # Divisions along y

# ----------------------------------------------------------------------------- Mesh

# Blocks
blocks = [m.mesh_block()]
blocks[0].set_coords (-1,                           # tag to be replicated to all elements
                      [(0,0), (L,0), (L,H), (0,H)], # vertices' coordinates
                      [(0,1), (1,2), (2,3), (3,0)], # edges
                      {(0,1):-10, (2,3):-20},       # edge tags
                      {},                           # face tags
                      [1 for i in range(nx)],       # weights x
                      [1 for i in range(ny)],       # weights y
                      [],                           # weights z
                      0, 1, 3, 0)                   # Origin, XPlus, YPlus, None

# Generate
print 'Mesh Generation: --------------------------------------------------------------'
ms = m.mesh_structured()
ne = ms.generate (blocks)
print ne, ' elements generated '
print

# ------------------------------------------------------------------------------ FEM

# Geometry
g = m.geom(2)

# Nodes brys
nbrys = [[L/2., 0.0, 0.0, 'ux', 0.0]] # x,y,z, key, val

# Edges brys
ebrys = [[-10, 'uy', 0.0], # [tag], [key], [val]
         [-20, 'fy',  -q]] # [tag], [key], [val]

# Elements attributes
eatts = [[-1, 'Quad4PStrain', 'LinElastic', 'E=%f nu=%f'%(E,nu), 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0']] # [tag], [type], [model], [prms], [inis]

# Set geometry: nodes, elements and boundaries
m.set_nodes_elems (ms, eatts, g)
m.set_brys        (ms, nbrys, ebrys, [], g) # [] => no face brys

# Solve
print 'Solution: ---------------------------------------------------------------------'
sol = m.solver('ForwardEuler')
sol.set_geom(g)
sol.solve()

# Output file
m.out_vtu(g, 'tpstrain02_py.vtu')
print 'File <tpstrain02_py.vtu> generated'

#----------------------------------------------------------------------------- Check

# Check
errors = 0.0

Sy = q
Ex = -nu*(1.0+nu)*Sy/E
Ey =  (1.0-nu*nu)*Sy/E
Sz = (E/(1.0+nu))*(nu/(1.0-2.0*nu))*(Ex+Ey)

# Stress and strains
for i in range(g.nelems()):
    errors += abs(g.ele(i).val("Ex" ) - (Ex))
    errors += abs(g.ele(i).val("Ey" ) - (Ey))
    errors += abs(g.ele(i).val("Exy") - (0.0))
    errors += abs(g.ele(i).val("Sx" ) - (0.0))
    errors += abs(g.ele(i).val("Sy" ) - (Sy ))
    errors += abs(g.ele(i).val("Sz" ) - (Sz ))
    errors += abs(g.ele(i).val("Sxy") - (0.0))

# Displacements
for i in range(g.nnodes()):
    # bottom nodes
    if (abs(g.nod(i).y())<1.0e-5):
        errors += abs(g.nod(i).val("uy") - (0.0))

    # fixed bottom (central) node
    if (abs(g.nod(i).y())<1.0e-5) and (abs(g.nod(i).x()-(L/2.0))<1.0e-5):
        errors += abs(g.nod(i).val("ux") - (0.0))

    # half-height nodes
    if (abs(g.nod(i).y()-(H/2.0))<1.0e-5):
        errors += abs(g.nod(i).val("uy") - (-0.5*H*Ey))

    # top nodes
    if (abs(g.nod(i).y()-(H))<1.0e-5):
        errors += abs(g.nod(i).val("uy") - (-H*Ey))

    # left nodes
    if (abs(g.nod(i).x())<1.0e-5):
        errors += abs(g.nod(i).val("ux") - (0.5*L*Ex))

    # right nodes
    if (abs(g.nod(i).x()-L)<1.0e-5):
        errors += abs(g.nod(i).val("ux") - (-0.5*L*Ex))

print 'Py:Errors = ', errors
