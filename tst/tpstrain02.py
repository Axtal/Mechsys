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
L     = 2.0   # Length
H     = 2.0   # Height
E     = 207.0 # Young
nu    = 0.3   # Poisson
q     = 1.0   # Load
ndivx = 2     # Divisions along x (must be even)
ndivy = 2     # Divisions along y

# ----------------------------------------------------------------------------- Mesh

# Blocks
blocks = [m.mesh_block()]
blocks[0].set ([[0.,  L, L, 0.,    L/2.,    L, L/2.,   0.] , # coordinates: x values
                [0., 0., H,  H,      0., H/2.,    H, H/2.]], # coordinates: y values
               [1 for i in range(ndivx)],                    # weights x
               [1 for i in range(ndivy)])                    # weights y

# Tags
blocks[0].set_etags ([0, 0, -10, -20])

# Generate
print 'Mesh Generation: --------------------------------------------------------------'
ms = m.mesh_struct()
ne = ms.generate (blocks)
print ne, ' elements generated '
print

# ------------------------------------------------------------------------------ FEM

# Problem dimension
m.dim (2)

# Nodes and elements
m.add_nodes_elems (ms, 'Quad4PStrain')

# Boundary conditions
m.set_node_brys (ms, [L/2.], [0.0], [0.0], ['ux'], [0.0])
m.set_face_brys (ms, [-10, -20], ['uy', 'fy'], [0.0, -1])

# Parameters and initial values
for i in range(m.n_elems()):
    m.elems(i).set_model('LinElastic', 'E=%f nu=%f'%(E,nu), 'Sx=0.0 Sy=0.0 Sz=0.0 Sxy=0.0')

# Solve
print 'Solution: ---------------------------------------------------------------------'
sol = m.solver('ForwardEuler')
sol.set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
sol.solve()

# Output file
m.write_vtu_equilib('tpstrain02_py.vtu')
print 'File <tpstrain02_py.vtu> generated'
