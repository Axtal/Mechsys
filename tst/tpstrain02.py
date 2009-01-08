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

import mechsys as ms

# Constants
L   = 1.0   # Length
H   = 1.0   # Height
E   = 207.0 # Young
nu  = 0.3   # Poisson
gam = 20   # Specific weight
q   = -1.0  # Load
nx  = 10    # Divisions along x (must be even)
ny  = 10    # Divisions along y

# ----------------------------------------------------------------------------- Mesh

#                      55        1
#                       Y+-------@
#                       |        |
#                     H |        |
#                       |    L   |
#                       O--------X+
#                       8        4

# Blocks
blocks = [ms.mesh_block()]
blocks[0].set    (-1,                                    # tag to be replicated to all elements
                  [(8,0,0), (4,L,0), (1,L,H), (55,0,H)], # vertices' coordinates
                  [(8,4), (4,1), (1,55), (55,8)],        # edges
                  [(8,4,-10), (1,55,-20)],               # edge tags
                  [],                                    # face tags
                  8, 4, 55, 0)                           # Origin, XPlus, YPlus, None
blocks[0].set_nx (nx)
blocks[0].set_ny (ny)

# Generate
print 'Mesh Generation: --------------------------------------------------------------'
mesh = ms.mesh_structured (False) # Is3D==False
mesh.set_blocks           (blocks)
mesh.generate             (True)

# ------------------------------------------------------------------------------ FEM

# Data and Solver
dat = ms.data   (2)
sol = ms.solver (dat, "tpstrain02_py")

# Elements attributes
eatts = [[-1, 'Quad4', 'PStrain', 'LinElastic', 'E=%f nu=%f'%(E,nu), 'ZERO', 'gam=%f'%(gam), True]]

# Set geometry: nodes and elements
dat.set_nodes_elems (mesh, eatts)

# Stage # 1 --------------------------------------
nbrys = [[L/2., 0.0, 0.0, 'ux', 0.0]]
ebrys = [[-10, 'uy', 0.0], [-20, 'fy',   q]]
dat.set_brys        (mesh, nbrys, ebrys, []) # [] => no face brys
sol.solve_with_info ()

#----------------------------------------------------------------------------- Check

err_eps = []
err_sig = []
err_dis = []

Sx  = 0.0
Sy  = q
Ex  = -nu*(1.0+nu)*Sy/E
Ey  =  (1.0-nu*nu)*Sy/E
Ez  = 0.0
Exy = 0.0
Sz  = (E/(1.0+nu))*(nu/(1.0-2.0*nu))*(Ex+Ey)
Sxy = 0.0

# Stress and strains
for i in range(dat.nelems()):
    for j in range(dat.ele(i).nnodes()):
        err_eps.append ( abs(dat.ele(i).val(j,"Ex" ) - Ex ) / (1.0+abs(Ex )) )
        err_eps.append ( abs(dat.ele(i).val(j,"Ey" ) - Ey ) / (1.0+abs(Ey )) )
        err_eps.append ( abs(dat.ele(i).val(j,"Ez" ) - Ez ) / (1.0+abs(Ez )) )
        err_eps.append ( abs(dat.ele(i).val(j,"Exy") - Exy) / (1.0+abs(Exy)) )
        err_sig.append ( abs(dat.ele(i).val(j,"Sx" ) - Sx ) / (1.0+abs(Sx )) )
        err_sig.append ( abs(dat.ele(i).val(j,"Sy" ) - Sy ) / (1.0+abs(Sy )) )
        err_sig.append ( abs(dat.ele(i).val(j,"Sz" ) - Sz ) / (1.0+abs(Sz )) )
        err_sig.append ( abs(dat.ele(i).val(j,"Sxy") - Sxy) / (1.0+abs(Sxy)) )

# Displacements
for i in range(dat.nnodes()):
    ux_correct = Ex*(dat.nod(i).x()-L/2.)
    uy_correct = Ey* dat.nod(i).y()
    err_dis.append ( abs(dat.nod(i).val("ux") - ux_correct) / (1.0+abs(ux_correct)) )
    err_dis.append ( abs(dat.nod(i).val("uy") - uy_correct) / (1.0+abs(uy_correct)) )

# Error summary
min_err_eps = min(err_eps)
min_err_sig = min(err_sig)
min_err_dis = min(err_dis)
max_err_eps = max(err_eps)
max_err_sig = max(err_sig)
max_err_dis = max(err_dis)
print
print 'Errors: -----------------------------------------------------------------------'
print "      Min            Max"
print "Eps   %8e   %8e" % (min_err_eps,max_err_eps)
print "Sig   %8e   %8e" % (min_err_sig,max_err_sig)
print "Dis   %8e   %8e" % (min_err_dis,max_err_dis)
print
