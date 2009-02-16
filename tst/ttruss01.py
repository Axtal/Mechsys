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

#      Small truss
#                                1.0 ^
#                                    |
#                                    |2
#                                    o---. 2.0
#                                  ,'|
#                                ,'  |
#                     E=200    ,'    |
#                     A=SQ2  ,'      | E=50
#                      [2] ,'        | A=1
#                        ,'          | [1]
#                      ,'            |
#                    ,'              |
#     y            ,'                |
#     |         0,'        [0]       |1
#     |         o--------------------o
#    (z)___x   /_\        E=100     /_\
#              ///        A=1       ___  
#    

import math
import mechsys as ms

#///////////////////////////////////////////////////////////////////////////////////////// Mesh /////

mesh = ms.mesh_generic (False) # False => 2D
mesh.set_nverts   (3)
mesh.set_nelems   (3)
mesh.set_vert     (0, True,  0.0,  0.0) # true => OnBry
mesh.set_vert     (1, True, 10.0,  0.0)
mesh.set_vert     (2, True, 10.0, 10.0)
mesh.set_elem     (0, -1, True, 3) # true => OnBry, 3 => VTK_LINE
mesh.set_elem     (1, -2, True, 3)
mesh.set_elem     (2, -3, True, 3)
mesh.set_elem_con (0, 0, 0);  mesh.set_elem_con(0, 1, 1)
mesh.set_elem_con (1, 0, 1);  mesh.set_elem_con(1, 1, 2)
mesh.set_elem_con (2, 0, 0);  mesh.set_elem_con(2, 1, 2)

#////////////////////////////////////////////////////////////////////////////////////////// FEM /////

# Data and Solver
dat = ms.data   (2) # 2D
sol = ms.solver (dat, "ttruss01_py")

# Elements attributes
eatts = [[-1, "Lin2", "Rod", "RodElastic", "E=100.0  A=1.0",               "ZERO", "gam=20", True],
         [-2, "Lin2", "Rod", "RodElastic", "E= 50.0  A=1.0",               "ZERO", "gam=20", True],
         [-3, "Lin2", "Rod", "RodElastic", "E=200.0  A=%f"%math.sqrt(2.0), "ZERO", "gam=20", True]]

# Set geometry: nodes and elements
dat.set_nodes_elems (mesh, eatts)

# Stage # 0 --------------------------------------------------------------
dat.nod(0).bry("ux", 0.0).bry("uy", -0.5) # Essential
dat.nod(1).               bry("uy",  0.4) # Essential
dat.nod(2).bry("fx", 2.0).bry("fy",  1.0) # Natural
sol.solve_with_info (1, 0.0, 0, 'First Stage') # nDiv,DTime,iStage,MoreInfo

#//////////////////////////////////////////////////////////////////////////////////////// Check /////

# Check
errors = 0.0

errors += abs(dat.nod(0).val('ux') - ( 0.0))
errors += abs(dat.nod(0).val('uy') - (-0.5))
errors += abs(dat.nod(1).val('ux') - ( 0.0))
errors += abs(dat.nod(1).val('uy') - ( 0.4))
errors += abs(dat.nod(2).val('ux') - (-0.5))
errors += abs(dat.nod(2).val('uy') - ( 0.2))

errors += abs(dat.nod(0).val('fx') - (-2.0))
errors += abs(dat.nod(0).val('fy') - (-2.0))
errors += abs(dat.nod(1).val('fx') - ( 0.0))
errors += abs(dat.nod(1).val('fy') - ( 1.0))
errors += abs(dat.nod(2).val('fx') - ( 2.0))
errors += abs(dat.nod(2).val('fy') - ( 1.0))

print '\nPy:Errors = ', errors
