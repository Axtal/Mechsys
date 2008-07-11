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
#             ##        A=1       ___  
#    

import mechsys as m

# 0) Geometry
g = m.geom(2) # 2D

# 1) Nodes
g.set_nnodes (3)
g.set_node   (0,  0.0,  0.0)
g.set_node   (1, 10.0,  0.0)
g.set_node   (2, 10.0, 10.0)

# 2) Elements
g.set_nelems (3)
g.set_elem   (0, "Rod", 1)
g.set_elem   (1, "Rod", 1)
g.set_elem   (2, "Rod", 1)

# 3) Set connectivity
g.ele(0).connect(0, g.nod(0)).connect(1, g.nod(1))
g.ele(1).connect(0, g.nod(1)).connect(1, g.nod(2))
g.ele(2).connect(0, g.nod(0)).connect(1, g.nod(2))

# 4) Boundary conditions (must be after set connectivity)
g.nod(0).bry("ux", 0.0).bry("uy", -0.5) # Essential
g.nod(1).               bry("uy",  0.4) # Essential
g.nod(2)                                # Essential
g.nod(2).bry("fx", 2.0).bry("fy",  1.0) # Natural

# 5) Parameters and initial values
g.ele(0).set_model("LinElastic", "E=100.0  A=1.0              ", "Sx=0.0")
g.ele(1).set_model("LinElastic", "E= 50.0  A=1.0              ", "Sx=0.0")
g.ele(2).set_model("LinElastic", "E=200.0  A=1.414213562373095", "Sx=0.0")

# 6) Solve
sol = m.solver('AutoME')
sol.set_geom(g).set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
sol.solve()

# Check
errors = 0.0

errors += abs(g.nod(0).val('ux') - ( 0.0))
errors += abs(g.nod(0).val('uy') - (-0.5))
errors += abs(g.nod(1).val('ux') - ( 0.0))
errors += abs(g.nod(1).val('uy') - ( 0.4))
errors += abs(g.nod(2).val('ux') - (-0.5))
errors += abs(g.nod(2).val('uy') - ( 0.2))

errors += abs(g.nod(0).val('fx') - (-2.0))
errors += abs(g.nod(0).val('fy') - (-2.0))
errors += abs(g.nod(1).val('fx') - ( 0.0))
errors += abs(g.nod(1).val('fy') - ( 1.0))
errors += abs(g.nod(2).val('fx') - ( 2.0))
errors += abs(g.nod(2).val('fy') - ( 1.0))

print 'Errors = ', errors
