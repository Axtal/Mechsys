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
# but WITHOUT ANY WARRANTY without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http:#www.gnu.org/licenses/>  #
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

import msysfem as m

# 0) Geometry type
m.dim(2) # 2D

# 1) Nodes
m.add_node( 0.0,  0.0) # 0
m.add_node(10.0,  0.0) # 1
m.add_node(10.0, 10.0) # 2

# 2) Elements
m.add_elem("Rod", 1) # 0
m.add_elem("Rod", 1) # 1
m.add_elem("Rod", 1) # 2

# 3) Set connectivity
m.elems(0).set_node(0, 0).set_node(1, 1)
m.elems(1).set_node(0, 1).set_node(1, 2)
m.elems(2).set_node(0, 0).set_node(1, 2)

# 4) Boundary conditions (must be after set connectivity)
m.nodes(0).bry("ux", 0.0).bry("uy", -0.5) # Essential
m.nodes(1).               bry("uy",  0.4) # Essential
m.nodes(2)                                # Essential
m.nodes(2).bry("fx", 2.0).bry("fy",  1.0) # Natural

# 5) Parameters and initial values
m.elems(0).set_model("LinElastic", "E=100.0  A=1.0              ", "Sx=0.0")
m.elems(1).set_model("LinElastic", "E= 50.0  A=1.0              ", "Sx=0.0")
m.elems(2).set_model("LinElastic", "E=200.0  A=1.414213562373095", "Sx=0.0")

# 6) Solve
sol = m.solver('AutoME')
sol.set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
sol.solve()

# Check
errors = 0.0

errors += abs(m.nodes(0).val('ux') - ( 0.0))
errors += abs(m.nodes(0).val('uy') - (-0.5))
errors += abs(m.nodes(1).val('ux') - ( 0.0))
errors += abs(m.nodes(1).val('uy') - ( 0.4))
errors += abs(m.nodes(2).val('ux') - (-0.5))
errors += abs(m.nodes(2).val('uy') - ( 0.2))

errors += abs(m.nodes(0).val('fx') - (-2.0))
errors += abs(m.nodes(0).val('fy') - (-2.0))
errors += abs(m.nodes(1).val('fx') - ( 0.0))
errors += abs(m.nodes(1).val('fy') - ( 1.0))
errors += abs(m.nodes(2).val('fx') - ( 2.0))
errors += abs(m.nodes(2).val('fy') - ( 1.0))

print 'Errors = ', errors
