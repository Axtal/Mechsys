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
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program. If not, see <http://www.gnu.org/licenses/>  #
########################################################################

#      F1        F2       F3    F1=F2=F3 = 1.0
#       ^         ^         ^
#       |         |         |
#
#     2 @---------@---------@ 6
#       |',       7         |
#       |  ',        e[1]   |
#       |    ',             |
#       |      ',  4        |
#     5 @        '@         @ 8
#       |          ',       |
#       |   e[0]     ',     |
#       |              ',   |
#       |         3      ', |
#     0 @---------@---------@ 1
#      /_\       /_\       /_\
#      ///       o o       o o
#
import msysfem as m

# 0) Geometry type
m.dim(2) # 2D

# 1) Nodes
m.add_node(0.0, 0.0) # 0
m.add_node(1.0, 0.0) # 1
m.add_node(0.0, 1.0) # 2
m.add_node(0.5, 0.0) # 3
m.add_node(0.5, 0.5) # 4
m.add_node(0.0, 0.5) # 5
m.add_node(1.0, 1.0) # 6
m.add_node(0.5, 1.0) # 7
m.add_node(1.0, 0.5) # 8

# 2) Elements
m.add_elem('Tri6PStrain', 1) # 0: 1=>Active
m.add_elem('Tri6PStrain', 1) # 1: 1=>Active

# 3) Set connectivity
m.elems(0).set_node(0,0).set_node(1,1).set_node(2,2).set_node(3,3).set_node(4,4).set_node(5,5)
m.elems(1).set_node(0,6).set_node(1,2).set_node(2,1).set_node(3,7).set_node(4,4).set_node(5,8)

# 4) Boundary conditions (must be after connectivity)
m.nodes(0).bry('ux', 0.0).bry('uy', 0.0)
m.nodes(1).bry('uy', 0.0)
m.nodes(3).bry('uy', 0.0)
m.nodes(2).bry('fy', 1.0)
m.nodes(7).bry('fy', 1.0)
m.nodes(6).bry('fy', 1.0)

# 5) Parameters and initial values
m.elems(0).set_model('LinElastic', 'E=10000.0 nu=0.25', 'Sx=0.0 Sy=0.0 Sxy=0.0')
m.elems(1).set_model('LinElastic', 'E=10000.0 nu=0.25', 'Sx=0.0 Sy=0.0 Sxy=0.0')

# 6) Solve
sol = m.solver('ForwardEuler')
sol.set_lin_sol('LA').set_num_div(1).set_delta_time(0.0)
sol.solve()

# Check
errors = 0.0;

errors += abs(m.elems(0).val(0, "Sx") - ( 1.56432140e-01))
errors += abs(m.elems(0).val(1, "Sx") - (-3.00686928e-01))
errors += abs(m.elems(0).val(2, "Sx") - ( 1.44254788e-01))
errors += abs(m.elems(0).val(3, "Sx") - (-3.19109076e-01))
errors += abs(m.elems(0).val(4, "Sx") - (-3.31286428e-01))
errors += abs(m.elems(0).val(5, "Sx") - ( 1.25832639e-01))

errors += abs(m.elems(0).val(0, "Sy") - (-2.05141549e-01))
errors += abs(m.elems(0).val(1, "Sy") - ( 1.15872190e+00))
errors += abs(m.elems(0).val(2, "Sy") - (-9.53580350e-01))
errors += abs(m.elems(0).val(3, "Sy") - (-2.22127394e+00))
errors += abs(m.elems(0).val(4, "Sy") - (-2.96971274e+00))
errors += abs(m.elems(0).val(5, "Sy") - (-4.33357619e+00))

errors += abs(m.elems(0).val(0, "Sxy") - (-1.56432140e-01))
errors += abs(m.elems(0).val(1, "Sxy") - (-6.74437968e-02))
errors += abs(m.elems(0).val(2, "Sxy") - ( 2.23875937e-01))
errors += abs(m.elems(0).val(3, "Sxy") - (-4.90216486e-02))
errors += abs(m.elems(0).val(4, "Sxy") - ( 3.31286428e-01))
errors += abs(m.elems(0).val(5, "Sxy") - ( 2.42298085e-01))

print 'Errors = ', errors
