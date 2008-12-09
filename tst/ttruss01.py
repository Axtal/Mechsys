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

# Data and Solver
dat = ms.data   (2) # 2D
sol = ms.solver (dat, "ttruss01_py")

# Nodes
dat.set_nnodes (3)
dat.set_node   (0,  0.0,  0.0)
dat.set_node   (1, 10.0,  0.0)
dat.set_node   (2, 10.0, 10.0)

# Elements
dat.set_nelems (3)
dat.set_elem   (0, "Rod2", True, -1)
dat.set_elem   (1, "Rod2", True, -1)
dat.set_elem   (2, "Rod2", True, -1)

# Set connectivity
dat.ele(0).connect(0, dat.nod(0)).connect(1, dat.nod(1))
dat.ele(1).connect(0, dat.nod(1)).connect(1, dat.nod(2))
dat.ele(2).connect(0, dat.nod(0)).connect(1, dat.nod(2))

# Parameters and initial values
dat.ele(0).set_model("LinElastic", "E=100.0  A=1.0",               "Sx=0.0")
dat.ele(1).set_model("LinElastic", "E= 50.0  A=1.0",               "Sx=0.0")
dat.ele(2).set_model("LinElastic", "E=200.0  A=%f"%math.sqrt(2.0), "Sx=0.0")

# State # 1
dat.nod(0).bry("ux", 0.0).bry("uy", -0.5) # Essential
dat.nod(1).               bry("uy",  0.4) # Essential
dat.nod(2).bry("fx", 2.0).bry("fy",  1.0) # Natural
sol.solve_with_info()

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
