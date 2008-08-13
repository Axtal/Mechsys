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
import timeit

########################################################################################### 2D #####

def test_2D():
    # Constants
    L = 1.0   # length
    H = 1.0   # height
    
    # Blocks
    blocks = [m.mesh_block()]
    blocks[0].set_coords (-1,                                       # tag to be replicated to all elements
                          False,                                    # is 3D
                          [(0,0), (L,0), (L,H), (0,H)],             # vertices
                          [(0,1), (1,2), (2,3), (3,0)],             # edges
                          {(0,1):-1, (1,2):-2, (2,3):-3, (3,0):-4}, # edge tags
                          {},                                       # face tags
                          [1,1],                                    # weights x
                          [1,1],                                    # weights y
                          [],                                       # weights z
                          0, 1, 3, 0)                               # Origin, XPlus, YPlus, None

    # Generate
    ms = m.mesh_structured()
    ne = ms.generate (blocks)
    print '2D => Generated ', ne, ' elements'

    # Output
    ms.write_vtu ('tmesh01_2D_py.vtu')
    print 'File <tmesh01_2D_py.vtu> created'

########################################################################################## Run #####
    
# 2D: Run and check speed
t2d = timeit.Timer('test_2D()', 'from __main__ import test_2D')
print t2d.timeit(number=1), ' seconds'
print ''
