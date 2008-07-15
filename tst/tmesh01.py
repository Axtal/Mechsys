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
    L  = 1.0   # length
    H  = 1.0   # height
    hL = L/2.0 # half-length
    hH = H/2.0 # half-height
    
    # Blocks
    blocks = [m.mesh_block()]
    blocks[0].set2d (-1,                                   # tag to be replicated to all elements
                     [[0.,  L, L, 0.,    hL,  L, hL, 0.] , # coordinates: x values
                      [0., 0., H,  H,    0., hH,  H, hH]], # coordinates: y values
                     [2,4,8,16,32],                        # weights x
                     [2,4,8,16,32,64])                     # weights y

    # Tags
    blocks[0].set_etags ([-10, -11, -12, -13])

    # Generate
    ms = m.mesh_structured()
    ne = ms.generate (blocks)
    print '2D => Generated ', ne, ' elements'

    # Output
    ms.write_vtu ('tmesh01_2D_py.vtu')
    print 'File <tmesh01_2D_py.vtu> created'

########################################################################################### 3D #####

def test_3D():
    # Constants
    L  = 1.0   # length
    H  = 1.0   # height
    D  = 1.0   # depth
    hL = L/2.0 # half-length
    hH = H/2.0 # half-height
    hD = D/2.0 # half-depth
    
    # Blocks
    blocks = [m.mesh_block()]
    blocks[0].set3d (-2,                                                                                            # tag to be replicated to all elements
                     [[0.,  L,  L, 0.,    0.,  L, L, 0.,    hL,  L, hL, 0.,    hL,  L, hL, 0.,    0.,  L,  L, 0.],  # coordinates: x values
                      [0., 0.,  H,  H,    0., 0., H,  H,    0., hH,  H, hH,    0., hH,  H, hH,    0., 0.,  H,  H],  # coordinates: y values
                      [0., 0., 0., 0.,     D,  D, D,  D,    0., 0., 0., 0.,     D,  D,  D,  D,    hD, hD, hD, hD]], # coordinates: z values
                     [2,4,8,16,32],                                                                                 # weights x
                     [2,4,8,16,32,64],                                                                              # weights y
                     [2,4,8,16,32,64,128])                                                                          # weights y

    # Tags
    blocks[0].set_etags ([-100, -101, -102, -103, -104, -105, -106, -107, -108, -109, -110, -111])
    blocks[0].set_ftags ([-10, -11, -12, -13, -14, -15])

    # Generate
    ms = m.mesh_structured()
    ne = ms.generate (blocks)
    print '3D => Generated ', ne, ' elements'

    # Output
    ms.write_vtu ('tmesh01_3D_py.vtu')
    print 'File <tmesh01_3D_py.vtu> created'
    
########################################################################################## Run #####
    
# 2D: Run and check speed
t2d = timeit.Timer('test_2D()', 'from __main__ import test_2D')
print t2d.timeit(number=1), ' seconds'
print ''

# 3D: Run and check speed
t3d = timeit.Timer('test_3D()', 'from __main__ import test_3D')
print t3d.timeit(number=1), ' seconds'
