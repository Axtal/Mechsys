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

import mechsys           as ms
import numpy             as np
import matplotlib.pyplot as plt

# Stress states
states = [(-2.0,-2.0,-2.0), (-1.08,-1.08,-3.84)] # tx, ty, tz

# Number of divisions
gndiv = 1 # global number of divisions for states
sndiv = 1 # solver number of divisions

# Parameters
lam  = 0.0778
kap  = 0.00824
phi  = 35.0
nu   = 0.0
vini = 1.699

k1   = 10.0
l1   = 0.1
b1   = 0.1
psi1 = 1.5

k2   = 10.0
l2   = 0.1
b2   = 0.1
psi2 = 2.0

k3   = 10.0
l3   = 0.1
b3   = 0.1
ev3  = -4.0/100

k4   = 10.0
l4   = 0.1
b4   = 0.1
ev4  = -4.0/100

# Mesh
#nxyz = 3  # Divisions along x
#iele = 13 # Index of element for output (Centre)
nxyz = 1
iele = 0

############################################################################## Mesh

#                     4----------------7
#                   ,'|              ,'|
#                 ,'  |            ,'  |
#               ,'    | -6    -1 ,'    |
#             ,'      |        ,'      |
#           5'===============6'        |
#           |         |      |    -4   |
#           |    -3   |      |         |
#           |         0- - - | -  - - -3
#           |       ,'       |       ,'
#           |     ,' -2      |     ,'
#           |   ,'        -5 |   ,'
#           | ,'             | ,'
#           1----------------2'

# Generate
mesh = ms.mesh_structured (True) # True => 3D
#mesh.set_o2               (True)
mesh.gen_box              (nxyz,nxyz,nxyz)

############################################################################### FEM

# Data and Solver
dat = ms.data   (3)
sol = ms.solver (dat, 'ttriax')

# Elements attributes
Sx = states[0][0]
Sy = states[0][1]
Sz = states[0][2]
#eatts = [[-1, 'Hex8', 'Equilib', 'PyEquilib,camclay.py', 'a0=%f a1=%f a2=%f a3=%f'%(lam,kap,phi,nu), 'Sx=%f Sy=%f Sz=%f v=%f'%(Sx,Sy,Sz,vini), 'gam=0', '', True]]
eatts = [[-1, 'Hex8', 'Equilib', 'LinElastic', 'E=%f nu=%f'%(1200,0.2), 'Sx=%f Sy=%f Sz=%f'%(Sx,Sy,Sz), 'gam=20', '', True]]
#eatts = [[-1, 'Hex8', 'Equilib', 'CamClay', 'lam=%f kap=%f phics=%f nu=%f'%(lam,kap,phi,nu), 'Sx=%f Sy=%f Sz=%f v=%f'%(Sx,Sy,Sz,vini), 'gam=2', '', True]]
#eatts = [[-1, 'Hex8', 'Equilib', 'Toto', 'k1=%f l1=%f b1=%f psi1=%f  k2=%f l2=%f b2=%f psi2=%f  k3=%f l3=%f b3=%f ev3=%f  k4=%f l4=%f b4=%f ev4=%f'%(k1,l1,b1,psi1, k2,l2,b2,psi2, k3,l3,b3,ev3, k4,l4,b4,ev4), 'Sx=%f Sy=%f Sz=%f'%(Sx,Sy,Sz), 'gam=0', '', True]]

# Set geometry: nodes and elements
dat.set_nodes_elems (mesh, eatts)

# Solve
P  = [dat.ele(iele).val('p')]
Q  = [dat.ele(iele).val('q')]
EV = [dat.ele(iele).val('Ev')]
ED = [dat.ele(iele).val('Ed')]
k  = 0
for i in range(1,len(states)):
    # increments
    dfx = (states[i][0]-Sx)/gndiv
    dfy = (states[i][1]-Sy)/gndiv
    dfz = (states[i][2]-Sz)/gndiv
    for j in range(gndiv):
        # solve
        fbrys = [[-1, 'ux', 0.0], [-3, 'uy', 0.0], [-5, 'uz', 0.0],
                 [-2, 'fx', dfx], [-4, 'fy', dfy], [-6, 'fz', dfz]]
        dat.set_brys (mesh, [], [], fbrys)
        sol.solve_with_info (sndiv,0.0,k,'')
        # results
        P .append(dat.ele(iele).val('p'))
        Q .append(dat.ele(iele).val('q'))
        EV.append(dat.ele(iele).val('Ev'))
        ED.append(dat.ele(iele).val('Ed'))
        # update state
        Sx += dfx
        Sy += dfy
        Sz += dfz
        k  += 1

############################################################################### Results

# Plot results
P  = np.array(P)
Q  = np.array(Q)
EV = np.array(EV)*100.0
ED = np.array(ED)*100.0
plt.rc ('text', usetex=True)
plt.rc ('font', family='serif')

# -q/p x Ed
plt.subplot (2,3,1)
plt.plot    (ED, -Q/P)
plt.xlabel  (r'$\varepsilon_d$ [\%]',fontsize=18)
plt.ylabel  (r'$-q/p$',fontsize=18)
plt.grid    ()

# -q/p x -Ev
plt.subplot (2,3,2)
plt.plot    (-EV, -Q/P)
plt.xlabel  (r'$-\varepsilon_v$ [\%]',fontsize=18)
plt.ylabel  (r'$-q/p$',fontsize=18)
plt.grid    ()

# q x -p
plt.subplot (2,3,3)
plt.plot    (-P, Q)
plt.xlabel  (r'-$p$',fontsize=18)
plt.ylabel  (r'$q$',fontsize=18)
plt.grid    ()

# Ev x Ed
plt.subplot (2,3,4)
plt.plot    (ED, EV)
plt.xlabel  (r'$\varepsilon_d$ [\%]',fontsize=18)
plt.ylabel  (r'$\varepsilon_v$ [\%]',fontsize=18)
plt.grid    ()

# Ev x -p
plt.subplot (2,3,5)
plt.plot    (-P, EV)
plt.xlabel  (r'-$p$',fontsize=18)
plt.ylabel  (r'$\varepsilon_v$ [\%]',fontsize=18)
plt.grid    ()

# Show plot
plt.show ()
