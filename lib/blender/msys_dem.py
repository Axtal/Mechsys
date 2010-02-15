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

from mechsys import *
from Blender import Mesh
from bpy     import data
from Blender.Mathutils import Vector

d = DEM_Domain()
d.AddCube (-1, (20,0,0), 1.0, 15, 1.0, 0.0, (1,0,0))
d.SetCamPos ((20,20,20))
P = []
d.GetParticles (P)
r = 1
scn = data.scenes.active

for p in P:
    # features
    V, E, F = p[0], p[1], p[2]

    # vertices
    for v in p[0]:
        m = Mesh.Primitives.UVsphere(16,16,r*2)
        o = scn.objects.new (m,'vert')
        o.setLocation(v[0],v[1],v[2])

    # edges
    for e in E:
        x0  = Vector(V[e[0]])
        x1  = Vector(V[e[1]])
        dx  = x1-x0
        mid = (x0+x1)/2.0
        m   = Mesh.Primitives.Tube(16,r*2,dx.length)
        o   = scn.objects.new (m,'edge')
        qua = dx.toTrackQuat ('z','x')
        o.setMatrix   (qua.toMatrix())
        o.setLocation (mid)

    # faces
    if len(F)>0:
        msh = data.meshes.new ('faces')
        msh.verts.extend      (V)
        obj = scn.objects.new (msh,'faces')
    for f in F:
        #ed0 = Vector(V[f[0][1]]) - Vector(V[f[0][0]])
        #ed1 = Vector(V[f[1][1]]) - Vector(V[f[1][0]])
        ed0 = Vector(V[f[1]]) - Vector(V[f[0]])
        ed1 = Vector(V[f[2]]) - Vector(V[f[1]])
        n   = ed0.cross(ed1)
        msh.faces.extend ([f])
