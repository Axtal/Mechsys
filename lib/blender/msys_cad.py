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

import math
import Blender
from   Blender.Mathutils import Vector
import bpy
import msys_dict as di


def export_points(filename):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    edm = Blender.Window.EditMode()
    if obj!=None and obj.type=='Mesh':
        if edm: Blender.Window.EditMode(0)
        key = Blender.sys.basename(Blender.sys.splitext(filename)[0])
        fil = open(filename, 'w')
        msh = obj.getData(mesh=1)
        lin = '  X      Y      Z\n'
        for v in msh.verts:
            lin = '%s %8s  %8s  %8s\n'%(lin,v.co[0],v.co[1],v.co[2])
        fil.write(lin)
        fil.close()
        if edm: Blender.Window.EditMode(1)
    else:
        if edm: Blender.Window.EditMode(1)
        raise Exception('Select a Mesh object before calling this function')


def add_point(x, y, z):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    edm = Blender.Window.EditMode()
    if obj==None or not edm:
        msh = bpy.data.meshes.new('points')
        obj = scn.objects.new(msh,'points')
    if obj!=None and obj.type=='Mesh':
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        msh.verts.extend(x,y,z)
        if edm: Blender.Window.EditMode(1)
        obj.select(1)
        Blender.Window.RedrawAll()
    else:
        if edm: Blender.Window.EditMode(1)
        raise Exception('Select a Mesh object before calling this function')


def add_points_from_file(filename):
    Blender.Window.WaitCursor(1)
    edm  = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0) # exit edm
    key  = Blender.sys.basename(Blender.sys.splitext(filename)[0])
    file = open(filename, 'r')
    msh  = bpy.data.meshes.new('points')
    scn  = bpy.data.scenes.active
    obj  = scn.objects.new(msh,key)
    for line in file.readlines():
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif words[0]=='x' or words[0]=='X': pass
        else:
            if len(words)>2: x, y, z = float(words[0]), float(words[1]), float(words[2])
            else:            x, y, z = float(words[0]), float(words[1]), 0
            msh.verts.extend(x,y,z)
    if edm: Blender.Window.EditMode(1) # enter edm
    Blender.Window.RedrawAll()
    Blender.Window.WaitCursor(0)


def add_spline_from_file(filename):
    Blender.Window.WaitCursor(1)
    edm  = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0) # exit edm
    key  = Blender.sys.basename(Blender.sys.splitext(filename)[0])
    file = open(filename, 'r')
    scn  = bpy.data.scenes.active
    spl  = bpy.data.curves.new(key,'Curve')
    obj  = scn.objects.new(spl,key)
    first = 1
    for line in file.readlines():
        words = line.split()
        if len(words)==0 or words[0].startswith('#'): pass
        elif words[0]=='x' or words[0]=='X': pass
        else:
            if len(words)>2: x, y, z = float(words[0]), float(words[1]), float(words[2])
            else:            x, y, z = float(words[0]), float(words[1]), 0
            point = Blender.BezTriple.New((x, y, z))
            if first==1:
                spl.appendNurb(point)
                path  = spl[0]
                first = 0
            else:
                path.append(point)
    for point in path:
        point.handleTypes = [Blender.BezTriple.HandleTypes.AUTO, Blender.BezTriple.HandleTypes.AUTO]
    spl.update()
    if edm: Blender.Window.EditMode(1) # enter edm
    Blender.Window.RedrawAll()
    Blender.Window.WaitCursor(0)


def arc_point(msh,sp_idx,ep_idx,cen,sp,ep,steps):
    spp = msh.verts[sp_idx]
    epp = msh.verts[ep_idx]
    dr1 = sp-cen
    dr2 = ep-cen
    axi = Blender.Mathutils.CrossVecs        (dr2,dr1)
    ang = Blender.Mathutils.AngleBetweenVecs (dr2,dr1)
    rot = Blender.Mathutils.RotationMatrix   (ang/steps, 3, 'r', axi)
    for i in range(steps-1):
        pt  = cen+rot*dr1
        dr1 = pt-cen
        msh.verts.extend(pt)
        if i==0: msh.edges.extend (spp, msh.verts[len(msh.verts)-1])
        else:    msh.edges.extend (msh.verts[len(msh.verts)-2], msh.verts[len(msh.verts)-1])
    msh.edges.extend (msh.verts[len(msh.verts)-1], epp)


def edge_intersect():
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        sel = di.get_selected_edges(msh)
        if len(sel)==2:
            # points (vertices)
            e1 = sel[0]
            e2 = sel[1]
            v1 = msh.edges[e1].v1
            v2 = msh.edges[e1].v2
            v3 = msh.edges[e2].v1
            v4 = msh.edges[e2].v2
            # intersection
            res = Blender.Mathutils.LineIntersect(v1.co,v2.co,v3.co,v4.co)
            if res==None: raise Exception('These edges are parallel (obj=%s)' % obj.name)
            else:
                i1, i2 = res
                if i1!=v1.co and i1!=v2.co and i1!=v3.co and i1!=v4.co and i2!=v1.co and i2!=v2.co and i2!=v3.co and i2!=v4.co:
                    e1v2 = msh.edges[e1].v2
                    e2v2 = msh.edges[e2].v2
                    if i1!=i2:
                        msh.verts.extend (i1)
                        msh.verts.extend (i2)
                        msh.edges[e1].v2 = msh.verts[-2]
                        msh.edges.extend (msh.verts[-2],e1v2)
                        msh.edges.extend (msh.verts[-2], msh.verts[-1])
                        msh.edges[e2].v2 = msh.verts[-1]
                        msh.edges.extend (msh.verts[-1],e2v2)
                    else:
                        msh.verts.extend (i1)
                        msh.edges[e1].v2 = msh.verts[-1]
                        msh.edges.extend (msh.verts[-1],e1v2)
                        msh.edges[e2].v2 = msh.verts[-1]
                        msh.edges.extend (msh.verts[-1],e2v2)
                Blender.Window.RedrawAll()
        else: raise Exception('Please, select exaclty two edges (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1)
    else: raise Exception('Please, select a Mesh object before calling this function')


def fillet(radius,steps):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        sel = di.get_selected_edges(msh)
        if len(sel)==2:
            # points (vertices)
            v1 = msh.edges[sel[0]].v1
            v2 = msh.edges[sel[0]].v2
            v3 = msh.edges[sel[1]].v1
            v4 = msh.edges[sel[1]].v2
            # intersection
            res = Blender.Mathutils.LineIntersect(v1.co,v2.co,v3.co,v4.co)
            if res==None: raise Exception('These edges are parallel (obj=%s)' % obj.name)
            else:
                i1, i2 = res
                if i1!=i2: raise Exception('These two edges do not intersect (the edges must be coplanar) (obj=%s)' % obj.name)
                else:
                    #                p1 ___o b
                    #           a _,--@'
                    #       _,-- @         p3
                    #   i1 @_---------------@
                    #        `-- @_   p2
                    #           c  '--@,____@ d
                    w1 = v1.co-i1
                    w2 = v2.co-i1
                    if w1.length>w2.length: a, b = v2, v1
                    else:                   a, b = v1, v2
                    w1 = v3.co-i1
                    w2 = v4.co-i1
                    if w1.length>w2.length: c, d, c_is_v4 = v4, v3, True
                    else:                   c, d, c_is_v4 = v3, v4, False
                    if radius>0.0:
                        if c==a:
                            msh.verts.extend(c.co)
                            if c_is_v4: msh.edges[sel[1]].v2 = msh.verts[-1]
                            else:       msh.edges[sel[1]].v1 = msh.verts[-1]
                            c = msh.verts[-1]
                        a.co = i1
                        c.co = i1
                        e1   = b.co-a.co; e1.normalize()
                        e2   = d.co-c.co; e2.normalize()
                        e1e2 = e1+e2
                        e3   = e1e2/e1e2.length
                        alp  = math.radians(Blender.Mathutils.AngleBetweenVecs(e1,e2))/2.0
                        l    = radius/math.sin(alp)
                        m    = radius/math.tan(alp)
                        p1   = i1+m*e1
                        p2   = i1+m*e2
                        p3   = i1+l*e3
                        a.co = p1
                        c.co = p2
                        arc_point(msh,c.index,a.index,p3,p2,p1,steps)
                    else:
                        a.co = i1
                        dcor = d.co
                        msh.verts.extend (dcor[0], dcor[1], dcor[2])
                        msh.edges.extend (a, d)
                        msh.edges.delete (sel[1])
                    Blender.Window.RedrawAll()
        else: raise Exception('Please, select exaclty two edges (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1)
    else: raise Exception('Please, select a Mesh object before calling this function')


def break_edge(at_mid=False):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        sel = di.get_selected_edges(msh)
        if len(sel)==1 or at_mid:
            for e in sel:
                # points (vertices)
                v1 = msh.edges[e].v1
                v2 = msh.edges[e].v2
                # break-point
                if at_mid: bp = 0.5*(v1.co+v2.co)
                else:      bp = Blender.Window.GetCursorPos()
                if bp!=v1.co and bp!=v2.co:
                    msh.verts.extend (bp[0],bp[1],bp[2])
                    new = msh.verts[-1]
                    v2  = msh.edges[e].v2
                    msh.edges[e].v2 = new
                    msh.edges.extend (new,v2)
        else: raise Exception('Please, select only one edge (obj=%s)' % obj.name)
        if edm: Blender.Window.EditMode(1)
    else: raise Exception('Please, select a Mesh object before calling this function')
