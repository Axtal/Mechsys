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
import bpy


def add_point(x, y, z):
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    if obj==None:
        msh = bpy.data.meshes.new('points')
        obj = scn.objects.new(msh,'points')
    if obj!=None and obj.type=='Mesh':
        edm = Blender.Window.EditMode()
        if edm: Blender.Window.EditMode(0)
        msh = obj.getData(mesh=1)
        msh.verts.extend(x,y,z)
        if edm: Blender.Window.EditMode(1)
        obj.select(1)
        Blender.Window.RedrawAll()
    else:
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


def distance(p1,p2):
    d = p1 - p2
    return math.sqrt(math.pow(d[0],2)+math.pow(d[1],2)+math.pow(d[2],2))


def closest(pt,pointlist):
    idx = 0
    min = distance(pt,pointlist[idx])
    for i in range(1,len(pointlist)):
        d = distance(pt,pointlist[i])
        if d<min:
            min = d
            idx = i
    return min, idx


def arc_point(msh,cen,sp,ep,steps):
    dr1 = sp-cen
    dr2 = ep-cen
    axi = Blender.Mathutils.CrossVecs(dr2,dr1)
    ang = Blender.Mathutils.AngleBetweenVecs(dr2,dr1)
    rot = Blender.Mathutils.RotationMatrix(ang/steps, 3, 'r', axi)
    for i in range(steps-1):
        pt  = cen+rot*dr1
        msh.verts.extend(pt)
        msh.edges.extend(msh.verts[len(msh.verts)-2],msh.verts[len(msh.verts)-1])
        dr1 = pt-cen


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
                # fillet
                if i1!=i2: raise Exception('These two edges do not intersect (obj=%s)' % obj.name)
                else:
                    #
                    #                   ep  _,--* p2
                    #                  _,-*'
                    #     a=|ep-p1|,--' /
                    #        _,--'     /radius  cen
                    #    p1 *_--------|----------*-  d=|c-p1|
                    #         `--,_    \
                    #              `--,_\
                    #                   `-*,_
                    #                   sp   `--* p3
                    di, i = closest(i1,[v1.co,v2.co])
                    dj, j = closest(i1,[v3.co,v4.co])
                    if i==0:
                        v1.co = i1
                        p1,p2 = v1,v2
                    else:
                        v2.co = i1
                        p1,p2 = v2,v1
                    if j==0:
                        msh.edges[sel[1]].v1 = p1
                        if v3!=p1: msh.verts.delete(v3)
                        p3 = msh.edges[sel[1]].v2
                    else:
                        msh.edges[sel[1]].v2 = p1
                        if v4!=p1: msh.verts.delete(v4)
                        p3 = msh.edges[sel[1]].v1
                    if radius>0.0:
                        # vectors along the edges
                        e1  = p2.co-p1.co; e1.normalize()
                        e2  = p3.co-p1.co; e2.normalize()
                        alp = Blender.Mathutils.AngleBetweenVecs(e1,e2)/2.0
                        a   = radius/math.tan(math.radians(alp))
                        d   = math.sqrt(math.pow(a,2)+math.pow(radius,2))
                        ep  = p1.co+a*e1
                        sp  = p1.co+a*e2
                        mid = Blender.Mathutils.MidpointVecs(sp,ep); mid-=p1.co; mid.normalize()
                        cen = p1.co+d*mid
                        # add points to the arch
                        msh.verts.extend(sp)
                        msh.edges.extend(p3,msh.verts[len(msh.verts)-1])
                        arc_point(msh,cen,sp,ep,steps)
                        msh.verts.extend(ep)
                        msh.edges.extend(msh.verts[-2],msh.verts[-1])
                        msh.edges.extend(msh.verts[-1],p2)
                        msh.verts.delete(p1)
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
