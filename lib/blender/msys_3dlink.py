#!BPY

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

import Blender
from   Blender import BGL, Draw, Window
from   Blender.Mathutils import Vector
import bpy
import math
import msys_dict as di

def sgn(val):
    if val<0.0: return -1.0
    else:       return  1.0

def draw_arrow_2d(x0,y0,dx,dy, h=-1.0, alp=15.0):
    L  = math.sqrt(dx*dx + dy*dy)
    if h<0: h = 0.1*L
    h  = h if h<=L else L
    s  = h*math.tan(alp*math.pi/360.0)
    xm = x0+(L-h)*dx/L
    ym = y0+(L-h)*dy/L
    BGL.glBegin    (BGL.GL_LINES)
    BGL.glVertex3f (x0,    y0,    0.0)
    BGL.glVertex3f (x0+dx, y0+dy, 0.0)
    BGL.glEnd      ()
    BGL.glBegin    (BGL.GL_POLYGON)
    BGL.glVertex3f (xm-s*dy/L, ym+s*dx/L, 0.0)
    BGL.glVertex3f (xm+s*dy/L, ym-s*dx/L, 0.0)
    BGL.glVertex3f (x0+dx,     y0+dy,     0.0)
    BGL.glEnd      ()

def draw_xfix_2d(x0,y0, h, dx):
    BGL.glBegin    (BGL.GL_LINES)
    BGL.glVertex3f (x0-dx/2.0, y0-h/2.0, 0.0)
    BGL.glVertex3f (x0-dx/2.0, y0+h/2.0, 0.0)
    BGL.glEnd      ()
    BGL.glBegin    (BGL.GL_LINES)
    BGL.glVertex3f (x0+dx/2.0, y0-h/2.0, 0.0)
    BGL.glVertex3f (x0+dx/2.0, y0+h/2.0, 0.0)
    BGL.glEnd      ()

def draw_yfix_2d(x0,y0, h, dy):
    BGL.glBegin    (BGL.GL_LINES)
    BGL.glVertex3f (x0-h/2.0, y0-dy/2.0, 0.0)
    BGL.glVertex3f (x0+h/2.0, y0-dy/2.0, 0.0)
    BGL.glEnd      ()
    BGL.glBegin    (BGL.GL_LINES)
    BGL.glVertex3f (x0-h/2.0, y0+dy/2.0, 0.0)
    BGL.glVertex3f (x0+h/2.0, y0+dy/2.0, 0.0)
    BGL.glEnd      ()


# Transformation matrix
if di.key('show_props'):
    # Buffer
    view_matrix = Window.GetPerspMatrix()
    view_buffer = [view_matrix[i][j] for i in xrange(4) for j in xrange(4)]
    view_buffer = BGL.Buffer(BGL.GL_FLOAT, 16, view_buffer)

    # Initialize
    BGL.glLoadIdentity ()
    BGL.glMatrixMode   (Blender.BGL.GL_PROJECTION)
    BGL.glPushMatrix   ()
    BGL.glLoadMatrixf  (view_buffer)


# Mesh properties
if di.key('show_props'):
    # Draw
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for obj in obs:
        if obj!=None and obj.type=='Mesh':

            # draw only if active layer corresponds to this object.Layer
            if Blender.Window.GetActiveLayer()==obj.Layer:

                # get mesh and transform to global coordinates
                msh = obj.getData(mesh=1)
                ori = [v for v in msh.verts] # create a copy in local coordinates
                msh.transform(obj.matrix)

                # show overlapping edges
                if obj.properties.has_key('over'):
                    for i, eid in enumerate(obj.properties['over']):
                        ed = msh.edges[eid]
                        if i%2==0: BGL.glColor3f (0.0, 0.757, 1.0)
                        else:      BGL.glColor3f (0.0, 0.0,   0.9)
                        BGL.glBegin    (BGL.GL_LINES)
                        BGL.glVertex3f (ed.v1.co[0], ed.v1.co[1], ed.v1.co[2])
                        BGL.glVertex3f (ed.v2.co[0], ed.v2.co[1], ed.v2.co[2])
                        BGL.glEnd      ()

                # draw vertices IDs
                if di.key('show_v_ids'):
                    BGL.glColor3f (1.0, 1.0, 0.0)
                    for v in msh.verts:
                        BGL.glRasterPos3f (v.co[0], v.co[1], v.co[2])
                        Draw.Text         (str(v.index))

                # draw nodes IDs
                if di.key('show_n_ids') and obj.properties.has_key('msh_name'):
                    msh_obj = bpy.data.objects[obj.properties['msh_name']]
                    msh_msh = msh_obj.getData(mesh=1)
                    msh_ori = [v for v in msh_msh.verts] # create a copy before transforming to global coordinates
                    msh_msh.transform(msh_obj.matrix)
                    BGL.glColor3f (1.0, 1.0, 0.0)
                    for v in msh_msh.verts:
                        BGL.glRasterPos3f (v.co[0], v.co[1], v.co[2])
                        Draw.Text         (str(v.index))
                    msh_msh.verts = msh_ori # restore mesh to local coordinates

                # draw edges IDs
                if di.key('show_e_ids'):
                    BGL.glColor3f (1.0, 1.0, 1.0)
                    for e in msh.edges:
                        mid = 0.5*(e.v1.co+e.v2.co)
                        BGL.glRasterPos3f (mid[0], mid[1], mid[2])
                        Draw.Text         (str(e.index))

                # draw vertex tags
                if di.key('show_vtags'):
                    if obj.properties.has_key('vtags'):
                        stg       = 'stg_'+str(di.key('fem_stage'))
                        BGL.glColor3f (0.551, 1.0, 0.370)
                        for k, v in obj.properties['vtags'].iteritems():
                            vid = int(k)
                            pos = msh.verts[vid].co
                            BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                            Draw.Text         (str(v))

                # draw edge tags
                if di.key('show_etags'):
                    if obj.properties.has_key('etags'):
                        stg       = 'stg_'+str(di.key('fem_stage'))
                        p0, p1, s = di.bounding_box(msh)
                        BGL.glColor3f (0.551, 1.0, 0.370)
                        for k, v in obj.properties['etags'].iteritems():
                            eid = int(k)
                            dP  = msh.edges[eid].v2.co-msh.edges[eid].v1.co
                            pos = msh.edges[eid].v1.co + 0.60*dP
                            BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                            Draw.Text         (str(v))
                            if obj.properties.has_key(stg):
                                if obj.properties[stg].has_key('ebrys'):
                                    pty = obj.properties['pty']
                                    for m, n in obj.properties[stg]['ebrys'].iteritems():
                                        if int(n[0])==v:
                                            val = sgn(n[2])*0.1*s
                                            pos = [msh.edges[eid].v1.co + c*dP for c in [0.0,0.25,0.5,0.75,1.0]]
                                            dfv = di.key('pty2Fdfv')[pty][n[1]]
                                            if dfv=='qy':
                                                L = 0.05*s # arrow length
                                                h = 0.4*L  # tip height
                                                if n[2]<0:
                                                    for p in pos: draw_arrow_2d (p[0],p[1]+L, 0.0,-L, h)
                                                else:
                                                    for p in pos: draw_arrow_2d (p[0],p[1], 0.0,L, h)
                                            if dfv=='qx':
                                                for p in pos: draw_arrow_2d (p[0],p[1], val,0.0, 0.05*s)
                                            if dfv=='ux':
                                                for p in pos: draw_xfix_2d  (p[0],p[1], 0.01*s, 0.005*s)
                                            if dfv=='uy':
                                                for p in pos: draw_yfix_2d  (p[0],p[1], 0.01*s, 0.005*s)
                                            break

                # draw face tags
                if di.key('show_ftags'):
                    if obj.properties.has_key('ftags'):
                        for k, v in obj.properties['ftags'].iteritems():
                            vids = [int(vid) for vid in k.split('_')]
                            cen  = msh.verts[vids[0]].co/len(vids)
                            for i in range(1,len(vids)): cen += msh.verts[vids[i]].co/len(vids)
                            BGL.glColor3f     (0.0, 0.0, 0.0)
                            BGL.glRasterPos3f (cen[0], cen[1], cen[2])
                            Draw.Text         (str(v))

                # draw block IDs
                if di.key('show_blks'):
                    if obj.properties.has_key('blks'):
                        BGL.glColor3f (0.930, 0.830, 0.810)
                        for k, v in obj.properties['blks'].iteritems():
                            neds = int(v[17])
                            cen  = 0.5*(msh.edges[int(v[18])].v1.co+msh.edges[int(v[18])].v2.co)/neds
                            for i in range(19,19+neds-1):
                                cen += 0.5*(msh.edges[int(v[i])].v1.co+msh.edges[int(v[i])].v2.co)/neds
                            BGL.glRasterPos3f (cen[0], cen[1], cen[2])
                            Draw.Text         ('blk:%d:%d'%(int(k),int(v[0])))

                # draw local axes
                if di.key('show_axes'):
                    if obj.properties.has_key('blks'):
                        for k, v in obj.properties['blks'].iteritems():
                            eids = [int(v[1]), int(v[2]), int(v[3])]
                            clrs = [(1.0,0.1,0.1), (0.1,1.0,0.1), (0.1,0.1,1.0)]
                            for i, eid in enumerate(eids):
                                if eid>=0:
                                    ed = msh.edges[eid]
                                    BGL.glColor3f  (clrs[i][0], clrs[i][1], clrs[i][2])
                                    BGL.glBegin    (BGL.GL_LINES)
                                    BGL.glVertex3f (ed.v1.co[0], ed.v1.co[1], ed.v1.co[2])
                                    BGL.glVertex3f (ed.v2.co[0], ed.v2.co[1], ed.v2.co[2])
                                    BGL.glEnd      ()
                            origin, xp, yp, zp = int(v[4]), int(v[5]), int(v[6]), int(v[7])
                            if origin>=0:
                                BGL.glColor3f     (0.0, 0.0, 0.0)
                                BGL.glRasterPos3f (msh.verts[origin].co[0], msh.verts[origin].co[1], msh.verts[origin].co[2])
                                Draw.Text         ('o')
                            if xp>=0:
                                BGL.glColor3f     (0.0, 0.0, 0.0)
                                BGL.glRasterPos3f (msh.verts[xp].co[0], msh.verts[xp].co[1], msh.verts[xp].co[2])
                                Draw.Text         ('+')
                            if yp>=0:
                                BGL.glColor3f     (0.0, 0.0, 0.0)
                                BGL.glRasterPos3f (msh.verts[yp].co[0], msh.verts[yp].co[1], msh.verts[yp].co[2])
                                Draw.Text         ('+')
                            if zp>=0:
                                BGL.glColor3f     (0.0, 0.0, 0.0)
                                BGL.glRasterPos3f (msh.verts[zp].co[0], msh.verts[zp].co[1], msh.verts[zp].co[2])
                                Draw.Text         ('+')

                # draw regions and holes
                if di.key('show_regs'):
                    if obj.properties.has_key('regs'):
                        BGL.glColor3f (0.054, 0.675, 1.0)
                        for k, v in obj.properties['regs'].iteritems():
                            BGL.glRasterPos3f (v[2],v[3],v[4])
                            Draw.Text         ('reg:%d:%d'%(int(k),int(v[0])))
                    if obj.properties.has_key('hols'):
                        BGL.glColor3f (0.054, 0.675, 1.0)
                        for k, v in obj.properties['hols'].iteritems():
                            BGL.glRasterPos3f (v[0],v[1],v[2])
                            Draw.Text         ('hol:%d'%int(k))

                # draw elements information
                if di.key('show_elems'):
                    if obj.properties.has_key('elems'):
                        BGL.glColor3f (0.0, 1.0, 0.0)
                        for k, v in obj.properties['elems'].iteritems():
                            BGL.glRasterPos3f (v[2], v[3], v[4])
                            Draw.Text         (k+'('+str(int(v[0]))+')')

                # draw linear elements
                if obj.properties.has_key('lines') and obj.properties.has_key('msh_type'):
                    msh_obj = di.get_msh_obj (obj, False)
                    if msh_obj!=None:
                        msh = msh_obj.getData(mesh=1)
                        ori = [v for v in msh.verts] # create a copy before transforming to global coordinates
                        msh.transform (obj.matrix)   # transform to global coordinates
                        BGL.glColor3f (0.483, 0.709, 0.257)
                        for k, v in obj.properties['lines'].iteritems():
                            pos = msh.verts[v[1]].co + 0.40*(msh.verts[v[2]].co-msh.verts[v[1]].co)
                            # text
                            BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                            Draw.Text         (str(v[0]))
                            # line
                            BGL.glBegin    (BGL.GL_LINES)
                            BGL.glVertex3f (msh.verts[v[1]].co[0], msh.verts[v[1]].co[1], msh.verts[v[1]].co[2])
                            BGL.glVertex3f (msh.verts[v[2]].co[0], msh.verts[v[2]].co[1], msh.verts[v[2]].co[2])
                            BGL.glEnd      ()

                # restore mesh to local coordinates
                msh.verts = ori
