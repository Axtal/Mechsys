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
import bpy
import msys_dict as di


# Transformation matrix
dict = di.load_dict()
if dict['show_props'] or dict['show_results']:
    # Buffer
    view_matrix = Window.GetPerspMatrix()
    view_buffer = [view_matrix[i][j] for i in xrange(4) for j in xrange(4)]
    view_buffer = BGL.Buffer(BGL.GL_FLOAT, 16, view_buffer)

    # Initialize
    BGL.glPushMatrix   ()
    BGL.glEnable       (Blender.BGL.GL_DEPTH_TEST)
    BGL.glLoadIdentity ()
    BGL.glMatrixMode   (Blender.BGL.GL_PROJECTION)
    BGL.glLoadMatrixf  (view_buffer)


# Mesh properties
if dict['show_props']:
    # Draw
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    edm = Blender.Window.EditMode()
    for obj in obs:
        if obj!=None and obj.type=='Mesh':

            # draw only if active layer corresponds to this object.Layer
            if Blender.Window.GetActiveLayer()==obj.Layer:

                # draw block tags
                btag = di.get_btag (obj)
                if btag<0:
                    loc = obj.getLocation()
                    BGL.glColor3f     (0.0, 1.0, 0.0)
                    BGL.glRasterPos3f (loc[0], loc[1], loc[2])
                    Draw.Text         ('btag=%d'%btag)

                # get mesh and transform to global coordinates
                msh = obj.getData(mesh=1)
                ori = msh.verts[:] # create a copy before transforming to global coordinates
                msh.transform(obj.matrix)

                # draw vertices IDs
                if dict['show_v_ids']:
                    for v in msh.verts:
                        BGL.glColor3f     (1.0, 1.0, 0.0)
                        BGL.glRasterPos3f (v.co[0], v.co[1], v.co[2])
                        Draw.Text         (str(v.index))

                # draw edges IDs
                if dict['show_e_ids']:
                    for e in msh.edges:
                        mid = 0.5*(e.v1.co+e.v2.co)
                        BGL.glColor3f     (1.0, 1.0, 1.0)
                        BGL.glRasterPos3f (mid[0], mid[1], mid[2])
                        Draw.Text         (str(e.index))

                # draw faces IDs
                if dict['show_f_ids']:
                    for f in msh.faces:
                        BGL.glColor3f     (1.0, 0.1, 0.2)
                        BGL.glRasterPos3f (f.cent[0], f.cent[1], f.cent[2])
                        Draw.Text         (str(f.index))

                # if there are properties (local axes, ndivs, tags, etc.)
                if len(obj.getAllProperties())>0:

                    # draw local axes
                    if dict['show_axes']:

                        # draw local x-axis
                        ix = di.get_local_axis (obj, 'x')
                        if ix>-1:
                            ex = msh.edges[ix]
                            BGL.glColor3f  (1.0, 0.1, 0.1)
                            BGL.glBegin    (BGL.GL_LINES)
                            BGL.glVertex3f (ex.v1.co[0], ex.v1.co[1], ex.v1.co[2])
                            BGL.glVertex3f (ex.v2.co[0], ex.v2.co[1], ex.v2.co[2])
                            BGL.glEnd      ()
                            # ndivs
                            nx = di.get_ndiv (obj, 'x')
                            if nx>-1:
                                pos = ex.v1.co + 0.40*(ex.v2.co-ex.v1.co)
                                BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                                Draw.Text         (str(nx))

                        # draw local y-axis
                        iy = di.get_local_axis (obj, 'y')
                        if iy>-1:
                            ey = msh.edges[iy]
                            BGL.glColor3f  (0.1, 1.0, 0.1)
                            BGL.glBegin    (BGL.GL_LINES)
                            BGL.glVertex3f (ey.v1.co[0], ey.v1.co[1], ey.v1.co[2])
                            BGL.glVertex3f (ey.v2.co[0], ey.v2.co[1], ey.v2.co[2])
                            BGL.glEnd      ()
                            # ndivs
                            ny = di.get_ndiv (obj, 'y')
                            if ny>-1:
                                pos = ey.v1.co + 0.40*(ey.v2.co-ey.v1.co)
                                BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                                Draw.Text         (str(ny))

                        # draw local z-axis
                        iz = di.get_local_axis (obj, 'z')
                        if iz>-1:
                            ez = msh.edges[iz]
                            BGL.glColor3f  (0.1, 0.1, 1.0)
                            BGL.glBegin    (BGL.GL_LINES)
                            BGL.glVertex3f (ez.v1.co[0], ez.v1.co[1], ez.v1.co[2])
                            BGL.glVertex3f (ez.v2.co[0], ez.v2.co[1], ez.v2.co[2])
                            BGL.glEnd      ()
                            # ndivs
                            nz = di.get_ndiv (obj, 'z')
                            if nz>-1:
                                pos = ez.v1.co + 0.40*(ez.v2.co-ez.v1.co)
                                BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                                Draw.Text         (str(nz))

                        # draw local system
                        origin, x_plus, y_plus, z_plus = di.get_local_system (obj)
                        if origin>-1:
                            BGL.glColor3f     (1.0, 0.1, 0.1)
                            BGL.glRasterPos3f (msh.verts[x_plus].co[0], msh.verts[x_plus].co[1], msh.verts[x_plus].co[2])
                            Draw.Text         ('+')
                            BGL.glColor3f     (0.1, 1.0, 0.1)
                            BGL.glRasterPos3f (msh.verts[y_plus].co[0], msh.verts[y_plus].co[1], msh.verts[y_plus].co[2])
                            Draw.Text         ('+')
                            if z_plus>-1:
                                BGL.glColor3f     (0.1, 0.1, 1.0)
                                BGL.glRasterPos3f (msh.verts[z_plus].co[0], msh.verts[z_plus].co[1], msh.verts[z_plus].co[2])
                                Draw.Text         ('+')

                    # draw regions
                    rgs = di.get_regs (obj)
                    for r in rgs:
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (float(r[2]), float(r[3]), float(r[4]))
                        Draw.Text         ('region_'+r[0])

                    # draw holes
                    hls = di.get_hols (obj)
                    for h in hls:
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (float(h[0]), float(h[1]), float(h[2]))
                        Draw.Text         ('hole')

                # draw edge tags
                if dict['show_etags']:
                    if obj.properties.has_key('etags'):
                        for k, v in obj.properties['etags'].iteritems():
                            eid = int(k)
                            pos = msh.edges[eid].v1.co + 0.60*(msh.edges[eid].v2.co-msh.edges[eid].v1.co)
                            BGL.glColor3f     (0.0, 0.0, 0.0)
                            BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                            Draw.Text         (str(v[0]))

                # draw face tags
                if dict['show_ftags']:
                    if obj.properties.has_key('ftags'):
                        BGL.glBlendFunc (BGL.GL_SRC_ALPHA, BGL.GL_ONE)
                        for k, v in obj.properties['ftags'].iteritems():
                            clr      = di.hex2rgb(v[1])
                            ids      = [int(id) for id in k.split('_')]
                            eds, vds = di.sort_edges_and_verts (msh, ids, msh.edges[ids[0]].v1.index) # will erase ids
                            cen      = msh.verts[vds[0]].co/len(eds)
                            for i in range(1,len(eds)):
                                cen += msh.verts[vds[i]].co/len(eds)
                            BGL.glColor4f     (clr[0], clr[1], clr[2], dict['ftags_opac'])
                            if dict['ftags_opac']<0.9:
                                BGL.glEnable  (BGL.GL_BLEND)
                                BGL.glDisable (BGL.GL_DEPTH_TEST)
                            BGL.glBegin       (BGL.GL_TRIANGLE_FAN)
                            BGL.glVertex3f    (cen[0], cen[1], cen[2])
                            for i in range(len(eds)):
                                BGL.glVertex3f(msh.verts[vds[i]].co[0], msh.verts[vds[i]].co[1], msh.verts[vds[i]].co[2])
                            BGL.glVertex3f    (msh.verts[vds[0]].co[0], msh.verts[vds[0]].co[1], msh.verts[vds[0]].co[2])
                            BGL.glEnd         ()
                            if dict['ftags_opac']<0.9:
                                BGL.glDisable (BGL.GL_BLEND)
                                BGL.glEnable  (BGL.GL_DEPTH_TEST)
                            BGL.glColor3f     (0.0, 0.0, 0.0)
                            BGL.glRasterPos3f (cen[0], cen[1], cen[2])
                            if dict['show_tags_txt']: Draw.Text (str(v[0]))

                # draw elements information
                if dict['show_elems']:
                    try:    nelems = obj.properties['nelems']
                    except: nelems = 0
                    if nelems>0:
                        for id in obj.properties['elems']['cons']:
                            x, y, z = di.get_cg (msh, obj.properties['elems']['cons'][id], obj.properties['elems']['vtks'][int(id)])
                            BGL.glColor3f     (0.0, 1.0, 0.0)
                            BGL.glRasterPos3f (x, y, z)
                            Draw.Text         (str(id)+'('+str(obj.properties['elems']['tags'][int(id)])+')')

                # Resore mesh to local coordinates
                msh.verts = ori


# Results (visualisation)
if dict['show_results']:
    # Draw
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    edm = Blender.Window.EditMode()
    for obj in obs:
        if obj!=None and obj.type=='Mesh':

            # draw only if active layer corresponds to this object.Layer
            if Blender.Window.GetActiveLayer()==obj.Layer:

                # get mesh and transform to global coordinates
                msh = obj.getData(mesh=1)
                ori = msh.verts[:] # create a copy before transforming to global coordinates
                msh.transform(obj.matrix)

                # Post-processing
                if dict['show_scalar']:
                    key = dict['scalar_key']
                    try:
                        vals = obj.properties['scalars'][key]
                        BGL.glColor3f (0.0, 0.0, 0.0)
                        for v in msh.verts:
                            BGL.glRasterPos3f (v.co[0], v.co[1], v.co[2])
                            Draw.Text         (str(vals[v.index]))
                    except: pass

                if dict['show_warp']:
                    ux = []
                    uy = []
                    uz = []
                    try:
                        ux = obj.properties['scalars']['ux']
                        uy = obj.properties['scalars']['uy']
                        uz = obj.properties['scalars']['uz']
                    except: pass
                    if len(ux)>0 or len(uy)>0 or len(uz)>0:
                        if len(ux)<1: ux = [0 for i in range(len(msh.verts))]
                        if len(uy)<1: uy = [0 for i in range(len(msh.verts))]
                        if len(uz)<1: uz = [0 for i in range(len(msh.verts))]
                        m = float(dict['warp_scale'])
                        BGL.glColor3f (0.85, 0.85, 0.85)
                        for e in msh.edges:
                            BGL.glBegin    (BGL.GL_LINES)
                            BGL.glVertex3f (e.v1.co[0]+m*ux[e.v1.index], e.v1.co[1]+m*uy[e.v1.index], e.v1.co[2]+m*uz[e.v1.index])
                            BGL.glVertex3f (e.v2.co[0]+m*ux[e.v2.index], e.v2.co[1]+m*uy[e.v2.index], e.v2.co[2]+m*uz[e.v2.index])
                            BGL.glEnd      ()

                # Resore mesh to local coordinates
                msh.verts = ori


# Restore transformation matrix
if dict['show_props'] or dict['show_results']:
    # Restore
    BGL.glPopMatrix ()
    BGL.glDisable   (Blender.BGL.GL_DEPTH_TEST)
