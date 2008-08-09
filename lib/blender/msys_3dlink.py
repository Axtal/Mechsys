#!BPY

import Blender
from   Blender import BGL, Draw, Window
import bpy
import msys_dict as di

dict = di.load_dict()
if dict['show_props']:
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

                    # draw edge tags
                    if dict['show_etags']:
                        for t in di.get_tags (obj, 'edge'):
                            if t[1]<0:
                                pos = msh.edges[t[0]].v1.co + 0.60*(msh.edges[t[0]].v2.co-msh.edges[t[0]].v1.co)
                                BGL.glColor3f     (0.0, 0.0, 0.0)
                                BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                                Draw.Text         (str(t[1]))

                    # draw regions
                    rgs = di.get_regs (obj)
                    for r in rgs:
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (float(r[1]), float(r[2]), float(r[3]))
                        Draw.Text         ('region_'+r[0])

                    # draw holes
                    hls = di.get_hols (obj)
                    for h in hls:
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (float(h[0]), float(h[1]), float(h[2]))
                        Draw.Text         ('hole')

                # draw face tags
                if dict['show_ftags']:
                    if obj.properties.has_key('ftags'):
                        BGL.glBlendFunc (BGL.GL_SRC_ALPHA, BGL.GL_ONE)
                        for eids in obj.properties['ftags']:
                            ids      = [int(id) for id in eids.split('_')]
                            eds, vds = di.sort_edges_and_verts (msh, ids, msh.edges[ids[0]].v1.index) # will erase ids
                            BGL.glColor4f (1.0, 1.0, 0.0, 0.3)
                            BGL.glEnable  (BGL.GL_BLEND)
                            BGL.glDisable (BGL.GL_DEPTH_TEST)
                            if len(eds)==8:
                                cen = (msh.verts[vds[1]].co + msh.verts[vds[5]].co)/2.0
                                BGL.glBegin    (BGL.GL_POLYGON)
                                BGL.glVertex3f (msh.verts[vds[7]].co[0], msh.verts[vds[7]].co[1], msh.verts[vds[7]].co[2])
                                BGL.glVertex3f (msh.verts[vds[0]].co[0], msh.verts[vds[0]].co[1], msh.verts[vds[0]].co[2])
                                BGL.glVertex3f (msh.verts[vds[1]].co[0], msh.verts[vds[1]].co[1], msh.verts[vds[1]].co[2])
                                BGL.glVertex3f (                 cen[0],                  cen[1],                  cen[2])
                                BGL.glEnd      ()
                                BGL.glBegin    (BGL.GL_POLYGON)
                                BGL.glVertex3f (msh.verts[vds[1]].co[0], msh.verts[vds[1]].co[1], msh.verts[vds[1]].co[2])
                                BGL.glVertex3f (msh.verts[vds[2]].co[0], msh.verts[vds[2]].co[1], msh.verts[vds[2]].co[2])
                                BGL.glVertex3f (msh.verts[vds[3]].co[0], msh.verts[vds[3]].co[1], msh.verts[vds[3]].co[2])
                                BGL.glVertex3f (                 cen[0],                  cen[1],                  cen[2])
                                BGL.glEnd      ()
                                BGL.glBegin    (BGL.GL_POLYGON)
                                BGL.glVertex3f (msh.verts[vds[3]].co[0], msh.verts[vds[3]].co[1], msh.verts[vds[3]].co[2])
                                BGL.glVertex3f (msh.verts[vds[4]].co[0], msh.verts[vds[4]].co[1], msh.verts[vds[4]].co[2])
                                BGL.glVertex3f (msh.verts[vds[5]].co[0], msh.verts[vds[5]].co[1], msh.verts[vds[5]].co[2])
                                BGL.glVertex3f (                 cen[0],                  cen[1],                  cen[2])
                                BGL.glEnd      ()
                                BGL.glBegin    (BGL.GL_POLYGON)
                                BGL.glVertex3f (msh.verts[vds[5]].co[0], msh.verts[vds[5]].co[1], msh.verts[vds[5]].co[2])
                                BGL.glVertex3f (msh.verts[vds[6]].co[0], msh.verts[vds[6]].co[1], msh.verts[vds[6]].co[2])
                                BGL.glVertex3f (msh.verts[vds[7]].co[0], msh.verts[vds[7]].co[1], msh.verts[vds[7]].co[2])
                                BGL.glVertex3f (                 cen[0],                  cen[1],                  cen[2])
                                BGL.glEnd      ()
                            else:
                                cen = (msh.verts[vds[0]].co + msh.verts[vds[2]].co)/2.0
                                BGL.glBegin    (BGL.GL_POLYGON)
                                BGL.glVertex3f (msh.verts[vds[0]].co[0], msh.verts[vds[0]].co[1], msh.verts[vds[0]].co[2])
                                BGL.glVertex3f (msh.verts[vds[1]].co[0], msh.verts[vds[1]].co[1], msh.verts[vds[1]].co[2])
                                BGL.glVertex3f (msh.verts[vds[2]].co[0], msh.verts[vds[2]].co[1], msh.verts[vds[2]].co[2])
                                BGL.glVertex3f (msh.verts[vds[3]].co[0], msh.verts[vds[3]].co[1], msh.verts[vds[3]].co[2])
                                BGL.glEnd      ()
                            BGL.glDisable     (BGL.GL_BLEND)
                            BGL.glEnable      (BGL.GL_DEPTH_TEST)
                            BGL.glColor3f     (0.0, 0.0, 0.0)
                            BGL.glRasterPos3f (cen[0], cen[1], cen[2])
                            Draw.Text         (str(obj.properties['ftags'][eids]))

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

    # Restore
    BGL.glPopMatrix ()
    BGL.glDisable   (Blender.BGL.GL_DEPTH_TEST)
