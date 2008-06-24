#!BPY

import Blender
from Blender import BGL, Draw, Window
import bpy
import msys_dict as di

dict = di.load_dict()
if dict['display_props']:
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
            msh = obj.getData(mesh=1)
            ori = msh.verts[:]
            msh.transform(obj.matrix)
            if dict['vertex_ids']:
                for v in msh.verts:
                    BGL.glColor3f     (1.0, 0.5, 0.0)
                    BGL.glRasterPos3f (v.co[0], v.co[1], v.co[2])
                    Draw.Text         ('%d'% v.index)
            if dict['edge_ids']:
                for e in msh.edges:
                    mid = 0.5*(e.v1.co+e.v2.co)
                    BGL.glColor3f     (1.0, 1.0, 1.0)
                    BGL.glRasterPos3f (mid[0], mid[1], mid[2])
                    Draw.Text         ('%d'% e.index)
            if dict['face_ids']:
                for f in msh.faces:
                    cen = f.cent
                    BGL.glColor3f     (1.0, 0.1, 0.2)
                    BGL.glRasterPos3f (cen[0], cen[1], cen[2])
                    Draw.Text         ('%d'% f.index)
            if len(obj.getAllProperties())>0:
                pro = di.get_all_props(obj)
                if pro['edge_x'].data>=0 and pro['edge_y'].data>=0:
                    ex  = msh.edges[pro['edge_x'].data]
                    ey  = msh.edges[pro['edge_y'].data]
                    BGL.glColor3f  (1.0, 0.0, 0.0)
                    BGL.glBegin    (BGL.GL_LINES)
                    BGL.glVertex3f (ex.v1.co[0], ex.v1.co[1], ex.v1.co[2])
                    BGL.glVertex3f (ex.v2.co[0], ex.v2.co[1], ex.v2.co[2])
                    BGL.glEnd      ()
                    BGL.glColor3f  (0.1, 1.0, 0.1)
                    BGL.glBegin    (BGL.GL_LINES)
                    BGL.glVertex3f (ey.v1.co[0], ey.v1.co[1], ey.v1.co[2])
                    BGL.glVertex3f (ey.v2.co[0], ey.v2.co[1], ey.v2.co[2])
                    BGL.glEnd      ()
                    if dict['disp_ndivs']:
                        ex_pos = ex.v1.co + 0.40*(ex.v2.co-ex.v1.co)
                        ey_pos = ey.v1.co + 0.40*(ey.v2.co-ey.v1.co)
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (ex_pos[0], ex_pos[1], ex_pos[2])
                        Draw.Text         ('%d'% pro['ndiv_x'].data)
                        BGL.glRasterPos3f (ey_pos[0], ey_pos[1], ey_pos[2])
                        Draw.Text         ('%d'% pro['ndiv_y'].data)
            msh.verts = ori

    # Restore
    BGL.glPopMatrix ()
    BGL.glDisable   (Blender.BGL.GL_DEPTH_TEST)
