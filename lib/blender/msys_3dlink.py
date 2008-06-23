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
            for v in msh.verts:
                BGL.glColor4f     (1.0,0.5,0.0,1.0)
                BGL.glRasterPos3f (v.co[0],v.co[1],v.co[2])
                Draw.Text         ('%d'% v.index)
            for e in msh.edges:
                mid = 0.5*(e.v1.co+e.v2.co)
                BGL.glColor3f     (1,1,1)
                BGL.glRasterPos3f (mid[0],mid[1],mid[2])
                Draw.Text         ('%d'% e.index)
            if len(obj.getAllProperties())>0:
                pro = di.get_all_props(obj)
                if pro['edge_x'].data>=0 and pro['edge_y'].data>=0:
                    xyz = obj.getLocation()
                    ex  = msh.edges[pro['edge_x'].data]
                    ey  = msh.edges[pro['edge_y'].data]
                    BGL.glColor3f  (1.0, 0.0, 0.0)
                    BGL.glBegin    (BGL.GL_LINES)
                    BGL.glVertex3f (xyz[0]+ex.v1.co[0], xyz[1]+ex.v1.co[1], xyz[2]+ex.v1.co[2])
                    BGL.glVertex3f (xyz[0]+ex.v2.co[0], xyz[1]+ex.v2.co[1], xyz[2]+ex.v2.co[2])
                    BGL.glEnd      ()
                    BGL.glColor3f  (0.1, 1.0, 0.1)
                    BGL.glBegin    (BGL.GL_LINES)
                    BGL.glVertex3f (xyz[0]+ey.v1.co[0], xyz[1]+ey.v1.co[1], xyz[2]+ey.v1.co[2])
                    BGL.glVertex3f (xyz[0]+ey.v2.co[0], xyz[1]+ey.v2.co[1], xyz[2]+ey.v2.co[2])
                    BGL.glEnd      ()

    # Restore
    BGL.glPopMatrix ()
    BGL.glDisable   (Blender.BGL.GL_DEPTH_TEST)
