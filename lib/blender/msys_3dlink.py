#!BPY

import Blender
from Blender import BGL, Draw, Window
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
            msh = obj.getData(mesh=1)
            ori = msh.verts[:] # create a copy before transforming to global coordinates
            msh.transform(obj.matrix)
            if dict['show_v_ids']:
                for v in msh.verts:
                    BGL.glColor3f     (1.0, 0.5, 0.0)
                    BGL.glRasterPos3f (v.co[0], v.co[1], v.co[2])
                    Draw.Text         ('%d'% v.index)
            if dict['show_e_ids']:
                for e in msh.edges:
                    mid = 0.5*(e.v1.co+e.v2.co)
                    BGL.glColor3f     (1.0, 1.0, 1.0)
                    BGL.glRasterPos3f (mid[0], mid[1], mid[2])
                    Draw.Text         ('%d'% e.index)
            if dict['show_f_ids']:
                for f in msh.faces:
                    BGL.glColor3f     (1.0, 0.1, 0.2)
                    BGL.glRasterPos3f (f.cent[0], f.cent[1], f.cent[2])
                    Draw.Text         ('%d'% f.index)
            if len(obj.getAllProperties())>0:
                # draw local axes
                p  = di.get_local_axes_props(obj)
                ex = msh.edges[p['edge_x'].data]
                ey = msh.edges[p['edge_y'].data]
                if p['edge_x'].data>=0 and p['edge_y'].data>=0:
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
                    if dict['show_ndivs']:
                        pp     = di.get_ndivs_props(obj)
                        ex_pos = ex.v1.co + 0.40*(ex.v2.co-ex.v1.co)
                        ey_pos = ey.v1.co + 0.40*(ey.v2.co-ey.v1.co)
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (ex_pos[0], ex_pos[1], ex_pos[2])
                        Draw.Text         ('%d'% pp['ndiv_x'].data)
                        BGL.glRasterPos3f (ey_pos[0], ey_pos[1], ey_pos[2])
                        Draw.Text         ('%d'% pp['ndiv_y'].data)
                # draw vertex boundary marks
                vp = di.get_v_bry_props(obj)
                if len(vp['v_bry_ids'].data)>0:
                    str_ids = vp['v_bry_ids'].data.split()
                    str_mks = vp['v_bry_mks'].data.split()
                    for i, si in enumerate(str_ids):
                        iv = int(si)
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (msh.verts[iv].co[0], msh.verts[iv].co[1], msh.verts[iv].co[2])
                        Draw.Text         (str_mks[i])
                # draw edge boundary marks
                ep = di.get_e_bry_props(obj)
                if len(ep['e_bry_ids'].data)>0:
                    str_ids = ep['e_bry_ids'].data.split()
                    str_mks = ep['e_bry_mks'].data.split()
                    for i, si in enumerate(str_ids):
                        ie  = int(si)
                        pos = msh.edges[ie].v1.co + 0.60*(msh.edges[ie].v2.co-msh.edges[ie].v1.co)
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                        Draw.Text         (str_mks[i])
                # draw face boundary marks
                fp = di.get_f_bry_props(obj)
                if len(fp['f_bry_ids'].data)>0:
                    str_ids = fp['f_bry_ids'].data.split()
                    str_mks = fp['f_bry_mks'].data.split()
                    for i, si in enumerate(str_ids):
                        ifa = int(si)
                        BGL.glColor3f     (0.0, 0.0, 0.0)
                        BGL.glRasterPos3f (msh.faces[ifa].cent[0], msh.faces[ifa].cent[1], msh.faces[ifa].cent[2])
                        Draw.Text         (str_mks[i])
            msh.verts = ori # restore copy

    # Restore
    BGL.glPopMatrix ()
    BGL.glDisable   (Blender.BGL.GL_DEPTH_TEST)
