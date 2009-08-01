#!BPY

from Blender import *
from bpy     import *

# Get dictionary
try:    d = Registry.GetKey('gui_dict')
except: d = {'show':False}

# Draw 
if d['show']:
    scn = data.scenes.active
    obj = scn.objects.active
    if obj!=None and obj.type=='Mesh':

        # transformation matrix
        view_matrix = Window.GetPerspMatrix()
        view_buffer = [view_matrix[i][j] for i in xrange(4) for j in xrange(4)]
        view_buffer = BGL.Buffer(BGL.GL_FLOAT, 16, view_buffer)
        BGL.glLoadIdentity ()
        BGL.glMatrixMode   (BGL.GL_PROJECTION)
        BGL.glPushMatrix   ()
        BGL.glLoadMatrixf  (view_buffer)

        # get mesh and transform to global coordinates
        msh = obj.getData(mesh=1)
        ori = [v for v in msh.verts] # create a copy before transforming to global coordinates
        msh.transform(obj.matrix)

        # draw vertices IDs
        if d['showvid']:
            BGL.glColor3f (1.0, 1.0, 0.0)
            for v in msh.verts:
                BGL.glRasterPos3f (v.co[0], v.co[1], v.co[2])
                Draw.Text         ('%d'%(v.index))

        # draw edges IDs
        if d['showeid']:
            BGL.glColor3f (1.0, 1.0, 1.0)
            for e in msh.edges:
                mid = 0.5*(e.v1.co+e.v2.co)
                BGL.glRasterPos3f (mid[0], mid[1], mid[2])
                Draw.Text         ('%d'%(e.index))

        # draw solids IDs
        if d['showsid']:
            BGL.glColor3f (0.7, 0.0, 0.0)
            for f in msh.faces:
                BGL.glRasterPos3f (f.cent[0], f.cent[1], f.cent[2])
                Draw.Text         ('%d'%(f.index))

        # draw vertices tags
        if d['showvta']:
            if obj.properties.has_key('vtags'):
                BGL.glColor3f (0.8, 0.8, 0.70)
                for k, v in obj.properties['vtags'].iteritems():
                    vid = int(k)
                    pos = msh.verts[vid].co
                    BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                    Draw.Text         ('    %d'%(v))

        # draw edges tags
        if d['showeta']:
            if obj.properties.has_key('etags'):
                BGL.glColor3f (0.551, 1.0, 0.370)
                for k, v in obj.properties['etags'].iteritems():
                    eid = int(k)
                    dP  = msh.edges[eid].v2.co-msh.edges[eid].v1.co
                    pos = msh.edges[eid].v1.co + 0.60*dP
                    BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                    Draw.Text         ('%d'%(v))

        # draw edges tags
        if d['showsta']:
            if obj.properties.has_key('stags'):
                BGL.glColor3f (0.9, 0.1, 0.3)
                for k, v in obj.properties['stags'].iteritems():
                    sid = int(k)
                    pos = msh.faces[sid].cent
                    BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                    Draw.Text         ('    %d'%(v))

        # restore mesh to local coordinates
        msh.verts = ori
