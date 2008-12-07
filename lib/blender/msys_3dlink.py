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
import msys_dict as di


# Transformation matrix
if di.key('show_props') or di.key('show_res') or di.key('show_reinfs'):
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
    edm = Blender.Window.EditMode()
    for obj in obs:
        if obj!=None and obj.type=='Mesh':

            # draw only if active layer corresponds to this object.Layer
            if Blender.Window.GetActiveLayer()==obj.Layer:

                # get mesh and transform to global coordinates
                msh = obj.getData(mesh=1)
                ori = msh.verts[:] # create a copy before transforming to global coordinates
                msh.transform(obj.matrix)

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
                    msh_ori = msh_msh.verts[:] # create a copy before transforming to global coordinates
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

                # draw edge tags
                if di.key('show_etags'):
                    if obj.properties.has_key('etags'):
                        BGL.glColor3f (0.551, 1.0, 0.370)
                        for k, v in obj.properties['etags'].iteritems():
                            eid = int(k)
                            pos = msh.edges[eid].v1.co + 0.60*(msh.edges[eid].v2.co-msh.edges[eid].v1.co)
                            BGL.glRasterPos3f (pos[0], pos[1], pos[2])
                            Draw.Text         (str(v[0]))

                # draw face tags
                if di.key('show_ftags'):
                    if obj.properties.has_key('ftags'):
                        BGL.glBlendFunc (BGL.GL_SRC_ALPHA, BGL.GL_ONE)
                        for k, v in obj.properties['ftags'].iteritems():
                            clr      = di.hex2rgb(v[1])
                            ids      = [int(id) for id in k.split('_')]
                            eds, vds = di.sort_edges_and_verts (msh, ids, msh.edges[ids[0]].v1.index) # will erase ids
                            cen      = msh.verts[vds[0]].co/len(eds)
                            for i in range(1,len(eds)):
                                cen += msh.verts[vds[i]].co/len(eds)
                            BGL.glColor4f     (clr[0], clr[1], clr[2], di.key('show_opac'))
                            if di.key('show_opac')<0.9:
                                BGL.glEnable  (BGL.GL_BLEND)
                                BGL.glDisable (BGL.GL_DEPTH_TEST)
                            BGL.glBegin       (BGL.GL_TRIANGLE_FAN)
                            BGL.glVertex3f    (cen[0], cen[1], cen[2])
                            for i in range(len(eds)):
                                BGL.glVertex3f(msh.verts[vds[i]].co[0], msh.verts[vds[i]].co[1], msh.verts[vds[i]].co[2])
                            BGL.glVertex3f    (msh.verts[vds[0]].co[0], msh.verts[vds[0]].co[1], msh.verts[vds[0]].co[2])
                            BGL.glEnd         ()
                            if di.key('show_opac')<0.9:
                                BGL.glDisable (BGL.GL_BLEND)
                                BGL.glEnable  (BGL.GL_DEPTH_TEST)
                            BGL.glColor3f     (0.0, 0.0, 0.0)
                            BGL.glRasterPos3f (cen[0], cen[1], cen[2])
                            Draw.Text         (str(v[0]))

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

                # Resore mesh to local coordinates
                msh.verts = ori


# draw reinforcement tags
if di.key('show_reinfs'):
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    for obj in obs:
        if obj!=None and obj.type=='Mesh':

            # draw only if active layer corresponds to this object.Layer
            if Blender.Window.GetActiveLayer()==obj.Layer:

                if obj.properties.has_key('reinfs'):
                    BGL.glColor3f (0.730, 0.681, 1.0)
                    for k, v in obj.properties['reinfs'].iteritems():
                        # text
                        xc = (v[1]+v[4])/2.0
                        yc = (v[2]+v[5])/2.0
                        zc = (v[3]+v[6])/2.0
                        BGL.glRasterPos3f (xc,yc,zc)
                        Draw.Text         (str(v[0]))
                        # line
                        BGL.glBegin    (BGL.GL_LINES)
                        BGL.glVertex3f (v[1],v[2],v[3])
                        BGL.glVertex3f (v[4],v[5],v[6])
                        BGL.glEnd      ()


def sgn(val):
    if val<0: return -1
    else:     return  1

# Results (visualisation)
if di.key('show_res'):
    # Draw
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    edm = Blender.Window.EditMode()
    for obj in obs:
        if obj!=None and obj.type=='Mesh':

            # find mesh type
            with_mesh = False
            mesh_type = obj.properties['mesh_type'] if obj.properties.has_key('mesh_type') else ''
            if mesh_type=='frame': with_mesh = True
            else: with_mesh = obj.properties.has_key('msh_name')

            # draw only if active layer corresponds to this object.Layer
            if Blender.Window.GetActiveLayer()==obj.Layer and obj.properties.has_key('res') and with_mesh:

                # get msh object
                if mesh_type=='frame': msh_obj = obj
                else: msh_obj = bpy.data.objects[obj.properties['msh_name']]

                # get mesh and transform to global coordinates
                msh = msh_obj.getData(mesh=1)
                ori = msh.verts[:] # create a copy before transforming to global coordinates
                msh.transform(msh_obj.matrix)

                # current stage
                s = str(di.key('res_stage'))

                # draw scalars text
                if di.key('res_show_scalar'):
                    lbl = obj.properties['res'][s]['idx2lbl'][str(di.key('res_lbl'))]
                    BGL.glColor3f (0.0, 0.0, 0.0)
                    for v in msh.verts:
                        BGL.glRasterPos3f (v.co[0], v.co[1], v.co[2])
                        Draw.Text         ('%g' % obj.properties['res'][s][lbl][v.index])

                # draw warped mesh
                if di.key('res_show_warp'):
                    ux = obj.properties['res'][s]['ux'] if obj.properties['res'][s].has_key('ux') else []
                    uy = obj.properties['res'][s]['uy'] if obj.properties['res'][s].has_key('uy') else []
                    uz = obj.properties['res'][s]['uz'] if obj.properties['res'][s].has_key('uz') else []
                    if len(ux)>0 or len(uy)>0 or len(uz)>0:
                        if len(ux)==0: ux = [0 for i in range(len(msh.verts))]
                        if len(uy)==0: uy = [0 for i in range(len(msh.verts))]
                        if len(uz)==0: uz = [0 for i in range(len(msh.verts))]
                        m = float(di.key('res_warp_scale'))
                        BGL.glColor3f (1.0, 1.0, 1.0)
                        for e in msh.edges:
                            BGL.glBegin    (BGL.GL_LINES)
                            BGL.glVertex3f (e.v1.co[0]+m*ux[e.v1.index], e.v1.co[1]+m*uy[e.v1.index], e.v1.co[2]+m*uz[e.v1.index])
                            BGL.glVertex3f (e.v2.co[0]+m*ux[e.v2.index], e.v2.co[1]+m*uy[e.v2.index], e.v2.co[2]+m*uz[e.v2.index])
                            BGL.glEnd      ()

                if di.key('res_show_extra') and len(obj.properties['res'][s]['extra'])>0:
                    clrs = [(0.276,0.276,1.0), (0.934, 0.643, 0.19), (0.69,0.81,0.57)]
                    exts = ['N', 'M', 'V']
                    idx  = di.key('res_ext')
                    ext  = exts[idx]
                    key  = 'max_'+ext
                    sca  = float(di.key('res_ext_scale'))
                    BGL.glColor3f (clrs[idx][0], clrs[idx][1], clrs[idx][2])
                    for ide in obj.properties['res'][s]['extra']:
                        maxv  = obj.properties['res'][s][key][0] if obj.properties['res'][s][key][0]>0 else 1.0
                        va    =        obj.properties['res'][s]['extra'][ide]['values']
                        co    =        obj.properties['res'][s]['extra'][ide]['coords']
                        no    = Vector(obj.properties['res'][s]['extra'][ide]['normal'])
                        m     = va[ext][0]
                        if ext=='N':
                            sf   = m*sca/maxv
                            epo1 = Vector([co['X'][0], co['Y'][0]]) + sf*no/2.0
                            epo2 = Vector([co['X'][0], co['Y'][0]]) - sf*no/2.0
                            if m<0.0: BGL.glColor3f (0.276, 0.276, 1.0)  # blue ==> + compression
                            else:     BGL.glColor3f (0.8,   0.0,   0.0)  # red  ==> - tension
                            for i, m in enumerate(va[ext]):
                                sf  = m*sca/maxv
                                sp  = Vector([co['X'][i], co['Y'][i]])
                                ep1 = sp + sf*no/2.0
                                ep2 = sp - sf*no/2.0
                                BGL.glBegin    (BGL.GL_LINES)
                                BGL.glVertex3f (sp [0], sp [1], 0.0)
                                BGL.glVertex3f (ep1[0], ep1[1], 0.0)
                                BGL.glEnd      ()
                                BGL.glBegin    (BGL.GL_LINES)
                                BGL.glVertex3f (sp [0], sp [1], 0.0)
                                BGL.glVertex3f (ep2[0], ep2[1], 0.0)
                                BGL.glEnd      ()
                                if (i>0):
                                    BGL.glBegin    (BGL.GL_LINES)
                                    BGL.glVertex3f (epo1[0], epo1[1], 0.0)
                                    BGL.glVertex3f (ep1 [0], ep1 [1], 0.0)
                                    BGL.glEnd      ()
                                    BGL.glBegin    (BGL.GL_LINES)
                                    BGL.glVertex3f (epo2[0], epo2[1], 0.0)
                                    BGL.glVertex3f (ep2 [0], ep2 [1], 0.0)
                                    BGL.glEnd      ()
                                epo1 = ep1
                                epo2 = ep2
                                if di.key('res_ext_txt'):
                                    BGL.glRasterPos3f (sp[0], sp[1], 0.0)
                                    Draw.Text ('%g' % m)
                        else:
                            c     = sgn(va['M'][1]) if ext=='V' else 1.0
                            sf    = m*sca/maxv
                            epold = Vector([co['X'][0], co['Y'][0]]) - sf*c*no # "-" ==> Moment are plotted to the tensioned side
                            nv    = len(va[ext]) # number of values
                            if ext=='V': BGL.glColor3f (0.69,  0.81,  0.57)
                            else:        BGL.glColor3f (0.934, 0.643, 0.19)
                            for i, m in enumerate(va[ext]):
                                if ext=='V':
                                    if   i==0:    c = -sgn(va['M'][1])
                                    elif i==nv-1: c = -sgn(va['M'][nv-2])
                                    else:         c = -sgn(va['M'][i])
                                else: c = 1.0
                                sf = m*sca/maxv
                                sp = Vector([co['X'][i], co['Y'][i]])
                                ep = sp - sf*c*no  # "-" ==> Moment are plotted to the tensioned side
                                BGL.glBegin    (BGL.GL_LINES)
                                BGL.glVertex3f (sp[0], sp[1], 0.0)
                                BGL.glVertex3f (ep[0], ep[1], 0.0)
                                BGL.glEnd      ()
                                if (i>0):
                                    BGL.glBegin    (BGL.GL_LINES)
                                    BGL.glVertex3f (epold[0], epold[1], 0.0)
                                    BGL.glVertex3f (ep[0], ep[1], 0.0)
                                    BGL.glEnd      ()
                                epold = ep
                                if di.key('res_ext_txt'):
                                    BGL.glRasterPos3f (ep[0], ep[1], 0.0)
                                    if ext=='V': Draw.Text ('%g' % m)
                                    else:        Draw.Text ('%g' % abs(m)) # show Moment as absolute values
                    res = obj.properties['res'][s][key]
                    m, e, x, y = res[0], res[1], res[2], res[3] # max, elem, x, y
                    BGL.glColor3f (0.0, 0.0, 0.0)
                    BGL.glRasterPos3f (x, y, 0.0)
                    Draw.Text ('%g at e%d' % (m,e))


                # Resore mesh to local coordinates
                msh.verts = ori
