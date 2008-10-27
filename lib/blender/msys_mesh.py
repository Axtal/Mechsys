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
import time
import Blender
import bpy
import mechsys as ms
import msys_dict as di


def print_timing(func):
    def wrapper(*arg):
        t1  = time.time ()
        res = func      (*arg)
        t2  = time.time ()
        print '[1;34mMechSys[0m: %s took [1;31m%f[0m [1;32mseconds[0m' % (func.func_name, (t2-t1))
        return res
    return wrapper


def get_list_etags(obj,msh):
    # edge tags
    etags = {}
    if obj.properties.has_key('etags'):
        for k, v in obj.properties['etags'].iteritems():
            eid = int(k)
            etags[(msh.edges[eid].v1.index, msh.edges[eid].v2.index)] = v[0]
    return etags

def get_list_ftags_fclrs(obj,msh):
    # face tags and colors
    ftags = {}
    fclrs = {}
    if obj.properties.has_key('ftags'):
        for k, v in obj.properties['ftags'].iteritems():
            # tags
            ids = [int(id) for id in k.split('_')]
            ftags[tuple(ids)] = v[0]
            # colors
            fclrs[v[0]] = v[1]
    return ftags, fclrs


# =========================================================================== Structured mesh

def gen_struct_mesh():
    # get objects
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    edm = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0)

    # generate blocks
    bks   = []
    fclrs = {} # face colors {tag:color}
    for obj in obs:
        if obj!=None and obj.type=='Mesh':
            # get mesh
            if len(obj.getAllProperties())==0: raise Exception('Please, assign all mesh properties to this object(%s) first' % obj.name)
            msh = obj.getData(mesh=1)

            # set block
            origin, x_plus, y_plus, z_plus = di.get_local_system (obj)
            if origin>-1:
                # vertices coordinates
                ori = msh.verts[:]                                       # create a copy in local coordinates
                msh.transform (obj.matrix)                               # transform mesh to global coordinates
                verts = [(v.co[0], v.co[1], v.co[2]) for v in msh.verts] # list of tuples
                msh.verts = ori                                          # restore local coordinates

                # number of divisions
                nx = di.get_ndiv (obj, 'x')
                ny = di.get_ndiv (obj, 'y')
                nz = di.get_ndiv (obj, 'z')
                if nx<1: nx = 1
                if ny<1: ny = 1
                if nz<1: nz = 1
                ax = float(di.get_acoef (obj, 'x'))
                ay = float(di.get_acoef (obj, 'y'))
                az = float(di.get_acoef (obj, 'z'))
                if di.get_nonlin (obj, 'x')==0: wx = [1.0+ax*float(i)  for i in range(nx)]
                else:                           wx = [float(i+1.0)**ax for i in range(nx)]
                if di.get_nonlin (obj, 'y')==0: wy = [1.0+ay*float(i)  for i in range(ny)]
                else:                           wy = [float(i+1.0)**ay for i in range(ny)]
                if di.get_nonlin (obj, 'z')==0: wz = [1.0+az*float(i)  for i in range(nz)]
                else:                           wz = [float(i+1.0)**az for i in range(nz)]

                # edges
                edges = [(ed.v1.index, ed.v2.index) for ed in msh.edges]

                # edge tags
                etags = get_list_etags(obj,msh)

                # face tags and colors
                ftags, fclrs = get_list_ftags_fclrs(obj,msh)

                # new block
                bks.append(ms.mesh_block());
                bks[-1].set_coords (di.get_btag(obj),               # tag to be replicated to all elements
                                    verts,                          # vertices' coordinates
                                    edges,                          # edges
                                    etags,                          # edge tags
                                    ftags,                          # face tags
                                    wx, wy, wz,                     # weigths x, y, and z
                                    origin, x_plus, y_plus, z_plus) # Origin, XPlus, YPlus, None
            else:
                raise Exception('Please, define local axes first (obj=%s)' % obj.name)

    # generate mesh and draw results
    if len(bks)>0:
        # set cursor
        Blender.Window.WaitCursor(1)

        # generate mesh and write VTU file for ParaView
        mms = ms.mesh_structured(1.0e-4)
        ne  = mms.generate (bks)
        fn  = Blender.sys.makename(ext='_MESH_'+obj.name+'.vtu')
        mms.write_vtu (fn)
        print '[1;34mMechSys[0m: [1;33m%d[0m elements generated' % ne
        print '[1;34mMechSys[0m: File <[1;33m%s[0m> created' % fn

        # draw generated mesh
        add_mesh (mms, fclrs)

        # restore cursor
        Blender.Window.WaitCursor(0)


# ========================================================================= Unstructured mesh

@print_timing
def gen_unstruct_mesh():
    # get objects
    scn = bpy.data.scenes.active
    obj = scn.objects.active
    edm = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0)

    # set input polygon
    if obj!=None and obj.type=='Mesh':
        Blender.Window.WaitCursor(1)
        msh = obj.getData(mesh=1)

        # Transform mesh to global coordinates
        ori = msh.verts[:] # create a copy in local coordinates
        msh.transform (obj.matrix)

        # Set polygon
        mu  = ms.mesh_unstructured()
        ets = get_list_etags(obj,msh)
        rgs = di.get_regs  (obj)
        hls = di.get_hols  (obj)
        #mu.set_3d(0)
        mu.set_poly_size (len(msh.verts), len(msh.edges), len(rgs), len(hls))
        for v in msh.verts:
            mu.set_poly_point (v.index, v.co[0], v.co[1], v.co[2])
        for e in msh.edges:
            key = (e.v1.index, e.v2.index)
            if key in ets: mu.set_poly_segment (e.index, e.v1.index, e.v2.index, ets[key])
            else:          mu.set_poly_segment (e.index, e.v1.index, e.v2.index)
        for i, r in enumerate(rgs):
            print r
            #                   i      Tag      MaxArea          X           Y            Z
            mu.set_poly_region (i, int(r[0]), float(r[1]), float(r[2]), float(r[3]), float(r[4]))
        for i, h in enumerate(hls):
            mu.set_poly_hole (i, float(h[0]), float(h[1]), float(h[2]))

        # Restore local coordinates
        msh.verts = ori

        # Generate
        maxarea  = di.get_maxarea  (obj)
        minangle = di.get_minangle (obj)
        mu.generate      (float(maxarea), float(minangle))
        add_mesh (mu, {})
        Blender.Window.WaitCursor(0)


# ====================================================================================== Draw

@print_timing
def set_elems(obj, nelems, elems):
    obj.properties['nelems'] = nelems
    obj.properties['elems']  = elems
    obj.properties['eatts']  = {}
    for t in obj.properties['elems']['tags']:
        tag = str(t)
        if not obj.properties['eatts'].has_key(tag):
            obj.properties['eatts'][tag] = '0 0 E=200_nu=0.2 Sx=0_Sy=0_Sz=0_Sxy=0' # ElemType Model Prms Inis


@print_timing
def set_etags(obj, msh, etags):
    obj.properties['etags'] = {}
    obj.properties['ebrys'] = {}
    for et in etags:
        edge_id = msh.findEdges (et[0], et[1])
        obj.properties['etags'][str(edge_id)] = [etags[et], 0] # tag, type
        if not obj.properties['ebrys'].has_key(str(etags[et])):
            obj.properties['ebrys'][str(etags[et])] = [0, 0.0] # DOFVar==ux,uy,..., Val


@print_timing
def set_ftags(obj, msh, ftags, fclrs):
    obj.properties['ftags'] = {}
    obj.properties['fbrys'] = {}
    for ft in ftags:
        eids = ''
        vids = ft.split('_')
        for i, pair in enumerate(vids):
            vs = pair.split(',')
            edge_id = msh.findEdges (int(vs[0]), int(vs[1]))
            if i>0: eids += '_'+str(edge_id)
            else:   eids +=     str(edge_id)
        print fclrs
        #obj.properties['ftags'][eids] = [ftags[ft], fclrs[ftags[ft]]] # tag, color
        #if not obj.properties['fbrys'].has_key(str(ftags[ft])):
            #obj.properties['fbrys'][str(ftags[ft])] = [0, 0.0, fclrs[ftags[ft]]] # DOFVar==uz..., Val, Clr


@print_timing
def add_mesh(mms, fclrs):
    # get vertices and edges
    verts = []
    edges = []
    mms.get_verts (verts)
    mms.get_edges (edges)

    # add new mesh to Blender
    key     = di.get_key()
    scn     = bpy.data.scenes.active
    new_msh = bpy.data.meshes.new      (key+'_structured')
    new_obj = scn.objects.new (new_msh, key+'_structured')
    new_msh.verts.extend (verts)
    new_msh.edges.extend (edges)
    new_obj.select       (1)
    print '[1;34mMechSys[0m: Mesh extended'

    # Vertices on boundary
    verts_bry = []
    mms.get_verts_bry (verts_bry)
    new_obj.properties['verts_bry'] = verts_bry;

    # 2D or 3D mesh ?
    new_obj.properties['is3d'] = mms.is_3d()

    # set elements
    elems  = {}
    nelems = mms.get_elems (elems)
    set_elems (new_obj, nelems, elems)

    # set etags
    etags = {}
    mms.get_etags (etags)
    set_etags (new_obj, new_msh, etags)

    # set ftags
    ftags = {}
    mms.get_ftags (ftags)
    set_ftags (new_obj, new_msh, ftags, fclrs)

    # redraw
    Blender.Window.QRedrawAll()
