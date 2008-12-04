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

# =========================================================================== Linear mesh

@print_timing
def gen_frame_mesh(gen_script=False,txt=None):
    Blender.Window.WaitCursor(1)

    # get selected object and mesh
    edm, obj, msh = di.get_msh()

    # 3D mesh?
    is3d = obj.properties['3dmesh']

    # transform vertices coordinates
    ori = msh.verts[:]         # create a copy in local coordinates
    msh.transform (obj.matrix) # transform mesh to global coordinates

    # MechSys mesh
    mesh = ms.mesh_generic(is3d)

    # etags
    etags = {}
    if obj.properties.has_key('etags'):
        for k, v in obj.properties['etags'].iteritems():
            eid = int(k)
            etags[(msh.edges[eid].v1.index, msh.edges[eid].v2.index)] = v[0]

    if gen_script:
        if txt==None:
            txt = Blender.Text.New(obj.name+'_lmesh')
            txt.write('import mechsys   as ms\n')
            txt.write('import msys_mesh as me\n\n')
        if is3d: txt.write('is3d     = True\n')
        else:    txt.write('is3d     = False\n')
        txt.write('onbry    = True\n')
        txt.write('VTK_LINE = 3\n\n')
        txt.write('mesh = ms.mesh_generic(is3d)\n')

        # vertices
        txt.write('mesh.set_nverts    (%d)\n'%len(msh.verts))
        for i, v in enumerate(msh.verts):
            if is3d: txt.write('mesh.set_vert      (%d, onbry, %g,%g,%g)\n'%(i, v.co[0], v.co[1], v.co[2]))
            else:    txt.write('mesh.set_vert      (%d, onbry, %g,%g)\n'   %(i, v.co[0], v.co[1]))

        # elements
        txt.write('mesh.set_nelems    (%d)\n'%len(msh.edges))
        for i, e in enumerate(msh.edges):
            key = (e.v1.index, e.v2.index)
            if not key in etags: raise Exception('All edges must have a edge tag')
            etag = etags[key]
            txt.write('mesh.set_elem      (%d,%d,onbry,VTK_LINE)\n'%(i,etag)) # 3 == VTK_LINE
            txt.write('mesh.set_elem_con  (%d,0,%d)\n'%(i, e.v1.index))       # 0 == local node index
            txt.write('mesh.set_elem_con  (%d,1,%d)\n'%(i, e.v2.index))       # 1 == local node index
            txt.write('mesh.set_elem_etag (%d,0,%d)\n'%(i, etags[key]))       # 0 == local edge index
    else:
        onbry    = True
        VTK_LINE = 3
        mesh     = ms.mesh_generic(is3d)

        # vertices
        mesh.set_nverts (len(msh.verts))
        for i, v in enumerate(msh.verts):
            if is3d: mesh.set_vert (i, onbry, v.co[0], v.co[1], v.co[2])
            else:    mesh.set_vert (i, onbry, v.co[0], v.co[1])

        # elements
        mesh.set_nelems (len(msh.edges))
        for i, e in enumerate(msh.edges):
            key = (e.v1.index, e.v2.index)
            if not key in etags: raise Exception('All edges must have a edge tag')
            etag = etags[key]
            mesh.set_elem      (i,etag,onbry,VTK_LINE) # 3 == VTK_LINE
            mesh.set_elem_con  (i, 0, e.v1.index)      # 0 == local node index
            mesh.set_elem_con  (i, 1, e.v2.index)      # 1 == local node index
            mesh.set_elem_etag (i, 0, etags[key])      # 0 == local edge index

    # Restore local coordinates
    msh.verts = ori
    if edm: Blender.Window.EditMode(1)

    # Add mesh
    if gen_script: txt.write('me.add_mesh(mesh, [], True)\n')
    else:          add_mesh (mesh, [], True)

    Blender.Window.WaitCursor(0)

    
# =========================================================================== Structured mesh

@print_timing
def gen_struct_mesh(gen_script=False,txt=None):
    Blender.Window.WaitCursor(1)

    # get selected objects
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    if len(obs)>2: raise Exception('Please, select at most two Mesh objects')
    if len(obs)==1 and obs[0].properties.has_key('rtags'): return # skip reinforcements

    # find main object (skip reinforcements)
    for o in obs:
        if o==None or o.type!='Mesh': raise Exception('Selected objects must be of Mesh type')
        if not o.properties.has_key('rtags'): obj = o

    # get mesh
    edm = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0)
    msh = obj.getData(mesh=1)

    # check for blocks
    if not obj.properties.has_key('blks'): raise Exception('Please, assign blocks first')

    # 3D mesh?
    is3d = obj.properties['3dmesh']

    # transform vertices coordinates
    ori = msh.verts[:]         # create a copy in local coordinates
    msh.transform (obj.matrix) # transform mesh to global coordinates

    # create text for script
    if gen_script:
        if txt==None:
            txt = Blender.Text.New(obj.name+'_smesh')
            txt.write('import mechsys   as ms\n')
            txt.write('import msys_mesh as me\n')
        txt.write('bks = []\n')

    # generate list with blocks and face colors
    bks   = []
    fclrs = {}
    for k, v in obj.properties['blks'].iteritems():
        # origin and local system
        origin, xp, yp, zp = int(v[4]), int(v[5]), int(v[6]), int(v[7])
        if origin<0:
            msh.verts = ori # restore local coordinates
            if edm: Blender.Window.EditMode(1)
            raise Exception('Please, assign local axes to this block (%d:%d)'%(int(k),int(v[0])))
        if is3d and zp<0:
            msh.verts = ori # restore local coordinates
            if edm: Blender.Window.EditMode(1)
            raise Exception('Please, define the Z-axis of this (3D) block (%d:%d)'%(int(k),int(v[0])))

        # divisions and weights
        nx,   ny,   nz   = int(v[ 8]), int(v[ 9]), int(v[10])
        ax,   ay,   az   =     v[11],      v[12],      v[13]
        linx, liny, linz = int(v[14]), int(v[15]), int(v[16])
        if linx: wx = [1.0+ax*float(i)  for i in range(nx)]
        else:    wx = [float(i+1.0)**ax for i in range(nx)]
        if liny: wy = [1.0+ay*float(i)  for i in range(ny)]
        else:    wy = [float(i+1.0)**ay for i in range(ny)]
        if linz: wz = [1.0+az*float(i)  for i in range(nz)]
        else:    wz = [float(i+1.0)**az for i in range(nz)]

        # vertices and edges
        verts = {}
        edges = []
        eids  = []
        for i in range(18,18+int(v[17])): eids.append(int(v[i]))
        for e in eids:
            v1 = msh.edges[e].v1.index
            v2 = msh.edges[e].v2.index
            if not verts.has_key(v1): verts[v1] = (msh.verts[v1].co[0], msh.verts[v1].co[1], msh.verts[v1].co[2])
            if not verts.has_key(v2): verts[v2] = (msh.verts[v2].co[0], msh.verts[v2].co[1], msh.verts[v2].co[2])
            edges.append((v1,v2))

        # etags
        etags = {}
        if obj.properties.has_key('etags'):
            for m, n in obj.properties['etags'].iteritems():
                e = int(m)
                if e in eids: etags[(msh.edges[e].v1.index, msh.edges[e].v2.index)] = n[0]

        # ftags
        ftags = {}
        if is3d and obj.properties.has_key('ftags'):
            for m, n in obj.properties['ftags'].iteritems():
                # tags
                eds = [int(eid) for eid in m.split('_')]
                face_is_in_block = True
                for i in eds:
                    if not i in eids:
                        face_is_in_block = False
                        break
                if face_is_in_block: ftags[tuple(eds)] = n[0]
                # colors
                if not fclrs.has_key(n[0]): fclrs[n[0]] = n[1]

        # new block
        if gen_script:
            txt.write('bks.append(ms.mesh_block())\n')
            txt.write('bks[-1].set_coords('+str(int(v[0]))+',\n')
            txt.write('                   '+str(verts)+',\n')
            txt.write('                   '+str(edges)+',\n')
            txt.write('                   '+str(etags)+',\n')
            txt.write('                   '+str(ftags)+',\n')
            txt.write('                   '+str(wx)+',\n')
            txt.write('                   '+str(wy)+',\n')
            txt.write('                   '+str(wz)+',\n')
            txt.write('                   '+str(origin)+',\n')
            txt.write('                   '+str(xp)+',\n')
            txt.write('                   '+str(yp)+',\n')
            txt.write('                   '+str(zp)+')\n')
        else:
            bks.append(ms.mesh_block())
            bks[-1].set_coords (int(v[0]),          # tag to be replicated to all elements
                                verts,              # vertices' coordinates
                                edges,              # edges
                                etags,              # edge tags
                                ftags,              # face tags
                                wx, wy, wz,         # weigths x, y, and z
                                origin, xp, yp, zp) # Origin, XPlus, YPlus, ZPlus

    # Restore local coordinates
    msh.verts = ori
    if edm: Blender.Window.EditMode(1)

    # generate mesh and draw results
    if gen_script:
        if is3d: txt.write('mesh = ms.mesh_structured(True)\n')
        else:    txt.write('mesh = ms.mesh_structured(False)\n')
        txt.write('face_colours = '+fclrs.__str__()+'\n')
        txt.write('mesh.set_tol    (1.0e-4)\n')
        txt.write('mesh.set_blocks (bks)\n')
        txt.write('mesh.generate   (True)\n')
        txt.write('me.add_mesh     (mesh, face_colours)\n')
    else:
        if len(bks)>0:
            mesh = ms.mesh_structured(is3d)
            mesh.set_tol    (1.0e-4)
            mesh.set_blocks (bks)
            mesh.generate   (True)
            add_mesh        (mesh, fclrs)

    Blender.Window.WaitCursor(0)


# ========================================================================= Unstructured mesh

@print_timing
def gen_unstruct_mesh(gen_script=False,txt=None):
    Blender.Window.WaitCursor(1)

    # get selected objects
    scn = bpy.data.scenes.active
    obs = scn.objects.selected
    if len(obs)>2: raise Exception('Please, select at most two Mesh objects')
    if len(obs)==1 and obs[0].properties.has_key('rtags'): return # skip reinforcements

    # find main object (skip reinforcements)
    for o in obs:
        if o==None or o.type!='Mesh': raise Exception('Selected objects must be of Mesh type')
        if not o.properties.has_key('rtags'): obj = o

    # get mesh
    edm = Blender.Window.EditMode()
    if edm: Blender.Window.EditMode(0)
    msh = obj.getData(mesh=1)

    # 3D mesh?
    is3d = obj.properties['3dmesh']

    # transform vertices coordinates
    ori = msh.verts[:]         # create a copy in local coordinates
    msh.transform (obj.matrix) # transform mesh to global coordinates

    # number of regions and holes
    nregs = len(obj.properties['regs']) if obj.properties.has_key('regs') else 0
    nhols = len(obj.properties['hols']) if obj.properties.has_key('hols') else 0

    # get etags
    etags = {}
    if obj.properties.has_key('etags'):
        for k, v in obj.properties['etags'].iteritems():
            eid = int(k)
            etags[(msh.edges[eid].v1.index, msh.edges[eid].v2.index)] = v[0]

    if gen_script:
        # unstructured mesh instance
        if txt==None:
            txt = Blender.Text.New(obj.name+'_umesh')
            txt.write('import mechsys   as ms\n')
            txt.write('import msys_mesh as me\n')
        if is3d: txt.write('mesh = ms.mesh_unstructured(True)\n')
        else:    txt.write('mesh = ms.mesh_unstructured(False)\n')

        # set polygon
        txt.write('mesh.set_poly_size    (%d,%d,%d,%d)'%(len(msh.verts), len(msh.edges), nregs, nhols)+'\n')

        # set vertices and edges
        for v in msh.verts:
            txt.write('mesh.set_poly_point   (%d,%f,%f,%f)'%(v.index, v.co[0], v.co[1], v.co[2])+'\n')
        for e in msh.edges:
            key = (e.v1.index, e.v2.index)
            if key in etags: txt.write('mesh.set_poly_segment (%d,%d,%d,%d)'%(e.index, e.v1.index, e.v2.index, etags[key])+'\n')
            else:            txt.write('mesh.set_poly_segment (%d,%d,%d)'   %(e.index, e.v1.index, e.v2.index)+'\n')

        # set regions and holes
        if obj.properties.has_key('regs'):
            for k, v in obj.properties['regs'].iteritems():
                txt.write('mesh.set_poly_region  (%d,%d,%f,%f,%f,%f)'%(int(k), int(v[0]), v[1], v[2], v[3], v[4])+'\n')
        if obj.properties.has_key('hols'):
            for k, v in obj.properties['hols'].iteritems():
                txt.write('mesh.set_poly_hole    (%d,%f,%f,%f)'%(int(k), v[0], v[1], v[2])+'\n')

    else:
        # unstructured mesh instance
        mesh = ms.mesh_unstructured(is3d)

        # set polygon
        mesh.set_poly_size (len(msh.verts), len(msh.edges), nregs, nhols)

        # set vertices and edges
        for v in msh.verts:
            mesh.set_poly_point (v.index, v.co[0], v.co[1], v.co[2])
        for e in msh.edges:
            key = (e.v1.index, e.v2.index)
            if key in etags: mesh.set_poly_segment (e.index, e.v1.index, e.v2.index, etags[key])
            else:            mesh.set_poly_segment (e.index, e.v1.index, e.v2.index)

        # set regions and holes
        if obj.properties.has_key('regs'):
            for k, v in obj.properties['regs'].iteritems():
                mesh.set_poly_region (int(k), int(v[0]), v[1], v[2], v[3], v[4])
        if obj.properties.has_key('hols'):
            for k, v in obj.properties['hols'].iteritems():
                mesh.set_poly_hole (int(k), v[0], v[1], v[2])

    # Restore local coordinates
    msh.verts = ori
    if edm: Blender.Window.EditMode(1)

    # generate mesh and draw results
    maxa = obj.properties['maxarea'] if obj.properties.has_key('maxarea') else -1.0
    mina = obj.properties['minang']  if obj.properties.has_key('minang')  else -1.0
    if gen_script:
        if maxa>0: txt.write('mesh.set_max_area_global  (%f)'%(maxa)+'\n')
        if mina>0: txt.write('mesh.set_min_angle_global (%f)'%(mina)+'\n')
        txt.write('mesh.generate             (True)\n')
        txt.write('me.add_mesh               (mesh)\n')
    else:
        if maxa>0: mesh.set_max_area_global  (maxa)
        if mina>0: mesh.set_min_angle_global (mina)
        mesh.generate (True)
        add_mesh      (mesh)

    Blender.Window.WaitCursor(0)


# ====================================================================================== Draw

@print_timing
def set_elems(obj, nelems, elems):
    stg                          = 'stg_'+str(di.key('fem_stage')) # stage
    obj.properties['nelems']     = nelems
    obj.properties['elems']      = elems
    obj.properties[stg]['eatts'] = {}
    temp                         = {}
    id                           = 0
    for i, tag in enumerate(obj.properties['elems']['tags']):
        if not temp.has_key(tag):
            temp[tag] = True
            obj.properties[stg]['eatts'][str(id)]    = di.new_eatt_props()
            obj.properties[stg]['eatts'][str(id)][0] = tag
            obj.properties[stg]['eatts'][str(id)][1] = di.key('vtk2ety')[obj.properties['elems']['vtks'][i]]
            # add new text for properties
            tid = di.props_push_new('texts', 'gam=20') # returns text_id  == tid
            obj.properties[stg]['eatts'][str(id)][3] = tid
            id += 1


@print_timing
def set_etags(obj, msh, etags):
    stg                          = 'stg_'+str(di.key('fem_stage')) # stage
    obj.properties['etags']      = {}
    obj.properties[stg]['ebrys'] = {}
    temp                         = {}
    id                           = 0
    for k, v in etags.iteritems():
        tag = v
        eid = msh.findEdges (k[0], k[1])
        obj.properties['etags'][str(eid)] = [tag, 0] # tag, type
        if not temp.has_key(tag):
            temp[tag] = True
            obj.properties[stg]['ebrys'][str(id)]    = di.new_ebry_props()
            obj.properties[stg]['ebrys'][str(id)][0] = tag
            id += 1


@print_timing
def set_ftags(obj, msh, ftags, fclrs):
    stg                          = 'stg_'+str(di.key('fem_stage')) # stage
    obj.properties['ftags']      = {}
    obj.properties[stg]['fbrys'] = {}
    temp                         = {}
    id                           = 0
    for k, v in ftags.iteritems():
        tag  = v
        vids = k.split('_')
        eids = ''
        for i, vert_pair in enumerate(vids):
            vs  = vert_pair.split(',')
            eid = msh.findEdges (int(vs[0]), int(vs[1]))
            if i>0: eids += '_'+str(eid)
            else:   eids +=     str(eid)
        obj.properties['ftags'][eids] = [tag, fclrs[tag]] # tag, colour
        if not temp.has_key(tag):
            temp[tag] = True
            obj.properties[stg]['fbrys'][str(id)]    = di.new_fbry_props()
            obj.properties[stg]['fbrys'][str(id)][0] = tag
            obj.properties[stg]['fbrys'][str(id)][3] = fclrs[tag]
            id += 1


@print_timing
def add_mesh(mesh, fclrs={}, frame=False):
    # get vertices and edges
    verts = []
    edges = []
    mesh.get_verts (verts)
    mesh.get_edges (edges)

    # add new mesh to Blender
    key     = di.get_file_key()
    scn     = bpy.data.scenes.active
    new_msh = bpy.data.meshes.new      (key+'_mesh')
    new_obj = scn.objects.new (new_msh, key+'_mesh')
    new_msh.verts.extend (verts)
    new_msh.edges.extend (edges)
    new_obj.select       (1)
    print '[1;34mMechSys[0m: Mesh extended'

    # Vertices on boundary
    verts_bry = []
    mesh.get_verts_bry (verts_bry)
    new_obj.properties['verts_bry'] = verts_bry;

    # 2D or 3D mesh ?
    new_obj.properties['3dmesh'] = mesh.is_3d()

    # Frame mesh ?
    new_obj.properties['frame'] = frame

    # new stage
    di.props_push_new_stage()

    # set elements
    elems  = {}
    nelems = mesh.get_elems (elems)
    set_elems (new_obj, nelems, elems)

    # set etags
    etags = {}
    mesh.get_etags (etags)
    set_etags (new_obj, new_msh, etags)

    # set ftags
    ftags = {}
    mesh.get_ftags (ftags)
    set_ftags (new_obj, new_msh, ftags, fclrs)

    # redraw
    Blender.Window.QRedrawAll()
