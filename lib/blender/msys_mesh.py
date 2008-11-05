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


# =========================================================================== Structured mesh

@print_timing
def gen_struct_mesh(gen_script=False,txt=None):
    # get selected object and mesh
    edm, obj, msh = di.get_msh()
    if not obj.properties.has_key('blks'):
        if edm: Blender.Window.EditMode(1)
        raise Exception('Please, assign blocks first')

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

    # generate mesh and draw results
    if gen_script:
        if is3d: txt.write('msm = ms.mesh_structured(True)\n')
        else:    txt.write('msm = ms.mesh_structured(False)\n')
        txt.write('nel = msm.generate(bks,1.0e-4)\n')
        txt.write('me.add_mesh(msm)\n')
    else:
        if len(bks)>0:
            Blender.Window.WaitCursor(1)
            msm = ms.mesh_structured(is3d)
            nel = msm.generate (bks,1.0e-4)
            print '[1;34mMechSys[0m: [1;33m%d[0m elements generated' % nel
            add_mesh (msm, fclrs)
            Blender.Window.WaitCursor(0)


# ========================================================================= Unstructured mesh

@print_timing
def gen_unstruct_mesh(gen_script=False,txt=None):
    # get selected object and mesh
    edm, obj, msh = di.get_msh()

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
        if is3d: txt.write('msm = ms.mesh_unstructured(True)\n')
        else:    txt.write('msm = ms.mesh_unstructured(False)\n')

        # set polygon
        txt.write('msm.set_poly_size    (%d,%d,%d,%d)'%(len(msh.verts), len(msh.edges), nregs, nhols)+'\n')

        # set vertices and edges
        for v in msh.verts:
            txt.write('msm.set_poly_point   (%d,%f,%f,%f)'%(v.index, v.co[0], v.co[1], v.co[2])+'\n')
        for e in msh.edges:
            key = (e.v1.index, e.v2.index)
            if key in etags: txt.write('msm.set_poly_segment (%d,%d,%d,%d)'%(e.index, e.v1.index, e.v2.index, etags[key])+'\n')
            else:            txt.write('msm.set_poly_segment (%d,%d,%d)'   %(e.index, e.v1.index, e.v2.index)+'\n')

        # set regions and holes
        if obj.properties.has_key('regs'):
            for k, v in obj.properties['regs'].iteritems():
                txt.write('msm.set_poly_region  (%d,%d,%f,%f,%f,%f)'%(int(k), int(v[0]), v[1], v[2], v[3], v[4])+'\n')
        if obj.properties.has_key('hols'):
            for k, v in obj.properties['hols'].iteritems():
                txt.write('msm.set_poly_hole    (%d,%f,%f,%f)'%(int(k), v[0], v[1], v[2])+'\n')

    else:
        # unstructured mesh instance
        msm = ms.mesh_unstructured(is3d)

        # set polygon
        msm.set_poly_size (len(msh.verts), len(msh.edges), nregs, nhols)

        # set vertices and edges
        for v in msh.verts:
            msm.set_poly_point (v.index, v.co[0], v.co[1], v.co[2])
        for e in msh.edges:
            key = (e.v1.index, e.v2.index)
            if key in etags: msm.set_poly_segment (e.index, e.v1.index, e.v2.index, etags[key])
            else:            msm.set_poly_segment (e.index, e.v1.index, e.v2.index)

        # set regions and holes
        if obj.properties.has_key('regs'):
            for k, v in obj.properties['regs'].iteritems():
                msm.set_poly_region (int(k), int(v[0]), v[1], v[2], v[3], v[4])
        if obj.properties.has_key('hols'):
            for k, v in obj.properties['hols'].iteritems():
                msm.set_poly_hole (int(k), v[0], v[1], v[2])

    # Restore local coordinates
    msh.verts = ori

    # generate mesh and draw results
    maxa = obj.properties['maxarea'] if obj.properties.has_key('maxarea') else -1.0
    mina = obj.properties['minang']  if obj.properties.has_key('minang')  else -1.0
    if gen_script:
        txt.write('nel = msm.generate   (%f,%f)'%(maxa, mina)+'\n')
        txt.write('me.add_mesh(msm)\n')
    else:
        Blender.Window.WaitCursor(1)
        ne = msm.generate (maxa, mina)
        print '[1;34mMechSys[0m: [1;33m%d[0m elements generated' % ne
        add_mesh (msm)
        Blender.Window.WaitCursor(0)


# ====================================================================================== Draw

@print_timing
def set_elems(obj, nelems, elems):
    obj.properties['nelems'] = nelems
    obj.properties['elems']  = elems
    obj.properties['eatts']  = {}
    temp                     = {}
    id                       = 0
    for i, tag in enumerate(obj.properties['elems']['tags']):
        if not temp.has_key(tag):
            temp[tag] = True
            obj.properties['eatts'][str(id)]    = di.new_eatt_props()
            obj.properties['eatts'][str(id)][0] = tag
            obj.properties['eatts'][str(id)][1] = di.key('vtk2ety')[obj.properties['elems']['vtks'][i]]
            id += 1


@print_timing
def set_etags(obj, msh, etags):
    obj.properties['etags'] = {}
    obj.properties['ebrys'] = {}
    temp                    = {}
    id                      = 0
    for k, v in etags.iteritems():
        tag = v
        eid = msh.findEdges (k[0], k[1])
        obj.properties['etags'][str(eid)] = [tag, 0] # tag, type
        if not temp.has_key(tag):
            temp[tag] = True
            obj.properties['ebrys'][str(id)]    = di.new_ebry_props()
            obj.properties['ebrys'][str(id)][0] = tag
            id += 1


@print_timing
def set_ftags(obj, msh, ftags, fclrs):
    obj.properties['ftags'] = {}
    obj.properties['fbrys'] = {}
    temp                    = {}
    id                      = 0
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
            obj.properties['fbrys'][str(id)]    = di.new_fbry_props()
            obj.properties['fbrys'][str(id)][0] = tag
            obj.properties['fbrys'][str(id)][3] = fclrs[tag]
            id += 1


@print_timing
def add_mesh(msm, fclrs={}):
    # get vertices and edges
    verts = []
    edges = []
    msm.get_verts (verts)
    msm.get_edges (edges)

    # add new mesh to Blender
    key     = di.get_file_key()
    scn     = bpy.data.scenes.active
    new_msh = bpy.data.meshes.new      (key+'_structured')
    new_obj = scn.objects.new (new_msh, key+'_structured')
    new_msh.verts.extend (verts)
    new_msh.edges.extend (edges)
    new_obj.select       (1)
    print '[1;34mMechSys[0m: Mesh extended'

    # Vertices on boundary
    verts_bry = []
    msm.get_verts_bry (verts_bry)
    new_obj.properties['verts_bry'] = verts_bry;

    # 2D or 3D mesh ?
    new_obj.properties['is3d'] = msm.is_3d()

    # set elements
    elems  = {}
    nelems = msm.get_elems (elems)
    set_elems (new_obj, nelems, elems)

    # set etags
    etags = {}
    msm.get_etags (etags)
    set_etags (new_obj, new_msh, etags)

    # set ftags
    ftags = {}
    msm.get_ftags (ftags)
    set_ftags (new_obj, new_msh, ftags, fclrs)

    # redraw
    Blender.Window.QRedrawAll()
