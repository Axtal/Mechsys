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
def gen_frame_mesh(txt=None):
    Blender.Window.WaitCursor(1)

    # get selected object and mesh
    edm, obj, msh = di.get_msh()

    # 3D mesh?
    is3d = obj.properties['is3d']

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

    if txt!=None:
        if is3d: txt.write ('mesh = ms.mesh_generic(True) # True=>3D\n')
        else:    txt.write ('mesh = ms.mesh_generic(False) # False=>2D\n')

        # vertices
        info = '# VertIdx, OnBry==True, X, Y, Z'
        txt.write('mesh.set_nverts    (%d) # set number of vertices\n'%len(msh.verts))
        for i, v in enumerate(msh.verts):
            if is3d: txt.write('mesh.set_vert      (%d, True, %g,%g,%g) %s\n'%(i, v.co[0], v.co[1], v.co[2], info))
            else:    txt.write('mesh.set_vert      (%d, True, %g,%g) %s\n'   %(i, v.co[0], v.co[1],          info))
            info = ''

        # elements
        inf1 = '# ElemIdx, ETag, OnBry==True, VTK_LINE==3'
        inf2 = '# ElemIdx, Vert1, Vert2'
        inf3 = '# ElemIdx, LocalEdgeID==0, Tag'
        txt.write('mesh.set_nelems    (%d) # set number of elements\n'%len(msh.edges))
        for i, e in enumerate(msh.edges):
            key = (e.v1.index, e.v2.index)
            if not key in etags: raise Exception('All edges must have a edge tag')
            etag = etags[key]
            txt.write('mesh.set_elem      (%d,%d,True,3) %s\n'%(i, etag,  inf1)) # 3 == VTK_LINE
            txt.write('mesh.set_elem_con  (%d,0,%d) %s\n'%(i, e.v1.index, inf2)) # 0 == local node index
            txt.write('mesh.set_elem_con  (%d,1,%d) %s\n'%(i, e.v2.index, inf2)) # 1 == local node index
            txt.write('mesh.set_elem_etag (%d,0,%d) %s\n'%(i, etags[key], inf3)) # 0 == local edge index
            inf1 = inf2 = inf3 = ''
    else:
        mesh = ms.mesh_generic(is3d)

        # vertices
        mesh.set_nverts (len(msh.verts))
        for i, v in enumerate(msh.verts):
            if is3d: mesh.set_vert (i, True, v.co[0], v.co[1], v.co[2])
            else:    mesh.set_vert (i, True, v.co[0], v.co[1])

        # elements
        mesh.set_nelems (len(msh.edges))
        for i, e in enumerate(msh.edges):
            key = (e.v1.index, e.v2.index)
            if not key in etags: raise Exception('All edges must have a edge tag')
            etag = etags[key]
            mesh.set_elem      (i,etag,True,3)    # 3 == VTK_LINE
            mesh.set_elem_con  (i, 0, e.v1.index) # 0 == local node index
            mesh.set_elem_con  (i, 1, e.v2.index) # 1 == local node index
            mesh.set_elem_etag (i, 0, etags[key]) # 0 == local edge index

    # Restore local coordinates
    msh.verts = ori
    if edm: Blender.Window.EditMode(1)

    # generate mesh
    Blender.Window.WaitCursor(0)
    if txt==None: return mesh

    
# =========================================================================== Structured mesh

@print_timing
def gen_struct_mesh(gen_script=False,txt=None):
    Blender.Window.WaitCursor(1)

    # get active object
    edm, obj, msh = di.get_msh()

    # check for blocks
    if not obj.properties.has_key('blks'): raise Exception('Please, assign blocks first')

    # 3D mesh?
    is3d = obj.properties['is3d']

    # transform vertices coordinates
    ori = msh.verts[:]         # create a copy in local coordinates
    msh.transform (obj.matrix) # transform mesh to global coordinates

    # create text for script
    if gen_script:
        if txt==None:
            txt = Blender.Text.New(obj.name+'_smesh')
            txt.write('import Blender, bpy\n')
            txt.write('import mechsys   as ms\n')
            txt.write('import msys_mesh as me\n')
        txt.write('bks = []\n')

    # generate list with blocks and face colors
    bks = []
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

        # new block
        if gen_script:
            txt.write('bks.append(ms.mesh_block())\n')
            txt.write('bks[-1].set_coords('+str(int(v[0]))+', # Tag\n')
            txt.write('                   '+str(verts)+', # Vertices\n')
            txt.write('                   '+str(edges)+', # Edges\n')
            txt.write('                   '+str(etags)+', # Edge Tags\n')
            txt.write('                   '+str(ftags)+', # Face Tags\n')
            txt.write('                   '+str(wx)+', # X-Weights\n')
            txt.write('                   '+str(wy)+', # Y-Weights\n')
            txt.write('                   '+str(wz)+', # Z-Weights\n')
            txt.write('                   '+str(origin)+', '+str(xp)+', '+str(yp)+', '+str(zp)+') # Origin, X-Plus, Y-Plus, Z-Plus\n')
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

    # generate mesh
    if gen_script:
        if is3d: txt.write('mesh = ms.mesh_structured(True) # True=>3D\n')
        else:    txt.write('mesh = ms.mesh_structured(False) # False=>2D\n')
        txt.write('mesh.set_tol    (1.0e-4)\n')
        txt.write('mesh.set_blocks (bks)\n')
        txt.write('mesh.generate   (True) # True=>WithInfo\n')
        Blender.Window.WaitCursor(0)
        return txt
    else:
        if len(bks)>0:
            mesh = ms.mesh_structured(is3d)
            mesh.set_tol    (1.0e-4)
            mesh.set_blocks (bks)
            mesh.generate   (True)
            print
            Blender.Window.WaitCursor(0)
            return mesh


# ========================================================================= Unstructured mesh

@print_timing
def gen_unstruct_mesh(gen_script=False,txt=None):
    Blender.Window.WaitCursor(1)

    # get active object
    edm, obj, msh = di.get_msh()

    # 3D mesh?
    is3d = obj.properties['is3d']

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

    # get rverts and rtags
    rverts = [] # vertices on reinforcements
    rtags  = {} # reinforcements connectivities and edge IDs
    if obj.properties.has_key('rtags'):
        for k, v in obj.properties['rtags'].iteritems():
            eid = int(k)
            rverts.append (msh.edges[eid].v1.index)
            rverts.append (msh.edges[eid].v2.index)
            rtags[(msh.edges[eid].v1.index, msh.edges[eid].v2.index)] = v[0]

    # number of vertices and segments
    nverts = len(msh.verts)-len(rverts)
    nedges = len(msh.edges)-len(rtags)

    if gen_script:
        # unstructured mesh instance
        if txt==None:
            txt = Blender.Text.New(obj.name+'_umesh')
            txt.write('import Blender, bpy\n')
            txt.write('import mechsys   as ms\n')
            txt.write('import msys_mesh as me\n')
        if is3d: txt.write('mesh = ms.mesh_unstructured(True) # True=>3D\n')
        else:    txt.write('mesh = ms.mesh_unstructured(False) # False=>2D\n')

        # set polygon
        txt.write('mesh.set_poly_size    (%d,%d,%d,%d)'%(nverts, nedges, nregs, nhols)+' # nVerts, nSegments, nRegs, nHols\n')

        # set vertices and edges
        info = ' # VertIdx, X, Y, Z'
        for v in msh.verts:
            if not v.index in rverts:
                txt.write('mesh.set_poly_point   (%d,%g,%g,%g)'%(v.index, v.co[0], v.co[1], v.co[2])+info+'\n')
                info = ''
        info = ' # SegmentIdx, VertIdx1, VertIdx2, Tag'
        for e in msh.edges:
            key = (e.v1.index, e.v2.index)
            if not key in rtags:
                if key in etags: txt.write('mesh.set_poly_segment (%d,%d,%d,%d)'%(e.index, e.v1.index, e.v2.index, etags[key])+info+'\n')
                else:            txt.write('mesh.set_poly_segment (%d,%d,%d)'   %(e.index, e.v1.index, e.v2.index)+info+'\n')
                info = ''

        # set regions and holes
        info = ' # RegIdx, Tag, MaxArea, X, Y, Z'
        if obj.properties.has_key('regs'):
            for k, v in obj.properties['regs'].iteritems():
                txt.write('mesh.set_poly_region  (%d,%d,%g,%g,%g,%g)'%(int(k), int(v[0]), v[1], v[2], v[3], v[4])+info+'\n')
                info = ''
        info = ' # HolIdx, X, Y, Z'
        if obj.properties.has_key('hols'):
            for k, v in obj.properties['hols'].iteritems():
                txt.write('mesh.set_poly_hole    (%d,%g,%g,%g)'%(int(k), v[0], v[1], v[2])+info+'\n')
                info = ''

    else:
        # unstructured mesh instance
        mesh = ms.mesh_unstructured(is3d)

        # set polygon
        mesh.set_poly_size (nverts, nedges, nregs, nhols)

        # set vertices and edges
        for v in msh.verts:
            if not v.index in rverts:
                mesh.set_poly_point (v.index, v.co[0], v.co[1], v.co[2])
        for e in msh.edges:
            key = (e.v1.index, e.v2.index)
            if not key in rtags:
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

    # generate mesh
    maxa = obj.properties['maxarea'] if obj.properties.has_key('maxarea') else -1.0
    mina = obj.properties['minang']  if obj.properties.has_key('minang')  else -1.0
    if gen_script:
        if maxa>0: txt.write('mesh.set_max_area_global  (%f)'%(maxa)+'\n')
        if mina>0: txt.write('mesh.set_min_angle_global (%f)'%(mina)+'\n')
        txt.write('mesh.generate         (True) # True=>WithInfo\n')
        Blender.Window.WaitCursor(0)
        return txt
    else:
        if maxa>0: mesh.set_max_area_global  (maxa)
        if mina>0: mesh.set_min_angle_global (mina)
        mesh.generate (True)
        print
        Blender.Window.WaitCursor(0)
        return mesh


# ====================================================================================== Draw

@print_timing
def set_elems(obj, elems):
    # delete old eatts from all stages
    sids = [k for k in obj.properties['stages']]
    di.props_del_all_fem (sids, 'eatts', False)
    # loop over element tags
    obj.properties['elems'] = elems
    temp                    = {}
    for k, v in elems.iteritems():
        tag = int(v[0])
        vtk = int(v[1])
        if not temp.has_key(tag):
            temp[tag] = True
            eatt      = di.new_eatt_props()
            eatt[0]   = tag
            eatt[1]   = di.key('vtk2ety')[vtk]
            di.props_push_new_fem (sids, 'eatts', eatt)

@print_timing
def set_ebrys(obj):
    # set edges boundaries
    if obj.properties.has_key('etags'):
        for sid in obj.properties['stages']:
            stg = 'stg_'+sid
            if not obj.properties[stg].has_key('ebrys'): obj.properties[stg]['ebrys'] = {} 
            temp = {}
            id   = 0
            for k, v in obj.properties['etags'].iteritems():
                tag = int(v[0])
                if not temp.has_key(tag):
                    temp[tag] = True
                    obj.properties[stg]['ebrys'][str(id)]    = di.new_ebry_props()
                    obj.properties[stg]['ebrys'][str(id)][0] = tag
                    id += 1

@print_timing
def set_fbrys(obj):
    # set faces boundaries
    if obj.properties['is3d'] and obj.properties.has_key('ftags'):
        for sid in obj.properties['stages']:
            stg = 'stg_'+sid
            if not obj.properties[stg].has_key('fbrys'): obj.properties[stg]['fbrys'] = {} 
            temp = {}
            id   = 0
            for k, v in obj.properties['ftags'].iteritems():
                tag = int(v[0])
                if not temp.has_key(tag):
                    temp[tag] = True
                    obj.properties[stg]['fbrys'][str(id)]    = di.new_fbry_props()
                    obj.properties[stg]['fbrys'][str(id)][0] = tag
                    obj.properties[stg]['fbrys'][str(id)][3] = 0
                    id += 1

@print_timing
def add_mesh(obj, mesh, mesh_type):
    # delete old mesh
    if obj.properties.has_key('msh_name'):
        msh_name = obj.properties['msh_name']
        old_msh  = bpy.data.objects[msh_name]
        scn      = bpy.data.scenes.active
        scn.objects.unlink (old_msh)
    else: msh_name = obj.name+'_msh'

    # add new object/mesh to Blender
    scn                        = bpy.data.scenes.active
    new_msh                    = bpy.data.meshes.new      (msh_name)
    new_obj                    = scn.objects.new (new_msh, msh_name)
    new_obj.restrictSelect     = True
    obj.properties['msh_name'] = new_obj.name
    verts                      = []
    edges                      = []
    mesh.get_verts       (verts)
    mesh.get_edges       (edges)
    new_msh.verts.extend (verts)
    new_msh.edges.extend (edges)
    new_obj.select       (0)
    obj.makeParent       ([new_obj])

    # new stage
    if not obj.properties.has_key('stages'): di.props_push_new_stage() # add to the current active object

    # set elements, edges brys, and faces brys
    elems = {}
    mesh.get_elems (elems)
    set_elems      (obj, elems)
    set_ebrys      (obj)
    set_fbrys      (obj)

    # set mesh type
    obj.properties['mesh_type'] = mesh_type

    # redraw
    Blender.Window.QRedrawAll()
